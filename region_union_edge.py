#!/usr/bin/env python3
"""Union the footprints of small-field FK5 DS9 regions into polygon outlines.

Dependencies: numpy, astropy, shapely >= 2.0. Python >= 3.10.
See region_union_README.md for geometry conventions and limitations.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
import re
import sys
import warnings

import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.validation import explain_validity


AREA_SHAPES = {"circle", "ellipse", "box", "polygon", "annulus",
               "panda", "epanda", "bpanda"}
NONAREA = {"point", "line", "vector", "segment", "text", "ruler", "compass"}
FRAMES = {"fk5", "j2000", "fk4", "b1950", "icrs", "galactic", "ecliptic",
          "image", "physical", "linear", "wcs", "detector", "amplifier"}
DEG = np.pi / 180.0


class RegionError(ValueError):
    """An unsupported or invalid input; never silently discard area regions."""


@dataclass
class RegionSpec:
    shape: str
    args: list[str]
    include: bool
    source: str
    file_index: int


def number(token: str) -> float:
    value = float(token)
    if not math.isfinite(value):
        raise RegionError(f"Non-finite value: {token}")
    return value


def angular(token: str) -> float:
    """DS9 angular size/angle in degrees; bare numbers mean degrees."""
    token = token.strip().lower()
    factors = {'"': 1 / 3600, "'": 1 / 60, "d": 1, "r": 180 / np.pi}
    if token.endswith(("p", "i")):
        raise RegionError("Pixel units need a FITS WCS and are not supported")
    if token and token[-1] in factors:
        return number(token[:-1]) * factors[token[-1]]
    return number(token)


def position(ra: str, dec: str) -> tuple[float, float]:
    def parse(token, is_ra):
        token = token.lower()
        if ":" in token:
            value = Angle(token, unit=u.hourangle if is_ra else u.deg).degree
        elif "h" in token or ("d" in token and not token.endswith("d")):
            value = Angle(token).degree
        else:
            value = angular(token)
        if not np.isfinite(value):
            raise RegionError("Non-finite sky coordinate")
        return float(value)
    x, y = parse(ra, True), parse(dec, False)
    if not -90 <= y <= 90:
        raise RegionError(f"Declination outside [-90, 90]: {y}")
    return x % 360, y


def count(token: str) -> int:
    value = number(token)
    if value < 1 or int(value) != value:
        raise RegionError(f"Subdivision count must be a positive integer: {token}")
    return int(value)


def property_flag(properties: str, name: str, default: bool) -> bool:
    # Do not mistake text={include=0} or callback strings for properties.
    clean = re.sub(r'\{[^}]*\}|"[^"]*"|\'[^\']*\'', "", properties)
    matches = re.findall(rf"(?:^|\s){name}\s*=\s*([01])(?:\s|$)", clean)
    return bool(int(matches[-1])) if matches else default


def read_regions(paths, *, ignore_nonarea=False) -> list[RegionSpec]:
    """Read ordinary DS9 records. Composite boolean expressions are rejected."""
    records = []
    for file_index, path in enumerate(paths):
        frame = None
        default_include, default_fixed = True, False
        for line_no, raw in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
            source = f"{path}:{line_no}"
            try:
                if re.match(r"\s*#\s*composite\b", raw, re.I):
                    raise RegionError("Dissolve composite regions into ordinary shapes first")
                body, _, properties = raw.partition("#")
                if not body.strip():
                    continue
                if "||" in body or "&&" in body:
                    raise RegionError("Composite boolean operators are not supported")
                statements = [s.strip() for s in body.split(";") if s.strip()]
                for index, statement in enumerate(statements):
                    local = properties if index == len(statements) - 1 else ""
                    lower = statement.lower()
                    if lower.startswith("global "):
                        default_include = property_flag(statement, "include", default_include)
                        default_fixed = property_flag(statement, "fixed", default_fixed)
                        continue
                    if lower in FRAMES or re.fullmatch(r"wcs[a-z0-9]", lower):
                        frame = lower
                        continue
                    match = re.fullmatch(r"([+-]?)\s*([a-zA-Z]+)\s*(.*)", statement)
                    if not match:
                        raise RegionError(f"Cannot parse record: {statement}")
                    sign, shape, args = match.groups()
                    shape = shape.lower()
                    if shape in NONAREA and ignore_nonarea:
                        warnings.warn(f"{source}: skipped non-area shape {shape}", stacklevel=2)
                        continue
                    if shape not in AREA_SHAPES:
                        raise RegionError(f"Unsupported shape '{shape}' (area shapes: "
                                          + ", ".join(sorted(AREA_SHAPES)) + ")")
                    if frame not in {"fk5", "j2000"}:
                        raise RegionError(f"Explicit fk5/j2000 required; current frame={frame!r}")
                    if property_flag(local, "fixed", default_fixed):
                        raise RegionError("fixed=1 is display-size dependent; not supported")
                    if args.startswith("(") and args.endswith(")"):
                        args = args[1:-1]
                    elif "(" in args or ")" in args:
                        raise RegionError("Unbalanced parentheses")
                    tokens = [t for t in re.split(r"[\s,]+", args.strip()) if t]
                    included = property_flag(local, "include", default_include)
                    if sign:
                        included = sign == "+"
                    records.append(RegionSpec(shape, tokens, included, source, file_index))
            except (ValueError, TypeError) as exc:
                raise RegionError(f"{source}: {exc}") from exc
    if not records:
        raise RegionError("No area regions found")
    return records


def unit_vector(ra, dec):
    ra, dec = np.asarray(ra) * DEG, np.asarray(dec) * DEG
    return np.stack((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra),
                     np.sin(dec)), axis=-1)


def basis(ra, dec):
    a, d = ra * DEG, dec * DEG
    return (unit_vector(ra, dec), np.array([-np.sin(a), np.cos(a), 0.0]),
            np.array([-np.sin(d) * np.cos(a), -np.sin(d) * np.sin(a), np.cos(d)]))


class TangentPlane:
    """A common gnomonic plane (east=x, north=y), measured in degree units."""
    def __init__(self, ra, dec, max_field_deg=5.0):
        self.ra, self.dec = ra, dec
        self.center, self.east, self.north = basis(ra, dec)
        self.max_field_deg = max_field_deg
        self.min_dot = np.cos(max_field_deg * DEG)

    def project(self, vectors):
        vectors = np.asarray(vectors)
        z = vectors @ self.center
        if np.any(z < self.min_dot):
            raise RegionError(f"Boundary extends more than {self.max_field_deg:g} deg from "
                              "the tangent origin; split this wide field into smaller fields")
        return np.stack((vectors @ self.east / z / DEG,
                         vectors @ self.north / z / DEG), axis=-1)

    def unproject(self, xy):
        xy = np.asarray(xy) * DEG
        vec = self.center + xy[..., 0, None] * self.east + xy[..., 1, None] * self.north
        vec /= np.linalg.norm(vec, axis=-1, keepdims=True)
        ra = np.arctan2(vec[..., 1], vec[..., 0]) / DEG % 360
        dec = np.arcsin(np.clip(vec[..., 2], -1, 1)) / DEG
        return np.column_stack((ra, dec))


def choose_plane(records, origin=None, max_field_deg=5.0):
    if origin is not None:
        ra, dec = position(*origin)
    else:
        centers = []
        for record in records:
            try:
                if record.shape == "polygon":
                    if len(record.args) < 6 or len(record.args) % 2:
                        raise RegionError("polygon needs at least three RA/Dec pairs")
                    coords = [position(*record.args[i:i + 2]) for i in range(0, len(record.args), 2)]
                    v = unit_vector(*np.array(coords).T).mean(axis=0)
                else:
                    if len(record.args) < 2:
                        raise RegionError("Missing center coordinates")
                    v = unit_vector(*position(*record.args[:2]))
                norm = np.linalg.norm(v)
                if norm < 1e-8:
                    raise RegionError("Cannot find a small-field center")
                centers.append(v / norm)
            except ValueError as exc:
                raise RegionError(f"{record.source}: {exc}") from exc
        v = np.mean(centers, axis=0)
        if np.linalg.norm(v) < 1e-8:
            raise RegionError("No stable tangent origin; split the sky coverage")
        v /= np.linalg.norm(v)
        ra = np.arctan2(v[1], v[0]) / DEG % 360
        dec = np.arcsin(v[2]) / DEG
    return TangentPlane(ra, dec, max_field_deg)


def checked_polygon(points):
    polygon = Polygon(points)
    if polygon.is_empty or polygon.area <= 0 or not polygon.is_valid:
        raise RegionError(f"Invalid/zero-area polygon: {explain_validity(polygon)}")
    return polygon


class ShapeBuilder:
    """Adaptive boundary sampling, followed by planar polygon operations."""
    def __init__(self, plane, tolerance_arcsec=0.2, max_vertices=200000):
        self.plane = plane
        self.tolerance = tolerance_arcsec / 3600
        self.max_vertices = max_vertices
        self.evaluations = 0

    def _sample(self, function, knots):
        points = []

        def evaluate(t):
            self.evaluations += 1
            if self.evaluations > self.max_vertices * 8:
                raise RegionError("Sampling limit exceeded; increase --tolerance-arcsec")
            return np.asarray(function(t))

        def refine(a, b, pa, pb, depth=0):
            mid = (a + b) / 2
            pm = evaluate(mid)
            chord = pb - pa
            length2 = chord @ chord
            fraction = np.clip((pm - pa) @ chord / length2, 0, 1) if length2 else 0
            error = np.linalg.norm(pm - (pa + fraction * chord))
            if error > self.tolerance:
                if depth >= 25:
                    raise RegionError("Adaptive sampling did not converge")
                refine(a, mid, pa, pm, depth + 1)
                refine(mid, b, pm, pb, depth + 1)
            else:
                points.append(pa)
                if len(points) > self.max_vertices:
                    raise RegionError("Too many vertices; increase --tolerance-arcsec")

        for a, b in zip(knots[:-1], knots[1:]):
            refine(a, b, evaluate(a), evaluate(b))
        points.append(evaluate(knots[-1]))
        return np.array(points)

    def _local_map(self, ra, dec, rotation):
        center, east, north = basis(ra, dec)
        c, s = np.cos(rotation * DEG), np.sin(rotation * DEG)

        def convert(point):
            x, y = point
            # DS9 convention verified on a north-up/east-left TAN image:
            # 0=west, 90=north; positive rotation is counterclockwise on screen.
            e, n = -(c * x - s * y) * DEG, (s * x + c * y) * DEG
            radius = np.hypot(e, n)
            if radius >= self.plane.max_field_deg * DEG:
                raise RegionError("Individual shape is too large for the small-field model")
            factor = np.sinc(radius / np.pi)
            vector = np.cos(radius) * center + factor * (e * east + n * north)
            return self.plane.project(vector)
        return convert

    def _loop(self, kind, axes, transform):
        a, b = axes
        if a <= 0 or b <= 0:
            raise RegionError("Outer radii/widths must be positive")
        if kind == "box":
            corners = np.array([[-a, -b], [a, -b], [a, b], [-a, b], [-a, -b]])
            points = []
            for p, q in zip(corners[:-1], corners[1:]):
                points.extend(self._sample(lambda t: transform(p + t * (q - p)), [0, 1])[:-1])
            return checked_polygon(points)
        points = self._sample(lambda t: transform([a * np.cos(t), b * np.sin(t)]),
                              np.linspace(0, 2 * np.pi, 73))
        return checked_polygon(points[:-1])

    def _sector(self, kind, axes, start, sweep, transform):
        """Ray-angle sector of an ellipse/box, with box corners inserted exactly."""
        a, b = axes
        if a <= 0 or b <= 0:
            raise RegionError("Sector outer dimensions must be positive")
        if sweep == 360:
            return self._loop(kind, axes, transform)
        start, end = start * DEG, (start + sweep) * DEG
        knots = list(np.linspace(start, end, max(2, int(np.ceil(sweep / 5)) + 1)))
        if kind == "box":
            for angle in [np.arctan2(y * b, x * a) for x, y in [(1, 1), (-1, 1), (-1, -1), (1, -1)]]:
                for k in range(-2, 4):
                    t = angle + k * 2 * np.pi
                    if start < t < end:
                        knots.append(t)

        def edge(t):
            c, s = np.cos(t), np.sin(t)
            if kind == "box":
                radius = 1 / max(abs(c) / a, abs(s) / b)
            else:
                radius = 1 / np.hypot(c / a, s / b)
            return transform([radius * c, radius * s])
        # Radial segments are great circles and straight in the common TAN plane.
        points = self._sample(edge, sorted(set(knots)))
        return checked_polygon(np.vstack((transform([0, 0]), points)))

    def build(self, record):
        self.evaluations = 0
        try:
            result = self._build(record)
            if result.is_empty or not result.is_valid or result.geom_type not in {"Polygon", "MultiPolygon"}:
                raise RegionError(f"Invalid shape result: {explain_validity(result)}")
            return result
        except (ValueError, IndexError, TypeError) as exc:
            raise RegionError(f"{record.source} ({record.shape}): {exc}") from exc

    def _build(self, record):
        shape, args = record.shape, list(record.args)
        n_tokens = [t for t in args if t.lower().startswith("n=")]
        if n_tokens:
            if shape not in {"annulus", "ellipse", "box"} or len(n_tokens) != 1:
                raise RegionError("n= is supported only for annulus/ellipse-annulus/box-annulus")
            count(n_tokens[0][2:])
            args.remove(n_tokens[0])
        if shape == "polygon":
            if len(args) < 6 or len(args) % 2:
                raise RegionError("polygon needs at least three RA/Dec pairs")
            coords = [position(*args[i:i + 2]) for i in range(0, len(args), 2)]
            # Input polygon edges use shortest great-circle segments.
            return checked_polygon(self.plane.project(unit_vector(*np.array(coords).T)))
        if len(args) < 2:
            raise RegionError("Missing center")
        ra, dec = position(*args[:2])
        p = args[2:]
        rotation, start, sweep = 0.0, 0.0, 360.0
        inner = None
        kind = "ellipse"
        if shape == "circle":
            if len(p) != 1:
                raise RegionError("circle expects x,y,radius")
            outer = (angular(p[0]),) * 2
        elif shape == "annulus":
            if len(p) < 2:
                raise RegionError("annulus expects x,y,inner,outer[,r3,...]")
            radii = np.array([angular(t) for t in p])
            if radii[0] < 0 or np.any(np.diff(radii) <= 0):
                raise RegionError("Annulus radii must be nonnegative and strictly increasing")
            inner, outer = (radii[0],) * 2, (radii[-1],) * 2
        elif shape in {"ellipse", "box"}:
            kind = shape
            if len(p) % 2:
                rotation = angular(p.pop())
            if len(p) < 2 or len(p) % 2:
                raise RegionError("Expected axis/size pairs followed by an optional angle")
            axes = np.array([angular(t) for t in p]).reshape(-1, 2)
            if n_tokens and len(axes) < 2:
                raise RegionError("n= requires inner and outer sizes")
            if shape == "box":
                axes /= 2  # DS9 box values are FULL widths/heights.
            self._validate_axes(axes)
            outer = axes[-1]
            inner = axes[0] if len(axes) > 1 else None
        elif shape in {"panda", "epanda", "bpanda"}:
            expected = (6,) if shape == "panda" else (8, 9)
            if len(p) not in expected:
                raise RegionError(f"Incorrect {shape} argument count; see README syntax table")
            start, end = angular(p[0]), angular(p[1])
            delta = end - start
            if abs(delta) > 360 + 1e-8:
                raise RegionError("Panda angle range must not exceed one revolution")
            sweep = delta % 360
            if abs(sweep) < 1e-10:
                sweep = 360.0
            start %= 360
            count(p[2])
            if shape == "panda":
                count(p[5])
                inner, outer = (angular(p[3]),) * 2, (angular(p[4]),) * 2
            else:
                count(p[7])
                rotation = angular(p[8]) if len(p) == 9 else 0
                inner = np.array([angular(t) for t in p[3:5]])
                outer = np.array([angular(t) for t in p[5:7]])
                if shape == "bpanda":
                    kind = "box"
                    inner, outer = inner / 2, outer / 2
            self._validate_axes(np.array([inner, outer]))
        else:
            raise RegionError(f"Unsupported shape: {shape}")

        transform = self._local_map(ra, dec, rotation)
        result = self._sector(kind, outer, start, sweep, transform)
        if inner is not None and max(inner) > 0:
            result = result.difference(self._sector(kind, inner, start, sweep, transform))
        return result

    @staticmethod
    def _validate_axes(axes):
        if np.any(axes < 0) or np.any(axes[-1] <= 0):
            raise RegionError("Dimensions must be nonnegative, with positive outer dimensions")
        if np.any((axes[:, 0] == 0) != (axes[:, 1] == 0)):
            raise RegionError("An inner size may be (0,0), but not have only one zero axis")
        if len(axes) > 1:
            if np.any(np.diff(axes, axis=0) <= 0):
                raise RegionError("Annulus dimensions must be strictly increasing")
            nonzero = axes[axes[:, 1] > 0]
            ratios = nonzero[:, 0] / nonzero[:, 1]
            if not np.allclose(ratios, ratios[-1], rtol=1e-4, atol=0):
                raise RegionError("DS9 annulus/panda axis ratios must be equal")


def polygon_parts(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == "Polygon":
        return [geometry]
    if geometry.geom_type in {"MultiPolygon", "GeometryCollection"}:
        return [p for g in geometry.geoms for p in polygon_parts(g)]
    return []


def union_regions(records, plane, *, tolerance_arcsec=0.2, exclude_scope="file"):
    builder = ShapeBuilder(plane, tolerance_arcsec)
    groups = {}
    for record in records:
        key = record.file_index if exclude_scope == "file" else 0
        included, excluded = groups.setdefault(key, ([], []))
        (included if record.include else excluded).append(builder.build(record))
    results = []
    for key, (included, excluded) in groups.items():
        if not included:
            raise RegionError(f"Input group {key} has exclusions but no included area; "
                              "use --exclude-scope global for a separate exclusion file")
        result = unary_union(included)
        if excluded:
            result = result.difference(unary_union(excluded))
        results.append(result)
    result = unary_union(results)
    if not polygon_parts(result):
        raise RegionError("Union is empty after exclusions")
    return result


def prepare_outline(geometry, *, keep_holes=False, smooth_arcsec=0, tolerance_arcsec=0.2):
    """Optional morphological closing rounds concave joins, changing the area."""
    if smooth_arcsec > 0:
        distance = smooth_arcsec / 3600
        ratio = min(1.0, tolerance_arcsec / smooth_arcsec)
        angle = 2 * math.acos(max(-1.0, 1 - ratio))
        resolution = max(16, math.ceil(np.pi / (2 * max(angle, 1e-5))))
        if resolution > 4096:
            raise RegionError("Smoothing resolution too large; increase --tolerance-arcsec")
        rounded = geometry.buffer(distance, quad_segs=resolution).buffer(-distance, quad_segs=resolution)
        geometry = unary_union([geometry, rounded])
    if not keep_holes:
        # Filling a hole can absorb a disconnected island; union again to remove it.
        geometry = unary_union([Polygon(p.exterior) for p in polygon_parts(geometry)])
    return geometry


def serialize_regions(geometry, plane, *, keep_holes=False, tolerance_arcsec=0.2,
                      smooth_arcsec=0, exclude_scope="file", color="cyan"):
    polygons = sorted(polygon_parts(geometry), key=lambda p: (-p.area, p.centroid.x, p.centroid.y))
    if not polygons:
        raise RegionError("No nonzero-area outline remains")
    holes = [Polygon(ring) for p in polygons for ring in p.interiors] if keep_holes else []
    # DS9 exclusions apply globally. A positive island inside a negative hole
    # cannot be represented by these exterior/-hole records without losing it.
    for hole in holes:
        if any(hole.contains(p.representative_point()) for p in polygons):
            raise RegionError("Nested island inside a hole cannot be serialized as DS9 "
                              "global exclusions; omit --keep-holes for outer boundaries")
    lines = ["# Region file format: DS9 version 4.1",
             "# Generated by region_union_fk5.py - FK5 J2000, small-field spherical model",
             f"# Boundary sampling target: {tolerance_arcsec:g} arcsec (not a WCS accuracy guarantee)",
             f"# Holes: {'kept' if keep_holes else 'filled for exterior-only output'}",
             f"# Exclusion scope: {exclude_scope}",
             f"# Morphological closing radius: {smooth_arcsec:g} arcsec (0 = no smoothing)",
             f"# Tangent origin: {plane.ra:.10f}, {plane.dec:.10f} deg",
             f'global color={color} width=2 fill=0 include=1', "fk5"]
    vertices = 0
    for sign, rings in [("", [p.exterior for p in polygons]), ("-", [h.exterior for h in holes])]:
        for ring in rings:
            coords = plane.unproject(np.array(ring.coords)[:-1])
            vertices += len(coords)
            text = ",".join(f"{ra:.10f},{dec:.10f}" for ra, dec in coords)
            lines.append(f"{sign}polygon({text})")
    return "\n".join(lines) + "\n", len(polygons), len(holes), vertices


def positive_float(token):
    value = float(token)
    if not np.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("Must be finite and positive")
    return value


def nonnegative_float(token):
    value = float(token)
    if not np.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError("Must be finite and nonnegative")
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="+", type=Path, help="Input DS9 .reg files, explicitly fk5/j2000")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output DS9 polygon .reg")
    parser.add_argument("--tolerance-arcsec", type=positive_float, default=0.2,
                        help="Boundary sampling target in arcsec (default: 0.2)")
    parser.add_argument("--smooth-arcsec", type=nonnegative_float, default=0,
                        help="Optional closing radius; rounds concave joins and can join gaps")
    parser.add_argument("--keep-holes", action="store_true", help="Write hole boundaries as -polygon exclusions")
    parser.add_argument("--exclude-scope", choices=["file", "global"], default="file",
                        help="Subtract exclusions within each file before union (default), or globally")
    parser.add_argument("--origin", nargs=2, metavar=("RA", "DEC"), help="Optional FK5 projection origin")
    parser.add_argument("--max-field-deg", type=positive_float, default=5.0,
                        help="Maximum distance from tangent origin, <=10 deg (default: 5)")
    parser.add_argument("--ignore-nonarea", action="store_true", help="Skip point/line/text/etc. with warnings")
    parser.add_argument("--color", choices=["cyan", "green", "red", "blue", "yellow", "magenta", "white"], default="cyan")
    parser.add_argument("--overwrite", action="store_true", help="Allow replacing an existing output")
    args = parser.parse_args(argv)
    try:
        if args.max_field_deg > 10:
            raise RegionError("--max-field-deg must be <=10; this is a small-field tool")
        if args.output.resolve() in {p.resolve() for p in args.inputs}:
            raise RegionError("Output must differ from all input paths")
        if args.output.exists() and not args.overwrite:
            raise RegionError("Output exists; choose a new filename or use --overwrite")
        records = read_regions(args.inputs, ignore_nonarea=args.ignore_nonarea)
        plane = choose_plane(records, args.origin, args.max_field_deg)
        geometry = union_regions(records, plane, tolerance_arcsec=args.tolerance_arcsec,
                                 exclude_scope=args.exclude_scope)
        initial_holes = sum(len(p.interiors) for p in polygon_parts(geometry))
        geometry = prepare_outline(geometry, keep_holes=args.keep_holes,
                                   smooth_arcsec=args.smooth_arcsec,
                                   tolerance_arcsec=args.tolerance_arcsec)
        text, components, holes, vertices = serialize_regions(
            geometry, plane, keep_holes=args.keep_holes,
            tolerance_arcsec=args.tolerance_arcsec, smooth_arcsec=args.smooth_arcsec,
            exclude_scope=args.exclude_scope, color=args.color)
        # Validate everything before opening the destination for writing.
        with args.output.open("w" if args.overwrite else "x", encoding="ascii") as handle:
            handle.write(text)
        print(f"Wrote {args.output}: {components} component(s), {holes} hole(s), {vertices} vertices")
        if initial_holes and not args.keep_holes:
            print(f"Exterior-only output: filled {initial_holes} interior hole(s)")
        if args.smooth_arcsec:
            print("Smoothing changes the selected area and may connect nearby components")
        return 0
    except (RegionError, OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
