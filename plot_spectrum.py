#!/usr/bin/env python3
"""XSPEC の QDP + PCO からスペクトルと残差を表示する (自動保存なし)。

依存: numpy, matplotlib。実行例: python3 plot_spectrum.py --help
対象: READ SERR 1 2、スペクトル1ブロック + delchi 1ブロック。
QDP/PCO の命令は実行せず、必要なデータと設定のみ読み取る。
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
import re
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter, FuncFormatter, MaxNLocator


# ===== よく変更する設定 =====
DEFAULT_QDP = Path(
    "/Users/mitsuhashitaiki/Documents/analysis/Centaurus_NE/Region_big_fit/"
    "xa000138000xtd_p030000010_beta_58_1apec_0p4keV.qdp"
)


@dataclass
class PlotStyle:
    data_color: str = "#303030"
    model_color: str = "#D43F3A"
    residual_color: str = "#303030"
    data_linewidth: float = 0.8       # 誤差棒の太さ (pt)
    model_linewidth: float = 2.0
    component_linewidth: float = 1.5
    residual_linewidth: float = 0.8
    model_linestyle: str = "-"
    component_linestyle: str = "--"
    component_alpha: float = 0.85
    marker: str = "."                # "none" でマーカーなし
    markersize: float = 2.0
    capsize: float = 0.0
    x_errorbars: bool = True
    show_components: bool = True
    show_legend: bool = True
    grid: bool = False
    figsize: tuple[float, float] = (9.0, 6.8)
    height_ratios: tuple[float, float] = (3.0, 1.0)
    font_family: str = "Times New Roman"
    residual_label_fontsize: float | None = 25.0  # 残差の縦軸ラベル (pt)。None は既定値
    legend_fontsize: float | None = 13.0          # 凡例 (pt)。None は既定値
    legend_ncol: int = 6                       # 枠外の凡例の列数
    axes_linewidth: float = 1.0
    xscale: str | None = None         # None: PCO に従う
    yscale: str | None = None
    xlim: tuple[float, float] | None = None
    ylim: tuple[float, float] | None = None
    residual_ylim: tuple[float, float] | None = None
    title: str | None = None          # "" でタイトルなし
    xlabel: str = "Energy (keV)"
    ylabel: str = r"Counts s$^{-1}$ keV$^{-1}$"
    residual_ylabel: str = "Residual"
    component_colors: tuple[str, ...] = (
        "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
        "#8C564B", "#7B61A8", "#BCBD22", "#666666",
    )
    # キーは成分番号 (1 = QDP の6列目)。ゼロ成分を省略しても番号は不変。
    # 例: {1: {"label": "ICM", "color": "royalblue", "linewidth": 1.8}}
    # 物理名は XSPEC のモデル式を確認してから設定してください。
    components: dict[int, dict] = field(default_factory=dict)


STYLE = PlotStyle()
STYLE.components = {1: {"label": "ICM", "color": "#0072B2"},
                    2: {"label": "NXB", "color": "#E69F00"},
                    3: {"label": r"$^{55}$Fe", "color": "#009E73"},
                    4: {"label": "CXB", "color": "#CC79A7"},
                    5: {"label": "GH", "color": "#56B4E9"},
                    6: {"label": "LHB", "color": "#8C564B"},
                    7: {"label": "0.8 keV", "color": "#7B61A8"},
                    8: {"label": "OI", "color": "#BCBD22"},
                    9: {"label": "OVII", "color": "#666666"}}

def read_qdp(path: Path) -> tuple[list[np.ndarray], str | None]:
    """NO のみの行を区切り、個々の NO は欠損値として扱う。"""
    blocks, rows = [], []
    serr, referenced_pco = None, None

    def finish():
        if rows:
            if len({len(row) for row in rows}) != 1:
                raise ValueError("同一 QDP ブロック内で列数が一致しません。")
            blocks.append(np.asarray(rows, dtype=float))
            rows.clear()

    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        line = line.split("!", 1)[0].strip()
        if not line:
            continue
        tokens = line.split()
        if line.upper().startswith("READ SERR"):
            serr = [int(t) for t in tokens[2:]]
            continue
        if line.upper().startswith("READ TERR"):
            raise ValueError("READ TERR (非対称誤差) はこのスクリプトの対象外です。")
        if line.startswith("@"):
            referenced_pco = line[1:].strip().strip('\"\'')
            continue
        if all(t.upper() == "NO" for t in tokens):
            finish()
            continue
        if not (re.match(r"^[+\-.\d]", line) or tokens[0].upper() == "NO"):
            continue
        try:
            rows.append([np.nan if t.upper() == "NO" else
                         float(t.replace("D", "E").replace("d", "e")) for t in tokens])
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: 数値行を読み取れません。") from exc
    finish()
    if serr != [1, 2]:
        raise ValueError(f"READ SERR 1 2 が必要です (検出: {serr})。")
    if len(blocks) not in (1, 2) or blocks[0].shape[1] < 5:
        raise ValueError("スペクトル1ブロックと任意の残差1ブロックを指定してください。")
    for block in blocks:
        if block.shape[1] < 4:
            raise ValueError("エネルギー・幅・値・誤差の4列が必要です。")
        if np.any(block[:, [1, 3]] < 0):
            raise ValueError("エネルギー幅または誤差に負の値があります。")
        if not np.all(np.isfinite(block[:, :2])) or np.any(np.diff(block[:, 0]) <= 0):
            raise ValueError("エネルギーと幅は有限値で、エネルギーは昇順である必要があります。")
    return blocks, referenced_pco


def read_pco(path: Path | None) -> dict:
    """タイトル、軸範囲、対数設定のみを読む (PCO 全体の再現器ではない)。"""
    meta = {"title": "", "xscale": "log", "yscale": "log"}
    if path is None:
        return meta
    window = 1
    for line in path.read_text().splitlines():
        tokens = line.split()
        if not tokens:
            continue
        upper = [t.upper() for t in tokens]
        if upper[0] == "WIN":
            window = int(tokens[1]) if tokens[1].isdigit() else 0
        elif upper[:2] == ["LAB", "T"] and window == 1:
            meta["title"] = " ".join(tokens[2:]).strip('\"')
        elif upper[0] == "R" and len(tokens) == 4:
            axis = re.fullmatch(r"([XY])(\d*)", upper[1])
            if axis and "NO" not in upper[2:]:
                target = int(axis[2]) if axis[2] else window
                limits = tuple(float(t.replace("D", "E")) for t in tokens[2:])
                key = "xlim" if axis[1] == "X" and target in (1, 2) else (
                    "ylim" if axis[1] == "Y" and target == 1 else (
                        "residual_ylim" if axis[1] == "Y" and target == 2 else None))
                if key:
                    meta[key] = limits
        elif upper[0] == "LOG" and len(tokens) >= 3:
            targets = [int(t) for t in tokens[3:] if t.isdigit()] or [window]
            if 1 in targets and upper[1] in ("X", "Y"):
                meta[upper[1].lower() + "scale"] = "log" if upper[2] == "ON" else "linear"
    return meta


def bin_curve(ax, x, dx, y, *, log_y, **kwargs):
    """実際の各ビン端で段状に描く。隙間・欠損値を線でつながない。"""
    left, right = x - dx, x + dx
    valid = np.isfinite(y) & ((y > 0) if log_y else True)
    px, py = [], []
    previous_right = None
    for lo, hi, value, ok in zip(left, right, y, valid):
        if not ok:
            px.append(np.nan); py.append(np.nan)
            previous_right = None
            continue
        if previous_right is not None and not np.isclose(lo, previous_right, rtol=2e-6, atol=1e-8):
            px.append(np.nan); py.append(np.nan)
        px.extend([lo, hi]); py.extend([value, value])
        previous_right = hi
    ax.plot(px, py, **kwargs)


def plot_spectrum(qdp: Path, pco: Path | None, style: PlotStyle):
    """図を表示し、Figure を返す。ファイルへの保存は行わない。"""
    if style.legend_ncol < 1:
        raise ValueError("legend_ncol は1以上を指定してください。")
    blocks, reference = read_qdp(qdp)
    if pco is None:
        candidates = [qdp.parent / reference] if reference else []
        candidates.append(qdp.with_suffix(".pco"))
        pco = next((p for p in candidates if p.is_file()), None)
        if pco is None:
            warnings.warn("PCO が見つからないため、標準設定を使用します。")
    meta = read_pco(pco)
    spectrum = blocks[0]
    x, dx, y, dy, model = spectrum[:, :5].T
    if len(blocks) == 2:
        residual = blocks[1]
        if residual.shape[0] != len(x) or not np.allclose(residual[:, :2], spectrum[:, :2]):
            raise ValueError("スペクトルと残差のエネルギービンが一致しません。")
        with np.errstate(divide="ignore", invalid="ignore"):
            expected = (y - model) / dy
        valid = np.isfinite(expected) & np.isfinite(residual[:, 2])
        if not valid.any() or not np.allclose(expected[valid], residual[valid, 2], atol=2e-4, rtol=2e-4):
            raise ValueError("第2ブロックが (data-model)/error と一致しません。残差の種類を確認してください。")
        rx, rdx, ry, rdy = residual[:, :4].T
    else:
        warnings.warn("残差ブロックがないため、(data-model)/error を計算します。")
        rx, rdx = x, dx
        ry = np.divide(y-model, dy, out=np.full_like(y, np.nan), where=dy > 0)
        rdy = np.ones_like(ry)

    xscale, yscale = style.xscale or meta["xscale"], style.yscale or meta["yscale"]
    xlim = style.xlim or meta.get("xlim") or (float(np.min(x-dx)), float(np.max(x+dx)))
    ylim = style.ylim or meta.get("ylim")
    rlim = style.residual_ylim  # 未指定時は表示中の残差から自動設定
    for limits, scale, label in [(xlim, xscale, "xlim"), (ylim, yscale, "ylim"), (rlim, "linear", "residual_ylim")]:
        if limits is not None and (not np.all(np.isfinite(limits)) or limits[0] >= limits[1] or (scale == "log" and limits[0] <= 0)):
            raise ValueError(f"{label}: 昇順の有限な範囲が必要です。対数軸は正の値を指定してください。")

    with plt.rc_context({"font.family": style.font_family,
                         "mathtext.fontset": "custom", "mathtext.rm": style.font_family,
                         "mathtext.it": style.font_family + ":italic",
                         "mathtext.bf": style.font_family + ":bold",
                         "axes.linewidth": style.axes_linewidth, "pdf.fonttype": 42,
                         "ps.fonttype": 42, "savefig.facecolor": "white"}):
        fig, (ax, res) = plt.subplots(2, 1, figsize=style.figsize, sharex=True,
                                     gridspec_kw={"height_ratios": style.height_ratios, "hspace": 0.04})
        ax.set_xscale(xscale); ax.set_yscale(yscale)
        ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        if style.show_components:
            for i, values in enumerate(spectrum[:, 5:].T, 1):
                if not np.any(np.isfinite(values) & (values != 0)):
                    continue
                options = dict(label=f"Component {i}", color=style.component_colors[(i-1) % len(style.component_colors)],
                               linewidth=style.component_linewidth, linestyle=style.component_linestyle,
                               alpha=style.component_alpha)
                options.update(style.components.get(i, {}))
                bin_curve(ax, x, dx, values, log_y=yscale == "log", zorder=1, **options)
        good = np.isfinite(y) & np.isfinite(dy)
        if yscale == "log":
            hidden = int(np.sum(good & (y <= 0)))
            if hidden:
                warnings.warn(f"対数軸で表示できない非正の観測値 {hidden} 点は上段のみ省略します。")
            good &= y > 0
        common = dict(fmt=style.marker, markersize=style.markersize, capsize=style.capsize)
        ax.errorbar(x[good], y[good], xerr=dx[good] if style.x_errorbars else None,
                    yerr=dy[good], color=style.data_color, elinewidth=style.data_linewidth,
                    label="Data", zorder=2, **common)
        bin_curve(ax, x, dx, model, log_y=yscale == "log", color=style.model_color,
                  linewidth=style.model_linewidth, linestyle=style.model_linestyle,
                  label="Total model", zorder=3)
        mask = np.isfinite(ry) & np.isfinite(rdy)
        res.errorbar(rx[mask], ry[mask], xerr=rdx[mask] if style.x_errorbars else None,
                     yerr=rdy[mask], color=style.residual_color,
                     elinewidth=style.residual_linewidth, **common)
        res.axhline(0, color=style.model_color, linewidth=1.0, linestyle="--", zorder=0)
        if rlim is not None:
            res.set_ylim(rlim)
        else:
            visible = mask & (rx + rdx >= xlim[0]) & (rx - rdx <= xlim[1])
            if visible.any():
                # 誤差棒の端まで含め、上下それぞれ5%の余白だけを加える。
                # 上下対称にはせず、外れた点も切り落とさない。
                lower = float(np.min(ry[visible] - rdy[visible]))
                upper = float(np.max(ry[visible] + rdy[visible]))
                span = upper - lower
                padding = 0.05 * span if span > 0 else max(abs(lower) * 0.05, 0.1)
                res.set_ylim(lower - padding, upper + padding)
            else:
                res.set_ylim(-1, 1)
        res.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
        ax.set_ylabel(style.ylabel)
        res.set_ylabel(style.residual_ylabel, fontsize=style.residual_label_fontsize)
        res.set_xlabel(style.xlabel)
        # タイトルと枠外の凡例を別々に配置し、自動レイアウトで余白を確保する。
        title = meta["title"] if style.title is None else style.title
        if title:
            fig.suptitle(title, fontsize=plt.rcParams["axes.titlesize"],
                         fontweight=plt.rcParams["axes.titleweight"])
        if xscale == "log":
            res.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            res.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            res.xaxis.set_minor_locator(LogLocator(base=10, subs=(3, 4, 6, 7, 8, 9)))
            res.xaxis.set_minor_formatter(NullFormatter())
        for axis in (ax, res):
            axis.tick_params(which="both", direction="in", top=True, right=True)
            axis.tick_params(which="major", length=6)
            axis.tick_params(which="minor", length=3)
            if style.grid:
                axis.grid(which="major", alpha=0.16, linewidth=0.6)
        if style.show_legend:
            ax.legend(fontsize=style.legend_fontsize, ncol=style.legend_ncol,
                      loc="lower center", bbox_to_anchor=(0.5, 1.02),
                      borderaxespad=0.0)
        # 表示のみ。保存する場合はウィンドウの保存ボタンを使用する。
        plt.show(block=True)
    print(f"Spectrum: {len(x)} bins; components: {spectrum.shape[1]-5}; PCO: {pco}")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("qdp", nargs="?", type=Path, default=DEFAULT_QDP, help="入力QDP (省略時は今回のファイル)")
    parser.add_argument("--pco", type=Path, help="省略時は QDP の @ 行、次に同名PCOを探索")
    for option in ("data-color", "model-color", "residual-color", "title"):
        parser.add_argument("--" + option)
    for option in ("data-linewidth", "model-linewidth", "component-linewidth", "residual-linewidth", "markersize",
                   "residual-label-fontsize", "legend-fontsize"):
        parser.add_argument("--" + option, type=float)
    for option in ("xlim", "ylim", "residual-ylim", "figsize"):
        parser.add_argument("--" + option, nargs=2, type=float)
    for option in ("xscale", "yscale"):
        parser.add_argument("--" + option, choices=("log", "linear"))
    parser.add_argument("--legend-ncol", type=int, help="枠外の凡例の列数 (既定値: 3)")
    parser.add_argument("--no-components", action="store_true")
    parser.add_argument("--no-xerr", action="store_true")
    args = parser.parse_args()
    for key in vars(STYLE):
        if hasattr(args, key) and getattr(args, key) is not None:
            setattr(STYLE, key, getattr(args, key))
    if args.no_components:
        STYLE.show_components = False
    if args.no_xerr:
        STYLE.x_errorbars = False
    plot_spectrum(args.qdp.expanduser(), args.pco.expanduser() if args.pco else None, STYLE)


if __name__ == "__main__":
    main()
