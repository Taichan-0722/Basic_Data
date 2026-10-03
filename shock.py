import numpy as np

v_shock = 1500.0    # Shock velocity (km/s)
kT1_keV = 4.0       # Upstream temperature (keV)
P1 = 1.20           # Upstream pressure (A.U.)
K1 = 8.54           # Upstream entropy (A.U.)

gamma = 5.0 / 3.0            # Specific heat ratio
mu = 0.61                    # Mean molecular weight
keV_to_J = 1.602176634e-16   # 1 keV (J)
m_p = 1.67262192369e-27      # Proton mass (kg)

kT1_J = kT1_keV * keV_to_J                         # Convert to J
v_sound_ms = np.sqrt(gamma * kT1_J / (mu * m_p))   # Sound speed (m/s)
v_sound_kms = v_sound_ms / 1000.0                  # Convert to km/s
M1 = v_shock / v_sound_kms                         # Mach number

density_ratio = ((gamma + 1.0) * M1**2 / ((gamma - 1.0) * M1**2 + 2.0))
pressure_ratio = ((2.0 * gamma * M1**2 - (gamma - 1.0)) / (gamma + 1.0))
temperature_ratio = ((2.0 * gamma * M1**2 - (gamma - 1.0)) * ((gamma - 1.0) * M1**2 + 2.0) / ((gamma + 1.0)**2 * M1**2))
entropy_ratio = (pressure_ratio / density_ratio**gamma)

kT2_keV = kT1_keV * temperature_ratio   # Downstream temperature (keV)
P2 = P1 * pressure_ratio                # Downstream pressure (A.U.)
K2 = K1 * entropy_ratio                 # Downstream entropy (A.U.)

print()
print("===== Shock properties =====")
print(f"Shock velocity            = {v_shock:.3f} km/s")
print(f"Upstream temperature      = {kT1_keV:.3f} keV")
print(f"Upstream pressure         = {P1:.3f} A.U.")
print(f"Upstream entropy          = {K1:.3f} A.U.")
print(f"Specific heat ratio       = {gamma:.3f}")
print(f"Mean molecular weight     = {mu:.3f}")
print()

print("===== Mach number =====")
print(f"Sound speed               = {v_sound_kms:.5f} km/s")
print(f"Mach number M1            = {M1:.5f}")
print()

print("===== Shock jump ratios =====")
print(f"Density ratio   rho2/rho1 = {density_ratio:.5f}")
print(f"Pressure ratio      P2/P1 = {pressure_ratio:.5f}")
print(f"Temperature ratio   T2/T1 = {temperature_ratio:.5f}")
print(f"Entropy ratio       K2/K1 = {entropy_ratio:.5f}")
print()

print("===== Downstream =====")
print(f"Downstream temperature    = {kT2_keV:.5f} keV")
print(f"Downstream pressure       = {P2:.5f} A.U.")
print(f"Downstream entropy        = {K2:.5f} A.U.")

if M1 <= 1.0:
    print()
    print("WARNING: M1 <= 1.")
    print("The input velocity does not correspond to a supersonic shock.")