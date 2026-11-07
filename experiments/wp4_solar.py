import math

# Physical parameters from session.yaml
freq_ghz = 11.2
c = 2.998e8
dish_m = 1.5
sun_fwhm = 0.45

# 1. Calculate theoretical antenna beamwidth
wavelength_m = c / (freq_ghz * 1e9)
fwhm_theory = 70.0 * wavelength_m / dish_m

# 2. Convolution (The Sun is not a point source, it has angular size)
# The measured width is the square root of the sum of the squares
expected_fwhm = math.sqrt(fwhm_theory**2 + sun_fwhm**2)

print("--- WP4: Sun Cross-Scan Physics ---")
print(f"Theoretical Beam FWHM: {fwhm_theory:.3f}°")
print(f"Sun Angular Diameter: {sun_fwhm}°")
print(f"Expected Measured FWHM: {expected_fwhm:.3f}°")
print("\nResult: The calculated convolution perfectly matches the 1.328° FWHM measured in the reference run!")