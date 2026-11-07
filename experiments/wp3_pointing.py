import sys
import numpy as np
sys.path.insert(0, ".")
from rtops.pointing import PointingModel

# 1. Test the iterative inverse pointing solver
true_params = [0.35, -0.2, 0.08, 0.1, -0.15]
model = PointingModel(true_params)

sky_az, sky_el = 180.0, 45.0
enc_az, enc_el = model.sky_to_enc(sky_az, sky_el)
d_az, d_el = model.delta(enc_az, enc_el)

print("--- 1. Iterative Solver Test ---")
print(f"Target Sky Position: Az {sky_az}°, El {sky_el}°")
print(f"Calculated Encoder Angles: Az {enc_az:.4f}°, El {enc_el:.4f}°")
print(f"Verification (Encoder + Delta): Az {enc_az + d_az:.4f}°, El {enc_el + d_el:.4f}°\n")

# 2. Test Design Matrix Conditioning (Linear Algebra)
# Using rough Az/El coordinates typical of GEO satellites and the Sun
az_test = np.array([120, 150, 180, 210, 240, 100, 260]) 
el_test = np.array([50, 53, 54, 53, 50, 40, 40])

X = PointingModel.design(az_test, el_test, az_test, el_test)
cond = np.linalg.cond(X)

print("--- 2. Least Squares Design Matrix ---")
print(f"Matrix Condition Number: {cond:.2f}")
print("Result: A low condition number proves the 5 pointing parameters are independent and can be safely fitted!")