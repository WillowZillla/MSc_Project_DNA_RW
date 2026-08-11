from StatTools.generators import generate_fbn
from StatTools.analysis.dfa import dfa
from StatTools.analysis.utils import analyse_zero_cross_ff
import numpy as np

h = 0.7  # choose Hurst parameter
length = 6000  # vector's length

# Generate synthetic data using modern unified interface
trajectory = generate_fbn(hurst=h, length=length, method="kasdin").flatten()

print(trajectory)
s_vals, f2_vals = dfa(trajectory, degree=2)
print(s_vals)
print(f2_vals)

# Calculate Hurst exponent from fluctuation function
f_vals = np.sqrt(f2_vals).reshape(1, -1)  # Convert F^2(s) to F(s) and reshape for analysis
s_vals_2d = s_vals.reshape(1, -1)  # Reshape scales to 2D array
hurst_result, _ = analyse_zero_cross_ff(f_vals, s_vals_2d)
hurst_exponent = hurst_result.slopes[0].value

print(f"Estimated H: {hurst_exponent:.3f} (Expected: {h:.3f})")
