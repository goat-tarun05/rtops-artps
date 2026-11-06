import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, ".")
from rtops.sky import refraction_deg

el = np.linspace(10, 90, 100)
ref = refraction_deg(el)

plt.figure(figsize=(8, 5))
plt.plot(el, ref, label="Saemundsson Refraction")
plt.xlabel("True Elevation (degrees)")
plt.ylabel("Refraction (degrees)")
plt.title("Simulated Atmospheric Refraction vs Elevation")
plt.grid(True)
plt.savefig("results/WP1/refraction_plot.png")
print("Plot saved to results/WP1/refraction_plot.png")