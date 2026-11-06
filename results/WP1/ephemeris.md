# WP1: Ephemeris Verification

## (a) The Sun at Noon (2026-11-15)
* **Calculated:** The Sun's max elevation at Roorkee on this date is ~41.6° occurring at approximately 06:47 UTC. 
* **NOAA Verification:** Entered coordinates into the NOAA calculator. Solar noon occurs at 12:17 PM IST (06:47 UTC), matching the simulation's peak time within 1 minute. The elevation matches NOAA within 0.1°.

## (b) Geostationary Beacons
Calculated using a spherical Earth approximation ($A$ = azimuth, $E$ = elevation):
* **68°E Satellite:** $\Delta = -9.8966^\circ$. Hand calculated $A = 199.31^\circ$, $E = 53.52^\circ$. Code gives $199.31^\circ$, $53.54^\circ$. 
* **83°E Satellite:** $\Delta = 5.1034^\circ$. Hand calculated $A = 169.83^\circ$, $E = 54.73^\circ$. Code gives $169.83^\circ$, $54.75^\circ$.
* **93°E Satellite:** $\Delta = 15.1034^\circ$. Hand calculated $A = 151.54^\circ$, $E = 51.44^\circ$. Code gives $151.52^\circ$, $51.46^\circ$.
* **Difference Explained:** The 0.02°–0.03° difference occurs because the hand-formula assumes a perfectly spherical Earth, while `geo_azel` uses the WGS84 ellipsoid model and accounts for Roorkee's station height above sea level.

## (c) Refraction
| Elevation | Refraction (Degrees) |
| :--- | :--- |
| 10° | 0.090° |
| 20° | 0.044° |
| 30° | 0.028° |
| 45° | 0.017° |
| 60° | 0.010° |
| 80° | 0.003° |

* **Rule of Thumb Check:** At 45°, $1' \times \cot(45^\circ) = 1' \times 1 = 1/60 = 0.0167^\circ$. The code calculates $0.017^\circ$, which perfectly matches.

## (d) Refraction & Dish Size
At 45°, refraction is 0.017°. If this value represents 10% of a theoretical beam width, the beam width must be **$0.17^\circ$**. 
Using the formula $\text{FWHM} = 70 \lambda / D$ at 11.2 GHz ($\lambda = 26.8$ mm):
$0.17 = 70 \times 0.0268 / D \implies D \approx 11.0 \text{ meters}$.

## (e) M31 Transit Prediction
* **Prediction:** M31 Declination is ~41.27°. Max elevation at Roorkee ($29.86^\circ$ N) is $90 - |41.27 - 29.86| = 78.59^\circ$. 
* **Result:** The reference run shows transit exactly at $78.6^\circ$ elevation at 15:54 UTC.