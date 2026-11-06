# rtops — a software-only radio-telescope control and data system (digital observatory)

A *digital twin* of a small radio telescope. Nothing here touches hardware: a virtual dish, mount, wind, receivers and sky generate the data,
and the control and reduction software treats them exactly as a real observatory would treat the real thing. Your own raw recordings can be
plugged in later through the adapters in `rtops/ingest.py` (Track R).

```
pip install astropy skyfield sgp4 scipy matplotlib pyyaml pytest
python -m rtops.run_session --out session            # the full 16-hour virtual day (~10 s), reduction, A/B comparison, report
python -m rtops.run_session --out session --no-ab    # same without the "no pointing calibration" comparison
open session/report.html
pytest -q                                            # starter tests
RTOPS_CONFIG=config/mini.yaml python -m rtops.run_session --out mini --no-ab   # 2-hour mini session
```

## What happens in a session
1. **Queue and scheduler** (`scheduler.py`): observation blocks (OBs) with priorities, windows and fixed times. The dispatcher fills gaps, inserts
   calibrations, picks a satellite pass that satisfies elevation / azimuth-rate / keyhole / weather constraints, and re-queues interrupted blocks.
2. **Telescope control** (`tcs.py`): `slew`, `track`, `scan`, `acquire` primitives combined into block types (calibration, pointing run, solar drift
   scan, satellite tracking, H I ON/OFF tracking, H I stare). Wind above the limit triggers stow and resume.
3. **Ephemeris** (`sky.py`): apparent az/el with refraction for sidereal targets, the Sun, geostationary beacons and satellites (SGP4).
4. **Pointing model** (`pointing.py`): TPOINT-style alt-az terms IA, IE, NPAE, AN, AW, inverted for every command and fitted from pointing runs.
5. **Mount** (`mount.py`): rate and acceleration limits, feed-forward tracking, actuator lag, encoder quantisation, wind-induced jitter, cable wrap (azimuth travel −90° to 450°: a pass that crosses north is tracked on the azimuth branch that keeps the whole pass inside this range, and the scheduler rejects passes that cannot fit).
6. **Receivers** (`receiver.py`): Ku total power and L-band spectra computed from where the beam REALLY points (command + hidden pointing error + wind).
7. **Archive** (`archive.py`): one FITS file per block (`session/archive/<date>/<block>.fits`) and an SQLite catalogue with calibration links.
8. **Reduction** (`reduce.py`, using `rtpipe/`): reads the archive only. Solar cross-scan, satellite light curve, H I line fit.
9. **Report** (`report.py`): `session/report.html`.

## Layout
| Path | Content |
|---|---|
| `rtops/` | the control and data system (modules above, plus `experiment.py`, `ingest.py`, `run_session.py`) |
| `rtpipe/` | science reduction library (solar, satellites, H I) used by `rtops/reduce.py` |
| `config/session.yaml` | the whole day: site, telescope, mount, receivers, hidden truth, weather, observation blocks |
| `config/mini.yaml`, `config/mini_solar.yaml` | short sessions for tests and quick experiments |
| `experiments/servo_lib.py` | helpers for the servo experiments (WP2) |
| `tests/` | starter tests |
| `reference/` | reference report and numbers from the unmodified configuration |

## Running an experiment (change a setting, run, compare)
```python
from rtops.experiment import run_variant
r = run_variant({"pointing_truth": [0.8, -0.5, 0.1, 0.2, -0.3], "session.seed": 11}, out="runs/big_errors")
pr = r["results"]["pointing_runs"][0]; print(pr["map_rms_true_before_deg"], pr["map_rms_true_after_deg"])
print(r["reduced"]["solar"]["sun_peak_k"], r["reduced"]["sat"]["mean_eff"], r["reduced"]["hi"]["tracked ON/OFF"]["vsys_helio_kms"])
```
Keys are dotted paths into `config/session.yaml`; list positions are numbers (`"observations.0.params.geo"`, `"weather.gusts.0.speed"`).
Use `reduce=False` when you only need the control-system results (faster). Everything is seeded: the same settings give the same numbers.

## Settings you will change
| Key | Meaning |
|---|---|
| `telescope.dish_m`, `efficiency` | dish size; sets the beams (Ku 70λ/D, L-band) and the effective area |
| `pointing_truth` | HIDDEN pointing errors (IA, IE, NPAE, AN, AW in degrees) the system must discover |
| `mount.*` | speeds, acceleration, loop gain `kp`, actuator lag `tau_act`, wind jitter, limits |
| `weather.gusts`, `stow_limit_ms`, `resume_limit_ms` | wind events and the safety thresholds |
| `observations` | the blocks: `type`, `priority`, `earliest`, `deadline`, `fixed_center`, `duration_min`, `params` |
| `satellites.catalogue_err_km` | accuracy of the published orbital elements |
| `receiver.*` | noise, gain drift, Sun flux, beacon level, M31 flux |
| `session.seed` | noise realisation (use for Monte Carlo) |

## Reference results (unmodified `config/session.yaml`, seed 7)
| Quantity | Value |
|---|---|
| Pointing: true sky error before → after PR1 | 0.384° → 0.016° |
| Pointing: verification run PR2 (model not refitted) | 0.016° |
| Servo following error, sidereal tracking | 0.0015° rms |
| True pointing error while tracking the satellite | 0.19° rms (catalogue error); efficiency 0.94 |
| Time on the satellite above half power: tracked vs fixed beam | 244 s vs 1.9 s |
| Sun cross-scan: beam FWHM · pointing offset · Sun height · T_sys | 1.328° · (+0.003°, −0.021°) · 872 K · 128 K |
| M31 tracked ON/OFF: v_sys · line S/N | −292.8 ± 18.7 km/s · 5.6 |
| M31 fixed-dish drift: v_sys · line S/N | −285.9 ± 19.3 km/s · 5.8 |
| Weather stow | 1 event, 21 min; the H I block is re-queued and completes |

## Plugging in real recordings (Track R)
`rtops.ingest.ingest_solar(arc, cfg, npz, meta, obid)` and `ingest_hi(...)` write the same kind of FITS block with `ORIGIN=REAL`. `reduce.py` then works unchanged.
Formats: solar `npz` with `unix, counts, hot, cold` and JSON `az_deg, el_deg, el_offset_set_deg, freq_ghz`; H I `npz` with `freq_hz, unix, spec, hot` (+ `is_on`)
and JSON `mode` (`switch` or `drift`), `target`, and `az_deg, el_deg` (drift) or `off` (switch). The round trip is tested in `tests/`.

## How v_sys and its error are obtained
The receding edge of M31 is blended with Milky Way hydrogen, so the data fix the LEFT edge (v_sys − V sin i) well and v_sys itself only through a prior on V sin i
(`hi.vrot_prior_kms`, 245 ± 25, deliberately 16 km/s above the injected 229). The error is the standard deviation of the posterior over v_sys built from the
profile likelihood, not a local curvature (which was found by Monte Carlo to be wildly optimistic at S/N of a few). Expect a +15 km/s systematic from the prior.

## What this does and does not prove
The simulator and the control software share conventions (Gaussian beams, the pointing-model form, refraction formula, SGP4), so recovery of the hidden
terms is optimistic; real telescopes add flexure and other terms. Satellite and beacon brightness are assumptions. All results are validation of the
software on simulated data, never observations.
