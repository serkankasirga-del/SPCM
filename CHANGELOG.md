# Changelog

## 0.11.0 (2026-10-15) - pre-release, in testing
- **Resolution from a sharp edge**: *Measure resolution* above the maps puts a green line
  on the reflection map; the derivative of the reflection along it is fitted with a
  Gaussian and its FWHM is the resolution. Live while dragging and during a scan;
  *Save PNG + CSV*
- **Line × parameter scans**: a line drawn on the camera (or a map) is scanned once per
  value of the bias, the laser power on the sample, or the chopping frequency, giving a
  distance × parameter image. Power and frequency are put back afterwards
- Fixes from a code review:
  - the window no longer talks to the stage or lock-ins while a scan, autofocus or
    calibration runs (could mix up replies and stop the scan)
  - *Set power now* respects *Never set SINE OUT above*
  - moves longer than 500 µm from a map double-click / *Go to start* ask first
  - saving *Hardware → Settings* no longer resets scan size, camera flip, etc.
  - stopping a fly scan during its speed test no longer leaves a broken .h5
  - focus (z) moves wait until z is really at the target (autofocus)
  - the window fits smaller screens (~1240 px wide); nothing cut off on the left
  - *Stop & discard* no longer shows an error if the scan just finished

## 0.10.0 (2026-10-15) - pre-release, in testing
- **Point sweeps** on one spot: power, bias, or bias × power (bias outer). Powers set with the
  2 mV SINE OUT step closest to the target; N readings per point with mean and spread;
  live linear and log-log plots with the fitted exponent α; bias–power map; .h5 + CSV + PNG
- **Line scans**: draw a line on a map or on the camera image, scan along it at any angle;
  with a bias list you get a distance × bias image
- **Camera calibration**: one button; afterwards double-click the camera to move the laser,
  and draw scan frames and lines on the live camera image
- Lists accept `1:100:*2` (multiply each step) and always include the stop value
- `tools/selftest.py`: tests everything in simulation; runs on GitHub for every change

## 0.9.0 (2026-10-14) - in testing
First version in the repository. Includes:
- Settings: your own `config.yaml` is created on first start from `config.default.yaml`
  and is never touched by updates (unzip or pull over the old folder safely)
- Camera: the Thorlabs camera package is found automatically in the ThorCam installation,
  so the camera now also works in SPCM.exe and without the pip step
- README: updating = Download ZIP, unzip over the old folder, double-click Install SPCM.bat
- Step and Fly scans with the Märzhäuser Tango; two SR830s (photocurrent/photovoltage + reflection)
- SR570 (with bias sweeps) / SR560 / lock-in-only signal conversion (× √2, peak amplitude)
- Origin (0, 0) at the current stage position; axes matching the joystick; rotated scans;
  draw a frame on a map and scan it
- Safety: no moves outside the scan area, arrival check, distance confirmation, stage log
- PM100USB power meter: live reading, fixed range, meter/sample ratio, set the laser power
  before a scan via the lock-in SINE OUT
- CS165CU camera: live view, flip, autofocus on the laser spot (Tango z)
- Tidy .h5 files (one group per map), CSV export, publication-style PNG per map
- Experiment notes, dark theme, number fields with SI prefixes
