Converting Point Forecasts to probabilistic wind power scenarios for a set of wind farms

- [ ] Data
  - [x] New sample of data with 35+ wind farms and less times
  - [ ] Left 6 wf for oos checks

- [ ] Figures
  - [ ] Regime behaviour
  - [ ] Seasonal patterns
  - [ ] Spatial correlation in errors

- [ ] Fit models
  - [x] Job files for each
  - [x] PR-NEN (1)
  - [x] PR-AR-NEN (2)
  - [x] ST-NEN (3)
  - [x] ST-PR-NEN (4)
  - [ ] ST2-PR-NEN (non-separable)
  - [ ] Regime switching

- [ ] Plot effects
  - [x] Effects used before
  - [ ] Spatial related effects
    - [x] Extract field at locations used in fitted next 24 hours
    - [ ] Extract field at unseen locations next 24 hours

- [ ] Generate scenarios
  - [x] Test/ fix code to get samples from spatial models
  - [x] Adapt code to get plots from all locations
  - [ ] Nonparametric scenarios for all locations
  - [ ] PR-NEN
  - [ ] PR-AR-NEN
  - [ ] ST-NEN
  - [x] ST-PR-NEN


- [ ] Scenario scoring
  - [ ] Ramp distribution
  - [ ] CRPS, Energy, Variogram
  - [ ] PIT diagrams
  - [ ] Calibration diagrams
  - [ ] Spatial OOS