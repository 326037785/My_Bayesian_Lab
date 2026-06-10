# Bayesian Filter

Multi-target tracking framework with GM-PHD, EKF/UKF/CKF, JPDA, MHT, and particle filters.

## Features

### Filters
- **KF** — Kalman Filter (linear Gaussian)
- **EKF** — Extended Kalman Filter (Jacobian linearization)
- **UKF** — Unscented Kalman Filter (sigma points)
- **CKF** — Cubature Kalman Filter (cubature points)
- **PF / APF / RBPF / UPF** — Particle filter variants

### Multi-Target Data Association
- **NN / KNN** — Nearest Neighbor
- **JPDA** — Joint Probabilistic Data Association
- **MHT** — Multiple Hypothesis Tracking
- **GM-PHD** — Gaussian Mixture PHD filter with:
  - Fixed birth (MATLAB-style peripheral strategy, dynamically scaled)
  - Adaptive birth (Vo early gating-based approach)
  - KF / EKF / UKF / CKF backends

### Measurement Models
- Linear (Cartesian position)
- Polar (range + bearing)
- Clutter: Poisson / Uniform / Gaussian

### Ground Truth
- CV / CA / CT / RW motion models
- Multi-target scenario generation

### Metrics
- RMSE / OSPA / GOSPA

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

```bash
# Linear tracking
python run_demo.py --demo linear

# Nonlinear tracking (polar measurements)
python run_demo.py --demo nonlinear --filter UKF

# Multi-target with PHD
python run_demo.py --demo multi --association phd

# Multi-target with JPDA
python run_demo.py --demo multi --association jpda

# Performance comparison
python run_demo.py --demo compare --scenario nonlinear
```

## API Usage

```python
from data_association import PHDFilter, GaussianComponent
from ground_truth import ScenarioManager
from measurements import MeasurementSimulator, MeasurementType, LinearMeasurementNoise, create_uniform_clutter
import numpy as np

# Create scenario
scenario = ScenarioManager(time_step=1.0, process_noise_std=0.1)
scenario.create_linear_scenario(n_targets=2, duration=50.0,
                                x_range=(-500, 500), y_range=(-500, 500))
scenario_data = scenario.generate_scenario(50.0)

# Generate measurements
sim = MeasurementSimulator(MeasurementType.LINEAR,
    LinearMeasurementNoise(1.0, 1.0),
    create_uniform_clutter(5.0),
    detection_probability=0.95)
sim.set_surveillance_region((-600, 600), (-600, 600))
measurements = sim.generate_scenario_measurements(scenario_data, 1.0)

# Configure PHD with scene-appropriate birth
birth = [
    GaussianComponent(weight=0.025,
        mean=np.array([x, 0., y, 0.]),
        covariance=np.diag([150., 15., 150., 15.]) ** 2)
    for x in [-200, 0, 200] for y in [-200, 0, 200]
]

phd = PHDFilter(
    birth_components=birth,
    surveillance_bounds=((-600, 600), (-600, 600)),
    detection_probability=0.95,
    clutter_rate=5.0,
    measurement_noise_matrix=np.diag([1., 1.]),
    filter_type='KF'
)

# Run PHD filter
for t in sorted(measurements.keys()):
    phd.predict(1.0)
    ms = measurements[t]
    if ms.size > 0:
        phd.update(ms.get_measurements_array())

# Extract estimates
states, n_targets = phd.extract_states()
print(f"Estimated {n_targets} targets")
```

## Project Structure

```
BayesianFilter/
  data_association/   PHD, JPDA, MHT, NN/KNN
  filter/             KF, EKF, UKF, CKF, PF, FilterBackend
  ground_truth/       Scenario generation, motion models
  measurements/       Measurement simulation, clutter
  demo/               Demo scripts
  utils/              Coordinate transforms, math utilities
  visualize/          OSPA, RMSE, plotting
```

## References

1. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis density filter. *IEEE TSP*, 54(11), 4091-4104.
2. Arulampalam, M. S., et al. (2002). A tutorial on particle filters. *IEEE TSP*, 50(2), 174-188.
3. Ristic, B., et al. (2012). Adaptive target birth intensity for PHD and CPHD filters. *IEEE TAES*, 48(2), 1656-1668.

## License

MIT
