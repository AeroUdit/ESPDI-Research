# ESPDI — Analytical Powered-Descent Initiation

**Author:** Udit Karthik

> An analytical powered-descent guidance framework for variable-mass rocket landing.

[![Research Paper](https://img.shields.io/badge/Research%20Paper-PDF-blue)](paper/ESPDI_Final_Research_Paper.pdf)
[![Simulation](https://img.shields.io/badge/Simulation-Python-green)](simulation/espdi_simulation_repaired.py)
[![Validation](https://img.shields.io/badge/Optimizer%20Validation-Python-orange)](validation/validate_optimizer_t_v3.py)

---

## Overview

Powered landing of a reusable rocket is a constrained guidance problem: the vehicle must convert an off-nominal position and velocity into a precise terminal state while respecting finite thrust and propellant.

Numerical optimal-control methods can solve this problem, but they require iterative computation.

**ESPDI (Enhanced Simplified Powered-Descent Initiation)** investigates whether the essential ignition decision can instead be obtained analytically for a defined class of powered-descent problems.

The primary contribution is **ESPDI-T**, a three-dimensional analytical translational method that calculates:

1. The required burn time.
2. The constant inertial thrust direction.
3. The corresponding ignition point.

A secondary extension, **ESPDI-A**, introduces a terminal thrust-direction constraint while retaining an analytical trajectory construction.

The numerical optimizer in this repository is an **independent benchmark**. It is not part of ESPDI.

---

## Research Question

> **Can a variable-mass analytical powered-descent initiation method reproduce the minimum-propellant translational solution of numerical optimal control without performing trajectory optimization online?**

A secondary question is:

> **Can the analytical framework be extended to impose a terminal thrust-direction requirement while retaining low computational complexity?**

---

## ESPDI-T

ESPDI-T solves the drag-free, constant-gravity, variable-mass point-mass landing problem.

The core calculation is:

### 1. Burn-time equation

The available variable-mass delta-v is

[
\Delta v_T =
I_{sp}g_0
\ln\left(\frac{m_0}{m_f}\right)
]

and the required velocity change is

[
\mathbf D(t_b)
==============

-\mathbf v_0-\mathbf g t_b.
]

Therefore the burn time is obtained from

[
I_{sp}g_0
\ln\left[
\frac{m_0}
{m_0-\dot m_p t_b}
\right]
=======

\left|
-\mathbf v_0-\mathbf g t_b
\right|.
]

This is a **scalar nonlinear equation**, not a trajectory optimization problem.

### 2. Constant thrust direction

Once the burn time is known,

[
\hat{\mathbf u}_{ESPDI}
=======================

\frac{
-\mathbf v_0-\mathbf g t_b
}{
\left|
-\mathbf v_0-\mathbf g t_b
\right|
}.
]

### 3. Ignition point

The variable-mass position integral gives the ignition-point location:

[
\mathbf r_I
===========

\mathbf r_L
-\mathbf v_0t_b
-K\hat{\mathbf u}_{ESPDI}
-\frac12\mathbf g t_b^2.
]

Thus the method directly converts the current state into a recoverable ignition boundary.

---

## ESPDI-A

ESPDI-A extends the analytical framework by imposing a terminal thrust-direction condition.

Instead of a constant thrust vector, the method uses a deterministic quadratic acceleration profile:

[
\mathbf a(\tau)
===============

\mathbf c_0
+
\mathbf c_1\tau
+
\mathbf c_2\tau^2,
\qquad
\tau=\frac{t}{t_f}.
]

The coefficients are obtained analytically from terminal position, velocity, and acceleration boundary conditions.

The required thrust is then

[
T_{req}(t)
==========

m(t)\left|
\mathbf a(t)-\mathbf g
\right|,
]

with variable mass

[
\dot m
======

-\frac{T_{req}}
{I_{sp}g_0}.
]

ESPDI-A therefore remains deterministic and does not perform online vector-valued trajectory optimization.

---

## Numerical Optimal-Control Benchmark

The optimizer solves the same underlying physical problem independently.

Its objective is minimum propellant:

[
J=
\frac{1}{I_{sp}g_0}
\int_0^{t_f}
|\mathbf T(t)|,dt.
]

Subject to

[
\dot{\mathbf r}=\mathbf v,
]

[
\dot{\mathbf v}
===============

\frac{\mathbf T}{m}
+
\mathbf g,
]

[
\dot m
======

-\frac{|\mathbf T|}
{I_{sp}g_0},
]

and the terminal constraints

[
\mathbf r(t_f)=\mathbf 0,
\qquad
\mathbf v(t_f)=\mathbf 0,
]

with

[
|\mathbf T|\le T_{max},
\qquad
m(t_f)\ge m_{dry}.
]

The optimizer is used **only as a benchmark** to determine how closely ESPDI reproduces the numerical optimum.

---

## Final Simulation

The simulation uses a three-panel 3D visualization:

| Left    | Center    | Right   |
| ------- | --------- | ------- |
| ESPDI-T | Optimizer | ESPDI-A |

Each method has its own physical trajectory.

The visualization includes:

* 3D trajectory
* moving rocket
* landing pad at ((0,0,0))
* ignition point
* position and velocity
* terminal velocity
* thrust magnitude
* fuel used
* guidance angle
* landing status

The simulation uses the same initial state for each comparison pair.

---

## Vehicle Model

| Parameter        |                       Value |
| ---------------- | --------------------------: |
| Initial mass     |                   12,000 kg |
| Dry mass         |                    4,500 kg |
| Maximum thrust   |                      180 kN |
| Specific impulse |                       320 s |
| Initial T/W      |                      1.5291 |
| Gravity          |                   9.81 m/s² |
| Drag             |                    Excluded |
| Model            | 3D variable-mass point mass |

The baseline model intentionally excludes aerodynamic drag and full rotational dynamics so that the primary analytical translational problem remains tractable and computationally lightweight.

---

## Final Monte Carlo Results

The final performance study contains:

* **500 Pair-T cases**
* **500 Pair-A cases**

### ESPDI-T vs Optimizer-T

Across all 500 paired translational cases:

* ESPDI-T success: **500/500**
* Optimizer-T success: **500/500**
* Mean ESPDI-T fuel: **874.040 kg**
* Mean Optimizer-T fuel: **874.040 kg**
* Mean fuel gap: **approximately (-5.28\times10^{-11}%)**
* Mean ESPDI-T runtime: **0.558 s**
* Mean Optimizer-T runtime: **13.316 s**

This corresponds to an average runtime ratio of approximately

[
\boxed{23.9\times}
]

in favor of ESPDI-T in the tested implementation.

The fuel difference is at numerical roundoff scale.

### ESPDI-A vs Optimizer-A

The Pair-A benchmark produced:

* ESPDI-A analytical cases: **500/500**
* Optimizer-A converged cases: **499/500**
* Mean ESPDI-A fuel: **1175.485 kg**
* Mean Optimizer-A fuel: **1131.205 kg**
* Mean fuel gap: **3.894%**
* Median fuel gap: **3.652%**
* Fuel-gap range: **0.221%–7.614%**
* Mean ESPDI-A runtime: **4.003 s**
* Mean Optimizer-A runtime: **24.083 s**

The results show that ESPDI-A retains a substantial computational advantage but does not generally reproduce the fully optimized terminal-direction solution exactly.

---

## Independent Optimizer-T Validation

A separate validation program was used to test whether the optimizer was simply reproducing ESPDI because ESPDI was supplied as its initial guess.

The validation:

* reproduced the main simulator's Pair-T state generator;
* used the same fixed seed;
* verified state agreement;
* tested an explicit ESPDI-seed control;
* tested non-ESPDI initializations;
* tested multiple discretization meshes.

The validation dataset contains **20 Pair-T states**, with two non-ESPDI starts at each of two meshes.

The best independent solution remained within approximately:

[
10^{-7}%
]

of ESPDI-T fuel usage across the tested validation states.

This supports the interpretation that the numerical optimizer independently rediscovers the same translational solution rather than merely being forced to reproduce ESPDI.

See:

[`validation/validate_optimizer_t_v3.py`](validation/validate_optimizer_t_v3.py)

and:

[`results/optimizer_t_validation_v3_results.csv`](results/optimizer_t_validation_v3_results.csv)

---

## Results and Figures

### ESPDI-T Fuel Gap

![ESPDI-T Fuel Gap](figures/fig_t_fuel_gap.png)

### ESPDI-A Fuel Gap

![ESPDI-A Fuel Gap](figures/fig_a_fuel_gap.png)

### Computation Time

![Computation Time](figures/fig_runtime.png)

### Independent Optimizer Validation

![Optimizer Validation](figures/fig_validation.png)

### ESPDI Architecture

![ESPDI Architecture](figures/fig_architecture.png)

---

## Repository Structure

```text
ESPDI/
│
├── README.md
│
├── simulation/
│   └── espdi_simulation_repaired.py
│
├── validation/
│   └── validate_optimizer_t_v3.py
│
├── results/
│   ├── espdi_monte_carlo_results.csv
│   └── optimizer_t_validation_v3_results.csv
│
├── figures/
│   ├── fig_architecture.png
│   ├── fig_t_fuel_gap.png
│   ├── fig_a_fuel_gap.png
│   ├── fig_runtime.png
│   └── fig_validation.png
│
└── paper/
    └── ESPDI_Final_Research_Paper.pdf
```

---

## Running the Simulation

Install the required Python packages:

```bash
pip install numpy scipy matplotlib
```

Then run:

```bash
python simulation/espdi_simulation_repaired.py
```

The simulator can run a single deterministic case or the Monte Carlo experiment by changing the settings at the top of the script.

---

## Running the Optimizer Validation

Run:

```bash
python validation/validate_optimizer_t_v3.py
```

The validation script performs independent optimizer tests and records:

* convergence status
* physical feasibility
* fuel usage
* terminal position error
* terminal velocity error
* final mass
* thrust-limit compliance
* optimization iterations
* runtime

---

## Research Paper

The complete mathematical derivation, methodology, experiments, validation, results, limitations, and discussion are available here:

**[Read the Full Research Paper](paper/ESPDI_Final_Research_Paper.pdf)**

---

## Limitations

ESPDI is evaluated only within the defined baseline model.

The current research does not include:

* aerodynamic drag
* aerodynamic moments
* full quaternion attitude dynamics
* time-varying inertia
* detailed center-of-mass migration
* engine ignition delay
* actuator lag
* gimbal-rate limits
* propellant slosh
* terrain constraints
* atmospheric uncertainty

ESPDI-T is therefore **not claimed to be a universal replacement for optimal control**.

Its claim is restricted to the defined drag-free, constant-gravity, variable-mass translational problem.

ESPDI-A is an analytical terminal-direction extension rather than a complete 6DOF attitude controller.

---

## Future Work

Potential extensions include:

* uncertainty and robustness analysis
* aerodynamic drag
* full 6DOF attitude dynamics
* variable inertia and center-of-mass migration
* actuator constraints
* terrain and glide-slope constraints
* additional powered-descent guidance benchmarks
* higher-fidelity propulsion and atmospheric models

---

## Citation

If you reference this project, please cite:

```text
Karthik, U. ESPDI: Analytical Powered-Descent Initiation for
Variable-Mass Rocket Landing. Research Project, 2026.
```

---

## Author

**Udit Karthik**

ESPDI is an independent research project investigating the mathematical structure and computational efficiency of analytical powered-descent guidance.
