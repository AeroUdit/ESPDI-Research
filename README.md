ESPDI — Analytical Powered-Descent Initiation

Author: Udit Karthik

An analytical powered-descent guidance framework for variable-mass rocket landing.

Research Paper · Simulation · Optimizer Validation

Overview

Powered landing of a reusable rocket is a constrained guidance problem. The vehicle must convert an off-nominal position and velocity into a precise terminal state while respecting finite thrust and propellant.

Numerical optimal-control methods can solve this problem, but they require iterative computation.

ESPDI (Enhanced Simplified Powered-Descent Initiation) investigates whether the essential ignition decision can instead be obtained analytically for a defined class of powered-descent problems.

The primary contribution is ESPDI-T, a three-dimensional analytical translational method that calculates:

the required burn time;

the constant inertial thrust direction;

the corresponding ignition point.

A secondary extension, ESPDI-A, introduces a terminal thrust-direction condition while retaining an analytical trajectory construction.

The numerical optimizer in this repository is an independent benchmark. It is not part of ESPDI.

Research Question

Can a variable-mass analytical powered-descent initiation method reproduce the minimum-propellant translational solution of numerical optimal control without performing trajectory optimization online?

A secondary question is whether the same analytical framework can be extended to impose a terminal thrust-direction requirement while retaining low computational complexity.

ESPDI-T

ESPDI-T solves the drag-free, constant-gravity, variable-mass point-mass landing problem.

1. Variable-mass model

For maximum thrust, the propellant mass-flow rate is

\dot{m}_p = \frac{T_{\max}}{I_{sp}g_0}

and the vehicle mass during the burn is

m(t) = m_0 - \dot{m}_p t

The final mass must satisfy

m_f \geq m_{\mathrm{dry}}

2. Required velocity change

For zero terminal velocity, the required thrust-generated velocity change is

\mathbf{D}(t_b) = -\mathbf{v}_0 - \mathbf{g}t_b

where t_b is the burn time.

3. Variable-mass delta-v

The available thrust-generated delta-v is

\Delta v_T =
I_{sp}g_0
\ln\left(\frac{m_0}{m_f}\right)

Therefore, ESPDI-T solves the scalar equation

I_{sp}g_0
\ln\left(
\frac{m_0}
{m_0-\dot{m}_p t_b}
\right)
=
\left\|
-\mathbf{v}_0-\mathbf{g}t_b
\right\|

No trajectory-level optimization is performed.

4. Constant thrust direction

Once the burn time is known,

\hat{\mathbf{u}}_{\mathrm{ESPDI}}
=
\frac{
-\mathbf{v}_0-\mathbf{g}t_b
}{
\left\|
-\mathbf{v}_0-\mathbf{g}t_b
\right\|
}

5. Ignition point

Define the variable-mass position coefficient

K =
I_{sp}g_0
\left[
\left(
t_b-\frac{m_0}{\dot{m}_p}
\right)
\ln\left(\frac{m_0}{m_f}\right)
+
t_b
\right]

The ignition point is then

\mathbf{r}_I =
\mathbf{r}_L
-
\mathbf{v}_0 t_b
-
K\hat{\mathbf{u}}_{\mathrm{ESPDI}}
-
\frac{1}{2}\mathbf{g}t_b^2

ESPDI-A

ESPDI-A extends the analytical framework by imposing a terminal thrust-direction condition.

The acceleration profile is represented as

\mathbf{a}(\tau)
=
\mathbf{c}_0
+
\mathbf{c}_1\tau
+
\mathbf{c}_2\tau^2,
\qquad
\tau=\frac{t}{t_f}

with terminal position and velocity constraints

\mathbf{r}(t_f)=\mathbf{0}

\mathbf{v}(t_f)=\mathbf{0}

Define

\Delta\mathbf{r}
=
-\mathbf{r}_0-\mathbf{v}_0 t_f

and

\Delta\mathbf{v}
=
-\mathbf{v}_0

The analytical coefficients are

\mathbf{c}_0
=
\frac{12\Delta\mathbf{r}}{t_f^2}
-
\frac{6\Delta\mathbf{v}}{t_f}
+
\mathbf{a}_f

\mathbf{c}_1
=
-\frac{48\Delta\mathbf{r}}{t_f^2}
+
\frac{30\Delta\mathbf{v}}{t_f}
-
6\mathbf{a}_f

\mathbf{c}_2
=
\frac{36\Delta\mathbf{r}}{t_f^2}
-
\frac{24\Delta\mathbf{v}}{t_f}
+
6\mathbf{a}_f

For an upright terminal thrust direction,

\mathbf{a}_{T,f}
=
\mathbf{a}_f-\mathbf{g}
=
\frac{T_f}{m_f}\hat{\mathbf{b}}_{3,f}

The required thrust magnitude is

T_{\mathrm{req}}(t)
=
m(t)
\left\|
\mathbf{a}(t)-\mathbf{g}
\right\|

and the mass equation becomes

\dot{m}
=
-\frac{
m(t)\left\|
\mathbf{a}(t)-\mathbf{g}
\right\|
}{
I_{sp}g_0
}

ESPDI-A determines the feasible final time with a scalar feasibility calculation rather than a vector-valued trajectory optimizer.

Numerical Optimal-Control Benchmark

The optimizer is an independent benchmark, not part of ESPDI.

It minimizes propellant:

J =
\frac{1}{I_{sp}g_0}
\int_0^{t_f}
\left\|
\mathbf{T}(t)
\right\|dt

subject to

\dot{\mathbf{r}}=\mathbf{v}

\dot{\mathbf{v}}
=
\frac{\mathbf{T}}{m}
+
\mathbf{g}

\dot{m}
=
-\frac{\|\mathbf{T}\|}
{I_{sp}g_0}

and the terminal constraints

\mathbf{r}(t_f)=\mathbf{0}

\mathbf{v}(t_f)=\mathbf{0}

with

\|\mathbf{T}(t)\| \leq T_{\max}

and

m(t_f)\geq m_{\mathrm{dry}}

The optimizer is tested from multiple initial guesses and multiple discretization meshes so that the comparison is not based only on an ESPDI-provided initial trajectory.

Final Simulation

The simulator uses a three-panel 3D visualization:

ESPDI-T

Optimizer

ESPDI-A

Analytical translation

Independent benchmark

Analytical extension

Each method has its own physical trajectory and telemetry.

The visualization includes:

3D trajectory

moving rocket

landing pad at (0, 0, 0)

trajectory trail

ignition point

position and velocity

terminal velocity

thrust magnitude

guidance angle

fuel used

landing status

Vehicle Model

Parameter

Value

Initial mass

12,000 kg

Dry mass

4,500 kg

Maximum thrust

180 kN

Specific impulse

320 s

Initial T/W

1.5291

Gravity

9.81 m/s²

Drag

Excluded

Model

3D variable-mass point mass

The baseline model intentionally excludes aerodynamic drag and full rotational dynamics so that the primary analytical translational problem remains tractable.

Final Monte Carlo Results

The final performance dataset contains:

500 Pair-T cases

500 Pair-A cases

ESPDI-T vs Optimizer-T

ESPDI-T success: 500/500

Optimizer-T success: 500/500

Mean ESPDI-T fuel: 874.040 kg

Mean Optimizer-T fuel: 874.040 kg

Mean paired fuel gap: -5.28 × 10⁻¹¹ %

Mean ESPDI-T runtime: 0.558 s

Mean Optimizer-T runtime: 13.316 s

The translational fuel difference is at numerical roundoff scale.

The mean computation-time ratio is approximately

\frac{13.316}{0.558} \approx 23.9

so ESPDI-T required about 23.9 times less mean computation time in this implementation.

ESPDI-A vs Optimizer-A

ESPDI-A analytical cases: 500/500

Optimizer-A converged cases: 499/500

Mean ESPDI-A fuel: 1175.485 kg

Mean Optimizer-A fuel: 1131.205 kg

Mean fuel gap: 3.894%

Median fuel gap: 3.652%

Fuel-gap range: 0.221%–7.614%

Mean ESPDI-A runtime: 4.003 s

Mean Optimizer-A runtime: 24.083 s

The single missing Optimizer-A result is retained as a numerical benchmark nonconvergence rather than being treated as a physical failure.

Independent Optimizer-T Validation

A separate validation program reproduces the Pair-T state-generation process using a fixed random seed.

The validation checks two things:

Same-seed control: the optimizer can reproduce the known ESPDI-T solution when initialized with the ESPDI trajectory.

Independent rediscovery: the optimizer can reach essentially the same solution from non-ESPDI initial thrust histories.

The final validation dataset contains:

20 Pair-T states

2 tested meshes: N = 8 and N = 12

2 independent non-ESPDI starts at each mesh

same-seed ESPDI control runs

The best independent optimizer solution across the 20 validation states differed from ESPDI-T by less than approximately

10^{-6}\%

in fuel at the tested meshes.

This supports the conclusion that the optimizer benchmark is not simply reproducing ESPDI because ESPDI was supplied as its initial trajectory.

Repository Structure

ESPDI/
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

How to Run the Simulation

Install the required packages:

pip install numpy scipy matplotlib

Then run:

python simulation/espdi_simulation_repaired.py

The simulator contains settings for single-case testing and Monte Carlo evaluation.

How to Run the Optimizer Validation

Run:

python validation/validate_optimizer_t_v3.py

The validation program records:

optimizer convergence

physical feasibility

propellant used

terminal position error

terminal velocity error

final mass

maximum thrust ratio

iteration count

runtime

Figures

ESPDI-T Fuel Gap



ESPDI-A Fuel Gap



Computation Time



Independent Optimizer Validation



ESPDI Architecture



Research Paper

Read the Full ESPDI Research Paper

Limitations

The current research does not include:

aerodynamic drag

aerodynamic moments

full quaternion attitude dynamics

time-varying inertia

detailed center-of-mass migration

engine ignition delay

actuator lag

gimbal-rate limits

propellant slosh

terrain constraints

atmospheric uncertainty

ESPDI-T is therefore not claimed to be a universal replacement for optimal control.

Its claim is restricted to the defined drag-free, constant-gravity, variable-mass translational problem.

ESPDI-A is a terminal-thrust-direction extension rather than a complete 6DOF attitude controller.

Future Work

Potential extensions include:

uncertainty and robustness analysis

aerodynamic drag

full 6DOF attitude dynamics

variable inertia and center-of-mass migration

actuator constraints

terrain and glide-slope constraints

higher-fidelity propulsion and atmosphere models

additional powered-descent guidance benchmarks

Citation

If you reference this project:

Karthik, U.
ESPDI: Analytical Powered-Descent Initiation for Variable-Mass Rocket Landing.
Research Project, 2026.

Author

Udit Karthik

ESPDI is an independent research project investigating the mathematical structure and computational efficiency of analytical powered-descent guidance.
