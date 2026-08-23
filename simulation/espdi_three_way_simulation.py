"""
ESPDI RESEARCH SIMULATOR — ROBUST VERSION
==========================================

Pair T:
    ESPDI-T and Optimizer-T start from the EXACT SAME analytical ESPDI-T
    ignition state.

Pair A:
    ESPDI-A and Optimizer-A start from a SEPARATE, identical Pair-A state.

Animation:
    LEFT   = ESPDI-T
    MIDDLE = Optimizer-T
    RIGHT  = ESPDI-A

The Optimizer-A calculation is also performed and reported, but the
three-panel animation is intentionally the clean research comparison:
ESPDI-T / Optimizer-T / ESPDI-A.

Model:
    3D point mass
    constant gravity
    no drag
    variable mass
    constant Isp
    thrust magnitude <= Tmax

Required packages:
    numpy
    scipy
    matplotlib
"""

from dataclasses import dataclass
from pathlib import Path
import csv
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.optimize import brentq, minimize
from scipy.integrate import solve_ivp


# ============================================================
# USER SETTINGS
# ============================================================

RANDOM_SEED = 20260818

# FIRST RUN: keep False.
RUN_MONTE_CARLO = True

# Recommended final benchmark after convergence testing:
MONTE_CARLO_CASES = 500

# Small first-pass benchmark. Increase after convergence testing.
N_OPT_T = 12
N_OPT_A = 12

# Multi-starts are intentionally small. The first guess is a
# known feasible ESPDI solution.
OPTIMIZER_STARTS = 2

SAVE_RESULTS = True

# ------------------------------------------------------------
# Convergence/constraint tolerances
# ------------------------------------------------------------

POSITION_TOL = 0.10       # m
VELOCITY_TOL = 0.05       # m/s
MASS_TOL = 1e-7           # kg

# Scaling used by the optimizer. Scaling is crucial for SLSQP.
POSITION_SCALE = 1000.0   # m
VELOCITY_SCALE = 100.0    # m/s
MASS_SCALE = 1000.0       # kg


# ============================================================
# VEHICLE
# ============================================================

G = 9.81
G0 = 9.80665

M0 = 12_000.0
M_DRY = 4_500.0
T_MAX = 180_000.0
ISP = 320.0

MDOT_MAX = T_MAX / (ISP * G0)

TARGET = np.zeros(3)

INITIAL_TW = T_MAX / (M0 * G)


# ============================================================
# DATA TYPES
# ============================================================

@dataclass
class Case:
    position: np.ndarray
    velocity: np.ndarray
    mass: float
    label: str


@dataclass
class Trajectory:
    name: str
    t: np.ndarray
    r: np.ndarray
    v: np.ndarray
    mass: np.ndarray
    thrust: np.ndarray
    u: np.ndarray
    angle_deg: np.ndarray
    final_r: np.ndarray
    final_v: np.ndarray
    prop_used: float
    ignition_time: float
    ignition_position: np.ndarray
    success: bool


@dataclass
class OptimizationInfo:
    converged: bool
    message: str
    iterations: int
    objective: float
    position_residual: float
    velocity_residual: float
    max_thrust_ratio: float
    final_mass: float


# ============================================================
# BASIC PHYSICS
# ============================================================

def gvec():
    return np.array([0.0, 0.0, -G])


def mass_after_burn(t, m_initial):
    return m_initial - MDOT_MAX * t


def max_propellant_time(m_initial):
    return max(
        0.0,
        (m_initial - M_DRY) / MDOT_MAX
    )


def thrust_delta_v(t_b, m_initial):
    mf = mass_after_burn(t_b, m_initial)
    if mf <= 0.0:
        return np.nan
    return ISP * G0 * np.log(m_initial / mf)


def required_dv_vector(t_b, velocity):
    return -velocity - gvec() * t_b


def angle_deg(u):
    horizontal = np.hypot(u[0], u[1])
    return float(np.degrees(np.arctan2(u[2], horizontal)))


# ============================================================
# ESPDI-T
# ============================================================

def solve_espdi_t_burn_time(velocity, m_initial):
    t_max = max_propellant_time(m_initial)
    if t_max <= 1e-8:
        return None

    def f(t):
        mf = mass_after_burn(t, m_initial)
        if mf < M_DRY:
            return np.nan

        required = np.linalg.norm(
            required_dv_vector(t, velocity)
        )
        available = thrust_delta_v(t, m_initial)
        return available - required

    grid = np.linspace(
        1e-6,
        t_max * (1.0 - 1e-10),
        1800
    )
    vals = np.array([f(t) for t in grid])

    for a, b, fa, fb in zip(
        grid[:-1],
        grid[1:],
        vals[:-1],
        vals[1:]
    ):
        if not (np.isfinite(fa) and np.isfinite(fb)):
            continue
        if fa == 0.0:
            return float(a)
        if fa * fb < 0.0:
            return float(
                brentq(
                    f,
                    a,
                    b,
                    xtol=1e-11,
                    rtol=1e-11
                )
            )
    return None


def espdi_t_solution(case):
    tb = solve_espdi_t_burn_time(
        case.velocity,
        case.mass
    )
    if tb is None:
        return None

    mf = mass_after_burn(tb, case.mass)
    D = required_dv_vector(tb, case.velocity)
    dnorm = np.linalg.norm(D)

    if dnorm < 1e-12:
        return None

    u = D / dnorm
    mdot = MDOT_MAX

    K = ISP * G0 * (
        (tb - case.mass / mdot)
        * np.log(case.mass / mf)
        + tb
    )

    r_ign = (
        TARGET
        - case.velocity * tb
        - K * u
        - 0.5 * gvec() * tb * tb
    )

    return tb, mf, u, r_ign


def propagate_espdi_t(case):
    sol = espdi_t_solution(case)
    if sol is None:
        return None

    tb, _, u, r_ign = sol

    mismatch = np.linalg.norm(
        case.position - r_ign
    )
    if mismatch > 1e-5:
        raise RuntimeError(
            "Pair-T state is not exactly its analytical ESPDI-T ignition point: "
            f"{mismatch:.6e} m."
        )

    def rhs(t, y):
        r = y[:3]
        v = y[3:6]
        m = max(y[6], M_DRY)

        a = T_MAX / m * u + gvec()
        return np.hstack((v, a, -MDOT_MAX))

    y0 = np.hstack(
        (case.position, case.velocity, case.mass)
    )

    soln = solve_ivp(
        rhs,
        [0.0, tb],
        y0,
        method="DOP853",
        rtol=1e-10,
        atol=1e-11,
        max_step=0.025
    )

    t = soln.t
    r = soln.y[:3].T
    v = soln.y[3:6].T
    mass = soln.y[6]

    u_hist = np.repeat(
        u[None, :],
        len(t),
        axis=0
    )

    thrust = np.full(
        len(t),
        T_MAX
    )

    angles = np.array([
        angle_deg(ui) for ui in u_hist
    ])

    return Trajectory(
        name="ESPDI-T",
        t=t,
        r=r,
        v=v,
        mass=mass,
        thrust=thrust,
        u=u_hist,
        angle_deg=angles,
        final_r=r[-1],
        final_v=v[-1],
        prop_used=case.mass - mass[-1],
        ignition_time=0.0,
        ignition_position=case.position.copy(),
        success=bool(
            np.linalg.norm(r[-1]) < 1e-4
            and np.linalg.norm(v[-1]) < 1e-6
            and mass[-1] >= M_DRY
        )
    )


# ============================================================
# ESPDI-A
# ============================================================

def quadratic_coefficients(tf, r0, v0, af):
    dr = -r0 - v0 * tf
    dv = -v0

    c0 = 12.0 * dr / tf**2 - 6.0 * dv / tf + af
    c1 = -48.0 * dr / tf**2 + 30.0 * dv / tf - 6.0 * af
    c2 = 36.0 * dr / tf**2 - 24.0 * dv / tf + 6.0 * af

    return c0, c1, c2


def espdi_a_profile(case, tf, n=900):
    t = np.linspace(0.0, tf, n)
    tau = t / tf

    mg = max(case.mass * 0.95, M_DRY + 1e-6)

    for _ in range(100):
        terminal_acc = T_MAX / mg
        af = gvec() + np.array(
            [0.0, 0.0, terminal_acc]
        )

        c0, c1, c2 = quadratic_coefficients(
            tf,
            case.position,
            case.velocity,
            af
        )

        a_net = (
            c0[None, :]
            + c1[None, :] * tau[:, None]
            + c2[None, :] * tau[:, None]**2
        )

        a_thrust = a_net - gvec()
        a_mag = np.linalg.norm(a_thrust, axis=1)

        # dm/dt = -m*a_thrust/(Isp*g0)
        lam = a_mag / (ISP * G0)

        # Accurate enough trapezoidal integration for the scalar mass ODE.
        integral = np.zeros_like(t)
        integral[1:] = np.cumsum(
            0.5 * (lam[1:] + lam[:-1])
            * np.diff(t)
        )

        mass = case.mass * np.exp(-integral)
        new_mg = mass[-1]

        if abs(new_mg - mg) < 1e-9:
            break

        mg = 0.5 * mg + 0.5 * new_mg

    thrust = mass * a_mag
    u = a_thrust / np.maximum(a_mag[:, None], 1e-12)

    # Exact analytical state reconstruction.
    tt = t[:, None]
    r = (
        case.position[None, :]
        + case.velocity[None, :] * tt
        + 0.5 * c0[None, :] * tt**2
        + c1[None, :] * tt**3 / (6.0 * tf)
        + c2[None, :] * tt**4 / (12.0 * tf**2)
    )

    v = (
        case.velocity[None, :]
        + c0[None, :] * tt
        + c1[None, :] * tt**2 / (2.0 * tf)
        + c2[None, :] * tt**3 / (3.0 * tf**2)
    )

    return t, r, v, mass, thrust, u


def espdi_a_margin(case, tf):
    _, _, _, mass, thrust, _ = espdi_a_profile(
        case,
        tf,
        n=500
    )
    return max(
        np.max(thrust) - T_MAX,
        100.0 * (M_DRY - np.min(mass))
    )


def solve_espdi_a(case):
    grid = np.linspace(2.0, 80.0, 240)
    vals = [espdi_a_margin(case, tf) for tf in grid]

    bracket = None
    for a, b, fa, fb in zip(
        grid[:-1],
        grid[1:],
        vals[:-1],
        vals[1:]
    ):
        if fa > 0.0 and fb <= 0.0:
            bracket = a, b
            break

    if bracket is None:
        return None

    tf = brentq(
        lambda x: espdi_a_margin(case, x),
        bracket[0],
        bracket[1],
        xtol=1e-7,
        rtol=1e-9
    )

    t, r, v, mass, thrust, u = espdi_a_profile(
        case,
        tf,
        n=1400
    )

    angles = np.array([angle_deg(ui) for ui in u])

    return Trajectory(
        name="ESPDI-A",
        t=t,
        r=r,
        v=v,
        mass=mass,
        thrust=thrust,
        u=u,
        angle_deg=angles,
        final_r=r[-1],
        final_v=v[-1],
        prop_used=case.mass - mass[-1],
        ignition_time=0.0,
        ignition_position=case.position.copy(),
        success=bool(
            np.linalg.norm(r[-1]) < 0.1
            and np.linalg.norm(v[-1]) < 0.01
            and np.max(thrust) <= T_MAX * 1.0001
            and np.min(mass) >= M_DRY
        )
    )


# ============================================================
# EXACT PROPAGATION OF A CONSTANT THRUST INTERVAL
# ============================================================

def exact_interval(r, v, m, Tvec, dt):
    """
    Exact interval propagation for constant Tvec over dt.

    This is the key repair for the optimizer.

    dm/dt = -T/(Isp*g0)
    dv/dt = Tvec/m + g
    dr/dt = v

    No midpoint approximation and no thrust clipping are used here.
    """
    T = np.linalg.norm(Tvec)

    if T < 1e-14:
        # Pure ballistic interval.
        r1 = r + v * dt + 0.5 * gvec() * dt**2
        v1 = v + gvec() * dt
        return r1, v1, m, T

    mdot = T / (ISP * G0)

    if m - mdot * dt <= M_DRY:
        # Return a mathematically bounded state; the optimizer's
        # inequality constraint will reject this candidate.
        dt_eff = max(
            0.0,
            (m - M_DRY) / max(mdot, 1e-15)
        )
        dt_use = min(dt, dt_eff)
    else:
        dt_use = dt

    m1 = m - mdot * dt_use

    u = Tvec / T

    # A(t) = ∫ T/m dt = Isp*g0*ln(m/m1)
    dv_thrust = ISP * G0 * np.log(m / m1)

    # Position coefficient
    K = ISP * G0 * (
        (dt_use - m / mdot)
        * np.log(m / m1)
        + dt_use
    )

    r1 = (
        r
        + v * dt_use
        + K * u
        + 0.5 * gvec() * dt_use**2
    )

    v1 = (
        v
        + dv_thrust * u
        + gvec() * dt_use
    )

    # If dt was beyond dry-mass reach, leave a clipped mass so the
    # constraint is violated rather than producing invalid logs.
    if dt_use < dt:
        m1 = M_DRY

    return r1, v1, m1, T


# ============================================================
# OPTIMIZER-T / OPTIMIZER-A
# ============================================================

def optimize_case(case, attitude_aware=False, n=N_OPT_T):
    """
    Direct-transcription nonlinear optimizer with scaled variables.

    Decision vector:
        x[0]        = normalized final time
        x[1:1+n]    = normalized Tx
        x[1+n:...]  = normalized Ty
        x[...]      = normalized Tz

    Normalization:
        tf = x_tf * t_scale
        T  = x_T * T_MAX

    Exact interval propagation is used.

    The initial guess is the analytical ESPDI solution, so the optimizer
    starts from a known feasible trajectory.
    """
    espdi_t = espdi_t_solution(case)
    if espdi_t is None:
        return None, OptimizationInfo(
            False,
            "ESPDI-T seed unavailable",
            0,
            np.nan,
            np.inf,
            np.inf,
            np.inf,
            np.nan
        )

    tf_seed = espdi_t[0]
    u_seed = espdi_t[2]

    # Use a sensible, fixed final-time scale.
    T_SCALE = T_MAX
    TF_SCALE = max(tf_seed, 10.0)

    if attitude_aware:
        a_seed = solve_espdi_a(case)
    else:
        a_seed = None

    if attitude_aware and a_seed is not None:
        tf_seed = a_seed.t[-1]
        TF_SCALE = max(tf_seed, 10.0)

        sample_t = np.linspace(
            0.0,
            tf_seed,
            n
        )

        Tmag = np.interp(
            sample_t,
            a_seed.t,
            a_seed.thrust
        )

        U = np.column_stack([
            np.interp(
                sample_t,
                a_seed.t,
                a_seed.u[:, j]
            )
            for j in range(3)
        ])

        Tvec = Tmag[:, None] * U
    else:
        Tvec = np.repeat(
            (T_MAX * u_seed)[None, :],
            n,
            axis=0
        )

    # IMPORTANT:
    # The constant-vector ESPDI-T solution is exactly representable
    # under this direct-thrust parameterization.
    z0 = np.concatenate([
        [tf_seed / TF_SCALE],
        Tvec[:, 0] / T_SCALE,
        Tvec[:, 1] / T_SCALE,
        Tvec[:, 2] / T_SCALE,
    ])

    def unpack(z):
        tf = z[0] * TF_SCALE
        tx = z[1:1+n] * T_SCALE
        ty = z[1+n:1+2*n] * T_SCALE
        tz = z[1+2*n:1+3*n] * T_SCALE
        return tf, tx, ty, tz

    def propagate(z, save=False):
        tf, tx, ty, tz = unpack(z)
        dt = tf / n

        r = case.position.copy()
        v = case.velocity.copy()
        m = case.mass

        if save:
            th = [0.0]
            rh = [r.copy()]
            vh = [v.copy()]
            mh = [m]
            Th = [0.0]
            uh = [np.array([0.0, 0.0, 1.0])]

        for k in range(n):
            Tvec_k = np.array([
                tx[k],
                ty[k],
                tz[k],
            ])

            r, v, m, Tmag = exact_interval(
                r,
                v,
                m,
                Tvec_k,
                dt
            )

            if save:
                th.append((k + 1) * dt)
                rh.append(r.copy())
                vh.append(v.copy())
                mh.append(m)
                Th.append(Tmag)

                if Tmag > 1e-12:
                    uh.append(Tvec_k / Tmag)
                else:
                    uh.append(
                        np.array([0.0, 0.0, 1.0])
                    )

        if save:
            return (
                np.asarray(th),
                np.asarray(rh),
                np.asarray(vh),
                np.asarray(mh),
                np.asarray(Th),
                np.asarray(uh),
            )

        return r, v, m

    def objective(z):
        tf, tx, ty, tz = unpack(z)

        dt = tf / n
        Tmag = np.sqrt(
            tx**2 + ty**2 + tz**2
        )

        # Propellant in kg.
        return (
            dt
            * np.sum(Tmag)
            / (ISP * G0)
        ) / M0

    def eq_constraints(z):
        rf, vf, _ = propagate(z)

        return np.hstack([
            (rf - TARGET) / POSITION_SCALE,
            vf / VELOCITY_SCALE,
        ])

    def ineq_constraints(z):
        tf, tx, ty, tz = unpack(z)

        # Dimensionless thrust-norm constraints.
        norm_sq = (
            tx**2 + ty**2 + tz**2
        ) / T_MAX**2

        thrust_margin = 1.0 - norm_sq

        _, _, mf = propagate(z)

        dry_margin = (
            mf - M_DRY
        ) / MASS_SCALE

        return np.hstack([
            thrust_margin,
            dry_margin,
        ])

    bounds = [
        (0.2, 5.0)
    ] + [
        (-1.0, 1.0)
    ] * (3 * n)

    if attitude_aware:
        # Final thrust vector must be vertical:
        # Tx_final = 0, Ty_final = 0, Tz_final >= 0.
        bounds[
            1 + n - 1
        ] = (0.0, 0.0)

        bounds[
            1 + 2*n - 1
        ] = (0.0, 0.0)

        bounds[
            1 + 3*n - 1
        ] = (0.0, 1.0)

    constraints = [
        {
            "type": "eq",
            "fun": eq_constraints
        },
        {
            "type": "ineq",
            "fun": ineq_constraints
        }
    ]

    # Exact feasible initial guess.
    guesses = [z0]

    rng = np.random.default_rng(20260818)

    # For extra starts, perturb the time and thrust vector slightly,
    # but keep the first guess untouched and feasible.
    for _ in range(
        max(0, OPTIMIZER_STARTS - 1)
    ):
        z = z0.copy()

        z[0] *= rng.uniform(
            0.99,
            1.05
        )

        z[1:] += (
            0.005
            * rng.standard_normal(
                3 * n
            )
        )

        z[1:] = np.clip(
            z[1:],
            -0.98,
            0.98
        )

        if attitude_aware:
            z[1+n-1] = 0.0
            z[1+2*n-1] = 0.0
            z[1+3*n-1] = max(
                0.01,
                z[1+3*n-1]
            )

        guesses.append(z)

    best = None
    best_result = None

    for guess in guesses:
        result = minimize(
            objective,
            guess,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={
                "ftol": 1e-9,
                "maxiter": 400,
                "disp": False,
            }
        )

        if result.success:
            if (
                best_result is None
                or result.fun < best_result.fun
            ):
                best = result.x
                best_result = result

    # If the first feasible ESPDI seed makes SLSQP fail to report success,
    # explicitly inspect whether the returned point is numerically feasible.
    if best_result is None:
        # Evaluate the exact initial feasible seed.
        eq0 = eq_constraints(z0)
        iq0 = ineq_constraints(z0)

        eq_norm = np.linalg.norm(eq0)
        min_ineq = np.min(iq0)

        return None, OptimizationInfo(
            False,
            "SLSQP failed to converge from a known feasible ESPDI seed.",
            0,
            np.nan,
            eq_norm,
            eq_norm * VELOCITY_SCALE,
            np.sqrt(max(0.0, 1.0 - min_ineq)),
            np.nan
        )

    t, r, v, mass, thrust, u = propagate(
        best,
        save=True
    )

    angles = np.array([
        angle_deg(ui) for ui in u
    ])

    final_r = r[-1]
    final_v = v[-1]

    info = OptimizationInfo(
        True,
        str(best_result.message),
        int(getattr(
            best_result,
            "nit",
            0
        )),
        float(best_result.fun * M0),
        float(np.linalg.norm(final_r)),
        float(np.linalg.norm(final_v)),
        float(np.max(thrust) / T_MAX),
        float(mass[-1])
    )

    traj = Trajectory(
        name=(
            "Optimizer-A"
            if attitude_aware
            else "Optimizer-T"
        ),
        t=t,
        r=r,
        v=v,
        mass=mass,
        thrust=thrust,
        u=u,
        angle_deg=angles,
        final_r=final_r,
        final_v=final_v,
        prop_used=case.mass - mass[-1],
        ignition_time=0.0,
        ignition_position=case.position.copy(),
        success=bool(
            np.linalg.norm(final_r) < POSITION_TOL
            and np.linalg.norm(final_v) < VELOCITY_TOL
            and mass[-1] >= M_DRY - MASS_TOL
            and np.max(thrust) <= T_MAX * (1.0 + 1e-8)
        )
    )

    return traj, info


# ============================================================
# PHYSICALLY CONSTRUCTED RANDOM CASES
# ============================================================

def make_pair_t(rng):
    """
    Pair-T state is EXACTLY an ESPDI-T ignition point.
    """
    for _ in range(2000):
        v = np.array([
            rng.uniform(-55.0, 55.0),
            rng.uniform(-55.0, 55.0),
            rng.uniform(-115.0, -50.0),
        ])

        dummy = Case(
            np.zeros(3),
            v,
            M0,
            "T"
        )

        sol = espdi_t_solution(dummy)
        if sol is None:
            continue

        tb, _, _, rI = sol

        if not (4.0 <= tb <= 40.0):
            continue

        if not (
            300.0 <= rI[2] <= 1800.0
        ):
            continue

        if (
            np.linalg.norm(rI[:2])
            > 1500.0
        ):
            continue

        return Case(
            position=rI,
            velocity=v,
            mass=M0,
            label="T"
        )

    raise RuntimeError(
        "Could not construct a valid Pair-T state."
    )


def make_pair_a(rng):
    """
    Separate Pair-A state.
    The generator only accepts states for which ESPDI-A is feasible.
    """
    for _ in range(2000):
        r = np.array([
            rng.uniform(-500.0, 500.0),
            rng.uniform(-500.0, 500.0),
            rng.uniform(700.0, 1800.0),
        ])

        v = np.array([
            rng.uniform(-45.0, 45.0),
            rng.uniform(-45.0, 45.0),
            rng.uniform(-100.0, -55.0),
        ])

        case = Case(
            position=r,
            velocity=v,
            mass=M0,
            label="A"
        )

        if solve_espdi_a(case) is not None:
            return case

    raise RuntimeError(
        "Could not construct a feasible Pair-A state."
    )


# ============================================================
# OUTPUT
# ============================================================

def print_trajectory(traj):
    if traj is None:
        print(
            "  NO CONVERGED NUMERICAL SOLUTION"
        )
        return

    print(
        f"\n{traj.name}:"
    )
    print(
        f"  burn time      = {traj.t[-1]:.9f} s"
    )
    print(
        f"  fuel used      = {traj.prop_used:.9f} kg"
    )
    print(
        "  final position = "
        + np.array2string(
            traj.final_r,
            precision=9
        )
        + " m"
    )
    print(
        "  final velocity = "
        + np.array2string(
            traj.final_v,
            precision=9
        )
        + " m/s"
    )
    print(
        f"  final |v|      = "
        f"{np.linalg.norm(traj.final_v):.9e} m/s"
    )
    print(
        f"  position error = "
        f"{np.linalg.norm(traj.final_r):.9e} m"
    )
    print(
        f"  max thrust     = "
        f"{np.max(traj.thrust):.6f} N"
    )
    print(
        f"  success        = {traj.success}"
    )


def print_optimizer_info(info):
    if info is None:
        return

    print(
        f"  optimizer status        = "
        f"{info.converged}"
    )
    print(
        f"  optimizer message      = "
        f"{info.message}"
    )
    print(
        f"  iterations             = "
        f"{info.iterations}"
    )
    print(
        f"  objective / fuel       = "
        f"{info.objective:.9f} kg"
    )
    print(
        f"  position residual      = "
        f"{info.position_residual:.9e} m"
    )
    print(
        f"  velocity residual      = "
        f"{info.velocity_residual:.9e} m/s"
    )
    print(
        f"  max thrust / Tmax      = "
        f"{info.max_thrust_ratio:.9f}"
    )
    print(
        f"  final mass             = "
        f"{info.final_mass:.9f} kg"
    )


# ============================================================
# ANIMATION
# ============================================================

def animate(case_t, espdi_t, optimizer_t, espdi_a):
    trajectories = [
        espdi_t,
        optimizer_t,
        espdi_a
    ]

    titles = [
        "ESPDI-T — Analytical Translation",
        "Optimizer-T — Independent Benchmark",
        "ESPDI-A — Analytical Extension"
    ]

    valid = [
        tr for tr in trajectories
        if tr is not None
    ]

    if not valid:
        return

    playback_end = max(
        tr.t[-1]
        for tr in valid
    )

    playback = np.linspace(
        0.0,
        playback_end,
        700
    )

    all_r = np.vstack([
        tr.r for tr in valid
    ])

    rmin = all_r.min(axis=0)
    rmax = all_r.max(axis=0)
    span = np.maximum(
        rmax - rmin,
        100.0
    )

    center = 0.5 * (rmin + rmax)
    half = 0.58 * np.max(span)

    fig = plt.figure(
        figsize=(19, 10),
        constrained_layout=True
    )

    gs = fig.add_gridspec(
        2,
        3,
        height_ratios=[3.8, 1.7]
    )

    axes = [
        fig.add_subplot(
            gs[0, i],
            projection="3d"
        )
        for i in range(3)
    ]

    info_axes = [
        fig.add_subplot(
            gs[1, i]
        )
        for i in range(3)
    ]

    for ax, title in zip(
        axes,
        titles
    ):
        ax.set_title(
            title,
            fontsize=12
        )

        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        ax.set_zlabel("Z (m)")

        ax.set_xlim(
            center[0] - half,
            center[0] + half
        )

        ax.set_ylim(
            center[1] - half,
            center[1] + half
        )

        ax.set_zlim(
            max(
                0.0,
                center[2] - half
            ),
            center[2] + half
        )

        ax.view_init(
            elev=25,
            azim=-58
        )

        ax.scatter(
            [0.0],
            [0.0],
            [0.0],
            marker="x",
            s=180,
            linewidths=3.5,
            color="red"
        )

        # Starting marker for this panel.
        tr = trajectories[
            len(axes) - len(axes)
        ]

    # Explicit initial markers.
    initial_cases = [
        case_t,
        (
            Case(
                optimizer_t.r[0],
                optimizer_t.v[0],
                optimizer_t.mass[0],
                "OT"
            )
            if optimizer_t is not None
            else None
        ),
        (
            Case(
                espdi_a.r[0],
                espdi_a.v[0],
                espdi_a.mass[0],
                "A"
            )
            if espdi_a is not None
            else None
        )
    ]

    for ax, case in zip(
        axes,
        initial_cases
    ):
        if case is not None:
            ax.scatter(
                [case.position[0]],
                [case.position[1]],
                [case.position[2]],
                marker="o",
                s=30
            )

    for ax in info_axes:
        ax.axis("off")

    artists = []
    ignition_markers = []
    texts = []

    for ax in axes:
        trail, = ax.plot(
            [],
            [],
            [],
            linewidth=2.0
        )

        rocket, = ax.plot(
            [],
            [],
            [],
            linewidth=5.0
        )

        nose, = ax.plot(
            [],
            [],
            [],
            marker="o",
            markersize=5
        )

        ignition = ax.scatter(
            [],
            [],
            [],
            marker="o",
            s=90,
            facecolors="none",
            linewidths=2.0
        )

        artists.append(
            (trail, rocket, nose)
        )

        ignition_markers.append(
            ignition
        )

    for ax in info_axes:
        texts.append(
            ax.text(
                0.01,
                0.98,
                "",
                transform=ax.transAxes,
                va="top",
                ha="left",
                family="monospace",
                fontsize=9,
                linespacing=1.35
            )
        )

    def sample_state(tr, t_now):
        local_t = min(
            t_now,
            tr.t[-1]
        )

        r = np.array([
            np.interp(
                local_t,
                tr.t,
                tr.r[:, j]
            )
            for j in range(3)
        ])

        v = np.array([
            np.interp(
                local_t,
                tr.t,
                tr.v[:, j]
            )
            for j in range(3)
        ])

        u = np.array([
            np.interp(
                local_t,
                tr.t,
                tr.u[:, j]
            )
            for j in range(3)
        ])

        mass = float(
            np.interp(
                local_t,
                tr.t,
                tr.mass
            )
        )

        thrust = float(
            np.interp(
                local_t,
                tr.t,
                tr.thrust
            )
        )

        angle = float(
            np.interp(
                local_t,
                tr.t,
                tr.angle_deg
            )
        )

        un = np.linalg.norm(u)
        if un > 1e-12:
            u = u / un

        return (
            local_t,
            r,
            v,
            u,
            mass,
            thrust,
            angle
        )

    def init():
        for trail, rocket, nose in artists:
            trail.set_data_3d(
                [],
                [],
                []
            )
            rocket.set_data_3d(
                [],
                [],
                []
            )
            nose.set_data_3d(
                [],
                [],
                []
            )

        for mark in ignition_markers:
            mark._offsets3d = (
                np.array([]),
                np.array([]),
                np.array([])
            )

        for text in texts:
            text.set_text("")

        return (
            [q for group in artists for q in group]
            + ignition_markers
            + texts
        )

    def update(frame):
        now = playback[frame]

        for i, tr in enumerate(
            trajectories
        ):
            if tr is None:
                continue

            (
                local_t,
                r,
                v,
                u,
                mass,
                thrust,
                angle
            ) = sample_state(
                tr,
                now
            )

            mask = (
                tr.t
                <= local_t + 1e-12
            )

            history = tr.r[mask]

            trail, rocket, nose = artists[i]

            trail.set_data_3d(
                history[:, 0],
                history[:, 1],
                history[:, 2]
            )

            length = max(
                30.0,
                0.035 * np.max(span)
            )

            p0 = r - 0.5 * length * u
            p1 = r + 0.5 * length * u

            rocket.set_data_3d(
                [p0[0], p1[0]],
                [p0[1], p1[1]],
                [p0[2], p1[2]]
            )

            nose.set_data_3d(
                [p1[0]],
                [p1[1]],
                [p1[2]]
            )

            if (
                tr.ignition_position is not None
            ):
                p = tr.ignition_position
                ignition_markers[i]._offsets3d = (
                    np.array([p[0]]),
                    np.array([p[1]]),
                    np.array([p[2]])
                )

            speed = np.linalg.norm(v)
            position_error = np.linalg.norm(r)
            fuel = (
                tr.mass[0] - mass
            )

            if (
                local_t >= tr.t[-1] - 1e-8
                and position_error < 0.1
                and speed < 0.05
            ):
                status = "LANDED"
            elif local_t < tr.t[-1]:
                status = "IN FLIGHT"
            else:
                status = "TERMINAL"

            texts[i].set_text(
                f"{tr.name}\n"
                f"time        {local_t:9.3f} s\n"
                f"x           {r[0]:10.3f} m\n"
                f"y           {r[1]:10.3f} m\n"
                f"z           {r[2]:10.3f} m\n"
                f"vx          {v[0]:10.3f} m/s\n"
                f"vy          {v[1]:10.3f} m/s\n"
                f"vz          {v[2]:10.3f} m/s\n"
                f"|v|         {speed:10.5f} m/s\n"
                f"angle       {angle:10.3f} deg\n"
                f"thrust      {thrust:10.1f} N\n"
                f"fuel used   {fuel:10.3f} kg\n"
                f"pos. error  {position_error:10.5f} m\n"
                f"status      {status}"
            )

        return (
            [q for group in artists for q in group]
            + ignition_markers
            + texts
        )

    animation = FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(playback),
        interval=25,
        blit=False,
        repeat=True
    )

    # Keep a strong reference.
    fig._espdi_animation = animation

    fig.suptitle(
        "ESPDI 3D Landing — Method-Specific Trajectories",
        fontsize=14
    )

    plt.show()


# ============================================================
# MONTE CARLO
# ============================================================

def monte_carlo(n_cases, seed):
    rng = np.random.default_rng(seed)

    rows = []

    for i in range(n_cases):
        print(
            f"\nMonte Carlo case "
            f"{i + 1}/{n_cases}"
        )

        try:
            # Pair T.
            case_t = make_pair_t(rng)

            t0 = time.perf_counter()
            et = propagate_espdi_t(case_t)
            rt_et = time.perf_counter() - t0

            t0 = time.perf_counter()
            ot, oi_t = optimize_case(
                case_t,
                attitude_aware=False,
                n=N_OPT_T
            )
            rt_ot = time.perf_counter() - t0

            # Pair A.
            case_a = make_pair_a(rng)

            t0 = time.perf_counter()
            ea = solve_espdi_a(case_a)
            rt_ea = time.perf_counter() - t0

            t0 = time.perf_counter()
            oa, oi_a = optimize_case(
                case_a,
                attitude_aware=True,
                n=N_OPT_A
            )
            rt_oa = time.perf_counter() - t0

            for case, tr, rt in [
                (case_t, et, rt_et),
                (case_t, ot, rt_ot),
                (case_a, ea, rt_ea),
                (case_a, oa, rt_oa)
            ]:
                if tr is None:
                    continue

                rows.append({
                    "case_id": i + 1,
                    "pair": case.label,
                    "method": tr.name,
                    "initial_x_m": case.position[0],
                    "initial_y_m": case.position[1],
                    "initial_z_m": case.position[2],
                    "initial_vx_mps": case.velocity[0],
                    "initial_vy_mps": case.velocity[1],
                    "initial_vz_mps": case.velocity[2],
                    "burn_time_s": tr.t[-1],
                    "propellant_kg": tr.prop_used,
                    "final_x_m": tr.final_r[0],
                    "final_y_m": tr.final_r[1],
                    "final_z_m": tr.final_r[2],
                    "final_vx_mps": tr.final_v[0],
                    "final_vy_mps": tr.final_v[1],
                    "final_vz_mps": tr.final_v[2],
                    "final_speed_mps": np.linalg.norm(tr.final_v),
                    "position_error_m": np.linalg.norm(tr.final_r),
                    "max_thrust_N": np.max(tr.thrust),
                    "runtime_s": rt,
                    "success": tr.success
                })

        except Exception as exc:
            print(
                "CASE FAILED:",
                exc
            )

    if rows and SAVE_RESULTS:
        path = Path(
            "espdi_monte_carlo_results.csv"
        )
        with path.open(
            "w",
            newline="",
            encoding="utf-8"
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=list(rows[0].keys())
            )
            writer.writeheader()
            writer.writerows(rows)

        print(
            "\nSaved Monte Carlo results:",
            path.resolve()
        )

    return rows


# ============================================================
# MAIN
# ============================================================

def main():
    rng = np.random.default_rng(
        RANDOM_SEED
    )

    print("=" * 90)
    print("ESPDI RESEARCH SIMULATION")
    print("=" * 90)
    print(
        f"m0={M0:.1f} kg | "
        f"dry={M_DRY:.1f} kg | "
        f"Tmax={T_MAX/1000:.1f} kN | "
        f"Isp={ISP:.1f} s | "
        f"T/W={INITIAL_TW:.4f}"
    )

    # --------------------------------------------------------
    # Pair T
    # --------------------------------------------------------
    case_t = make_pair_t(rng)

    print("\nPAIR T — SAME INITIAL STATE")
    print("r0 =", case_t.position)
    print("v0 =", case_t.velocity)

    t0 = time.perf_counter()
    espdi_t = propagate_espdi_t(case_t)
    rt_et = time.perf_counter() - t0

    print("\nESPDI-T:")
    print_trajectory(espdi_t)
    print(f"runtime = {rt_et:.6f} s")

    t0 = time.perf_counter()
    optimizer_t, opt_info_t = optimize_case(
        case_t,
        attitude_aware=False,
        n=N_OPT_T
    )
    rt_ot = time.perf_counter() - t0

    print("\nOPTIMIZER-T:")
    print_trajectory(optimizer_t)
    print_optimizer_info(opt_info_t)
    print(f"runtime = {rt_ot:.6f} s")

    if (
        espdi_t is not None
        and optimizer_t is not None
        and optimizer_t.success
    ):
        gap_t = 100.0 * (
            espdi_t.prop_used
            - optimizer_t.prop_used
        ) / optimizer_t.prop_used

        print(
            f"\nESPDI-T fuel gap vs "
            f"Optimizer-T = {gap_t:.6f}%"
        )

    # --------------------------------------------------------
    # Pair A
    # --------------------------------------------------------
    case_a = make_pair_a(rng)

    print("\nPAIR A — SAME INITIAL STATE")
    print("r0 =", case_a.position)
    print("v0 =", case_a.velocity)

    t0 = time.perf_counter()
    espdi_a = solve_espdi_a(case_a)
    rt_ea = time.perf_counter() - t0

    print("\nESPDI-A:")
    print_trajectory(espdi_a)
    print(f"runtime = {rt_ea:.6f} s")

    t0 = time.perf_counter()
    optimizer_a, opt_info_a = optimize_case(
        case_a,
        attitude_aware=True,
        n=N_OPT_A
    )
    rt_oa = time.perf_counter() - t0

    print("\nOPTIMIZER-A:")
    print_trajectory(optimizer_a)
    print_optimizer_info(opt_info_a)
    print(f"runtime = {rt_oa:.6f} s")

    if (
        espdi_a is not None
        and optimizer_a is not None
        and optimizer_a.success
    ):
        gap_a = 100.0 * (
            espdi_a.prop_used
            - optimizer_a.prop_used
        ) / optimizer_a.prop_used

        print(
            f"\nESPDI-A fuel gap vs "
            f"Optimizer-A = {gap_a:.6f}%"
        )

    # --------------------------------------------------------
    # Save single case
    # --------------------------------------------------------
    if SAVE_RESULTS:
        rows = []

        for case, tr, rt in [
            (case_t, espdi_t, rt_et),
            (case_t, optimizer_t, rt_ot),
            (case_a, espdi_a, rt_ea),
            (case_a, optimizer_a, rt_oa)
        ]:
            if tr is None:
                continue

            rows.append({
                "pair": case.label,
                "method": tr.name,
                "initial_x_m": case.position[0],
                "initial_y_m": case.position[1],
                "initial_z_m": case.position[2],
                "initial_vx_mps": case.velocity[0],
                "initial_vy_mps": case.velocity[1],
                "initial_vz_mps": case.velocity[2],
                "burn_time_s": tr.t[-1],
                "propellant_kg": tr.prop_used,
                "final_x_m": tr.final_r[0],
                "final_y_m": tr.final_r[1],
                "final_z_m": tr.final_r[2],
                "final_vx_mps": tr.final_v[0],
                "final_vy_mps": tr.final_v[1],
                "final_vz_mps": tr.final_v[2],
                "final_speed_mps": np.linalg.norm(tr.final_v),
                "position_error_m": np.linalg.norm(tr.final_r),
                "max_thrust_N": np.max(tr.thrust),
                "runtime_s": rt,
                "success": tr.success
            })

        out = Path(
            "espdi_single_case_results.csv"
        )

        with out.open(
            "w",
            newline="",
            encoding="utf-8"
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=list(rows[0].keys())
            )
            writer.writeheader()
            writer.writerows(rows)

        print(
            "\nSaved single-case results:",
            out.resolve()
        )

    # --------------------------------------------------------
    # Proper animated landing comparison
    # --------------------------------------------------------
    animate(
        case_t,
        espdi_t,
        optimizer_t,
        espdi_a
    )

    # --------------------------------------------------------
    # Monte Carlo
    # --------------------------------------------------------
    if RUN_MONTE_CARLO:
        monte_carlo(
            MONTE_CARLO_CASES,
            RANDOM_SEED
        )


if __name__ == "__main__":
    main()