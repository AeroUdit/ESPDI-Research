"""
INDEPENDENT OPTIMIZER-T VALIDATION v3
=====================================

This script is intentionally standalone: it does NOT import the main
simulator. It reproduces the main simulator's Pair-T random-state generator
so that the same RANDOM_SEED produces the same Pair-T initial condition.

It performs two separate checks:

1. SAME-STATE / ESPDI-SEED CONTROL
   The optimizer is initialized with the ESPDI-T solution. This verifies that
   the benchmark optimizer can reproduce the known analytical solution.

2. INDEPENDENT REDISCOVERY
   The optimizer is initialized with deliberately non-ESPDI thrust histories.
   The best feasible solution is compared with ESPDI-T.

It also checks mesh sensitivity.

IMPORTANT:
    Same RANDOM_SEED + same generator = same Pair-T state as the main
    simulator, because the generator below is copied from the main simulator's
    Pair-T generator.

Required packages:
    numpy
    scipy
"""

from dataclasses import dataclass
from pathlib import Path
import csv
import time

import numpy as np
from scipy.optimize import brentq, minimize


# ============================================================
# SETTINGS
# ============================================================

RANDOM_SEED = 20260818

# Smoke-test settings. Change ONLY after the smoke test works.
N_STATES = 20
INDEPENDENT_STARTS = 2
MESHES = (8, 12)
FINAL_MESH = None

# Final validation recommendation after smoke test:
# N_STATES = 10
# INDEPENDENT_STARTS = 8
# MESHES = (8, 12, 16)
# FINAL_MESH = 20

SAVE_RESULTS = True

POSITION_TOL = 0.10       # m
VELOCITY_TOL = 0.05       # m/s

POSITION_SCALE = 1000.0
VELOCITY_SCALE = 100.0
MASS_SCALE = 1000.0

G = 9.81
G0 = 9.80665
M0 = 12_000.0
M_DRY = 4_500.0
T_MAX = 180_000.0
ISP = 320.0
MDOT_MAX = T_MAX / (ISP * G0)
TARGET = np.zeros(3)


@dataclass
class Case:
    r0: np.ndarray
    v0: np.ndarray
    m0: float
    label: str


@dataclass
class OptResult:
    success: bool
    feasible: bool
    fuel: float
    tf: float
    position_error: float
    velocity_error: float
    final_mass: float
    max_thrust_ratio: float
    iterations: int
    runtime: float
    message: str
    start_type: str
    mesh: int
    z: np.ndarray | None


# ============================================================
# COMMON PHYSICS
# ============================================================

def gvec():
    return np.array([0.0, 0.0, -G])


def mass_after(t, m):
    return m - MDOT_MAX * t


def t_prop_max(m):
    return max(0.0, (m - M_DRY) / MDOT_MAX)


def required_dv(t_b, v0):
    return -v0 - gvec() * t_b


def thrust_dv(t_b, m0):
    mf = mass_after(t_b, m0)
    if mf <= 0.0:
        return np.nan
    return ISP * G0 * np.log(m0 / mf)


def solve_espdi_t(v0, m0):
    tmax = t_prop_max(m0)
    if tmax <= 1e-8:
        return None

    def f(t):
        mf = mass_after(t, m0)
        if mf < M_DRY:
            return np.nan
        return thrust_dv(t, m0) - np.linalg.norm(required_dv(t, v0))

    grid = np.linspace(1e-6, tmax * (1.0 - 1e-10), 1600)
    vals = np.array([f(t) for t in grid])

    for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if not (np.isfinite(fa) and np.isfinite(fb)):
            continue
        if fa == 0.0:
            tb = float(a)
            break
        if fa * fb < 0.0:
            tb = float(brentq(f, a, b, xtol=1e-11, rtol=1e-11))
            break
    else:
        return None

    D = required_dv(tb, v0)
    u = D / np.linalg.norm(D)
    return tb, u


def espdi_ignition_position(v0, m0):
    sol = solve_espdi_t(v0, m0)
    if sol is None:
        return None

    tb, u = sol
    mf = mass_after(tb, m0)
    mdot = MDOT_MAX

    K = ISP * G0 * (
        (tb - m0 / mdot) * np.log(m0 / mf) + tb
    )

    r_i = (
        TARGET
        - v0 * tb
        - K * u
        - 0.5 * gvec() * tb**2
    )

    return tb, u, r_i


# ============================================================
# IMPORTANT: THIS IS THE MAIN SIMULATOR'S PAIR-T GENERATOR
# ============================================================

def make_pair_t(rng, label="T"):
    """
    Verbatim mathematical structure used by the main simulator's Pair-T
    generator: sample velocity, solve ESPDI-T, and accept the analytical
    ignition point if it is inside the demonstration envelope.
    """
    for _ in range(2000):
        v = np.array([
            rng.uniform(-55.0, 55.0),
            rng.uniform(-55.0, 55.0),
            rng.uniform(-115.0, -50.0),
        ])

        solved = espdi_ignition_position(v, M0)
        if solved is None:
            continue

        tb, _, r_i = solved

        if not (4.0 <= tb <= 40.0):
            continue

        if not (300.0 <= r_i[2] <= 1800.0):
            continue

        if np.linalg.norm(r_i[:2]) > 1500.0:
            continue

        return Case(
            r0=r_i,
            v0=v,
            m0=M0,
            label=label,
        )

    raise RuntimeError("Could not construct a valid Pair-T state.")


# ============================================================
# EXACT CONSTANT-THRUST INTERVAL PROPAGATION
# ============================================================

def exact_interval(r, v, m, Tvec, dt):
    T = float(np.linalg.norm(Tvec))

    if T < 1e-14:
        r1 = r + v * dt + 0.5 * gvec() * dt**2
        v1 = v + gvec() * dt
        return r1, v1, m, T

    mdot = T / (ISP * G0)
    m1 = m - mdot * dt

    if m1 <= M_DRY:
        return (
            np.full(3, 1e6),
            np.full(3, 1e4),
            M_DRY - 1000.0,
            T,
        )

    u = Tvec / T
    L = np.log(m / m1)
    dv_thrust = ISP * G0 * L
    K = ISP * G0 * ((dt - m / mdot) * L + dt)

    r1 = r + v * dt + K * u + 0.5 * gvec() * dt**2
    v1 = v + dv_thrust * u + gvec() * dt

    return r1, v1, m1, T


def propagate(case, z, n):
    seed = solve_espdi_t(case.v0, case.m0)
    if seed is None:
        raise RuntimeError("ESPDI seed unavailable.")

    tf_scale = max(10.0, seed[0])
    tf = z[0] * tf_scale

    tx = z[1:1+n] * T_MAX
    ty = z[1+n:1+2*n] * T_MAX
    tz = z[1+2*n:1+3*n] * T_MAX

    dt = tf / n
    r = case.r0.copy()
    v = case.v0.copy()
    m = case.m0

    for k in range(n):
        Tvec = np.array([tx[k], ty[k], tz[k]])
        r, v, m, _ = exact_interval(r, v, m, Tvec, dt)

    return r, v, m


# ============================================================
# INITIAL GUESSES
# ============================================================

def seed_espdi(case, n):
    sol = solve_espdi_t(case.v0, case.m0)
    if sol is None:
        return None

    tb, u = sol
    z = np.concatenate([
        [tb / max(10.0, tb)],
        np.full(n, u[0]),
        np.full(n, u[1]),
        np.full(n, u[2]),
    ])
    return z


def seed_vertical(case, n):
    sol = solve_espdi_t(case.v0, case.m0)
    if sol is None:
        return None

    tb, _ = sol
    return np.concatenate([
        [tb / max(10.0, tb)],
        np.zeros(n),
        np.zeros(n),
        np.full(n, 0.95),
    ])


def seed_horizontal(case, n):
    sol = solve_espdi_t(case.v0, case.m0)
    if sol is None:
        return None

    tb, _ = sol
    return np.concatenate([
        [1.15 * tb / max(10.0, tb)],
        np.full(n, 0.55),
        np.zeros(n),
        np.full(n, 0.30),
    ])


def seed_wrong_fixed(case, n):
    sol = solve_espdi_t(case.v0, case.m0)
    if sol is None:
        return None

    tb, _ = sol
    # Constant thrust direction deliberately different from ESPDI.
    u = np.array([0.60, 0.20, 0.775])
    u = u / np.linalg.norm(u)

    return np.concatenate([
        [1.10 * tb / max(10.0, tb)],
        np.full(n, u[0]),
        np.full(n, u[1]),
        np.full(n, u[2]),
    ])


def seed_random(case, n, rng):
    z = seed_espdi(case, n)
    if z is None:
        return None

    # Deliberately erase the ESPDI thrust-direction information.
    z[1:1+n] = np.clip(
        rng.uniform(-0.75, 0.75, n),
        -0.95,
        0.95,
    )
    z[1+n:1+2*n] = np.clip(
        rng.uniform(-0.75, 0.75, n),
        -0.95,
        0.95,
    )
    z[1+2*n:1+3*n] = np.clip(
        rng.uniform(0.15, 0.95, n),
        -0.95,
        0.95,
    )

    return z


def make_non_espdi_starts(case, n, count, rng):
    generators = [
        seed_vertical,
        seed_horizontal,
        seed_wrong_fixed,
    ]

    starts = []

    # At least the requested count; cycle through deterministic seeds,
    # then use random seeds.
    for i in range(count):
        if i < len(generators):
            z = generators[i](case, n)
        else:
            z = seed_random(case, n, rng)

        if z is not None:
            starts.append(z)

    return starts


# ============================================================
# OPTIMIZER
# ============================================================

def optimize_from_seed(case, z0, n, start_type):
    seed = solve_espdi_t(case.v0, case.m0)
    if seed is None:
        return OptResult(
            False, False, np.nan, np.nan,
            np.inf, np.inf, np.nan, np.inf,
            0, 0.0, "ESPDI seed unavailable",
            start_type, n, z0,
        )

    tf_scale = max(10.0, seed[0])

    def unpack(z):
        tf = z[0] * tf_scale
        tx = z[1:1+n] * T_MAX
        ty = z[1+n:1+2*n] * T_MAX
        tz = z[1+2*n:1+3*n] * T_MAX
        return tf, tx, ty, tz

    def objective(z):
        tf, tx, ty, tz = unpack(z)
        dt = tf / n
        mag = np.sqrt(tx**2 + ty**2 + tz**2)
        return (
            dt * np.sum(mag) / (ISP * G0)
        ) / M0

    def equality(z):
        try:
            rf, vf, _ = propagate(case, z, n)
        except Exception:
            rf = np.full(3, 1e6)
            vf = np.full(3, 1e4)
        return np.hstack([
            (rf - TARGET) / POSITION_SCALE,
            vf / VELOCITY_SCALE,
        ])

    def inequality(z):
        tf, tx, ty, tz = unpack(z)
        ratio_sq = (
            tx**2 + ty**2 + tz**2
        ) / T_MAX**2
        try:
            _, _, mf = propagate(case, z, n)
        except Exception:
            mf = M_DRY - 1000.0
        return np.hstack([
            1.0 - ratio_sq,
            (mf - M_DRY) / MASS_SCALE,
        ])

    bounds = (
        [(0.20, 5.0)]
        + [(-1.0, 1.0)] * (3*n)
    )

    constraints = [
        {"type": "eq", "fun": equality},
        {"type": "ineq", "fun": inequality},
    ]

    t0 = time.perf_counter()

    try:
        result = minimize(
            objective,
            z0,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={
                "ftol": 1e-8,
                "maxiter": 500,
                "disp": False,
            },
        )
    except Exception as exc:
        runtime = time.perf_counter() - t0
        return OptResult(
            False, False, np.nan, np.nan,
            np.inf, np.inf, np.nan, np.inf,
            0, runtime, str(exc), start_type, n, z0,
        )

    runtime = time.perf_counter() - t0

    z = result.x
    tf, tx, ty, tz = unpack(z)

    rf, vf, mf = propagate(
        case,
        z,
        n,
    )

    mag = np.sqrt(
        tx**2 + ty**2 + tz**2
    )

    position_error = float(
        np.linalg.norm(rf - TARGET)
    )
    velocity_error = float(
        np.linalg.norm(vf)
    )

    max_ratio = float(
        np.max(mag) / T_MAX
    )

    feasible = bool(
        position_error <= POSITION_TOL
        and velocity_error <= VELOCITY_TOL
        and mf >= M_DRY
        and max_ratio <= 1.0 + 1e-7
    )

    fuel = float(
        objective(z) * M0
    )

    return OptResult(
        bool(result.success),
        feasible,
        fuel,
        float(tf),
        position_error,
        velocity_error,
        float(mf),
        max_ratio,
        int(getattr(result, "nit", 0)),
        runtime,
        str(result.message),
        start_type,
        n,
        z,
    )


# ============================================================
# STATE MATCHING CHECK
# ============================================================

def print_reference_state(case, state_index):
    print(
        f"Reference Pair-T state {state_index}:"
    )
    print(
        "  r0 =",
        np.array2string(case.r0, precision=12)
    )
    print(
        "  v0 =",
        np.array2string(case.v0, precision=12)
    )
    print(
        "  m0 =",
        f"{case.m0:.12f} kg"
    )


# ============================================================
# ONE STATE VALIDATION
# ============================================================

def validate_case(case, case_index, rng):
    print("\n" + "=" * 78)
    print(
        f"CASE {case_index}"
    )
    print_reference_state(case, case_index)

    espdi = solve_espdi_t(case.v0, case.m0)
    if espdi is None:
        raise RuntimeError("ESPDI-T failed for generated state.")

    tb, u = espdi
    espdi_fuel = MDOT_MAX * tb

    print(
        f"ESPDI burn = {tb:.9f} s"
    )
    print(
        f"ESPDI fuel = {espdi_fuel:.9f} kg"
    )

    all_results = []

    for n in MESHES:
        print(
            f"\n  MESH N={n}"
        )

        # A. Same-seed control.
        z_espdi = seed_espdi(case, n)
        print("    same-seed control...")
        control = optimize_from_seed(
            case,
            z_espdi,
            n,
            "espdi"
        )
        all_results.append(control)
        print(
            f"      success={control.success} | "
            f"feasible={control.feasible} | "
            f"fuel={control.fuel:.9f} | "
            f"poserr={control.position_error:.3e} | "
            f"verr={control.velocity_error:.3e} | "
            f"time={control.runtime:.2f}s"
        )

        # B. Independent/non-ESPDI starts.
        starts = make_non_espdi_starts(
            case,
            n,
            INDEPENDENT_STARTS,
            rng,
        )

        for j, z in enumerate(starts, start=1):
            name = (
                "independent_"
                + str(j)
            )
            print(
                f"    {name}...",
                flush=True
            )
            result = optimize_from_seed(
                case,
                z,
                n,
                name
            )
            all_results.append(result)
            print(
                f"      success={result.success} | "
                f"feasible={result.feasible} | "
                f"fuel={result.fuel:.9f} | "
                f"poserr={result.position_error:.3e} | "
                f"verr={result.velocity_error:.3e} | "
                f"time={result.runtime:.2f}s"
            )

    # Summaries for this case.
    feasible_independent = [
        r for r in all_results
        if r.start_type != "espdi"
        and r.feasible
    ]

    best_independent = (
        min(
            feasible_independent,
            key=lambda r: r.fuel
        )
        if feasible_independent
        else None
    )

    if best_independent is None:
        print(
            "\n  INDEPENDENT REDISCOVERY: FAILED"
        )
    else:
        best_gap = 100.0 * (
            best_independent.fuel - espdi_fuel
        ) / espdi_fuel
        print(
            "\n  BEST INDEPENDENT FUEL = "
            f"{best_independent.fuel:.9f} kg"
        )
        print(
            "  BEST INDEPENDENT GAP  = "
            f"{best_gap:+.6f}%"
        )
        print(
            "  INDEPENDENT REDISCOVERY = "
            f"{abs(best_gap) <= 0.05}"
        )

    return espdi_fuel, all_results


# ============================================================
# MAIN
# ============================================================

def main():
    rng = np.random.default_rng(
        RANDOM_SEED
    )

    print("=" * 78)
    print(
        "INDEPENDENT OPTIMIZER-T VALIDATION v3"
    )
    print("=" * 78)
    print(
        "Standalone validation: no import of the main simulator."
    )
    print(
        "Pair-T generator matches the main simulator's Pair-T generator."
    )
    print(
        f"Seed={RANDOM_SEED} | "
        f"Cases={N_STATES} | "
        f"Meshes={MESHES} | "
        f"Independent starts/mesh={INDEPENDENT_STARTS}"
    )

    rows = []

    for i in range(N_STATES):
        case = make_pair_t(
            rng,
            label="T"
        )

        espdi_fuel, results = validate_case(
            case,
            i + 1,
            rng,
        )

        for result in results:
            rows.append({
                "case": i + 1,
                "start_type": result.start_type,
                "mesh": result.mesh,
                "success": result.success,
                "feasible": result.feasible,
                "espdi_fuel_kg": espdi_fuel,
                "optimizer_fuel_kg": result.fuel,
                "fuel_gap_percent": (
                    100.0 * (result.fuel - espdi_fuel)
                    / espdi_fuel
                    if np.isfinite(result.fuel)
                    else np.nan
                ),
                "tf_s": result.tf,
                "position_error_m": result.position_error,
                "velocity_error_mps": result.velocity_error,
                "final_mass_kg": result.final_mass,
                "max_thrust_ratio": result.max_thrust_ratio,
                "iterations": result.iterations,
                "runtime_s": result.runtime,
                "message": result.message,
            })

    # --------------------------------------------------------
    # Write raw results.
    # --------------------------------------------------------
    out_raw = Path(
        "optimizer_t_validation_v3_results.csv"
    )

    if SAVE_RESULTS and rows:
        with out_raw.open(
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

    # --------------------------------------------------------
    # Summary.
    # --------------------------------------------------------
    print("\n" + "=" * 78)
    print("VALIDATION SUMMARY")
    print("=" * 78)

    # Group by mesh and start type.
    for mesh in MESHES:
        for start in [
            "espdi",
            *[
                "independent_" + str(i)
                for i in range(1, INDEPENDENT_STARTS + 1)
            ]
        ]:
            subset = [
                r for r in rows
                if r["mesh"] == mesh
                and r["start_type"] == start
                and r["feasible"]
            ]

            if not subset:
                continue

            gaps = np.array([
                r["fuel_gap_percent"]
                for r in subset
            ])

            print(
                f"N={mesh:2d} | "
                f"start={start:16s} | "
                f"feasible={len(subset)}/{N_STATES} | "
                f"mean gap={np.mean(gaps):+.6e}% | "
                f"max |gap|={np.max(np.abs(gaps)):.6e}%"
            )

    # Global best independent result by case.
    print("\nBEST INDEPENDENT REDISCOVERY BY CASE")
    for case_id in range(1, N_STATES + 1):
        subset = [
            r for r in rows
            if r["case"] == case_id
            and r["start_type"] != "espdi"
            and r["feasible"]
        ]
        if not subset:
            print(
                f"Case {case_id}: NO FEASIBLE INDEPENDENT SOLUTION"
            )
            continue

        best = min(
            subset,
            key=lambda r: r["optimizer_fuel_kg"]
        )

        print(
            f"Case {case_id}: "
            f"best gap={best['fuel_gap_percent']:+.6f}% | "
            f"mesh={best['mesh']} | "
            f"start={best['start_type']}"
        )

    if SAVE_RESULTS:
        print(
            "\nSaved raw results to:",
            out_raw.resolve()
        )


if __name__ == "__main__":
    main()