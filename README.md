## ESPDI-T

ESPDI-T solves the drag-free, constant-gravity, variable-mass point-mass landing problem.

The core calculation is:

### 1. Variable-mass model

The maximum-thrust propellant flow rate is

$$
\dot{m}_p = \frac{T_{\max}}{I_{sp}g_0}
$$

and the rocket mass during the burn is

$$
m(t)=m_0-\dot{m}_p t.
$$

The burn must satisfy the dry-mass constraint

$$
m_f \geq m_{\mathrm{dry}}.
$$

### 2. Required velocity change

To reach the landing point with zero terminal velocity,

$$
\mathbf{D}(t_b)
=
-\mathbf{v}_0-\mathbf{g}t_b
$$

where:

- $\mathbf{v}_0$ is the initial velocity
- $\mathbf{g}$ is the gravitational acceleration vector
- $t_b$ is the burn time

### 3. Variable-mass delta-v

The thrust-generated velocity change is

$$
\Delta v_T
=
I_{sp}g_0
\ln\left(\frac{m_0}{m_f}\right).
$$

Therefore the burn time is determined by the scalar equation

$$
I_{sp}g_0
\ln\left(
\frac{m_0}
{m_0-\dot{m}_p t_b}
\right)
=
\left\|
-\mathbf{v}_0-\mathbf{g}t_b
\right\|.
$$

No trajectory-level optimization is required.

### 4. Thrust direction

Once $t_b$ is known, the constant inertial thrust direction is

$$
\hat{\mathbf{u}}_{\mathrm{ESPDI}}
=
\frac{
-\mathbf{v}_0-\mathbf{g}t_b
}{
\left\|
-\mathbf{v}_0-\mathbf{g}t_b
\right\|
}.
$$

### 5. Ignition point

The variable-mass position integral gives

$$
K
=
I_{sp}g_0
\left[
\left(
t_b-\frac{m_0}{\dot{m}_p}
\right)
\ln\left(\frac{m_0}{m_f}\right)
+
t_b
\right].
$$

The ESPDI-T ignition point is then

$$
\mathbf{r}_I
=
\mathbf{r}_L
-
\mathbf{v}_0t_b
-
K\hat{\mathbf{u}}_{\mathrm{ESPDI}}
-
\frac{1}{2}\mathbf{g}t_b^2.
$$
