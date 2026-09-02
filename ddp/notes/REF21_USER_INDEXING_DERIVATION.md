# Ref. [21] in the user's period indexing

Source: S. El Kazdadi, J. Carpentier, and J. Ponce, "Equality Constrained
Differential Dynamic Programming," ICRA 2021, Ref. [21] of the FilterDDP
paper.  Official record: <https://inria.hal.science/hal-03184203v2>.

This note translates and derives the complete mathematical development of
Ref. [21].  It is deliberately not a blind symbol replacement: all derivative
orientations and index shifts are checked under the convention below.

## 1. Approved convention and exact index map

The user's periods are numbered `1,2,...,T`:

```text
x_0 --u_1--> x_1 --u_2--> ... --u_t--> x_t --...-- u_T--> x_T
```

- `x_{t-1}` is the state at the **start** of period `t`.
- `u_t` is the complete vector of period-`t` decision variables.
- `x_t` is the state at the **end** of period `t`.
- `t` is a period label, not a state-node label.
- `tau` is used as a summation index.

Ref. [21] uses nodes `i=0,...,N` and controls `u_i` on the transition from
`x_i` to `x_{i+1}`.  The exact map is

```text
paper i       <-> user t-1
paper i+1     <-> user t
paper N       <-> user T
paper u_i     <-> user u_t
paper f_i     <-> user f_t
paper h_i     <-> user h_t
paper V_i(x_i)<-> user V_t(x_{t-1})
paper V_N(x_N)<-> user V_{T+1}(x_T)
```

The last line is the important shift: the terminal value function is
`V_{T+1}`, not `V_T`.

All gradients below are **column vectors**.  Ref. [21] writes some gradients
as rows; the formulas here include the corresponding transposes explicitly.
For a vector map `f=(f_1,...,f_{n_x})`, the contraction

```text
V_x^+ . f_yy := sum_a (V_x^+)_a * nabla_yy^2 f_a
```

is a matrix.  A superscript `+` on a value derivative means "the derivative
of the next value function, `V_{t+1}`, evaluated at `x_t=f_t(...)`."

## 2. Equations (1)-(6): dynamics, objective, and Bellman recursion

### Ref. [21] Eq. (1): dynamics

For every period `t=1,...,T`,

```math
x_t=f_t(x_{t-1},u_t),
```

with fixed initial state `x_0`.

### Eqs. (2)-(3): finite-horizon objective

For a complete decision sequence `u_{1:T}`,

```math
J_1(x_0,u_{1:T})
=\sum_{\tau=1}^{T}\ell_\tau(x_{\tau-1},u_\tau)+\ell_F(x_T),
```

and the trajectory-optimization problem is

```math
u_{1:T}^\star=\arg\min_{u_{1:T}}J_1(x_0,u_{1:T})
```

subject to the dynamics above.  More generally, the cost remaining at the
start of period `t` is

```math
J_t(x_{t-1},u_{t:T})
=\sum_{\tau=t}^{T}\ell_\tau(x_{\tau-1},u_\tau)+\ell_F(x_T).
```

The value function is a **function**, not one already-known scalar value:

```math
V_t(x_{t-1})=\min_{u_{t:T}}J_t(x_{t-1},u_{t:T}).
```

### Eq. (4): Bellman recursion

The terminal boundary is

```math
V_{T+1}(x_T)=\ell_F(x_T).
```

For `t=T,T-1,...,1`,

```math
V_t(x_{t-1})
=\min_{u_t}\left[
\ell_t(x_{t-1},u_t)
+V_{t+1}\!\left(f_t(x_{t-1},u_t)\right)
\right].
```

For example, one step before the end,

```math
V_T(x_{T-1})
=\min_{u_T}\left[
\ell_T(x_{T-1},u_T)+\ell_F(f_T(x_{T-1},u_T))
\right],
```

and recursively,

```math
V_{T-1}(x_{T-2})
=\min_{u_{T-1}}\left[
\ell_{T-1}(x_{T-2},u_{T-1})
+V_T(f_{T-1}(x_{T-2},u_{T-1}))
\right].
```

Thus knowledge of `V_{t+1}` means knowledge of the function to insert into the
period-`t` minimization.  It does not mean that its numerical value is fixed
before `x_t` is known.

### Eqs. (5)-(6): one-period function to minimize

To avoid confusing Ref. [21]'s `Q` with OPF reactive power, this note denotes
its local Bellman function by `mathcal Q`:

```math
\mathcal Q_t(x_{t-1},u_t)
:=\ell_t(x_{t-1},u_t)
+V_{t+1}(f_t(x_{t-1},u_t)).
```

The period-`t` policy is therefore

```math
u_t^\star(x_{t-1})=\arg\min_{u_t}\mathcal Q_t(x_{t-1},u_t),
```

and `V_t(x_{t-1})=mathcal Q_t(x_{t-1},u_t^star(x_{t-1}))`.

## 3. Equations (7)-(12): ordinary unconstrained DDP

Let `(bar x_{t-1},bar u_t)` be the current nominal trajectory and define

```math
\delta x_{t-1}=x_{t-1}-\bar x_{t-1},\qquad
\delta u_t=u_t-\bar u_t.
```

### Eq. (7): quadratic value model

Around `bar x_{t-1}`,

```math
V_t(\bar x_{t-1}+\delta x_{t-1})
\approx v_t+V_{x,t}^{\mathsf T}\delta x_{t-1}
+\tfrac12\delta x_{t-1}^{\mathsf T}V_{xx,t}\delta x_{t-1}.
```

Here `v_t=V_t(bar x_{t-1})`, `V_{x,t}=nabla V_t(bar x_{t-1})`, and
`V_{xx,t}=nabla^2 V_t(bar x_{t-1})`.

The corresponding local model is

```math
\begin{aligned}
\mathcal Q_t(\bar x+\delta x,\bar u+\delta u)
\approx{}&q_t+\mathcal Q_{x,t}^{\mathsf T}\delta x
+\mathcal Q_{u,t}^{\mathsf T}\delta u\\
&+\tfrac12\delta x^{\mathsf T}\mathcal Q_{xx,t}\delta x
+\delta u^{\mathsf T}\mathcal Q_{ux,t}\delta x
+\tfrac12\delta u^{\mathsf T}\mathcal Q_{uu,t}\delta u.
\end{aligned}
```

### Eq. (8): terminal derivatives

Because `V_{T+1}=ell_F`,

```math
v_{T+1}=\ell_F(\bar x_T),\qquad
V_{x,T+1}=\nabla\ell_F(\bar x_T),\qquad
V_{xx,T+1}=\nabla^2\ell_F(\bar x_T).
```

### Eq. (9): derivatives of the local Bellman function

Evaluate all stage derivatives at `(bar x_{t-1},bar u_t)` and all next-value
derivatives at `bar x_t=f_t(bar x_{t-1},bar u_t)`.  The chain rule gives

```math
\mathcal Q_{x,t}=\ell_{x,t}+f_{x,t}^{\mathsf T}V_{x,t+1},
```

```math
\mathcal Q_{u,t}=\ell_{u,t}+f_{u,t}^{\mathsf T}V_{x,t+1},
```

```math
\mathcal Q_{xx,t}
=\ell_{xx,t}+f_{x,t}^{\mathsf T}V_{xx,t+1}f_{x,t}
+V_{x,t+1}\mathbin{\cdot}f_{xx,t},
```

```math
\mathcal Q_{ux,t}
=\ell_{ux,t}+f_{u,t}^{\mathsf T}V_{xx,t+1}f_{x,t}
+V_{x,t+1}\mathbin{\cdot}f_{ux,t},
```

```math
\mathcal Q_{uu,t}
=\ell_{uu,t}+f_{u,t}^{\mathsf T}V_{xx,t+1}f_{u,t}
+V_{x,t+1}\mathbin{\cdot}f_{uu,t},
```

```math
q_t=\ell_t(\bar x_{t-1},\bar u_t)+v_{t+1}.
```

The last terms in the three Hessians are the second-order dynamics terms.
Dropping them produces iLQR/Gauss-Newton-style dynamics curvature, not the full
DDP equations stated by Ref. [21].

### Eq. (10): feedforward and feedback steps

Stationarity of the quadratic model with respect to `delta u_t` gives

```math
0=\mathcal Q_{u,t}+\mathcal Q_{uu,t}\delta u_t
+\mathcal Q_{ux,t}\delta x_{t-1}.
```

Therefore

```math
\delta u_t=k_t+K_t\delta x_{t-1},
```

where

```math
k_t=-\mathcal Q_{uu,t}^{-1}\mathcal Q_{u,t},\qquad
K_t=-\mathcal Q_{uu,t}^{-1}\mathcal Q_{ux,t}.
```

`k_t` changes the nominal period-`t` decision.  `K_t` says how that decision
should change if the actual entering state differs from `bar x_{t-1}`.

### Eq. (11): backward value recursion

Substituting the minimizing policy into the quadratic model yields

```math
V_{xx,t}=\mathcal Q_{xx,t}-K_t^{\mathsf T}\mathcal Q_{uu,t}K_t,
```

```math
V_{x,t}=\mathcal Q_{x,t}-K_t^{\mathsf T}\mathcal Q_{uu,t}k_t,
```

```math
v_t=q_t-\tfrac12 k_t^{\mathsf T}\mathcal Q_{uu,t}k_t.
```

Equivalent expanded formulas are obtained by substituting the definitions of
`k_t` and `K_t`.  This is the backward transfer of the optimized future cost's
level, slope, and curvature from `x_t` to the preceding state `x_{t-1}`.

### Eq. (12): forward rollout

With line-search fraction `alpha in (0,1]`, start from the fixed state
`x_0^+=x_0` and roll forward:

```math
u_t^+=\bar u_t+\alpha k_t+K_t(x_{t-1}^+-\bar x_{t-1}),
```

```math
x_t^+=f_t(x_{t-1}^+,u_t^+),\qquad t=1,\ldots,T.
```

## 4. Equations (13)-(20): augmented-Lagrangian background

These equations are not yet the DDP algorithm.  Ref. [21] first explains the
generic augmented-Lagrangian mechanism that it will place inside DDP.

### Eqs. (13)-(18): equality-constrained quadratic program

For

```math
\min_z\ \tfrac12 z^{\mathsf T}Hz+g^{\mathsf T}z
\quad\text{subject to}\quad Az=b,
```

the augmented Lagrangian is

```math
\mathscr L_\rho(z,\lambda)
=\tfrac12 z^{\mathsf T}Hz+g^{\mathsf T}z
+\lambda^{\mathsf T}(Az-b)+\tfrac\rho2\|Az-b\|^2.
```

Minimizing it at fixed `lambda,rho` gives

```math
z^+=-(H+\rho A^{\mathsf T}A)^{-1}
\left(g+A^{\mathsf T}\lambda-\rho A^{\mathsf T}b\right),
```

followed by

```math
\lambda^+=\lambda+\rho(Az^+-b).
```

Here this note uses `rho` for Ref. [21]'s penalty `mu`, avoiding confusion with
the user's dynamics dual `mu[t]`.  The numerically preferable saddle-point
form is

```math
\begin{bmatrix}
H&A^{\mathsf T}\\
A&-\rho^{-1}I
\end{bmatrix}
\begin{bmatrix}z^+\\\lambda^+\end{bmatrix}
=
\begin{bmatrix}-g\\b-\rho^{-1}\lambda\end{bmatrix}.
```

This system is algebraically equivalent but avoids directly forming the badly
conditioned `H+rho A^T A` for large `rho`.

### Eqs. (19)-(20): nonlinear equality constraints

For

```math
\min_z F(z)\quad\text{subject to}\quad c(z)=0,
```

the augmented Lagrangian is

```math
\mathscr L_\rho(z,\lambda)
=F(z)+\lambda^{\mathsf T}c(z)+\tfrac\rho2\|c(z)\|^2.
```

Ref. [21] uses a BCL-style outer update: approximately minimize this function
at fixed `(lambda,rho)`; update `lambda` when both stationarity and feasibility
are good enough; otherwise increase `rho` when stationarity is good but
feasibility is not.  Its adaptive suggestion estimates a constraint-reduction
factor `phi` and, for desired exponent `beta in (0,1)`, chooses

```math
\rho_{n+1}=\phi^{1/(1-\beta)}.
```

The penalty `rho` is not an inequality dual variable.  It controls how strongly
constraint violation is penalized.

## 5. Equations (21)-(24): DDP with globally constant multipliers

### Eq. (21): period equalities

Add

```math
h_t(x_{t-1},u_t)=0,\qquad t=1,\ldots,T,
```

where `h_t` contains period-`t` equality constraints.

### Eq. (22): horizon augmented Lagrangian

With one multiplier vector `lambda_t` per period,

```math
\mathscr L_\rho(x_{0:T},u_{1:T},\lambda_{1:T})
=\sum_{t=1}^{T}\left[
\ell_t(x_{t-1},u_t)
+\lambda_t^{\mathsf T}h_t(x_{t-1},u_t)
+\tfrac\rho2\|h_t(x_{t-1},u_t)\|^2
\right]+\ell_F(x_T).
```

During an inner DDP solve, `lambda_t` and `rho` are held fixed.

### Eq. (23): constrained local derivatives

Let

```math
r_t:=\lambda_t+\rho h_t.
```

The augmented terms modify the ordinary derivatives as follows:

```math
\mathcal Q_{x,t}^{c}
=\mathcal Q_{x,t}+h_{x,t}^{\mathsf T}r_t,
```

```math
\mathcal Q_{u,t}^{c}
=\mathcal Q_{u,t}+h_{u,t}^{\mathsf T}r_t,
```

```math
\mathcal Q_{xx,t}^{c}
=\mathcal Q_{xx,t}
+r_t\mathbin{\cdot}h_{xx,t}
+\rho h_{x,t}^{\mathsf T}h_{x,t},
```

```math
\mathcal Q_{ux,t}^{c}
=\mathcal Q_{ux,t}
+r_t\mathbin{\cdot}h_{ux,t}
+\rho h_{u,t}^{\mathsf T}h_{x,t},
```

```math
\mathcal Q_{uu,t}^{c}
=\mathcal Q_{uu,t}
+r_t\mathbin{\cdot}h_{uu,t}
+\rho h_{u,t}^{\mathsf T}h_{u,t},
```

```math
q_t^c=q_t+\lambda_t^{\mathsf T}h_t+\tfrac\rho2\|h_t\|^2.
```

These formulas come directly from differentiating
`lambda_t^T h_t+(rho/2)h_t^T h_t`.  The DDP equations (10)-(12) then use the
superscript-`c` quantities.

### Eq. (24): outer multiplier update

After the inner augmented-Lagrangian DDP problem is sufficiently minimized,

```math
\lambda_t^+=\lambda_t+\rho h_t(x_{t-1}^+,u_t^+).
```

This global-multiplier method improves the nominal trajectory but its feedback
law does not itself guarantee that an off-nominal state remains on the equality
constraint manifold.

## 6. Equations (25)-(28): locally affine multiplier functions

Ref. [21]'s second method lets each equality multiplier change with the state
entering its period.

### Eq. (25): affine multiplier model

Near the nominal entering state `bar x_{t-1}`,

```math
\lambda_t(x_{t-1})
=\bar\lambda_t+\Lambda_t(x_{t-1}-\bar x_{t-1}),
```

where `Lambda_t=partial lambda_t/partial x_{t-1}`.

### Eq. (26): augmented Lagrangian with multiplier policies

```math
\mathscr L_\rho
=\sum_{t=1}^{T}\left[
\ell_t(x_{t-1},u_t)
+\lambda_t(x_{t-1})^{\mathsf T}h_t(x_{t-1},u_t)
+\tfrac\rho2\|h_t(x_{t-1},u_t)\|^2
\right]+\ell_F(x_T).
```

### Eq. (27): derivative corrections

At the nominal point, again let `r_t=bar lambda_t+rho h_t`.  Direct
differentiation gives

```math
\mathcal Q_{x,t}^{c}
=\mathcal Q_{x,t}+h_{x,t}^{\mathsf T}r_t+\Lambda_t^{\mathsf T}h_t,
```

```math
\mathcal Q_{u,t}^{c}
=\mathcal Q_{u,t}+h_{u,t}^{\mathsf T}r_t,
```

```math
\mathcal Q_{xx,t}^{c}
=\mathcal Q_{xx,t}
+r_t\mathbin{\cdot}h_{xx,t}
+\rho h_{x,t}^{\mathsf T}h_{x,t}
+h_{x,t}^{\mathsf T}\Lambda_t
+\Lambda_t^{\mathsf T}h_{x,t},
```

```math
\mathcal Q_{ux,t}^{c}
=\mathcal Q_{ux,t}
+r_t\mathbin{\cdot}h_{ux,t}
+\rho h_{u,t}^{\mathsf T}h_{x,t}
+h_{u,t}^{\mathsf T}\Lambda_t,
```

```math
\mathcal Q_{uu,t}^{c}
=\mathcal Q_{uu,t}
+r_t\mathbin{\cdot}h_{uu,t}
+\rho h_{u,t}^{\mathsf T}h_{u,t},
```

```math
q_t^c=q_t+\bar\lambda_t^{\mathsf T}h_t
+\tfrac\rho2\|h_t\|^2.
```

The first-gradient term is `Lambda_t^T h_t`.  Some searchable text copies
render it as `Lambda_t^T h_{x,t}`, which is dimensionally impossible.  The
formula above follows from the product rule applied to
`lambda_t(x)^T h_t(x,u)`.

When the nominal state is shifted from `bar x_{t-1}` to `x_{t-1}^+`, the same
affine multiplier function is recentered by

```math
\bar\lambda_t\leftarrow
\bar\lambda_t+\Lambda_t(x_{t-1}^+-\bar x_{t-1}),
```

while `Lambda_t` itself is unchanged by this recentering.

### Eq. (28): multiplier-policy update

After an inner DDP minimization,

```math
\bar\lambda_t^+=\bar\lambda_t+\rho h_t,
```

```math
\Lambda_t^+=\Lambda_t
+\rho\left(h_{x,t}+h_{u,t}K_t\right).
```

The second equation is the chain rule along the feedback policy:

```math
\frac{d}{d x_{t-1}}h_t(x_{t-1},u_t(x_{t-1}))
=h_{x,t}+h_{u,t}K_t.
```

This is why the second method can correct both the control and the equality
multiplier when the entering state is perturbed.

## 7. Complete backward/forward workflow in this indexing

For each outer augmented-Lagrangian iteration:

1. Hold `lambda_{1:T}` (and `Lambda_{1:T}` in the affine version) and `rho`
   fixed.
2. Set the terminal derivatives from `V_{T+1}(x_T)=ell_F(x_T)`.
3. For `t=T,T-1,...,1`:
   - build the ordinary derivatives of `mathcal Q_t` from Eq. (9);
   - add either the constant-multiplier Eq. (23) or affine-multiplier Eq. (27)
     corrections;
   - solve for `k_t,K_t` using the corrected `mathcal Q_u^c`,
     `mathcal Q_uu^c`, and `mathcal Q_ux^c`;
   - update `v_t,V_{x,t},V_{xx,t}` using Eq. (11).
4. Starting from fixed `x_0`, roll `t=1,...,T` forward using Eq. (12) and a
   line search.
5. Repeat DDP backward/forward passes until the augmented Lagrangian is
   sufficiently stationary.
6. If equality feasibility is also adequate, update the multipliers using
   Eq. (24), or both parts of Eq. (28).  If stationarity is adequate but
   feasibility is not, increase `rho` according to the BCL rule.
7. Stop when both stationarity and equality feasibility meet their tolerances.

## 8. Relation to FilterDDP

Ref. [21] is an augmented-Lagrangian equality-constrained DDP method.
FilterDDP cites it as prior work but does **not** simply implement these outer
updates.  FilterDDP instead solves a constrained saddle-point system inside
each backward stage and globalizes the iteration with a filter line search.
The unconstrained Bellman/DDP backbone in Sections 2-3 above is shared; the
constraint treatment in Sections 4-6 is specific to Ref. [21].
