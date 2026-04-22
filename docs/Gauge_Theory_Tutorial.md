# A Full Tutorial on Gauge Theory
### From First Principles to Non-Abelian Fields and Transformer Connections

**Prerequisites assumed:** Group theory (Lie groups, Lie algebras, representations),
calculus of variations (action principles, Euler-Lagrange equations), basic Riemannian
geometry (manifolds, tangent bundles, covariant derivatives, curvature tensors).

---

## Table of Contents

1. [Motivation: What Problem Does Gauge Theory Solve?](#1-motivation)
2. [The Simplest Case: U(1) and Electromagnetism](#2-u1-electromagnetism)
3. [The Gauge Principle: From Global to Local Symmetry](#3-the-gauge-principle)
4. [Fiber Bundles: The Geometric Foundation](#4-fiber-bundles)
5. [Connections on Principal Bundles](#5-connections)
6. [Covariant Derivatives and Parallel Transport](#6-covariant-derivatives)
7. [Curvature: The Field Strength Tensor](#7-curvature)
8. [Non-Abelian Gauge Theories](#8-non-abelian)
9. [The Yang-Mills Action and Equations of Motion](#9-yang-mills)
10. [Gauge Transformations: Active and Passive](#10-gauge-transformations)
11. [Holonomy, Wilson Loops, and Topological Invariants](#11-holonomy)
12. [The BRST Formalism and Gauge Fixing](#12-brst)
13. [Instantons and Topological Charge](#13-instantons)
14. [The Standard Model as a Gauge Theory](#14-standard-model)
15. [Gauge Theory in Condensed Matter: Berry Phase](#15-condensed-matter)
16. [Gauge Theory and Transformers: The Connection](#16-transformers)
    - [16.5 Statement Decoder](#165-statement-decoder)
17. [Summary: The Unified Picture](#17-summary)

---

## 1. Motivation: What Problem Does Gauge Theory Solve?

### 1.1 The Redundancy Problem

Consider a physical system described by a field $\varphi(x)$. Suppose two field
configurations $\varphi(x)$ and $\varphi'(x)$ produce identical observable predictions
for every possible measurement. Then $\varphi$ and $\varphi'$ are physically equivalent —
they describe the same physical state. The map $\varphi \to \varphi'$ is called a
**gauge transformation**, and the associated symmetry is called **gauge invariance**
or **gauge symmetry**.

Gauge symmetry is fundamentally different from ordinary physical symmetry:

```
Ordinary symmetry (e.g., rotation):
  Two distinct states are related by symmetry
  Rotating a hydrogen atom gives a different (but equivalent) orientation
  The states ARE physically distinct — they just have the same energy

Gauge symmetry:
  Two field configurations describe THE SAME physical state
  It is a redundancy in the mathematical description
  Not a physical symmetry — a mathematical overcounting
```

The insight of gauge theory is that this redundancy, far from being a nuisance, is
enormously productive: **requiring that physics be invariant under local gauge
transformations completely determines the form of the forces between particles.**

### 1.2 The Historical Path

The development of gauge theory proceeded in three stages:

**Stage 1 (Weyl, 1918):** Attempted to unify electromagnetism and gravity by making
the metric scale (gauge) a local degree of freedom. Failed physically but introduced
the word "gauge" (from German *Eich* = scale).

**Stage 2 (Weyl, 1929; London, 1927):** Recognized that the correct gauge transformation
in quantum mechanics is a phase transformation $\psi \to e^{i\alpha(x)} \psi$, not a
scale transformation. This gives electromagnetism as a $U(1)$ gauge theory.

**Stage 3 (Yang & Mills, 1954):** Generalized from the abelian group $U(1)$ to
non-abelian groups $SU(N)$. This gave the framework for the strong and weak forces.
The resulting Yang-Mills theories are the foundation of the Standard Model.

---

## 2. U(1) and Electromagnetism

### 2.1 The Free Schrödinger Field

Begin with a free complex scalar field $\psi(x, t) \in \mathbb{C}$ in quantum mechanics.
The Lagrangian density (or action for the Schrödinger equation) has a global $U(1)$
symmetry:

$$
\psi(x) \to e^{i\alpha} \psi(x), \qquad \alpha \in \mathbb{R} \text{ constant (global)}
$$

The key word is **global**: the same phase $\alpha$ is applied everywhere in spacetime
simultaneously. The Lagrangian

$$
\mathcal{L} = i\psi^{\ast} \partial_t \psi - \frac{1}{2m}|\nabla \psi|^2 - V|\psi|^2
$$

is invariant under this global phase rotation because $\psi^{\ast}\psi$ and
$|\nabla \psi|^2$ are both invariant (the phases cancel).

### 2.2 Promoting to a Local Symmetry

Now ask: what if we want invariance under a **local** $U(1)$ transformation?

$$
\psi(x) \to e^{i\alpha(x)} \psi(x), \qquad \alpha(x) \text{ position-dependent}
$$

Under this transformation, the gradient term transforms as

$$
\nabla \psi \to \nabla(e^{i\alpha} \psi) = e^{i\alpha}(\nabla \psi + i(\nabla \alpha)\psi).
$$

The extra term $i(\nabla \alpha)\psi$ breaks the invariance. The Lagrangian is **not**
invariant under local $U(1)$.

**The fix:** Replace the ordinary derivative $\partial_\mu$ with a **covariant
derivative**:

$$
D_\mu = \partial_\mu - iq A_\mu(x),
$$

where $A_\mu(x)$ is a new field (the gauge field) that transforms simultaneously with
$\psi$:

$$
\begin{aligned}
\psi(x) &\to e^{iq\alpha(x)} \psi(x), \\
A_\mu(x) &\to A_\mu(x) + \partial_\mu \alpha(x).
\end{aligned}
$$

Under these simultaneous transformations:

$$
\begin{aligned}
D_\mu \psi
&\to \partial_\mu(e^{iq\alpha} \psi) - iq(A_\mu + \partial_\mu \alpha)(e^{iq\alpha} \psi) \\
&= e^{iq\alpha}(\partial_\mu \psi + iq(\partial_\mu \alpha)\psi - iq A_\mu \psi - iq(\partial_\mu \alpha)\psi) \\
&= e^{iq\alpha} D_\mu \psi.
\end{aligned}
$$

The extra $iq(\partial_\mu \alpha)\psi$ terms cancel exactly. The covariant derivative
transforms covariantly: $D_\mu \psi \to e^{iq\alpha} D_\mu \psi$.

### 2.3 The Gauge Field IS the Electromagnetic Potential

The field $A_\mu(x) = (\varphi, \mathbf{A})$ that we were forced to introduce to restore
local invariance is precisely the electromagnetic four-potential:

- $A_0 = \varphi$: the electric scalar potential
- $A_i$: the magnetic vector potential

The electromagnetic fields are

$$
\begin{aligned}
E_i &= -\partial_i A_0 - \partial_0 A_i, \\
B_i &= \epsilon_{ijk} \partial_j A_k.
\end{aligned}
$$

In covariant form, the field strength tensor

$$
F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu
$$

is gauge invariant:

$$
\begin{aligned}
F_{\mu\nu} &\to \partial_\mu(A_\nu + \partial_\nu \alpha) - \partial_\nu(A_\mu + \partial_\mu \alpha) \\
&= F_{\mu\nu} + \partial_\mu \partial_\nu \alpha - \partial_\nu \partial_\mu \alpha \\
&= F_{\mu\nu} \qquad \text{(mixed partials commute).}
\end{aligned}
$$

### 2.4 The Maxwell Action from Gauge Invariance

The simplest Lorentz-invariant, gauge-invariant kinetic term for $A_\mu$ is

$$
S_{\text{Maxwell}} = -\frac{1}{4} \int F_{\mu\nu} F^{\mu\nu} d^4x.
$$

Adding the matter coupling (minimal coupling):

$$
S_{\text{matter}} = \int \psi^{\ast}(i D_t - H_0)\psi  d^4x
$$

gives the full quantum electrodynamics (QED) action. Varying with respect to $A_\mu$
gives Maxwell's equations:

$$
\begin{aligned}
\partial_\nu F^{\nu\mu} &= j^\mu && \text{(sourced Maxwell equations)}, \\
\partial_{[\mu} F_{\nu\rho]} &= 0 && \text{(Bianchi identity — automatic from } F = dA\text{).}
\end{aligned}
$$

**The profound conclusion:** Requiring local $U(1)$ invariance of the free matter
Lagrangian uniquely determines the form of electromagnetism. The photon field $A_\mu$
is not put in by hand — it is forced by the gauge principle.

---

## 3. The Gauge Principle

### 3.1 The Principle Stated

The gauge principle is the organizing idea of modern physics:

> **Gauge Principle:** If a physical theory has a global symmetry $G$, promote it to
> a local symmetry. The gauge fields required to restore local invariance are the
> force-carrying fields (gauge bosons).

| Global symmetry | Local gauge symmetry | Gauge field | Force |
|---|---|---|---|
| $U(1)$ phase | Local $U(1)$ | Photon $A_\mu$ | Electromagnetism |
| $SU(2)$ isospin | Local $SU(2)$ | $W^a_\mu$ (3 fields) | Weak force |
| $SU(3)$ color | Local $SU(3)$ | Gluons $G^a_\mu$ (8 fields) | Strong force |
| Diffeomorphism | Local Poincaré | Metric $g_{\mu\nu}$ | Gravity |

### 3.2 Counting Gauge Bosons

For a gauge group $G$ with Lie algebra $\mathfrak{g}$, the number of gauge bosons equals

$$
\dim(G) = \dim(\mathfrak{g}) = \text{number of generators of } G.
$$

| Group | dim | Gauge bosons | Physics |
|---|---|---|---|
| $U(1)$ | 1 | 1 photon | QED |
| $SU(2)$ | 3 | 3 weak bosons $W^{+}, W^{-}, Z$ | Weak interaction |
| $SU(3)$ | 8 | 8 gluons | QCD |
| $U(1) \times SU(2) \times SU(3)$ | 12 | 12 Standard Model bosons | Full SM |

### 3.3 Why Local Symmetry Is More Restrictive Than Global

Global symmetry: transforming all fields simultaneously by the same group element
leaves the physics unchanged. This is a statement about the structure of the theory.

Local symmetry: transforming fields by **different** group elements at each spacetime
point leaves the physics unchanged. This is a much stronger statement — it means
the theory has no physical information in the relative phases/orientations at
different spacetime points. It also implies:

1. The gauge boson is massless (a massive photon would break gauge invariance).
2. The coupling structure is uniquely determined (minimal coupling).
3. The self-interactions of gauge bosons (in the non-abelian case) are fixed.
4. There are conservation laws (Noether's theorem applied locally).

---

## 4. Fiber Bundles: The Geometric Foundation

### 4.1 The Bundle Picture

The geometric framework for gauge theory is the theory of fiber bundles. A fiber
bundle makes precise the idea of "a space with additional structure attached at
each point."

**Definition (Fiber Bundle).** A fiber bundle is a quadruple $(E, B, \pi, F)$ where:

- $E$ is the **total space**.
- $B$ is the **base space** (spacetime in physics).
- $\pi: E \to B$ is the **projection** (surjective).
- $F$ is the **typical fiber** (the space attached at each point).
- Local trivialization: for each $x \in B$, there exists an open set $U \ni x$ and
  a homeomorphism $\varphi: \pi^{-1}(U) \to U \times F$ such that
  $\pi = \mathrm{proj}_1 \circ \varphi$.

The fiber over a point $x \in B$ is the preimage $\pi^{-1}(x) \cong F$. Schematically:

$$
E \xrightarrow{\pi} B,
$$

so each point $b \in B$ has a copy of $F$ sitting above it: $\pi^{-1}(b) \cong F$.

### 4.2 The Three Bundles of Gauge Theory

Gauge theory involves three related bundles.

**Principal Bundle $P(B, G)$:**

- Fiber = the gauge group $G$ itself.
- Total space = collection of all gauge frames at each point.
- Right action of $G$ on $P$: for $p \in P$, $g \in G$, $R_g(p) = pg$.
- Transition functions are $G$-valued.

**Associated Vector Bundle $E = P \times_G V$:**

- Fiber = a representation $V$ of $G$ (e.g., $\mathbb{R}^n$, $\mathbb{C}^n$).
- Sections of $E$ are the matter fields (quarks, electrons, etc.).
- The gauge group acts on $V$ via the representation $\rho: G \to GL(V)$.

**Adjoint Bundle $\mathrm{ad}(P) = P \times_G \mathfrak{g}$:**

- Fiber = the Lie algebra $\mathfrak{g}$ of $G$.
- $G$ acts on $\mathfrak{g}$ by the adjoint representation:
  $\mathrm{Ad}_g(X) = g X g^{-1}$.
- Sections of $\mathrm{ad}(P)$ are the gauge field strengths.

### 4.3 Sections and Matter Fields

A **section** of a bundle $E \to B$ is a smooth map $s: B \to E$ such that
$\pi \circ s = \mathrm{id}_B$. Locally, a section looks like a function $B \to F$.

In physics:

$$
\begin{aligned}
\text{Matter fields } \psi &= \text{sections of associated vector bundle } E, \\
\text{Gauge potentials } A_\mu &= \text{connection 1-form on principal bundle } P, \\
\text{Field strengths } F_{\mu\nu} &= \text{curvature 2-form of the connection.}
\end{aligned}
$$

### 4.4 Transition Functions and Gauge Transformations

A bundle may not be globally trivial — it may be impossible to choose a consistent
global section. Instead, it is covered by local trivializations $\{U_\alpha\}$ with
transition functions $g_{\alpha\beta}: U_\alpha \cap U_\beta \to G$ satisfying the
cocycle condition

$$
g_{\alpha\beta}  g_{\beta\gamma}  g_{\gamma\alpha} = e \quad \text{on triple overlaps } U_\alpha \cap U_\beta \cap U_\gamma.
$$

A **gauge transformation** is precisely a change of local trivialization — a smooth
map $g: B \to G$ (or more precisely, a vertical automorphism of $P$).

Under a gauge transformation $g: U_\alpha \to G$:

$$
\begin{aligned}
\psi &\to \rho(g) \psi && \text{(matter field transforms in representation } \rho\text{)}, \\
A_\mu &\to g A_\mu g^{-1} + g \partial_\mu g^{-1} && \text{(gauge potential transforms as a connection)}, \\
F_{\mu\nu} &\to g F_{\mu\nu} g^{-1} && \text{(field strength transforms covariantly).}
\end{aligned}
$$

---

## 5. Connections on Principal Bundles

### 5.1 The Connection 1-Form

A **connection** on a principal bundle $P(B, G)$ is a $\mathfrak{g}$-valued 1-form
$\omega \in \Omega^1(P; \mathfrak{g})$ satisfying two conditions:

1. **Right equivariance:** $R_g^{\ast} \omega = \mathrm{Ad}_{g^{-1}} \omega$ for all
   $g \in G$.
2. **Vertical normalization:** $\omega(\hat V) = V$ for all $V \in \mathfrak{g}$,
   where $\hat V$ is the fundamental vector field generated by $V$.

The connection splits the tangent space $T_p P$ at each $p \in P$ into

$$
T_p P = V_p P \oplus H_p P,
$$

where $V_p P = \ker(d\pi_p)$ is the tangent space to the fiber (vertical), and
$H_p P = \ker(\omega_p)$ is the horizontal subspace (determined by the connection).

### 5.2 Local Representatives: The Gauge Potential

Given a local section $s: U \to P$ (a choice of gauge), the **local connection form**
(gauge potential) is the pullback

$$
A = s^{\ast}\omega \in \Omega^1(U; \mathfrak{g}).
$$

In physics notation $A = A_\mu  dx^\mu$, where $A_\mu(x) \in \mathfrak{g}$ is a
Lie-algebra-valued field.

- For $G = U(1)$: $A_\mu$ is a real-valued function (the electromagnetic potential).
- For $G = SU(N)$: $A_\mu = A_\mu^a T_a$ where $T_a$ are the generators of $SU(N)$,
  and $a = 1, \ldots, N^2 - 1$.

Under a gauge transformation $g: U \to G$ (change of section):

$$
A_\mu \to g A_\mu g^{-1} + g \partial_\mu g^{-1}.
$$

For $U(1)$ with $g = e^{i\alpha}$: $A_\mu \to A_\mu + \partial_\mu \alpha$ (recovers
the earlier formula).

### 5.3 The Covariant Derivative from the Connection

Given a connection $A$ on $P$ and a matter field $\psi$ in representation $\rho$ of
$G$, the **covariant derivative** is

$$
D_\mu \psi = \partial_\mu \psi + \rho(A_\mu) \psi.
$$

For the fundamental representation of $SU(N)$, $\rho(A_\mu) = A_\mu^a T_a^{\mathrm{fund}}$, so

$$
D_\mu \psi = \partial_\mu \psi + A_\mu^a T_a \psi \qquad \text{(matrix acting on } \psi \in \mathbb{C}^N\text{).}
$$

The covariant derivative measures how much $\psi$ deviates from parallel transport:
$D_\mu \psi = 0$ means $\psi$ is parallel-transported in the direction $\mu$.

**Properties:**

1. Gauge covariance: $D_\mu(g\psi) = g (D_\mu \psi)$ for gauge transformation $g$.
2. Leibniz rule: $D_\mu(\psi_1 \psi_2) = (D_\mu \psi_1)\psi_2 + \psi_1 (D_\mu \psi_2)$.
3. Reduces to $\partial_\mu$ when $A = 0$ (in the trivial connection).
4. The commutator $[D_\mu, D_\nu] = F_{\mu\nu}$ (the curvature — see §7).

---

## 6. Covariant Derivatives and Parallel Transport

### 6.1 Parallel Transport

Given a curve $\gamma: [0,1] \to B$ and a vector $v_0$ in the fiber over $\gamma(0)$,
the **parallel transport** of $v_0$ along $\gamma$ is the unique section $v(t)$ of
$E$ along $\gamma$ satisfying

$$
D_{\dot\gamma} v = 0 \quad \text{with} \quad v(0) = v_0,
$$

where $\dot\gamma = d\gamma / dt$ is the tangent vector to $\gamma$.

In local coordinates, with $\gamma(t)$ having components $x^\mu(t)$:

$$
\frac{d v^i}{dt} + A_\mu^{ij}(\gamma(t))  \dot x^\mu v^j = 0.
$$

This is a linear ODE — it always has a unique solution. The parallel transport map

$$
\tau_\gamma: E_{\gamma(0)} \to E_{\gamma(1)}
$$

is an isomorphism of the fibers (a group element in $G$).

### 6.2 Why Parallel Transport Depends on Path

In a curved bundle (non-zero curvature), parallel transport is path-dependent.
Transporting a vector around a closed loop $\gamma$ returns a vector that is NOT in
general equal to the original. The transformation

$$
\mathrm{Hol}_\gamma: E_{\gamma(0)} \to E_{\gamma(0)}
$$

is an element of $G$ called the **holonomy** of the connection around $\gamma$.
Non-trivial holonomy = non-zero curvature = physical force.

### 6.3 Analogy with Riemannian Geometry

The analogy with the Levi-Civita connection on a Riemannian manifold is exact:

| Riemannian geometry | Gauge theory |
|---|---|
| Tangent bundle $TM$ | Associated vector bundle $E$ |
| Levi-Civita connection $\Gamma^k_{ij}$ | Gauge potential $A_\mu$ |
| Covariant derivative $\nabla_\mu V^k = \partial_\mu V^k + \Gamma^k_{\mu j} V^j$ | $D_\mu \psi = \partial_\mu \psi + A_\mu \psi$ |
| Riemannian curvature $R^k{}_{l\mu\nu}$ | Field strength $F_{\mu\nu}$ |
| Geodesic: $\nabla_{\dot\gamma} \dot\gamma = 0$ | Parallel section: $D_\mu \psi = 0$ |
| Holonomy group $\subseteq O(n)$ | Gauge holonomy $\in G$ |
| Metric compatibility $\nabla g = 0$ | Gauge compatibility $D^\dagger D = \ldots$ |

The difference: in Riemannian geometry, the connection is determined by the metric
(Levi-Civita theorem). In gauge theory, the connection $A$ is an independent
dynamical field with its own equations of motion.

---

## 7. Curvature: The Field Strength Tensor

### 7.1 Definition via the Commutator of Covariant Derivatives

The **curvature** of a connection is the $\mathfrak{g}$-valued 2-form defined by

$$
F_{\mu\nu} = [D_\mu, D_\nu] \qquad \text{(commutator of covariant derivatives).}
$$

Expanding explicitly:

$$
\begin{aligned}
{[D_\mu, D_\nu]} \psi
&= [\partial_\mu + A_\mu, \partial_\nu + A_\nu] \psi \\
&= (\partial_\mu \partial_\nu - \partial_\nu \partial_\mu)\psi
   + (\partial_\mu A_\nu - \partial_\nu A_\mu)\psi
   + A_\mu A_\nu \psi - A_\nu A_\mu \psi \\
&= 0 + (\partial_\mu A_\nu - \partial_\nu A_\mu)\psi + [A_\mu, A_\nu] \psi.
\end{aligned}
$$

Therefore

$$
F_{\mu\nu} = \underbrace{\partial_\mu A_\nu - \partial_\nu A_\mu}_{\text{abelian part (same as Maxwell)}} + \underbrace{[A_\mu, A_\nu]}_{\text{non-abelian part (zero for } U(1)\text{)}}.
$$

### 7.2 The Critical Difference Between Abelian and Non-Abelian

For **$U(1)$** (electromagnetism):

- $A_\mu$ is a real number (or imaginary scalar).
- $[A_\mu, A_\nu] = A_\mu A_\nu - A_\nu A_\mu = 0$ (numbers commute).
- $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$ (the Maxwell tensor).
- $F_{\mu\nu}$ is gauge invariant: $F_{\mu\nu} \to F_{\mu\nu}$.

For **$SU(N)$** (non-abelian):

- $A_\mu = A_\mu^a T_a$ is a matrix (Lie-algebra-valued).
- $[A_\mu, A_\nu] = A_\mu^a A_\nu^b [T_a, T_b] = f^{abc} A_\mu^a A_\nu^b T_c \neq 0$.
- $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu]$
  (non-linear in $A$!).
- $F_{\mu\nu}$ is only covariant: $F_{\mu\nu} \to g F_{\mu\nu} g^{-1}$.

The non-abelian term $[A_\mu, A_\nu]$ is the source of self-interactions of gauge
bosons. In QED, photons do not interact with each other (superposition principle
holds exactly). In QCD, gluons carry color charge and interact with each other —
which is why $[A_\mu, A_\nu] \neq 0$ matters fundamentally.

### 7.3 The Bianchi Identity

The curvature satisfies the **Bianchi identity**

$$
D_{[\mu} F_{\nu\rho]} = 0,
$$

or equivalently

$$
D_\mu F_{\nu\rho} + D_\nu F_{\rho\mu} + D_\rho F_{\mu\nu} = 0,
$$

where $D_\mu F_{\nu\rho} = \partial_\mu F_{\nu\rho} + [A_\mu, F_{\nu\rho}]$ is the
covariant derivative of $F$.

For $U(1)$: this reduces to $\partial_{[\mu} F_{\nu\rho]} = 0$, which gives the
source-free Maxwell equations $\nabla \cdot B = 0$ and
$\nabla \times E + \partial B / \partial t = 0$.

### 7.4 Obstruction to a Scalar Potential

Here is the precise statement of why non-zero curvature obstructs conservative
dynamics.

**Theorem (Gauge Obstruction).** The force field $F(x) = -A_\mu(x)  dx^\mu$ is the
gradient of a scalar function (i.e., $A_\mu = \partial_\mu V$ for some $V$) if and
only if $F_{\mu\nu} = 0$ everywhere.

**Proof.** ($\Rightarrow$) If $A_\mu = \partial_\mu V$, then

$$
F_{\mu\nu} = \partial_\mu \partial_\nu V - \partial_\nu \partial_\mu V = 0 \qquad \text{(Clairaut's theorem).}
$$

($\Leftarrow$) If $F_{\mu\nu} = 0$ everywhere, the connection is **flat**. By the
Poincaré lemma (or more precisely, the de Rham cohomology vanishing theorem on
contractible spaces), a flat connection is locally pure gauge:
$A_\mu = g^{-1} \partial_\mu g$ for some smooth $G$-valued function $g$. In the
abelian case with $G = U(1)$ and $g = e^{iV}$, this gives $A_\mu = \partial_\mu V$.
$\blacksquare$

**Corollary.** For multi-head attention with $[A_\mu^{(h)}, A_\mu^{(h')}] \neq 0$,
the curvature $F_{\mu\nu} \neq 0$, and therefore no scalar potential exists that can
reproduce the attention force field. This is an exact mathematical obstruction, not
an approximation failure.

**Why the obstruction holds regardless of capacity.**

This point deserves emphasis because it is easy to misread the experimental
result ($R^2 = 0.04$–$0.20$ for GPT-2 middle layers) as a capacity problem — as
if a sufficiently large or deep $V_\psi$ would eventually succeed.

It cannot. The reason is Clairaut's theorem, which is not a statement about
approximation but an algebraic identity. For ANY smooth function
$V: \mathbb{R}^d \to \mathbb{R}$, however complex,

$$
\partial_i(-\partial_j V) - \partial_j(-\partial_i V) = -\partial_i \partial_j V + \partial_j \partial_i V = 0
$$

because mixed partial derivatives commute (Schwarz / Clairaut theorem).

So the curvature of $-\nabla V$ is **exactly zero** — not approximately zero,
not zero in the limit of large capacity — zero as a mathematical fact holding
for every smooth scalar $V$ without exception.

The attention force field $F$ has non-zero curvature (from
$[A^{(h)}, A^{(h')}] \neq 0$):

$$
\begin{aligned}
\mathrm{Curvature}(F) &= [A, A] \neq 0 && \text{(from cross-head commutators)}, \\
\mathrm{Curvature}(-\nabla V) &= 0 && \text{(Clairaut, for any } V\text{).}
\end{aligned}
$$

Therefore $F \neq -\nabla V$ for any $V$, regardless of its complexity or capacity.

The analogy: asking whether a large enough $V_\psi$ can overcome this is like
asking whether a large enough continuous map $f: [0,1] \to \mathbb{R}$ can be
surjective onto $\mathbb{R}^2$ — the answer is no, and adding more parameters to
$f$ does not change the topological obstruction.

The experimental $R^2$ failure is the **quantitative measurement** of a structural
impossibility, not evidence of a capacity bottleneck.

---

## 8. Non-Abelian Gauge Theories

### 8.1 The Yang-Mills Construction (1954)

Yang and Mills generalized the $U(1)$ gauge theory to $SU(2)$ by:

1. Identifying the matter field as an $SU(2)$ doublet $\psi = (\psi_1, \psi_2)^{\top}$.
2. Requiring invariance under local $SU(2)$: $\psi \to U(x) \psi$, $U(x) \in SU(2)$.
3. Introducing the gauge potential $A_\mu = A_\mu^a \tau_a / 2$ ($\tau_a$ = Pauli
   matrices).
4. Defining $D_\mu \psi = (\partial_\mu + ig A_\mu)\psi$.
5. Writing the gauge-invariant field strength
   $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu + ig[A_\mu, A_\nu]$.

The gauge transformation law. Under $U(x) = e^{i\alpha^a(x) \tau_a / 2} \in SU(2)$:

$$
\begin{aligned}
\psi &\to U \psi, \\
A_\mu &\to U A_\mu U^\dagger + \frac{i}{g}(\partial_\mu U) U^\dagger, \\
F_{\mu\nu} &\to U F_{\mu\nu} U^\dagger.
\end{aligned}
$$

### 8.2 The Lie Algebra Structure Constants

For a general gauge group $G$ with generators $T_a$ satisfying

$$
[T_a, T_b] = i f^{abc} T_c
$$

(where $f^{abc}$ are the structure constants of the Lie algebra $\mathfrak{g}$), the
field strength is

$$
F_{\mu\nu}^c = \partial_\mu A_\nu^c - \partial_\nu A_\mu^c - g f^{abc} A_\mu^a A_\nu^b.
$$

The term $f^{abc} A_\mu^a A_\nu^b$ is the source of three- and four-point
self-interactions of the gauge bosons.

- For $SU(2)$: $f^{abc} = \epsilon^{abc}$ (Levi-Civita symbol).
- For $SU(3)$: $f^{abc}$ are the $SU(3)$ structure constants
  ($8 \times 8 \times 8$ tensor).

### 8.3 Why Non-Abelian Theories Are Self-Interacting

The gauge bosons themselves carry the "charge" of the gauge group. In QED:

- The photon is electrically neutral (charge 0).
- Therefore photons do not couple to each other directly.
- Maxwell's equations are linear.

In QCD ($SU(3)$ gauge theory):

- Gluons carry color charge (they are in the adjoint representation of $SU(3)$).
- Therefore gluons couple to each other.
- The Yang-Mills equations are non-linear (3-gluon and 4-gluon vertices).

This non-linearity is responsible for:

- **Asymptotic freedom**: the coupling weakens at high energies
  (Gross, Politzer, Wilczek 1973).
- **Confinement**: quarks cannot be isolated at low energies
  (still not rigorously proved!).
- **Glueballs**: bound states of gluons with no quark content.

### 8.4 Representations and Matter Content

The matter content of a gauge theory is specified by choosing a representation of
$G$ for each matter field:

- **Fundamental representation:** quarks in QCD (dimension $N$ for $SU(N)$).
- **Anti-fundamental representation:** antiquarks.
- **Adjoint representation:** gauge bosons themselves (dimension $N^2 - 1$ for
  $SU(N)$).
- **Singlet (trivial representation):** no gauge coupling.

The covariant derivative in representation $R$:

$$
D_\mu \psi = \partial_\mu \psi + i g A_\mu^a T_a^R \psi,
$$

where $T_a^R$ are the generators in representation $R$.

---

## 9. The Yang-Mills Action and Equations of Motion

### 9.1 The Yang-Mills Action

The gauge-invariant kinetic term for the gauge field is

$$
S_{\text{YM}} = -\frac{1}{2 g^2} \int \mathrm{tr}(F_{\mu\nu} F^{\mu\nu})  d^4x
  = -\frac{1}{4 g^2} \int F_{\mu\nu}^a F^{\mu\nu a}  d^4x.
$$

The trace is over the Lie algebra indices (using the Killing form
$\mathrm{tr}(T_a T_b) = \delta_{ab} / 2$). This action is:

- Gauge invariant: $F_{\mu\nu} \to U F_{\mu\nu} U^\dagger$, so
  $\mathrm{tr}(F^2) \to \mathrm{tr}(U F^2 U^\dagger) = \mathrm{tr}(F^2)$.
- Lorentz invariant (indices contracted with the Minkowski metric).
- Second order in the gauge field (from the kinetic term
  $\partial_\mu A_\nu - \partial_\nu A_\mu$).
- Contains cubic and quartic self-interaction terms from $[A_\mu, A_\nu]^2$.

### 9.2 The Full Action with Matter

$$
\begin{aligned}
S &= S_{\text{YM}} + S_{\text{matter}}, \\
S_{\text{matter}} &= \int \bar\psi (i\gamma^\mu D_\mu - m)\psi  d^4x
  \qquad \text{(Dirac fermions in representation } R\text{)} \\
&= \int \bar\psi (i\gamma^\mu \partial_\mu - m)\psi  d^4x
  + g \int \bar\psi \gamma^\mu A_\mu^a T_a^R \psi  d^4x
  \qquad \text{(minimal coupling).}
\end{aligned}
$$

### 9.3 The Yang-Mills Equations

Varying $S$ with respect to $A_\mu^a$ gives the **Yang-Mills equations**

$$
D_\nu F^{\nu\mu a} = j^{\mu a} \qquad \text{(sourced Yang-Mills equations)},
$$

where
$D_\nu F^{\nu\mu a} = \partial_\nu F^{\nu\mu a} + f^{abc} A_\nu^b F^{\nu\mu c}$
(the covariant divergence), and
$j^{\mu a} = g \bar\psi \gamma^\mu T_a^R \psi$ is the matter current.

Together with the Bianchi identity $D_{[\mu} F_{\nu\rho]} = 0$, these are the
complete equations of motion for the Yang-Mills gauge field.

**Comparison with Maxwell:**

$$
\begin{aligned}
\text{Maxwell:} \quad & \partial_\nu F^{\nu\mu} = j^\mu && \text{(linear, abelian)}, \\
\text{Yang-Mills:} \quad & D_\nu F^{\nu\mu} = j^\mu && \text{(non-linear, non-abelian).}
\end{aligned}
$$

The covariant derivative on the Yang-Mills side includes the $[A, F]$ term.

---

## 10. Gauge Transformations: Active and Passive

### 10.1 Two Perspectives

**Passive gauge transformation:** A change of local coordinates in the fiber bundle
(change of basis in the internal space at each point). The physics does not change;
only our description does. This is the "redundancy" perspective.

**Active gauge transformation:** A genuine transformation of the fields that maps
one physical configuration to another equivalent one. In the Hamiltonian formalism,
this generates constraints (Gauss's law in electromagnetism).

### 10.2 Small vs. Large Gauge Transformations

**Small gauge transformations:** Continuously connected to the identity. In the
path integral formulation, these correspond to the gauge redundancy that must be
fixed by Faddeev-Popov procedure or BRST (see §12).

**Large gauge transformations:** Not continuously connected to the identity; they
are topologically non-trivial. Parametrized by the homotopy group $\pi_n(G)$ for
appropriate $n$.

For $G = SU(2)$ on $S^3$ (compactified $\mathbb{R}^3$):
$\pi_3(SU(2)) = \pi_3(S^3) = \mathbb{Z}$. The integer winding number classifies
large gauge transformations. This is related to instantons (see §13).

### 10.3 The Gauge Orbit and Physical States

The space of all gauge field configurations $\mathcal{A}$ is acted on by the group
$\mathcal{G}$ of gauge transformations. The **gauge orbit** of a configuration $A$
is

$$
\mathcal{O}(A) = \{A^g : g \in \mathcal{G}\} = \{g A g^{-1} + g \partial g^{-1} : g \in \mathcal{G}\}.
$$

All configurations in the same gauge orbit represent the same physical state. The
physical configuration space is

$$
\mathcal{A} / \mathcal{G} = \{\text{gauge orbits}\}.
$$

This quotient space is generically singular (Gribov copies, Singer-Atiyah theorem),
which creates subtleties in the path integral quantization.

---

## 11. Holonomy, Wilson Loops, and Topological Invariants

### 11.1 Wilson Lines and Loops

Given a path $\gamma: [0,1] \to B$, the **Wilson line** is the parallel transport
operator

$$
W(\gamma) = \mathcal{P} \exp\left( \int_\gamma A_\mu  dx^\mu \right),
$$

where $\mathcal{P}$ denotes **path ordering** (necessary because $A_\mu$ at different
points do not commute in the non-abelian case):

$$
\mathcal{P} \exp\left( \int_0^1 A(t)  dt \right)
  = \lim_{N \to \infty} [1 + A(t_1) \Delta t][1 + A(t_2) \Delta t] \cdots [1 + A(t_N) \Delta t].
$$

Under a gauge transformation $g$:
$W(\gamma) \to g(\gamma(0))  W(\gamma)  g(\gamma(1))^{-1}$.

For a **closed loop** $\gamma$ with $\gamma(0) = \gamma(1) = x_0$, the **Wilson
loop** is

$$
W(C) = \mathrm{tr}\left( \mathcal{P} \exp\left( \oint_C A_\mu  dx^\mu \right) \right) \qquad \text{(trace makes it gauge invariant).}
$$

Wilson loops are the fundamental gauge-invariant observables of Yang-Mills theory.

### 11.2 The Area Law and Confinement

In QCD, the expectation value of a Wilson loop $C$ of area $A$ and perimeter $P$
satisfies

$$
\begin{aligned}
\langle W(C) \rangle &\sim e^{-\sigma A} && \text{(area law — confining phase)}, \\
\langle W(C) \rangle &\sim e^{-\kappa P} && \text{(perimeter law — deconfined phase).}
\end{aligned}
$$

The **area law** signals confinement: the potential between a quark and antiquark
grows linearly with distance; $\sigma$ is the string tension.

### 11.3 The Aharonov-Bohm Effect

The cleanest physical manifestation of holonomy: an electron moving in a region
with $B = 0$ but $A \neq 0$ (outside a solenoid) acquires a phase

$$
\varphi_{\text{AB}} = \frac{q}{\hbar} \oint_C A_\mu  dx^\mu
  = \frac{q}{\hbar} \iint B \cdot dS
  = \frac{q}{\hbar} \Phi_B.
$$

This phase is physical (observable via interference) even though $E = B = 0$ along
the electron's path. It demonstrates that $A_\mu$ (not just $F_{\mu\nu}$) is the
fundamental physical quantity — the holonomy is physical even when the curvature
locally vanishes.

### 11.4 Chern Classes: Global Topological Invariants

The **Chern classes** are topological invariants built from the curvature that do
not depend on the specific connection.

**First Chern class (for $U(1)$ bundles):**

$$
c_1 = \frac{i}{2\pi} \int_M F \qquad (F = F_{\mu\nu}  dx^\mu \wedge dx^\nu, \text{ integrated over the base}).
$$

For a $U(1)$ bundle over $S^2$, $c_1 \in \mathbb{Z}$ is the magnetic monopole charge.

**Second Chern class (for $SU(2)$ bundles):**

$$
c_2 = \frac{1}{8 \pi^2} \int_M \mathrm{tr}(F \wedge F)
  = \frac{1}{8 \pi^2} \int_M \mathrm{tr}(F_{\mu\nu} F^{\mu\nu})  d^4x.
$$

This is the **instanton number** (topological charge). For $SU(2)$ on $S^4$:
$c_2 \in \mathbb{Z}$.

The Chern classes classify principal bundles up to isomorphism — they are the
topological quantum numbers of gauge field configurations.

---

## 12. The BRST Formalism and Gauge Fixing

### 12.1 The Problem with Path Integral Quantization

Naively quantizing a gauge theory via the path integral

$$
Z = \int \mathcal{D}A  e^{i S[A]}
$$

is problematic: the integral overcounts physically equivalent configurations
(gauge copies). The integrand is constant along gauge orbits, so the integral
over $\mathcal{A}$ diverges by the volume of the gauge group $\mathcal{G}$.

### 12.2 The Faddeev-Popov Procedure

**Step 1.** Insert
$1 = \Delta_{\text{FP}}[A] \int \mathcal{D}g  \delta(G(A^g))$ into the path
integral, where $G(A) = 0$ is the gauge-fixing condition (e.g., Lorenz gauge
$\partial^\mu A_\mu = 0$) and $\Delta_{\text{FP}}$ is the Faddeev-Popov determinant

$$
\Delta_{\text{FP}}[A] = \det\left(\frac{\partial G}{\partial \alpha}\right)|_{G = 0} \qquad \text{(functional determinant).}
$$

**Step 2.** Write $\Delta_{\text{FP}}$ as a Grassmann (fermionic) path integral:

$$
\Delta_{\text{FP}}[A] = \int \mathcal{D}\bar c  \mathcal{D}c  \exp\left( i \int \bar c  \frac{\partial G}{\partial \alpha}  c  dx \right),
$$

where $c$, $\bar c$ are **Faddeev-Popov ghosts** — anticommuting scalar fields with
the wrong spin-statistics relation (necessary for consistency of the path integral,
not physical particles).

**Step 3.** The gauge-fixed path integral:

$$
Z = \int \mathcal{D}A  \mathcal{D}c  \mathcal{D}\bar c  \exp\left( i (S_{\text{YM}} + S_{\text{gauge-fix}} + S_{\text{ghost}}) \right).
$$

### 12.3 BRST Symmetry

The gauge-fixed action has a residual fermionic symmetry called **BRST symmetry**
(Becchi, Rouet, Stora, Tyutin 1975), generated by the **BRST charge** $Q$. The BRST
transformations are

$$
\begin{aligned}
\delta_B A_\mu^a &= D_\mu^{ab} c^b  \varepsilon && (\varepsilon: \text{Grassmann parameter}), \\
\delta_B c^a &= -\frac{1}{2} f^{abc} c^b c^c  \varepsilon, \\
\delta_B \bar c^a &= b^a  \varepsilon, \\
\delta_B b^a &= 0.
\end{aligned}
$$

Key property: BRST is **nilpotent**: $Q^2 = 0$ (i.e., $\delta_B^2 = 0$ on all
fields).

**Physical states** are those annihilated by $Q$: $Q |\text{phys}\rangle = 0$ (modulo
$Q$-exact states). This is the cohomological characterization of the physical
Hilbert space:

$$
\mathcal{H}_{\text{phys}} = \ker(Q) / \mathrm{im}(Q) = H^0(Q).
$$

The BRST formalism is essential for:

- Proving unitarity and renormalizability of Yang-Mills theories
  (Becchi, Rouet, Stora 1975).
- Defining the path integral rigorously.
- String theory (the worldsheet has a BRST structure).

---

## 13. Instantons and Topological Charge

### 13.1 What Is an Instanton?

An **instanton** is a solution to the Euclidean Yang-Mills equations with finite
action that is localized in (Euclidean) spacetime. It is a tunneling amplitude
between different topological sectors of the gauge field.

For $SU(2)$ Yang-Mills in Euclidean $\mathbb{R}^4$ (compactified to $S^4$):

$$
\begin{aligned}
\text{Self-dual instanton:} \quad & F_{\mu\nu} = +\ast F_{\mu\nu} && \text{(BPS equation)}, \\
\text{Anti-self-dual:} \quad & F_{\mu\nu} = -\ast F_{\mu\nu},
\end{aligned}
$$

where $\ast F_{\mu\nu} = \epsilon_{\mu\nu\rho\sigma} F^{\rho\sigma} / 2$ is the
Hodge dual.

The BPST instanton (Belavin, Polyakov, Schwarz, Tyupkin 1975) with instanton
number $k = 1$:

$$
A_\mu^a = \frac{2}{g}  \eta_{\mu\nu}^a  \frac{x_\nu}{x^2 + \rho^2},
$$

where $\eta_{\mu\nu}^a$ are the 't Hooft symbols and $\rho$ is the instanton size.

### 13.2 The Instanton Number

The instanton is classified by its **topological charge** (second Chern number):

$$
k = c_2 = \frac{1}{8 \pi^2} \int \mathrm{tr}(F_{\mu\nu}  \ast F^{\mu\nu})  d^4x \in \mathbb{Z}.
$$

For the BPST instanton: $k = 1$ (one unit of topological charge).

The Yang-Mills action satisfies the bound

$$
S_{\text{YM}} \geq \frac{8 \pi^2 |k|}{g^2} \qquad \text{(Bogomolny bound)},
$$

with equality for (anti-)self-dual configurations. This is why instantons have
minimum action for their topological class.

### 13.3 Physical Consequences

**QCD vacuum.** The QCD vacuum is a superposition of all topological sectors:

$$
|\theta\rangle = \sum_n e^{i n \theta} |n\rangle,
$$

where $|n\rangle$ is the state with $n$ units of topological charge. The parameter
$\theta$ is the **QCD $\theta$-term** — it would be observable (CP violation) but is
measured to be $|\theta| < 10^{-10}$ (the strong CP problem).

**Axial anomaly.** Instantons are responsible for the axial $U(1)$ anomaly in QCD
— the would-be ninth Goldstone boson of chiral symmetry breaking (the $\eta'$) gets
a mass from instanton effects.

---

## 14. The Standard Model as a Gauge Theory

### 14.1 The Gauge Group

The Standard Model gauge group is

$$
G_{\text{SM}} = U(1)_Y \times SU(2)_L \times SU(3)_c.
$$

- **$U(1)_Y$:** hypercharge (1 generator, 1 gauge boson $B_\mu$).
- **$SU(2)_L$:** weak isospin (3 generators, 3 gauge bosons $W^a_\mu$).
- **$SU(3)_c$:** color (8 generators, 8 gluons $G^a_\mu$).

Total: 12 gauge bosons.

After electroweak symmetry breaking (Higgs mechanism):

- $W^a_\mu$ and $B_\mu$ mix to give $W^{\pm}$, $Z$, and $\gamma$.
- The photon $\gamma$ remains massless (gauge boson of unbroken $U(1)_{\text{EM}}$).
- $W^{\pm}$ and $Z$ get masses via the Higgs mechanism.

### 14.2 Matter Content and Representations

$$
\begin{aligned}
\text{Quarks:} \quad & (3, 2, 1/6)_L + (\bar 3, 1, -2/3)_R + (\bar 3, 1, 1/3)_R, \\
\text{Leptons:} \quad & (1, 2, -1/2)_L + (1, 1, -1)_R, \\
\text{Higgs:} \quad & (1, 2, 1/2).
\end{aligned}
$$

Format: ($SU(3)_c$ rep, $SU(2)_L$ rep, $U(1)_Y$ charge).

### 14.3 The Full Standard Model Lagrangian

$$
\begin{aligned}
\mathcal{L}_{\text{SM}} &= \mathcal{L}_{\text{gauge}} + \mathcal{L}_{\text{fermion}} + \mathcal{L}_{\text{Higgs}} + \mathcal{L}_{\text{Yukawa}}, \\[4pt]
\mathcal{L}_{\text{gauge}} &= -\frac{1}{4} B_{\mu\nu} B^{\mu\nu}
  - \frac{1}{4} W^a_{\mu\nu} W^{a\mu\nu}
  - \frac{1}{4} G^a_{\mu\nu} G^{a\mu\nu}, \\[4pt]
\mathcal{L}_{\text{fermion}} &= \sum_{\text{generations}} [
  \bar q_L  i\gamma^\mu D_\mu  q_L
  + \bar u_R  i\gamma^\mu D_\mu  u_R
  + \bar d_R  i\gamma^\mu D_\mu  d_R \\
  &\qquad\qquad + \bar\ell_L  i\gamma^\mu D_\mu  \ell_L
  + \bar e_R  i\gamma^\mu D_\mu  e_R ], \\[4pt]
\mathcal{L}_{\text{Higgs}} &= |D_\mu H|^2 - \lambda (H^\dagger H - v^2 / 2)^2, \\[4pt]
\mathcal{L}_{\text{Yukawa}} &= -y_u  \bar q_L \tilde H u_R
  - y_d  \bar q_L H d_R
  - y_e  \bar\ell_L H e_R + \text{h.c.}
\end{aligned}
$$

Everything is determined by:

1. The gauge group $G_{\text{SM}}$.
2. The matter representations.
3. The requirement of renormalizability.

The Standard Model has $\sim 19$ free parameters (masses, mixing angles, coupling
constants) but its structure is completely fixed by gauge invariance.

---

## 15. Gauge Theory in Condensed Matter: Berry Phase

### 15.1 The Adiabatic Connection

In condensed matter physics, gauge theory appears through the **Berry phase** — a
geometric phase acquired by quantum states under adiabatic evolution.

Consider a quantum system with Hamiltonian $H(R)$ depending on parameters
$R \in M$ (parameter space). For a normalized eigenstate $|n(R)\rangle$:

$$
H(R) |n(R)\rangle = E_n(R) |n(R)\rangle.
$$

The **Berry connection** (gauge potential on parameter space $M$):

$$
A^n_\mu(R) = i \langle n(R)| \frac{\partial}{\partial R^\mu} |n(R)\rangle \in i\mathbb{R} \qquad \text{(imaginary, for } U(1) \text{ bundle).}
$$

The **Berry curvature** (field strength):

$$
F^n_{\mu\nu}(R) = \partial_\mu A^n_\nu - \partial_\nu A^n_\mu
  = i \left( \langle \partial_\mu n | \partial_\nu n \rangle - \langle \partial_\nu n | \partial_\mu n \rangle \right).
$$

### 15.2 The Berry Phase as Holonomy

For a closed loop $C$ in parameter space,

$$
\gamma_n(C) = \oint_C A^n_\mu  dR^\mu = i \oint_C \langle n(R)| d |n(R)\rangle.
$$

This is the **Berry phase** — a gauge-theoretic holonomy in the space of quantum
states. It is observable (via interference) and non-trivial when $F \neq 0$.

### 15.3 The TKNN Invariant and Quantum Hall Effect

The **integer quantum Hall effect** is characterized by the quantization of the
Hall conductance $\sigma_{xy} = (e^2 / h) \times n$, where $n$ is an integer. This
integer is the **Chern number** (TKNN invariant, Thouless, Kohmoto, Nightingale,
den Nijs 1982):

$$
n = \frac{1}{2\pi} \int_{\text{BZ}} F^n_{k_x k_y}  dk_x  dk_y \in \mathbb{Z},
$$

integrated over the Brillouin zone (a torus $T^2$). The quantization is exact
because it is a topological invariant — it cannot change continuously.

This is the condensed matter realization of the first Chern number, and it is a
profound example of gauge theory governing macroscopic observable phenomena.

---

## 16. Gauge Theory and Transformers: The Connection

### 16.1 The Hidden State Space as a Base Manifold

In a transformer, the hidden state space $\mathbb{R}^d$ plays the role of the base
manifold $B$. At each hidden state $h \in \mathbb{R}^d$, the force exerted on the
hidden state by the attention mechanism defines a connection.

The gauge group structure enters through the multi-head structure:

$$
\begin{aligned}
\text{Head } h \text{ acts via:} \quad & A^{(h)}_\mu(x) = \mathrm{softmax}_\mu\left( \frac{Q_h x \cdot K_h x_j}{\sqrt{d_h}} \right) \cdot V_h x_j, \\[4pt]
\text{Each head:} \quad & A^{(h)} \in \mathrm{End}(\mathbb{R}^{d_h}) \cong \mathfrak{u}(d_h) \quad \text{(Lie algebra of } U(d_h)\text{)}, \\[4pt]
\text{After } W_O \text{ mixing:} \quad & A = W_O ( A^{(1)} \oplus \cdots \oplus A^{(H)} ) W_O^\dagger \in \mathrm{End}(\mathbb{R}^d) \cong \mathfrak{u}(d).
\end{aligned}
$$

### 16.2 Why Multi-Head Creates Non-Abelian Structure

**Single head ($H = 1$).** Before and after projection,

$$
A_\mu(x) \in \mathfrak{u}(d_h) \longrightarrow A_\mu(x) \in \mathfrak{u}(d) \text{ (but rank } d_h\text{)},
$$

and the curvature is

$$
F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu].
$$

- If $A_\mu = (\partial_\mu V) \cdot I$ (scalar times identity), then
  $[A_\mu, A_\nu] = 0$ — abelian, potentially conservative.
- If $A_\mu \neq (\partial_\mu V) \cdot I$ but rank-1, $[A_\mu, A_\nu]$ depends on
  the specific matrices.

**Multi-head ($H > 1$).** Summing over heads,

$$
\begin{aligned}
\text{Total force:} \quad & F = \sum_h f^{(h)}, \\
\text{Connection:} \quad & A = \sum_h A^{(h)}, \\
\text{Curvature:} \quad & F_{\mu\nu} = \sum_h F_{\mu\nu}^{(h)} + \sum_{h \neq h'} [A_\mu^{(h)}, A_\nu^{(h')}].
\end{aligned}
$$

The second sum is the **cross-head commutator** term. It is generically non-zero
because:

- $A^{(h)}$ acts in the $(Q_h, K_h)$ subspace.
- $A^{(h')}$ acts in the $(Q_{h'}, K_{h'})$ subspace.
- The $W_O$ mixing makes these subspaces overlap.

Therefore $[A^{(h)}, A^{(h')}] \neq 0$ generically.

**Why $[A^{(h)}, A^{(h')}] \neq 0$: the explicit argument.**

Before $W_O$ mixing, each head $h$ produces a matrix $A^{(h)}$ that acts only
within its own $d_h$-dimensional subspace $S_h \subset \mathbb{R}^d$. Two matrices
acting in **orthogonal** subspaces always commute:

$$
\text{If } S_h \perp S_{h'}: \quad A^{(h)} A^{(h')} = A^{(h')} A^{(h)} = 0 \Longrightarrow [A^{(h)}, A^{(h')}] = 0.
$$

But the output projection $W_O \in \mathbb{R}^{d \times H d_h}$ **recombines all
heads into the shared $d$-dimensional space**. After this mixing, $A^{(h)}$ and
$A^{(h')}$ both act on all of $\mathbb{R}^d$, but in directions determined by the
columns of $W_O$. Unless those directions happen to be orthogonal — which requires
special structure in $W_O$ that trained transformers do not have — the subspaces
overlap and

$$
[A^{(h)}, A^{(h')}] = A^{(h)} A^{(h')} - A^{(h')} A^{(h)} \neq 0 \qquad \text{(generically).}
$$

The commutator is the matrix that measures how much the two heads "interfere"
when their forces are composed in different orders. When this is non-zero, the
gauge group of the full attention mechanism is **non-abelian** (§1.4 of the Lie
Groups Tutorial), and non-abelian means non-zero curvature, which means no scalar
potential (§7.4).

### 16.3 The Obstruction Theorem for Multi-Head Attention

**Theorem.** For a transformer with $H \geq 2$ attention heads and generic weight
matrices $W_Q^h$, $W_K^h$, $W_V^h$, $W_O$, the attention force field

$$
F(h) = \sum_h \sum_j \mathrm{softmax}\left( \frac{Q_h h \cdot K_h h_j}{\sqrt{d_h}} \right) \cdot V_h h_j
$$

is NOT the gradient of any scalar potential $V: \mathbb{R}^d \to \mathbb{R}$.

**Proof sketch.** The force Jacobian

$$
\left( \frac{\partial F}{\partial h} \right)_{ij}
  = \sum_h \sum_j \left[ \mathrm{softmax} \cdot \frac{\partial (V_h h_j)}{\partial h_i}
  + \left( \frac{\partial  \mathrm{softmax}}{\partial h_i} \right) \cdot V_h h_j \right].
$$

The antisymmetric part of the Jacobian

$$
\Omega_{ij} = \frac{1}{2}\left[ \frac{\partial F}{\partial h} - \left(\frac{\partial F}{\partial h}\right)^{\top} \right]_{ij}
$$

contains terms from the cross-head interactions proportional to
$[A^{(h)}, A^{(h')}]_{ij}$. For generic $W$ matrices, these are non-zero. A non-zero
antisymmetric Jacobian means $F \neq -\nabla V$ for any $V$. $\blacksquare$

### 16.4 The Van Nierop Connection

Van Nierop (2024) showed that transformers are invariant under

$$
\text{For each head } h: \quad K_h \to U_h K_h, \quad Q_h \to U_h Q_h \quad \text{for } U_h \in GL(d_h).
$$

This is the **gauge invariance of the parameterization** — different weight
matrices can produce the same function. The gauge group of parameterization
redundancy is

$$
G_{\text{param}} = GL(d_1) \times GL(d_2) \times \cdots \times GL(d_H).
$$

Your paper establishes the **complementary** result: the **dynamical** gauge
structure of the force field on hidden states has curvature given by cross-head
commutators.

These are dual aspects:

- Van Nierop: redundancy in weight space (passive gauge).
- Your paper: curvature of hidden-state force field (active gauge / dynamics).

### 16.5 Statement Decoder

The statement

> *"cross-head commutators $[A^{(h)}, A^{(h')}]$ generate non-abelian curvature
> that obstructs any scalar potential on hidden-state space, regardless of
> capacity"*

is assembled from concepts developed across this tutorial. Here is the complete
term-by-term reading guide, with every word traced to its definition.

**"cross-head"** → Two different attention heads $h \neq h'$ (§16.1). Before $W_O$
mixing, each head acts in its own $d_h$-dimensional subspace. After $W_O$ mixing,
they share the full $d$-dimensional space and can interfere.

**"commutators $[A^{(h)}, A^{(h')}]$"** → The matrix commutator (§4.1a of the
Lie Groups Tutorial) of the two connection 1-forms:

$$
[A^{(h)}, A^{(h')}] = A^{(h)} A^{(h')} - A^{(h')} A^{(h)} \in \mathrm{End}(\mathbb{R}^d).
$$

$A^{(h)}$ is the gauge potential contributed by head $h$ after $W_O$ projection
(§16.1). This commutator is zero if and only if the two matrices commute — which
they generically do not after $W_O$ mixes the head subspaces (§16.2).

**"generate"** → Appear explicitly in the curvature formula (§7.1):

$$
F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu + \underbrace{[A_\mu, A_\nu]}_{=\sum_{h \neq h'} [A_\mu^{(h)}, A_\nu^{(h')}]}.
$$

When the cross-head commutators are non-zero, this term is non-zero, making
$F_{\mu\nu} \neq 0$.

**"non-abelian"** → Coming from a non-abelian gauge group (§1.4 of Lie Groups
Tutorial, §8 of this tutorial). "Non-abelian curvature" specifically means the
curvature arising from the $[A, A]$ commutator term — the term that is absent in
electromagnetism ($U(1)$ is abelian, $[A, A] = 0$) and present in Yang-Mills
theories ($SU(N)$ is non-abelian, $[A, A] \neq 0$).

**"curvature"** → The field strength
$F_{\mu\nu} = \partial A - \partial A + [A, A]$, measuring the failure of parallel
transport to commute around infinitesimal loops (§7.1). Non-zero curvature = the
connection is not flat = physics depends on the path.

**"obstructs"** → Prevents existence, in the mathematical sense of making it
impossible rather than merely difficult. Specifically: $F_{\mu\nu} \neq 0$ is an
algebraic obstruction to writing $F = -\nabla V$ (§7.4, Gauge Obstruction Theorem).

**"any scalar potential"** → Any smooth function $V: \mathbb{R}^d \to \mathbb{R}$,
however complex, however many parameters, however deeply nonlinear. "Any" is
unrestricted.

**"on hidden-state space"** → On the base manifold $\mathbb{R}^d$ where hidden
states live (§16.1). The connection $A$ and curvature $F$ are defined on this
space.

**"regardless of capacity"** → The obstruction is Clairaut's theorem (§7.4):
every smooth $V$ satisfies $\partial_i \partial_j V = \partial_j \partial_i V$, so
$\mathrm{Curvature}(-\nabla V) = 0$ exactly. Since $F_{\mu\nu} \neq 0$ and
$\mathrm{Curvature}(-\nabla V) = 0$, we have $F \neq -\nabla V$ for any $V$
whatsoever. This is an identity of calculus, not a statement about approximation
quality. Increasing $V_\psi$ capacity cannot change a non-zero number into zero.

**The four-step logical chain:**

1. $W_O$ mixing $\Longrightarrow$ head subspaces overlap.
2. Overlapping subspaces $\Longrightarrow$ $[A^{(h)}, A^{(h')}] \neq 0$ (matrices
   don't commute).
3. $[A, A] \neq 0 \Longrightarrow F_{\mu\nu} \neq 0$ (non-zero curvature).
4. $F_{\mu\nu} \neq 0 \Longrightarrow F \neq -\nabla V$ for any $V$ (Clairaut
   obstruction).

Every step is an algebraic identity or a theorem of differential geometry. None
is an approximation. The experimental shared-$V_\psi$ failure
($R^2 = 0.04$–$0.20$ for GPT-2 middle layers) measures the quantitative signature
of this structural impossibility.

---

## 17. Summary: The Unified Picture

### 17.1 The Core Ideas

1. **Gauge principle.** Requiring local symmetry $\Longrightarrow$ forces are
   mandatory (not optional). The form of all interactions is determined by the
   gauge group.

2. **Fiber bundle geometry.** Physical states are sections of associated bundles;
   forces are connections on principal bundles; field strengths are curvatures.

3. **Abelian vs. non-abelian.**
   - Abelian $U(1)$: photon, no self-interaction, linear Maxwell equations.
   - Non-abelian $SU(N)$: gauge bosons self-interact, non-linear field equations.
   - Key difference: $[A_\mu, A_\nu] \neq 0$.

4. **Topology matters.** Flat connections (zero curvature) are locally pure gauge,
   globally may have holonomy. Non-trivial topology gives Chern classes,
   instantons, Berry phase. Physical consequences: confinement, quantum Hall
   effect, Aharonov-Bohm effect.

5. **Obstruction theorem.** Non-zero curvature $\iff$ no scalar potential exists.
   Multi-head attention generates non-zero curvature $\Longrightarrow$ no scalar
   $V(h)$ can represent multi-head attention dynamics.

### 17.2 The Logical Chain From Gauge Principle to Physics

$$
\begin{aligned}
& \text{Global symmetry } G && \text{(e.g., } U(1) \text{ phase invariance)} \\
& \quad \downarrow \text{ promote to local} \\
& \text{Local symmetry } G(x) && \text{(different group element at each point)} \\
& \quad \downarrow \text{ requires} \\
& \text{Connection } A_\mu(x) && \text{(gauge potential / force carrier)} \\
& \quad \downarrow \text{ has} \\
& \text{Curvature } F_{\mu\nu} = \partial A - \partial A + [A, A] && \text{(field strength)} \\
& \quad \downarrow \text{ satisfies} \\
& \text{Yang-Mills equations } D_\nu F^{\nu\mu} = j^\mu && \text{(equations of motion)} \\
& \quad \downarrow \text{ classified by} \\
& \text{Topological invariants (Chern classes)} && \text{(winding numbers, instanton charge)} \\
& \quad \downarrow \text{ observed as} \\
& \text{Physical forces, Berry phases,} && \text{(electromagnetism, weak force,} \\
& \text{holonomy effects, confinement} && \text{ strong force, quantum Hall, \ldots)}
\end{aligned}
$$

### 17.3 Key Formulas Reference

$$
\begin{aligned}
\text{Covariant derivative:} \quad & D_\mu = \partial_\mu + A_\mu \quad (A_\mu \in \mathfrak{g}), \\
\text{Field strength:} \quad & F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu], \\
\text{Gauge transformation:} \quad & A_\mu \to g A_\mu g^{-1} + g \partial_\mu g^{-1} \quad (g: B \to G), \\
& F_{\mu\nu} \to g F_{\mu\nu} g^{-1}, \\
\text{Matter coupling:} \quad & D_\mu \psi = \partial_\mu \psi + \rho(A_\mu) \psi \quad (\rho: G \to GL(V)), \\
\text{Yang-Mills action:} \quad & S = -\frac{1}{4 g^2} \int F^a_{\mu\nu} F^{a\mu\nu}  d^4x, \\
\text{YM equations of motion:} \quad & D_\nu F^{\nu\mu} = j^\mu, \\
\text{Bianchi identity:} \quad & D_{[\mu} F_{\nu\rho]} = 0, \\
\text{Wilson loop:} \quad & W(C) = \mathrm{tr} \mathcal{P} \exp\left( \oint_C A_\mu  dx^\mu \right), \\
\text{Chern number:} \quad & c_2 = \frac{1}{8 \pi^2} \int \mathrm{tr}(F \wedge F) \in \mathbb{Z}, \\
\text{Holonomy:} \quad & \mathrm{Hol}_\gamma = \mathcal{P} \exp\left( \int_\gamma A \right) \in G, \\
\text{Berry phase:} \quad & \gamma = \oint A^n_\mu  dR^\mu \in \mathbb{R} / (2\pi \mathbb{Z}).
\end{aligned}
$$

### 17.4 Further Reading

**Textbooks:**

- **Nakahara, "Geometry, Topology and Physics"** — the best mathematical reference,
  covers fiber bundles, connections, characteristic classes, and instantons. Highly
  recommended given your Riemannian geometry background.

- **Bleecker, "Gauge Theory and Variational Principles"** — rigorous mathematical
  treatment using the calculus of variations you already know.

- **Hamilton, "Mathematical Gauge Theory"** — modern, self-contained, uses principal
  bundles from the start.

- **Peskin & Schroeder, "Introduction to Quantum Field Theory"** — standard physics
  reference for Yang-Mills, Faddeev-Popov, BRST.

- **Bertlmann, "Anomalies in Quantum Field Theory"** — for the topological aspects
  (instantons, anomalies, Chern-Simons theory).

**Papers:**

- Yang & Mills (1954): "Conservation of Isotopic Spin and Isotopic Gauge Invariance"
- Atiyah, Hitchin, Singer (1978): "Self-duality in four-dimensional Riemannian geometry"
- BPST (1975): "Pseudoparticle solutions of the Yang-Mills equations"
- Berry (1984): "Quantal Phase Factors Accompanying Adiabatic Changes"
- Thouless et al. (1982): "Quantized Hall Conductance in a Two-Dimensional Periodic Potential"
- Van Nierop (2024): "Transformer models are gauge invariant" (arXiv:2412.14543)

---

*Tutorial version: April 2026.*
*Pitched at readers with group theory, calculus of variations, and basic Riemannian geometry.*
