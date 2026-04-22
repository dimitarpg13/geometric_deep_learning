# Lie Groups for Gauge Theory
### A Graduate Tutorial from Matrix Groups to the Structure of U(1), SU(2), SU(3)

**Prerequisites assumed:** Linear algebra (eigenvalues, matrix exponential, determinants),
multivariable calculus, basic topology (open sets, continuity, compactness),
complex numbers and complex analysis basics.

**Goal:** By the end of this tutorial you will understand $U(1)$ as a Lie group from
every angle — as a manifold, as a matrix group, as an abstract group, as a topological
space, and as the gauge group of electromagnetism. You will understand how $U(1)$ relates
to $SU(2)$ and $SU(3)$ and why the structure of these groups completely determines the
form of the fundamental forces.

---

## Table of Contents

1. [What Is a Lie Group? The Three Faces](#1-what-is-a-lie-group)
   - [1.4 Abelian and Non-Abelian Lie Groups](#14-abelian-and-non-abelian-lie-groups)
2. [The Simplest Examples: $\mathbb{R}$, $S^1$, and $U(1)$](#2-simplest-examples)
3. [Matrix Lie Groups: GL, SL, O, U, S](#3-matrix-lie-groups)
4. [The Lie Algebra: Tangent Space at the Identity](#4-lie-algebra)
   - [4.0 Unpacking "Tangent Space at the Identity"](#40-unpacking-tangent-space-at-the-identity)
   - [4.1a The Lie Bracket in Detail: Bilinear, Antisymmetric, Jacobi](#41a-the-lie-bracket-in-detail)
5. [The Exponential Map](#5-exponential-map)
6. [The Adjoint Representation](#6-adjoint-representation)
7. [$U(1)$ in Full Detail](#7-u1-full-detail)
8. [$SU(2)$: The Three-Sphere and Rotations](#8-su2)
9. [$SO(3)$ vs $SU(2)$: Covering Groups and Spinors](#9-so3-vs-su2)
10. [$SU(3)$ and the Gell-Mann Matrices](#10-su3)
11. [Representations of Lie Groups](#11-representations)
12. [The Peter-Weyl Theorem and Harmonic Analysis](#12-peter-weyl)
13. [Roots, Weights, and the Classification of Simple Lie Algebras](#13-classification)
14. [Homomorphisms, Subgroups, and Quotients](#14-homomorphisms)
15. [The Topology of Lie Groups](#15-topology)
16. [Lie Groups in Gauge Theory: The Full Picture](#16-gauge-theory-connection)
17. [Reference: Key Facts About $U(1)$, $SU(2)$, $SU(3)$](#17-reference)

---

## 1. What Is a Lie Group? The Three Faces

### 1.1 The Three Structures

A **Lie group** is an object with three simultaneous structures that are mutually compatible:

**Face 1: A GROUP**
- Set $G$ with binary operation $\cdot : G \times G \to G$
- Associativity: $(ab)c = a(bc)$
- Identity element: $e \in G$ with $ge = eg = g$
- Inverses: for each $g \in G$, $\exists g^{-1}$ with $g g^{-1} = e$

**Face 2: A SMOOTH MANIFOLD**
- $G$ is a smooth manifold of dimension $n$
- Has an atlas of coordinate charts
- Transition maps are smooth ($C^{\infty}$)

**Face 3: COMPATIBILITY**
- The group operations are smooth maps:
  - $\mu: G \times G \to G$, $(g, h) \mapsto gh$ (smooth)
  - $\iota: G \to G$, $g \mapsto g^{-1}$ (smooth)

The compatibility condition is what makes Lie groups so powerful: the algebraic
structure (group) and the geometric structure (manifold) are woven together.
Every theorem about one structure applies to the other.

### 1.2 Why This Structure Matters for Physics

In physics, we use Lie groups as symmetry groups. The manifold structure means:
- We can talk about **continuous** symmetries (infinitesimal transformations)
- We can differentiate transformations $\to$ Lie algebra
- We can integrate $\to$ action of the group on fields

The group structure means:
- Compositions of symmetries are symmetries
- Every symmetry has an inverse
- The identity transformation is included

Without the smooth structure, we would only have discrete symmetries (like crystal
symmetries). The smoothness is what gives us conservation laws via Noether's theorem.

### 1.3 The Lie Group–Lie Algebra Correspondence

The most important theorem about Lie groups:

Every Lie group $G$ has an associated Lie algebra $\mathfrak{g} = T_e G$
(the tangent space at the identity).

This correspondence is:
- **Local:** $\mathfrak{g}$ determines $G$ near the identity
- **Global (for simply connected $G$):** the exponential map $\exp: \mathfrak{g} \to G$
  gives a local diffeomorphism near $0 \in \mathfrak{g}$

Simply connected Lie groups $\leftrightarrow$ Lie algebras (1:1) — **Lie's Third Theorem**.

The Lie algebra is a vector space with a bilinear antisymmetric bracket $[\cdot, \cdot]$
satisfying the Jacobi identity — each of these terms is defined precisely in
§4.1a. It is often much easier to work with (it is linear!) and then
"exponentiate" back to the group.

### 1.4 Abelian and Non-Abelian Lie Groups

These two terms appear throughout gauge theory and are worth defining precisely
at the outset, because the entire difference between electromagnetism (simple,
linear) and the strong/weak forces (complex, nonlinear) reduces to them.

**Definition (Abelian Lie group).** A Lie group $G$ is **abelian** (or commutative)
if group multiplication commutes:

$$
gh = hg \quad \text{for all } g, h \in G
$$

**Definition (Non-abelian Lie group).** A Lie group $G$ is **non-abelian**
(or non-commutative) if there EXIST $g, h \in G$ such that:

$$
gh \neq hg \quad \text{(at least one pair of elements that do not commute)}
$$

The word "non-abelian" means the group does not have the abelian property —
it is not required to fail commutativity everywhere, only somewhere.

**The Lie algebra test.** Because the Lie bracket $[X, Y]$ measures the
infinitesimal failure of group elements to commute (see §4.1a), a connected
Lie group is abelian if and only if its Lie algebra satisfies:

$$
G \text{ abelian} \iff [X, Y] = 0 \text{ for all } X, Y \in \mathfrak{g}
$$

This means: checking commutativity of the group reduces to checking whether
all brackets vanish in the Lie algebra — a purely linear algebra problem.

**Standard examples:**

| Group | Abelian? | Why |
|---|---|---|
| $U(1) = \{e^{i\theta}\}$ | **Yes** | $e^{i\alpha} e^{i\beta} = e^{i(\alpha+\beta)} = e^{i\beta} e^{i\alpha}$ |
| $(\mathbb{R}, +)$ | **Yes** | $x + y = y + x$ |
| $SU(2)$ | **No** | $\sigma_1 \sigma_2 = i\sigma_3 \neq -i\sigma_3 = \sigma_2 \sigma_1$ (Pauli matrices don't commute) |
| $SU(3)$ | **No** | Most generators don't commute |
| $GL(n)$, $n \geq 2$ | **No** | Matrix multiplication is not commutative |

**Why this matters for gauge theory:**

$G = U(1)$ (abelian):

$$
\begin{aligned}
[A_\mu, A_\nu] &= 0 \quad \text{for all } \mu, \nu \\
F_{\mu\nu} &= \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu]
           = \partial_\mu A_\nu - \partial_\nu A_\mu \quad \text{(linear in } A\text{)}
\end{aligned}
$$

- Photon does not self-interact
- Maxwell equations are linear
- Scalar potential CAN exist (and does: the electromagnetic potential)

$G = SU(2)$ or $SU(3)$ (non-abelian):

$$
\begin{aligned}
[A_\mu, A_\nu] &\neq 0 \quad \text{generically} \\
F_{\mu\nu} &= \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu] \quad \text{(nonlinear in } A\text{)}
\end{aligned}
$$

- Gauge bosons self-interact ($W$ bosons carry weak charge, gluons carry color)
- Yang-Mills equations are nonlinear
- No scalar potential can exist that reproduces $F_{\mu\nu}$

The last bullet is the **non-abelian gauge obstruction**: as soon as $[A_\mu, A_\nu] \neq 0$,
the curvature $F_{\mu\nu} \neq 0$, and by the obstruction theorem (§7.4 in the Gauge Theory
Tutorial), no scalar potential $V$ — of any form, any complexity — can satisfy
$F = -\nabla V$. This is an algebraic identity, not an approximation limit.

---

## 2. The Simplest Examples: $\mathbb{R}$, $S^1$, and $U(1)$

### 2.1 The Real Line $(\mathbb{R}, +)$

The real line with addition is the simplest Lie group:

- **Manifold:** $\mathbb{R}^1$ (a 1-dimensional smooth manifold, the real line)
- **Group op:** $(x, y) \mapsto x + y$
- **Identity:** $0$
- **Inverse:** $x \mapsto -x$
- **Lie algebra:** $\mathbb{R}$ (with trivial bracket $[x, y] = 0$)

The exponential map $\exp: \mathbb{R} \to \mathbb{R}$ sends the Lie algebra to itself.
(For additive groups, $\exp$ is literally the identity map on $\mathbb{R}$ when viewed abstractly,
though numerically $\exp(t) = e^t$ takes $\mathbb{R}$ to $\mathbb{R}_+$ as a map of manifolds.)

### 2.2 The Circle $S^1$ as a Lie Group

The unit circle in $\mathbb{R}^2$:

$$
S^1 = \{(x, y) \in \mathbb{R}^2 : x^2 + y^2 = 1\}
$$

with multiplication defined by angle addition:

$$
(x_1, y_1) \cdot (x_2, y_2) = (x_1 x_2 - y_1 y_2,\ x_1 y_2 + y_1 x_2)
$$

(this is just complex multiplication of unit complex numbers.)

Written in terms of angles $\theta \in [0, 2\pi)$:

$$
\theta_1 \cdot \theta_2 = \theta_1 + \theta_2 \pmod{2\pi}
$$

So $S^1$ with this multiplication is the group $\mathbb{R}/2\pi\mathbb{Z}$ (the real line modulo
the integers scaled by $2\pi$).

### 2.3 $U(1)$: The Complex Definition

The **unitary group $U(1)$** is the group of $1 \times 1$ unitary matrices — complex numbers
of modulus 1:

$$
U(1) = \{z \in \mathbb{C} : |z| = 1\} = \{e^{i\theta} : \theta \in \mathbb{R}\}
$$

with multiplication = complex multiplication.

**Claim: $U(1) \cong S^1$ as Lie groups.**

The isomorphism is: $e^{i\theta} \leftrightarrow (\cos\theta, \sin\theta)$. Under this identification:
- The group operation (complex multiplication) corresponds to angle addition
- The manifold structure (circle) is the same
- The smooth structure is the same

So $U(1)$, $S^1$, and $\mathbb{R}/2\pi\mathbb{Z}$ are three names for the same Lie group.

### 2.4 Why $U(1)$ Is Called a "Circle Group"

Topologically, $U(1) \cong S^1$. This has profound consequences:

$$
\begin{aligned}
\pi_0(U(1)) &= 0 \quad \text{(connected: one piece)} \\
\pi_1(U(1)) &= \mathbb{Z} \quad \text{(fundamental group: loops wind around the circle by integer amounts)} \\
\pi_2(U(1)) &= 0 \quad \text{(no non-trivial 2-spheres)}
\end{aligned}
$$

The fact that $\pi_1(U(1)) = \mathbb{Z}$ means there are topologically distinct ways to map a
loop into $U(1)$ — classified by the winding number. This is the origin of
**magnetic monopoles** in physics (magnetic charge is quantized because it is
measured by this $\mathbb{Z}$).

---

## 3. Matrix Lie Groups: GL, SL, O, U, S

### 3.1 The General Linear Group $GL(n, \mathbb{F})$

The group of invertible $n \times n$ matrices over a field $\mathbb{F}$ ($\mathbb{R}$ or $\mathbb{C}$):

$$
\begin{aligned}
GL(n, \mathbb{R}) &= \{A \in \mathrm{Mat}(n \times n, \mathbb{R}) : \det A \neq 0\} \\
GL(n, \mathbb{C}) &= \{A \in \mathrm{Mat}(n \times n, \mathbb{C}) : \det A \neq 0\}
\end{aligned}
$$

- Dimension (as a manifold): $n^2$ over $\mathbb{R}$, $2n^2$ over $\mathbb{R}$ (= $n^2$ over $\mathbb{C}$)
- $GL(n)$ is an open subset of $\mathrm{Mat}(n \times n)$ — the condition $\det A \neq 0$ is open
- Not compact (matrices can have arbitrarily large or small entries)
- Not connected: $GL(n, \mathbb{R})$ has two connected components ($\det > 0$ and $\det < 0$)

### 3.2 The Special Linear Group $SL(n)$

$$
\begin{aligned}
SL(n, \mathbb{R}) &= \{A \in GL(n, \mathbb{R}) : \det A = 1\} \\
SL(n, \mathbb{C}) &= \{A \in GL(n, \mathbb{C}) : \det A = 1\}
\end{aligned}
$$

- Dimension: $n^2 - 1$ (one constraint $\det = 1$ reduces dimension by 1)
- $SL(n)$ is a closed subgroup of $GL(n)$ ($\det = 1$ is a closed condition)
- $SL(n)$ is connected

### 3.3 The Orthogonal and Special Orthogonal Groups

$$
\begin{aligned}
O(n) &= \{A \in GL(n, \mathbb{R}) : A^T A = I\} \quad \text{(orthogonal matrices)} \\
SO(n) &= \{A \in O(n) : \det A = 1\} \quad \text{(special orthogonal = rotation matrices)}
\end{aligned}
$$

- $O(n)$ has two connected components: $\det = +1$ (rotations) and $\det = -1$ (reflections)
- $SO(n) = O(n)$ component containing identity = proper rotations
- Dimension: $n(n-1)/2$
- Compact: entries bounded by $|A_{ij}| \leq 1$ from $A^T A = I$

Physical examples:

$$
\begin{aligned}
SO(2) &\cong U(1) = \text{rotations in 2D} \quad (\dim = 1) \\
SO(3) &= \text{rotations in 3D} \quad (\dim = 3) \\
SO(3, 1) &= \text{Lorentz group} \quad (\dim = 6, \text{ but non-compact due to boosts})
\end{aligned}
$$

### 3.4 The Unitary and Special Unitary Groups

$$
\begin{aligned}
U(n) &= \{A \in GL(n, \mathbb{C}) : A^{\dagger} A = I\} \quad \text{(unitary matrices, } A^{\dagger} = \bar{A}^T\text{)} \\
SU(n) &= \{A \in U(n) : \det A = 1\} \quad \text{(special unitary)}
\end{aligned}
$$

- Dimension: $n^2$ (as a real manifold; $n^2$ real constraints from $U^{\dagger} U = I$ reduce $2n^2$ by $n^2$)
  More precisely: $\dim U(n) = n^2$, $\dim SU(n) = n^2 - 1$
- Both compact: unitary matrices have $|\text{eigenvalues}| = 1$, so entries bounded
- Both connected

**Explicit matrix forms:**

$$
U(1) = \{e^{i\theta} : \theta \in [0, 2\pi)\} \quad \text{— } 1 \times 1 \text{ unitary matrices, } \dim = 1
$$

$$
SU(2) = \left\{ \begin{pmatrix} a & -\bar{b} \\ b & \bar{a} \end{pmatrix} : a, b \in \mathbb{C},\ |a|^2 + |b|^2 = 1 \right\} \quad \text{— } \dim = 3
$$

$$
U(2) = \{A \in GL(2, \mathbb{C}) : A^{\dagger} A = I\} \quad \text{— } \dim = 4
$$

$$
SU(3) = \{A \in GL(3, \mathbb{C}) : A^{\dagger} A = I,\ \det A = 1\} \quad \text{— } \dim = 8
$$

### 3.5 The Relationship Between These Groups

$$
SU(n) \subset U(n) \subset GL(n, \mathbb{C})
$$

$$
U(n) \cong SU(n) \times U(1) / \mathbb{Z}_n \quad \text{(almost a direct product)}
$$

$$
U(n) \cong (SU(n) \times U(1)) / \{(e^{2\pi i k/n} I,\ e^{-2\pi i k/n}) : k = 0, \ldots, n-1\}
$$

The exact sequence:

$$
1 \to SU(n) \to U(n) \xrightarrow{\det} U(1) \to 1
$$

The determinant map $U(n) \to U(1)$ is surjective with kernel $SU(n)$, so $U(n)/SU(n) \cong U(1)$.
This means $U(n)$ is "$SU(n)$ times $U(1)$" with a discrete identification.

---

## 4. The Lie Algebra: Tangent Space at the Identity

### 4.0 Unpacking "Tangent Space at the Identity"

Before giving the formal definition, it is worth building up what this phrase means
from scratch — because it contains three separate ideas packed into four words.

---

#### Step 1: What Is a Tangent Vector?

Start with a smooth manifold $M$ (think of a surface in $\mathbb{R}^3$, or a more abstract space
like a Lie group). A **tangent vector** at a point $p \in M$ is the velocity vector of
a smooth curve passing through $p$.

Concretely: take any smooth curve $\gamma: (-\epsilon, \epsilon) \to M$ with $\gamma(0) = p$.
The tangent vector to $\gamma$ at $p$ is its velocity:

$$
\dot{\gamma}(0) = \left.\frac{d\gamma}{dt}\right|_{t=0} \in T_p M
$$

Two curves that pass through $p$ with the same velocity define the same tangent vector.
So a tangent vector is really an **equivalence class of curves** through $p$,
where two curves are equivalent if they have the same first-order behavior at $p$.

**In $\mathbb{R}^n$:** Every tangent vector at every point is literally a vector in $\mathbb{R}^n$.
The tangent space $T_p \mathbb{R}^n = \mathbb{R}^n$ for all $p$.

**On a sphere $S^2$:** At a point $p \in S^2$, the tangent space $T_p S^2$ is the plane
tangent to the sphere at $p$ — a 2-dimensional flat vector space sitting in $\mathbb{R}^3$.
Tangent vectors point along the surface; they do not point into or out of the sphere.

```
         T_p S² (tangent plane)
           ─────────
         /     ↑     \
        /    v ∈ T_p  \
       /               \
      |        p        |   ← the sphere S²
       \               /
        \             /
```

---

#### Step 2: The Tangent Space $T_p M$

The **tangent space** at $p$ is the collection of all tangent vectors at $p$:

$$
T_p M = \{\dot{\gamma}(0) : \gamma \text{ smooth curve in } M \text{ with } \gamma(0) = p\}
$$

Key facts:
- $T_p M$ is a **vector space** of dimension equal to $\dim(M)$
- Each point has its own tangent space — $T_p M$ and $T_q M$ are different spaces for $p \neq q$
- Tangent vectors at $p$ "live at $p$" — you cannot directly compare tangent vectors
  at different points without a connection (which is exactly what gauge theory provides)

**In coordinates:** If $(x^1, \ldots, x^n)$ are local coordinates near $p$, then
$T_p M$ has basis $\{\partial/\partial x^1|_p,\ \ldots,\ \partial/\partial x^n|_p\}$.
A general tangent vector is $v = v^i \partial/\partial x^i|_p$ for some real numbers $v^1, \ldots, v^n$.

**For a matrix Lie group $G \subset GL(n)$:**
$G$ is a subset of the space $\mathrm{Mat}(n \times n) \cong \mathbb{R}^{n^2}$. The tangent space at any point
$g \in G$ consists of the velocity vectors of smooth curves in $G$ passing through $g$:

$$
T_g G = \{\dot{\gamma}(0) : \gamma: (-\epsilon, \epsilon) \to G \text{ smooth},\ \gamma(0) = g\} \subset \mathrm{Mat}(n \times n)
$$

These are matrices — specifically the matrices that are "tangent to $G$ at $g$."

---

#### Step 3: Why the Identity Is Special

A Lie group $G$ has a distinguished point: the **identity element** $e$ (the matrix $I$
for matrix groups). The tangent space at the identity:

$$
T_e G = \{\dot{\gamma}(0) : \gamma \text{ smooth curve in } G,\ \gamma(0) = e\}
$$

is special for one reason: the group structure lets you **move any tangent space
back to $T_e G$** via left or right multiplication.

For any $g \in G$, left multiplication $L_g: G \to G$ defined by $L_g(h) = gh$ is a
diffeomorphism. Its derivative maps tangent spaces:

$$
(dL_g)_e : T_e G \to T_g G
$$

This map is a linear isomorphism — it identifies $T_e G$ with $T_g G$ for every $g$.
So all tangent spaces are isomorphic to $T_e G$, and $T_e G$ is the canonical
"home base" for all tangent vectors.

**Physical interpretation:** Every infinitesimal transformation in the group can
be "brought home" to the identity. The Lie algebra $T_e G$ collects all infinitesimal
group transformations in one place.

---

#### Step 4: Computing $T_e G$ for Matrix Groups

For a matrix Lie group $G \subset GL(n)$, a curve through the identity means:

$$
\gamma: (-\epsilon, \epsilon) \to G \subset GL(n) \quad \text{with } \gamma(0) = I
$$

The tangent vector at the identity is:

$$
\dot{\gamma}(0) = \left.\frac{d\gamma}{dt}\right|_{t=0} \in \mathrm{Mat}(n \times n)
$$

So $T_e G$ consists of all matrices $X$ that are the velocity of some smooth curve in $G$
starting at $I$.

**For $U(1)$:** Curves through the identity in $U(1) = \{e^{i\theta}\}$ look like $\gamma(t) = e^{i\alpha(t)}$
with $\alpha(0) = 0$. The velocity is:

$$
\dot{\gamma}(0) = i \dot{\alpha}(0) \in i\mathbb{R}
$$

Any purely imaginary number is achievable (choose $\alpha(t) = ct$ for any $c \in \mathbb{R}$).
So $T_e U(1) = i\mathbb{R}$ — the purely imaginary numbers. $\checkmark$

**For $SU(2)$:** Curves through $I$ in $SU(2) \subset \mathrm{Mat}(2 \times 2, \mathbb{C})$. Write $\gamma(t) = I + tX + O(t^2)$.
For $\gamma(t) \in SU(2)$ we need:
- $\gamma(t)^{\dagger} \gamma(t) = I$: differentiating $\to X^{\dagger} + X = 0$ (skew-Hermitian)
- $\det \gamma(t) = 1$: differentiating $\to \mathrm{tr} X = 0$ (traceless)

So $T_e SU(2) = \{X \in \mathrm{Mat}(2 \times 2, \mathbb{C}) : X^{\dagger} = -X,\ \mathrm{tr} X = 0\} = \mathfrak{su}(2)$. $\checkmark$

This computation — differentiate the defining equations of $G$ at the identity — is the
general recipe for finding the Lie algebra of any matrix Lie group.

---

#### Step 5: The Lie Algebra IS $T_e G$, as a Vector Space

The tangent space $T_e G$ is a vector space of dimension $n = \dim(G)$:

$$
\dim(T_e G) = \dim(G) \text{ as a manifold}
$$

Examples:

$$
\begin{aligned}
T_e U(1)   &= i\mathbb{R}                  & \dim &= 1 \\
T_e SU(2)  &= \mathfrak{su}(2)             & \dim &= 3 \\
T_e SU(3)  &= \mathfrak{su}(3)             & \dim &= 8 \\
T_e GL(n)  &= \mathrm{Mat}(n \times n)     & \dim &= n^2
\end{aligned}
$$

The Lie algebra $\mathfrak{g}$ is $T_e G$ **as a vector space**, plus the additional structure of
the Lie bracket $[\cdot, \cdot]$ coming from the group multiplication.

The bracket does not come from the tangent space structure alone — it comes from
differentiating the group commutator $g h g^{-1} h^{-1}$ twice at the identity. Concretely:

If $\gamma_1(t) = e^{tX}$ and $\gamma_2(t) = e^{tY}$ are curves through $I$, then:

$$
e^{tX} e^{tY} e^{-tX} e^{-tY} = I + t^2 [X, Y] + O(t^3)
$$

The bracket $[X, Y] = XY - YX$ measures the second-order failure to commute.
It is invisible at first order (the linear term vanishes), and it first appears at
second order — which is why we differentiate the group commutator **twice**
to extract the Lie algebra structure.

---

#### The Complete Picture: Three Things at Once

- $T_e G$ **as a set** $=$ equivalence classes of smooth curves through $e$
- $T_e G$ **as a vector space** $=$ infinitesimal directions you can move away from $e$ in $G$
- $T_e G$ **as a Lie algebra** $=$ infinitesimal generators of $G$, with bracket $[X, Y]$ encoding the non-commutativity of group multiplication

The phrase "tangent space at the identity" refers to the first two.
The full Lie algebra $\mathfrak{g} = T_e G$ adds the third — the bracket — on top.

---

#### A Visual Summary for $U(1)$

```
The group U(1) = S¹:

          1 ∈ U(1)
          |
    ─────●─────────  ←  T_e U(1) = iℝ (vertical line through 1)
         |↑
         | γ̇(0) = iα̇(0)   (tangent vector = purely imaginary number)
         |
  γ(t) = e^{iα(t)}   (a curve in U(1) starting at e = 1)
```

The circle $U(1)$ is 1-dimensional.
Its tangent space at $e = 1$ is a 1-dimensional real vector space: $i\mathbb{R}$.
The "vertical line" is the Lie algebra $\mathfrak{u}(1) = i\mathbb{R}$.
The exponential map wraps this line around the circle: $e^{i\theta} \mapsto U(1)$.

---

### 4.1 Definition

The **Lie algebra** $\mathfrak{g}$ of a Lie group $G$ is:

$$
\mathfrak{g} = T_e G \quad \text{(tangent space at the identity element } e\text{)}
$$

As a vector space, $\mathfrak{g}$ has dimension equal to the dimension of $G$ (as a manifold).

But $\mathfrak{g}$ has more structure than just a vector space: it has a **Lie bracket** $[\cdot, \cdot]$:

$$
[\cdot, \cdot]: \mathfrak{g} \times \mathfrak{g} \to \mathfrak{g}
$$

satisfying:

1. **Bilinearity:** $[\alpha X + \beta Y, Z] = \alpha[X, Z] + \beta[Y, Z]$
2. **Antisymmetry:** $[X, Y] = -[Y, X]$
3. **Jacobi identity:** $[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$

### 4.1a The Lie Bracket in Detail: Bilinear, Antisymmetric, Jacobi

The three axioms listed above — bilinearity, antisymmetry, and the Jacobi identity —
deserve a careful unpacking. They are not arbitrary: each encodes a precise geometric
or physical requirement.

#### What the bracket IS

The bracket $[\cdot, \cdot]$ is a rule that takes **two vectors in $\mathfrak{g}$ and produces a third**:

$$
[\cdot, \cdot] : \mathfrak{g} \times \mathfrak{g} \to \mathfrak{g}
$$

$$
[X, Y] = Z \quad \text{where } X, Y, Z \in \mathfrak{g}
$$

For matrix Lie algebras — the case you encounter in all gauge theory — it is the
**matrix commutator**:

$$
[X, Y] = XY - YX
$$

The three axioms are then properties of this commutator.

#### Bilinear

"Bilinear" means **linear in each slot separately**, with the other held fixed.

**Linear in the first slot:**

$$
[\alpha X + \beta Y, Z] = \alpha [X, Z] + \beta [Y, Z]
$$

**Linear in the second slot:**

$$
[X, \alpha Y + \beta Z] = \alpha [X, Y] + \beta [X, Z]
$$

Together: scaling or summing generators works consistently in either argument.
This is the same structure as the determinant (bilinear in rows), the inner product,
or matrix multiplication — "bi-linear" simply means linear in both inputs.

**Verification for the matrix commutator:**

$$
\begin{aligned}
[\alpha X + \beta Y, Z] &= (\alpha X + \beta Y) Z - Z (\alpha X + \beta Y) \\
                        &= \alpha XZ + \beta YZ - \alpha ZX - \beta ZY \\
                        &= \alpha (XZ - ZX) + \beta (YZ - ZY) \\
                        &= \alpha [X, Z] + \beta [Y, Z] \quad \checkmark
\end{aligned}
$$

**Why bilinearity is required:** The bracket must respect the vector space structure
of $\mathfrak{g}$. If $X$ generates a rotation by angle $\theta$, then $2X$ should generate a rotation by
$2\theta$, and the bracket should reflect this scaling consistently.

#### Antisymmetric

"Antisymmetric" means **swapping the two arguments negates the result**:

$$
[X, Y] = -[Y, X] \quad \text{for all } X, Y \in \mathfrak{g}
$$

**Immediate consequences:**

*The bracket of anything with itself is zero:*

$$
[X, X] = -[X, X] \quad \Rightarrow \quad 2[X, X] = 0 \quad \Rightarrow \quad [X, X] = 0
$$

*Over $\mathbb{R}$ or $\mathbb{C}$ (characteristic $\neq 2$): "$[X, X] = 0$ for all $X$" and "$[X, Y] = -[Y, X]$"
are equivalent conditions.*

**Verification for the matrix commutator:**

$$
[Y, X] = YX - XY = -(XY - YX) = -[X, Y] \quad \checkmark
$$

**Why antisymmetry is required — the geometric meaning:**

Antisymmetry encodes the fact that the bracket measures **non-commutativity** —
how much two group transformations fail to commute. The Baker-Campbell-Hausdorff
formula makes this precise:

$$
e^{tX} e^{tY} = e^{tY} e^{tX} \cdot e^{t^2 [X, Y] + O(t^3)}
$$

The bracket $[X, Y]$ is the leading correction when you reverse the order of two
transformations. Swapping $X \leftrightarrow Y$ reverses which ordering is "first" and which is
"second", so $[Y, X] = -[X, Y]$. When $[X, Y] = 0$, the two transformations commute
exactly — there is no correction.

This is the mathematical explanation of a crucial physical fact:

- $U(1)$: $[i\alpha, i\beta] = 0 \to$ phase rotations always commute $\to$ photons do not interact with each other $\to$ Maxwell equations are linear.
- $SU(2)$: $[T_a, T_b] \neq 0 \to$ weak isospin rotations do not commute $\to$ $W$ bosons carry weak charge and self-interact $\to$ Yang-Mills equations are nonlinear.

#### The Jacobi Identity

Together with bilinearity and antisymmetry, the bracket must satisfy:

$$
[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0
$$

This is the condition that makes $[\cdot, \cdot]$ a **Lie bracket** rather than an arbitrary
bilinear antisymmetric map.

**Intuition:** Define $\mathrm{ad}_X(Y) = [X, Y]$ (the operator "bracket with $X$ on the left").
The Jacobi identity says:

$$
\mathrm{ad}_X([Y, Z]) = [\mathrm{ad}_X(Y), Z] + [Y, \mathrm{ad}_X(Z)]
$$

This is a **Leibniz rule** (product rule) for $\mathrm{ad}_X$ acting on the bracket.
It says: "differentiating the bracket is like differentiating a product."
This is necessary for the exponential map to be well-behaved and for the Lie algebra
to faithfully represent the group multiplication near the identity.

**Verification for matrix commutators:**

$$
\begin{aligned}
[X, [Y, Z]] &= X(YZ - ZY) - (YZ - ZY) X = XYZ - XZY - YZX + ZYX \\
[Y, [Z, X]] &= YZX - YXZ - ZXY + XZY \\
[Z, [X, Y]] &= ZXY - ZYX - XYZ + YXZ
\end{aligned}
$$

Sum: each of the 12 terms appears exactly twice with opposite signs $\to 0$. $\checkmark$

#### The Complete Definition

A **Lie algebra** is a vector space $\mathfrak{g}$ over a field $\mathbb{F}$ with a bracket
$[\cdot, \cdot]: \mathfrak{g} \times \mathfrak{g} \to \mathfrak{g}$ satisfying:

- **L1. Bilinearity:**

$$
[\alpha X + \beta Y, Z] = \alpha [X, Z] + \beta [Y, Z], \qquad
[X, \alpha Y + \beta Z] = \alpha [X, Y] + \beta [X, Z]
$$

- **L2. Antisymmetry:**

$$
[X, Y] = -[Y, X]
$$

- **L3. Jacobi identity:**

$$
[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0
$$

#### A Concrete 3D Example: $\mathfrak{su}(2)$

The Lie algebra $\mathfrak{su}(2)$ has basis $\{T_1, T_2, T_3\}$ with brackets:

$$
[T_1, T_2] = iT_3, \qquad [T_2, T_3] = iT_1, \qquad [T_3, T_1] = iT_2
$$

Compactly: $[T_a, T_b] = i \epsilon_{abc} T_c$ where $\epsilon$ is the Levi-Civita symbol.

Check **bilinearity**:

$$
\begin{aligned}
[2T_1 + 3T_2, T_3] &= 2[T_1, T_3] + 3[T_2, T_3] \\
                   &= 2(-iT_2) + 3(iT_1) \\
                   &= 3iT_1 - 2iT_2 \quad \checkmark \quad \text{(a linear combination of basis vectors)}
\end{aligned}
$$

Check **antisymmetry**:

$$
[T_2, T_1] = -[T_1, T_2] = -iT_3 \quad \checkmark
$$

Check **Jacobi**:

$$
\begin{aligned}
& [T_1, [T_2, T_3]] + [T_2, [T_3, T_1]] + [T_3, [T_1, T_2]] \\
&\quad = [T_1, iT_1] + [T_2, iT_2] + [T_3, iT_3] \\
&\quad = i[T_1, T_1] + i[T_2, T_2] + i[T_3, T_3] \\
&\quad = i \cdot 0 + i \cdot 0 + i \cdot 0 = 0 \quad \checkmark \quad \text{(bracket of anything with itself is zero)}
\end{aligned}
$$

#### Summary Table

| Word | Meaning | Formula |
|---|---|---|
| **Bracket** | Binary operation producing a third element | $[X, Y] \in \mathfrak{g}$ |
| **Bilinear** | Linear in each argument separately | $[\alpha X + \beta Y, Z] = \alpha [X, Z] + \beta [Y, Z]$ |
| **Antisymmetric** | Swap arguments $\to$ negate result | $[X, Y] = -[Y, X]$ |
| **$[X, X] = 0$** | Consequence of antisymmetry | $[X, X] = 0$ for all $X$ |
| **Jacobi identity** | Cyclic sum of nested brackets $= 0$ | $[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$ |
| **For matrices** | All three hold for the commutator | $[X, Y] = XY - YX$ |
| **For $U(1)$** | All brackets vanish (abelian) | $[X, Y] = 0$ for all $X, Y$ |

---

### 4.2 Computing the Lie Algebra for Matrix Groups

For matrix Lie groups, the Lie algebra has a concrete description:

**Theorem.** For a matrix Lie group $G \subset GL(n)$, the Lie algebra $\mathfrak{g}$ consists of all
matrices $X$ such that $e^{tX} \in G$ for all $t \in \mathbb{R}$.

The Lie bracket is simply the matrix commutator:

$$
[X, Y] = XY - YX
$$

This is automatic from the matrix structure and is the reason matrix Lie groups
are particularly tractable.

### 4.3 The Lie Algebras of the Classical Groups

**$\mathfrak{gl}(n) = T_I GL(n)$:**
All $n \times n$ matrices (the identity is open in $GL(n)$, so the tangent space is the full
matrix space):

$$
\mathfrak{gl}(n, \mathbb{R}) = \mathrm{Mat}(n \times n, \mathbb{R}) \quad \text{with bracket } [X, Y] = XY - YX
$$

$$
\mathfrak{gl}(n, \mathbb{C}) = \mathrm{Mat}(n \times n, \mathbb{C})
$$

**$\mathfrak{sl}(n) = T_I SL(n)$:**
The condition $\det(e^{tX}) = 1$ differentiates to $\mathrm{tr}(X) = 0$:

$$
\mathfrak{sl}(n) = \{X \in \mathrm{Mat}(n \times n) : \mathrm{tr} X = 0\} \quad \text{(traceless matrices)}
$$

Proof: $\left.\frac{d}{dt}\right|_{t=0} \det(e^{tX}) = \mathrm{tr}(X) \cdot \det(I) = \mathrm{tr}(X)$.

**$\mathfrak{o}(n) = T_I O(n)$:**
The condition $(e^{tX})^T e^{tX} = I$ differentiates to $X^T + X = 0$:

$$
\mathfrak{o}(n) = \{X \in \mathrm{Mat}(n \times n, \mathbb{R}) : X^T = -X\} \quad \text{(antisymmetric matrices)}
$$

Dimension: $n(n-1)/2$ $\checkmark$

**$\mathfrak{u}(n) = T_I U(n)$:**
The condition $(e^{tX})^{\dagger} e^{tX} = I$ differentiates to $X^{\dagger} + X = 0$:

$$
\mathfrak{u}(n) = \{X \in \mathrm{Mat}(n \times n, \mathbb{C}) : X^{\dagger} = -X\} \quad \text{(skew-Hermitian matrices)}
$$

**$\mathfrak{su}(n) = T_I SU(n)$:**
Both traceless AND skew-Hermitian:

$$
\mathfrak{su}(n) = \{X \in \mathrm{Mat}(n \times n, \mathbb{C}) : X^{\dagger} = -X,\ \mathrm{tr} X = 0\}
$$

Dimension: $n^2 - 1$ $\checkmark$

### 4.4 The Lie Algebra of $U(1)$

$$
U(1) = \{e^{i\theta} : \theta \in \mathbb{R}\}
$$

Lie algebra $\mathfrak{u}(1)$:
- $\left.\frac{d}{dt}\right|_{t=0} e^{itX} = iX$ must be in $T_e U(1)$
- $e^{tX} \in U(1)$ for all $t$ means $|e^{tX}| = 1$
- $\to X$ must be purely imaginary: $X = i\theta$ for some $\theta \in \mathbb{R}$

$$
\mathfrak{u}(1) = i\mathbb{R} = \{i\theta : \theta \in \mathbb{R}\}
$$

As a vector space, $\mathfrak{u}(1) \cong \mathbb{R}$ (one-dimensional).
The Lie bracket: $[i\theta_1, i\theta_2] = i^2 \theta_1 \theta_2 - i^2 \theta_2 \theta_1 = 0$.

**$U(1)$ has an abelian Lie algebra — all brackets vanish.**

This is the fundamental reason electromagnetism is simpler than the weak and strong
forces: the gauge group $U(1)$ is abelian, so $[A_\mu, A_\nu] = 0$, and the photon
self-interaction term vanishes from the field strength tensor.

---

## 5. The Exponential Map

### 5.1 Definition and Convergence

The **exponential map** $\exp: \mathfrak{g} \to G$ is defined by:

$$
\exp(X) = e^X = \sum_{n=0}^{\infty} \frac{X^n}{n!} = I + X + \frac{X^2}{2!} + \frac{X^3}{3!} + \cdots
$$

For matrix Lie groups, this is literally the matrix exponential, which converges
absolutely for all matrices $X$.

### 5.2 Key Properties

1. $\exp(0) = I$ (identity at origin)
2. $\exp((s+t) X) = \exp(sX) \exp(tX)$ (one-parameter subgroup)
3. $\left.\frac{d}{dt}\right|_{t=0} \exp(tX) = X$ ($X$ is the velocity at the identity)
4. $\exp(-X) = (\exp X)^{-1}$ (inverse via negation)
5. $\det(\exp X) = e^{\mathrm{tr} X}$ (Liouville formula)
6. If $XY = YX$: $\exp(X + Y) = \exp(X) \exp(Y)$ (Baker-Campbell-Hausdorff reduces here)
7. In general: $\exp(X) \exp(Y) \neq \exp(X + Y)$ for non-commuting $X, Y$

Property 2 means $t \mapsto \exp(tX)$ is a **one-parameter subgroup** of $G$ — a smooth
homomorphism $\mathbb{R} \to G$. Every one-parameter subgroup arises this way.

### 5.3 The Baker-Campbell-Hausdorff Formula

For non-commuting $X, Y$:

$$
\exp(X) \exp(Y) = \exp\left( X + Y + \frac{1}{2}[X, Y] + \frac{1}{12}([X, [X, Y]] - [Y, [X, Y]]) + \cdots \right)
$$

The full series involves only nested Lie brackets. This shows that the group
multiplication near the identity is completely determined by the Lie bracket.

**Consequence for gauge theory:** The non-abelian structure $[A_\mu, A_\nu] \neq 0$ is
precisely the BCH correction when composing gauge transformations. It is the
reason why Yang-Mills theory has self-interacting gauge bosons.

### 5.4 When Is $\exp$ Surjective?

- **For compact, connected Lie groups** (e.g., $U(1)$, $SU(n)$, $SO(n)$): $\exp$ is surjective.
  Every group element can be written as $e^X$ for some $X \in \mathfrak{g}$.

- **For non-compact groups** (e.g., $GL(n)$, $SL(n)$): $\exp$ is not surjective in general.
  Example: in $SL(2, \mathbb{R})$, the matrix $\begin{pmatrix} -1 & 1 \\ 0 & -1 \end{pmatrix}$ is not in the image of $\exp$.

- **For $U(1)$ specifically:**
  $\exp(i\theta) = e^{i\theta}$, and every $z = e^{i\varphi} \in U(1)$ is hit by $i\varphi \in \mathfrak{u}(1)$.
  The map $\exp: \mathfrak{u}(1) = i\mathbb{R} \to U(1)$ is surjective but not injective:
  $\exp(i\theta) = \exp(i(\theta + 2\pi n))$ for any integer $n$.
  The kernel is $2\pi i \mathbb{Z} \subset i\mathbb{R}$.

### 5.5 The Local Diffeomorphism Property

Near the identity, $\exp$ is a local diffeomorphism: there exists a neighborhood
$U$ of $0$ in $\mathfrak{g}$ and a neighborhood $V$ of $e$ in $G$ such that $\exp: U \to V$ is a
diffeomorphism with smooth inverse $\log: V \to U$.

This is proved by the inverse function theorem: $(d \exp)_0 = \mathrm{id}$.

This is why **working in the Lie algebra is the same as working near the identity
in the Lie group** — they are locally diffeomorphic.

---

## 6. The Adjoint Representation

### 6.1 The Adjoint Action of $G$ on $\mathfrak{g}$

For $g \in G$, the conjugation map $C_g: G \to G$ defined by $C_g(h) = g h g^{-1}$ is a Lie
group automorphism (it preserves the group structure). Its derivative at the
identity gives a linear map on $\mathfrak{g}$:

$$
\mathrm{Ad}_g : \mathfrak{g} \to \mathfrak{g}, \qquad \mathrm{Ad}_g(X) = (dC_g)_e(X)
$$

For matrix groups: $\mathrm{Ad}_g(X) = g X g^{-1}$.

**Key property:** $\mathrm{Ad}_g$ is a Lie algebra automorphism:

$$
\mathrm{Ad}_g([X, Y]) = [\mathrm{Ad}_g(X),\ \mathrm{Ad}_g(Y)]
$$

The map $\mathrm{Ad}: G \to GL(\mathfrak{g})$ defined by $g \mapsto \mathrm{Ad}_g$ is a **representation** of $G$
(a Lie group homomorphism from $G$ to $GL$ of some vector space) — specifically
the **adjoint representation**.

### 6.2 The Adjoint Representation of the Lie Algebra

Differentiating $\mathrm{Ad}: G \to GL(\mathfrak{g})$ at the identity gives the **adjoint representation
of the Lie algebra**:

$$
\mathrm{ad}: \mathfrak{g} \to \mathfrak{gl}(\mathfrak{g}) \quad \text{(linear maps on } \mathfrak{g}\text{)}
$$

$$
\mathrm{ad}_X(Y) = [X, Y]
$$

This is just the Lie bracket! The Jacobi identity becomes:

$$
\mathrm{ad}_X \circ \mathrm{ad}_Y - \mathrm{ad}_Y \circ \mathrm{ad}_X = \mathrm{ad}_{[X, Y]} \quad \text{(Jacobi = representation property)}
$$

### 6.3 Why This Matters for Gauge Theory

In Yang-Mills theory, the field strength $F_{\mu\nu}$ transforms in the adjoint representation:

$$
F_{\mu\nu} \to g F_{\mu\nu} g^{-1} = \mathrm{Ad}_g(F_{\mu\nu})
$$

This is exactly the adjoint action. The gauge bosons (which are sections of the
adjoint bundle) transform under $\mathrm{Ad}_g$, not under the fundamental representation.

For $U(1)$: $\mathrm{Ad}_g = \mathrm{id}$ (since $U(1)$ is abelian: $g X g^{-1} = X$ for all $g, X$).
This is why the photon is neutral — the $U(1)$ adjoint representation is trivial.

For $SU(N)$: $\mathrm{Ad}_g$ is non-trivial. Gluons and $W/Z$ bosons carry their own charge.

### 6.4 The Killing Form

The **Killing form** is the symmetric bilinear form on $\mathfrak{g}$:

$$
B(X, Y) = \mathrm{tr}(\mathrm{ad}_X \circ \mathrm{ad}_Y) \quad \text{(trace in the adjoint representation)}
$$

For matrix Lie algebras, this simplifies to:

$$
B(X, Y) = 2N \cdot \mathrm{tr}(XY) \quad \text{(for } \mathfrak{su}(N) \text{ in the fundamental representation)}
$$

The Killing form is:
- **Negative definite** for compact semisimple Lie algebras (like $\mathfrak{su}(N)$)
- **Non-degenerate** for semisimple Lie algebras (Cartan's criterion)
- Used to raise/lower Lie algebra indices (like the metric tensor)
- Appears in the Yang-Mills action: $S_{YM} = -\frac{1}{4g^2} \int \mathrm{tr}(F_{\mu\nu} F^{\mu\nu}) d^4 x$

---

## 7. $U(1)$ in Full Detail

### 7.1 Five Equivalent Descriptions

- **Description 1 (Topological):** $S^1 =$ unit circle in $\mathbb{R}^2$ or $\mathbb{C}$
- **Description 2 (Algebraic):** $\mathbb{R}/\mathbb{Z}$ or $\mathbb{R}/2\pi\mathbb{Z}$ (quotient of $\mathbb{R}$ by integer lattice)
- **Description 3 (Matrix):** $\{e^{i\theta} : \theta \in \mathbb{R}\} \subset \mathbb{C} \cong GL(1, \mathbb{C})$
- **Description 4 (Representation):** smallest compact connected group whose representations classify $U(1)$ charges (integers $\mathbb{Z}$)
- **Description 5 (Gauge):** gauge group of electromagnetism

All five are isomorphic as Lie groups.

### 7.2 The Manifold Structure of $U(1)$

$U(1)$ is a 1-dimensional smooth manifold — a circle. The single coordinate chart
(with overlap) is:

**Chart 1:** $U_1 = U(1) \setminus \{-1\}$
- $\varphi_1 : U_1 \to (-\pi, \pi)$
- $\varphi_1(e^{i\theta}) = \theta$ for $\theta \in (-\pi, \pi)$

**Chart 2:** $U_2 = U(1) \setminus \{+1\}$
- $\varphi_2 : U_2 \to (0, 2\pi)$
- $\varphi_2(e^{i\theta}) = \theta$ for $\theta \in (0, 2\pi)$

**Transition map** on $U_1 \cap U_2 = U(1) \setminus \{\pm 1\}$:
- Two components: upper semicircle $\theta \in (0, \pi)$ and lower semicircle $\theta \in (-\pi, 0)$
- On each: $\varphi_2 \circ \varphi_1^{-1}$ is the identity (just the same coordinate)

$U(1)$ requires two charts because $S^1$ is not diffeomorphic to an open interval
(it is compact, while open intervals are not).

### 7.3 The Lie Algebra $\mathfrak{u}(1)$

$$
\mathfrak{u}(1) = T_{1} U(1) = i\mathbb{R} = \{i\theta : \theta \in \mathbb{R}\}
$$

- As a vector space: 1-dimensional over $\mathbb{R}$
- Generator: $T = i$ (or sometimes written as $T = -i$ by convention)
- Bracket: $[i\alpha, i\beta] = (i\alpha)(i\beta) - (i\beta)(i\alpha) = 0$ (abelian!)

The Lie algebra is abelian — the only 1-dimensional Lie algebra over $\mathbb{R}$.

**The standard physics convention** uses Hermitian generators:

Instead of working with $\mathfrak{u}(1) = i\mathbb{R}$, physicists write the generator as $T = 1$
(real, Hermitian), and write group elements as $g = e^{i\theta T} = e^{i\theta}$.
So the gauge potential $A_\mu = A_\mu \cdot T$ is real-valued.

This is why the gauge potential in electromagnetism is a real-valued 1-form
(not imaginary-valued), and why the field strength $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$ is real.

### 7.4 Representations of $U(1)$

The irreducible representations (irreps) of $U(1)$ are all 1-dimensional
(since $U(1)$ is abelian, by Schur's lemma). They are classified by integers $n \in \mathbb{Z}$:

$$
\rho_n : U(1) \to GL(1, \mathbb{C}) = \mathbb{C}^{\ast}, \qquad \rho_n(e^{i\theta}) = e^{in\theta}
$$

The integer $n$ is the **charge** of the representation. For electromagnetism:
- Electron: charge $n = -1$ (or $+1$ depending on sign convention)
- Positron: charge $n = +1$
- Photon: charge $n = 0$ (neutral, in the adjoint)
- $W$ boson: charge $n = \pm 1$ (but under $SU(2) \times U(1)$, not just $U(1)$)

**Why charges are integers:** The map $n \mapsto \rho_n$ is a group homomorphism from $\mathbb{Z}$
to the set of 1-d representations of $U(1)$. For $\rho_n$ to be a well-defined
representation (single-valued), we need:

$$
\rho_n(e^{i\theta}) = e^{in\theta} \quad \text{with } e^{in(\theta + 2\pi)} = e^{in\theta}
$$

This requires $n \in \mathbb{Z}$. The quantization of electric charge is a direct consequence
of the compactness of $U(1)$ ($\pi_1(U(1)) = \mathbb{Z}$ implies charge quantization).

### 7.5 $U(1)$ as a Quotient

The most illuminating way to understand $U(1)$ is as a quotient:

$$
U(1) \cong \mathbb{R}/2\pi\mathbb{Z}
$$

The quotient map $q: \mathbb{R} \to U(1)$ is:

$$
q(\theta) = e^{i\theta}
$$

This is a covering map:
- $q$ is a local homeomorphism (local diffeomorphism)
- Each point $e^{i\theta} \in U(1)$ has preimage $q^{-1}(e^{i\theta}) = \{\theta + 2\pi n : n \in \mathbb{Z}\}$
- The fiber is countably infinite: $\mathbb{Z}$

This is the **universal cover** of $U(1)$: $\mathbb{R}$ is simply connected, and $U(1) = \mathbb{R}/2\pi\mathbb{Z}$.

The covering structure explains:
- Why $\pi_1(U(1)) = \mathbb{Z}$ (loops in $U(1)$ lift to paths in $\mathbb{R}$; the winding number is the
  endpoint displacement divided by $2\pi$)
- Why charges are quantized (representations of $U(1)$ come from representations of
  $\mathbb{R}$ that are periodic with period $2\pi$ — i.e., $e^{in\theta}$ for $n \in \mathbb{Z}$)
- Why monopole charges are quantized (Dirac quantization condition)

### 7.6 The Hopf Fibration: $U(1)$ Inside $S^3$

The **Hopf fibration** is a fiber bundle:

$$
U(1) \to S^3 \to S^2
$$

$S^3$ is the 3-sphere $\{(z_1, z_2) \in \mathbb{C}^2 : |z_1|^2 + |z_2|^2 = 1\}$.
The bundle map $\pi: S^3 \to S^2$ sends $(z_1, z_2) \mapsto$ the point on $S^2$ given by the
ratio $z_1/z_2 \in \mathbb{C} \cup \{\infty\} \cong S^2$ (Riemann sphere).
The fiber over each point is a copy of $U(1) = \{e^{i\theta}(z_1, z_2)\}$.

Physical significance: the Hopf fibration is the geometric structure of a magnetic
monopole of strength $1/2$. The non-trivial topology of the bundle (it is not
$S^2 \times S^1$) is the reason the monopole exists and its charge is quantized.

---

## 8. $SU(2)$: The Three-Sphere and Rotations

### 8.1 The Matrix Realization

$$
SU(2) = \left\{ \begin{pmatrix} a & -\bar{b} \\ b & \bar{a} \end{pmatrix} : a, b \in \mathbb{C},\ |a|^2 + |b|^2 = 1 \right\}
$$

Parametrize with $a = x_0 + ix_3$, $b = x_2 + ix_1$, where $x_0^2 + x_1^2 + x_2^2 + x_3^2 = 1$:

$$
SU(2) = \{x_0 I + i(x_1 \sigma_1 + x_2 \sigma_2 + x_3 \sigma_3) : x_0^2 + |\mathbf{x}|^2 = 1\}
$$

where $\sigma_1, \sigma_2, \sigma_3$ are the Pauli matrices:

$$
\sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \qquad
\sigma_2 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \qquad
\sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$

**As a manifold:** $SU(2) \cong S^3$ (the 3-sphere in $\mathbb{R}^4 = \mathbb{C}^2$).
This is the key geometric fact about $SU(2)$:

$$
SU(2) \text{ is the unit 3-sphere: } x_0^2 + x_1^2 + x_2^2 + x_3^2 = 1
$$

### 8.2 The Lie Algebra $\mathfrak{su}(2)$

$$
\mathfrak{su}(2) = \{X \in \mathrm{Mat}(2 \times 2, \mathbb{C}) : X^{\dagger} = -X,\ \mathrm{tr} X = 0\}
      = \left\{ \begin{pmatrix} ia & -\bar{b} + ic \\ -\bar{b} - ic & -ia \end{pmatrix} : a, b, c \in \mathbb{R} \right\}
$$

Standard basis (using $i/2$ times Pauli matrices):

$$
T_1 = i\sigma_1/2 = \begin{pmatrix} 0 & i/2 \\ i/2 & 0 \end{pmatrix}, \qquad
T_2 = i\sigma_2/2 = \begin{pmatrix} 0 & 1/2 \\ -1/2 & 0 \end{pmatrix}, \qquad
T_3 = i\sigma_3/2 = \begin{pmatrix} i/2 & 0 \\ 0 & -i/2 \end{pmatrix}
$$

Or in physics convention with Hermitian generators (dropping the $i$):

$$
L_a = \sigma_a / 2 \quad (a = 1, 2, 3)
$$

Lie bracket: $[L_a, L_b] = i \epsilon_{abc} L_c$ ($\epsilon_{abc}$ is the Levi-Civita symbol).

This is the **angular momentum algebra** from quantum mechanics! $SU(2)$ is the
quantum mechanical rotation group.

### 8.3 The Structure Constants of $SU(2)$

The Lie bracket $[L_a, L_b] = i \epsilon_{abc} L_c$ means the **structure constants** are:

$$
f^{abc} = \epsilon^{abc} \quad \text{(Levi-Civita symbol)}
$$

$$
f^{123} = f^{231} = f^{312} = 1, \qquad f^{132} = f^{213} = f^{321} = -1, \qquad \text{all others} = 0
$$

In gauge theory, the structure constants determine:
- The 3-gluon vertex (in QCD with $SU(3)$, analogous structure)
- The self-interaction of the $W$ bosons (in $SU(2)$ weak theory)
- The non-abelian part of $F_{\mu\nu}$: $F_{\mu\nu}^c = \partial_\mu A_\nu^c - \partial_\nu A_\mu^c - g f^{abc} A_\mu^a A_\nu^b$

### 8.4 The Exponential Map for $SU(2)$

For a general element $X = \mathbf{n} \cdot \mathbf{L} = n_a L_a$ (unit vector $\mathbf{n}$, angle $\theta$):

$$
\exp(\theta \mathbf{n} \cdot \mathbf{L}) = \cos(\theta/2) I + 2i \sin(\theta/2) \mathbf{n} \cdot \mathbf{L}
                                    = \cos(\theta/2) I + i \sin(\theta/2)(n_1 \sigma_1 + n_2 \sigma_2 + n_3 \sigma_3)
$$

This is a rotation by angle $\theta$ about axis $\mathbf{n}$ in the spin-$1/2$ representation.

**Key property:** $\exp$ is surjective for $SU(2)$ (since $SU(2) \cong S^3$ is simply connected
and compact). Every element of $SU(2)$ is $e^X$ for some $X \in \mathfrak{su}(2)$.

### 8.5 $SU(2)$ Is Simply Connected

$$
\begin{aligned}
\pi_0(SU(2)) &= 0 \quad \text{(connected)} \\
\pi_1(SU(2)) &= 0 \quad \text{(simply connected — no non-trivial loops!)} \\
\pi_2(SU(2)) &= 0 \\
\pi_3(SU(2)) &= \mathbb{Z} \quad \text{(non-trivial 3-spheres — instantons!)}
\end{aligned}
$$

Contrast with:

$$
\pi_1(U(1)) = \mathbb{Z} \quad (U(1) \text{ is NOT simply connected}), \qquad
\pi_1(SO(3)) = \mathbb{Z}_2 \quad (SO(3) \text{ is NOT simply connected})
$$

The simple connectivity of $SU(2)$ is why:
- $SU(2)$ is the universal cover of $SO(3)$
- The exponential map $\exp: \mathfrak{su}(2) \to SU(2)$ determines $SU(2)$ completely
- There are no theta-term complications from $\pi_1$ (but there are instantons from $\pi_3$)

---

## 9. $SO(3)$ vs $SU(2)$: Covering Groups and Spinors

### 9.1 The 2:1 Covering Map

There is a surjective Lie group homomorphism:

$$
\varphi : SU(2) \to SO(3), \qquad \varphi(U) \mathbf{v} = U \mathbf{v} U^{\dagger} \quad \text{for } \mathbf{v} = v^a \sigma_a \in \mathfrak{su}(2) \cong \mathbb{R}^3
$$

This map is 2:1: both $U$ and $-U$ map to the same rotation $R \in SO(3)$.

$$
\ker(\varphi) = \{I, -I\} = \mathbb{Z}_2, \qquad SO(3) = SU(2)/\mathbb{Z}_2
$$

### 9.2 Why $\pm I$ Map to the Same Rotation

$$
\varphi(U) \mathbf{v} = U \mathbf{v} U^{\dagger}
$$

$$
\varphi(-U) \mathbf{v} = (-U) \mathbf{v} (-U)^{\dagger} = U \mathbf{v} U^{\dagger} \quad \text{(the two minus signs cancel!)}
$$

So $\varphi(U) = \varphi(-U)$ — they give the same rotation.

### 9.3 The Physical Consequence: Spinors

A spin-$1/2$ particle (electron) transforms under $SU(2)$, not $SO(3)$:

$$
\text{Under rotation by } 2\pi: \quad U = \exp(2\pi \cdot i\sigma_3/2) = \exp(i\pi \sigma_3) = -I
$$

So a $2\pi$ rotation sends a spin-$1/2$ state to its negative:

$$
|\psi\rangle \to -|\psi\rangle \quad \text{under } 2\pi \text{ rotation}
$$

This is observable! In neutron interferometry, the beam acquires a phase of $-1$
under $360^{\circ}$ rotation and returns to its original phase after $720^{\circ}$. This is the
spinor property, and it comes directly from the 2:1 covering $SU(2) \to SO(3)$.

### 9.4 The General Pattern: Covering Groups

For every Lie group $G$, there is a unique **universal cover** $\tilde{G}$:
- $\tilde{G}$ is simply connected ($\pi_1(\tilde{G}) = 0$)
- $G = \tilde{G}/\Gamma$ for a discrete normal subgroup $\Gamma \cong \pi_1(G)$

| Group $G$ | Universal cover $\tilde{G}$ | Quotient |
|---|---|---|
| $U(1) = S^1$ | $\mathbb{R}$ | $\mathbb{R}/\mathbb{Z}$ |
| $SO(2)$ | $\mathbb{R}$ | $\mathbb{R}/\mathbb{Z}$ |
| $SO(3)$ | $SU(2)$ | $SU(2)/\mathbb{Z}_2$ |
| $SO(n),\ n \geq 3$ | $\mathrm{Spin}(n)$ | $\mathrm{Spin}(n)/\mathbb{Z}_2$ |
| $SO(n, 1)$, Lorentz | $SL(2, \mathbb{C})$ | $SL(2, \mathbb{C})/\mathbb{Z}_2$ |

In gauge theory, the relevant group is always the simply connected cover,
because the path integral and representations are better behaved.

---

## 10. $SU(3)$ and the Gell-Mann Matrices

### 10.1 The Structure of $SU(3)$

$$
SU(3) = \{A \in GL(3, \mathbb{C}) : A^{\dagger} A = I,\ \det A = 1\}
$$

- Dimension: $8$ (as a real manifold)
- Rank: $2$ (dimension of maximal torus = dimension of maximal abelian subgroup)

$SU(3)$ is the gauge group of QCD (quantum chromodynamics). The 8 gauge bosons
are the gluons, corresponding to the 8 generators of $\mathfrak{su}(3)$.

### 10.2 The Gell-Mann Matrices

The standard basis for $\mathfrak{su}(3)$ uses the **Gell-Mann matrices** $\lambda_1, \ldots, \lambda_8$
(analogous to Pauli matrices for $SU(2)$):

$$
\lambda_1 = \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}, \quad
\lambda_2 = \begin{pmatrix} 0 & -i & 0 \\ i & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}, \quad
\lambda_3 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

$$
\lambda_4 = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 1 & 0 & 0 \end{pmatrix}, \quad
\lambda_5 = \begin{pmatrix} 0 & 0 & -i \\ 0 & 0 & 0 \\ i & 0 & 0 \end{pmatrix}, \quad
\lambda_6 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}
$$

$$
\lambda_7 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & -i \\ 0 & i & 0 \end{pmatrix}, \qquad
\lambda_8 = \frac{1}{\sqrt{3}} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & -2 \end{pmatrix}
$$

Properties:

$$
\mathrm{tr}(\lambda_a) = 0 \quad \text{(traceless)}
$$

$$
\mathrm{tr}(\lambda_a \lambda_b) = 2\delta_{ab} \quad \text{(orthonormality)}
$$

$$
[\lambda_a, \lambda_b] = 2i f^{abc} \lambda_c \quad \text{(Lie bracket, } f^{abc} = SU(3) \text{ structure constants)}
$$

Generators in physics convention: $T_a = \lambda_a / 2$.

### 10.3 The Structure of $SU(3)$: Cartan Subalgebra

The **Cartan subalgebra** of $\mathfrak{su}(3)$ is the maximal abelian subalgebra:

$$
\text{Cartan subalgebra} = \mathrm{span}\{T_3, T_8\} = \mathrm{span}\{\lambda_3/2,\ \lambda_8/2\}
$$

These are the simultaneously diagonalizable generators — they commute: $[T_3, T_8] = 0$.

In QCD: $T_3$ and $T_8$ correspond to the two independently conserved color charges.
In the quark model: $T_3$ is the isospin component, $T_8$ is related to hypercharge.

### 10.4 Roots and the Root System

The 6 remaining generators ($\lambda_1, \ldots, \lambda_7$ minus $\lambda_3, \lambda_8$) can be organized into
**raising and lowering operators** — ladder operators that shift between
eigenstates of the Cartan subalgebra.

For $SU(3)$, the 6 roots form the $A_2$ root system:

```
         ●  (raising by α₁+α₂)
        / \
       /   \
      ●     ●  (raising by α₂, α₁+α₂)
     / \   / \
    /   \ /   \
   ●     ×     ●   (×  = origin, Cartan subalgebra)
    \   / \   /
     \ /   \ /
      ●     ●  (lowering by α₂, α₁+α₂)
       \   /
        \ /
         ●  (lowering by α₁+α₂)
```

This hexagonal root diagram is what Gell-Mann used to predict the existence
of the $\Omega^-$ baryon (the missing corner of the decuplet) before it was discovered.

---

## 11. Representations of Lie Groups

### 11.1 Definition

A **representation** of a Lie group $G$ is a smooth group homomorphism:

$$
\rho: G \to GL(V)
$$

where $V$ is a finite-dimensional vector space (the **representation space**) and
$GL(V)$ is the group of invertible linear maps on $V$.

The **dimension** of the representation is $\dim V$.
The representation is **faithful** if $\ker(\rho) = \{e\}$ (injective).

### 11.2 The Induced Representation of the Lie Algebra

Given a representation $\rho: G \to GL(V)$, differentiating at the identity gives:

$$
d\rho: \mathfrak{g} \to \mathfrak{gl}(V) = \mathrm{End}(V), \qquad d\rho(X) = \left.\frac{d}{dt}\right|_{t=0} \rho(\exp(tX))
$$

This is a Lie algebra homomorphism: $d\rho([X, Y]) = [d\rho(X),\ d\rho(Y)]$.

For matrix groups: $\rho(\exp X) = \exp(d\rho(X))$ (exponential intertwines the two).

### 11.3 Irreducible Representations

A representation $\rho: G \to GL(V)$ is **irreducible** (an irrep) if there is no
proper non-zero subspace $W \subset V$ with $\rho(g) W \subseteq W$ for all $g \in G$.

**Schur's Lemma:** If $\rho_1, \rho_2$ are irreps, then any intertwining map $T: V_1 \to V_2$
(with $\rho_2(g) T = T \rho_1(g)$) is either $0$ or an isomorphism.

Consequence for abelian groups (like $U(1)$):
All irreps are 1-dimensional. (Because any 1-d subspace is invariant for
abelian $G$, so the only irreducible representation on $V$ is when $\dim V = 1$.)

### 11.4 The Irreps of $SU(2)$

The irreps of $SU(2)$ are labeled by a half-integer $j = 0, 1/2, 1, 3/2, 2, \ldots$:

$$
\rho_j : SU(2) \to GL(\mathbb{C}^{2j+1})
$$

Basis: $|j, m\rangle$ for $m = -j,\ -j+1,\ \ldots,\ j-1,\ j$.

Action:

$$
T_3 |j, m\rangle = m |j, m\rangle
$$

$$
T_{\pm} |j, m\rangle = \sqrt{j(j+1) - m(m \pm 1)}\ |j, m \pm 1\rangle \qquad (T_{\pm} = T_1 \pm i T_2)
$$

In physics:
- $j = 0$: singlet (spin-0, scalar particle)
- $j = 1/2$: doublet (spin-$1/2$, electron, quarks, neutrinos)
- $j = 1$: triplet (spin-1, vector particles like $W$ bosons)
- $j = 3/2$: quadruplet (spin-$3/2$, delta baryons)

**Dimension of irrep $j$:** $2j + 1$.

**Tensor product:** $\rho_{j_1} \otimes \rho_{j_2} = \bigoplus_{j = |j_1 - j_2|}^{j_1 + j_2} \rho_j$ (Clebsch-Gordan decomposition).

### 11.5 The Irreps of $SU(3)$

The irreps of $SU(3)$ are labeled by two non-negative integers $(p, q)$:

$$
\text{Dimension} = (p+1)(q+1)(p+q+2)/2
$$

Examples:
- $(0, 0)$: dimension $1$ — singlet (colorless)
- $(1, 0)$: dimension $3$ — fundamental representation (quarks)
- $(0, 1)$: dimension $\bar{3}$ — anti-fundamental (antiquarks)
- $(1, 1)$: dimension $8$ — adjoint representation (gluons)
- $(3, 0)$: dimension $10$ — decuplet (baryon resonances)
- $(2, 1)$: dimension $15$

The famous quark model prediction: hadrons are color singlets ($p = q = 0$),
built from quarks in the $3$ and $\bar{3}$ and their bound states.

---

## 12. The Peter-Weyl Theorem and Harmonic Analysis

### 12.1 Functions on a Compact Lie Group

For a compact Lie group $G$ with Haar measure $dg$ (the unique translation-invariant
measure with $\int_G dg = 1$), the space $L^2(G)$ of square-integrable functions has
a beautiful decomposition.

**Peter-Weyl Theorem:** The matrix coefficients of irreducible representations:

$$
\rho_{ij} : G \to \mathbb{C}, \qquad g \mapsto \langle e_i,\ \rho_j(g) e_j \rangle \quad (i, j = 1, \ldots, \dim \rho)
$$

form an orthonormal basis for $L^2(G)$:

$$
L^2(G) = \bigoplus_{\rho \text{ irreps}} V_\rho \otimes V_\rho^{\ast}
$$

where $V_\rho$ is the representation space of irrep $\rho$.

### 12.2 Application to $U(1)$: Fourier Series

For $G = U(1)$:

- Irreps: $\rho_n(e^{i\theta}) = e^{in\theta}$ for $n \in \mathbb{Z}$
- Matrix coefficients: $e^{in\theta}$
- Peter-Weyl: $L^2(U(1)) = L^2(S^1) = \bigoplus_{n \in \mathbb{Z}} \mathbb{C} e^{in\theta}$

This is **Fourier series**: every square-integrable function on the circle
has a Fourier expansion $f(\theta) = \sum_{n \in \mathbb{Z}} c_n e^{in\theta}$.

The Peter-Weyl theorem is the group-theoretic generalization of Fourier analysis
to any compact Lie group.

### 12.3 Application to $SU(2)$: Spherical Harmonics

For $G = SU(2) \cong S^3$:

- Irreps: $\rho_j$ for $j = 0, 1/2, 1, 3/2, \ldots$
- Peter-Weyl gives an orthonormal basis for $L^2(S^3)$ labeled by $(j, m, m')$
- Restricting to functions invariant under $U(1) \subset SU(2)$ gives **spherical harmonics**
  on $S^2 = SU(2)/U(1)$

The spherical harmonics $Y_l^m(\theta, \varphi)$ are precisely the matrix coefficients of
$SU(2)$ representations restricted to the appropriate quotient. This is why
spherical harmonics have the angular momentum quantum numbers they do.

---

## 13. Roots, Weights, and the Classification of Simple Lie Algebras

### 13.1 The Cartan Classification

Every simple Lie algebra over $\mathbb{C}$ belongs to one of:

Classical algebras:
- $A_n\ (n \geq 1)$: $\mathfrak{sl}(n+1, \mathbb{C}) \to SU(n+1)$ (compact real form)
- $B_n\ (n \geq 2)$: $\mathfrak{so}(2n+1, \mathbb{C}) \to SO(2n+1)$
- $C_n\ (n \geq 3)$: $\mathfrak{sp}(2n, \mathbb{C}) \to Sp(n)$ (symplectic)
- $D_n\ (n \geq 4)$: $\mathfrak{so}(2n, \mathbb{C}) \to SO(2n)$

Exceptional algebras:
- $G_2$ (dim $14$)
- $F_4$ (dim $52$)
- $E_6$ (dim $78$)
- $E_7$ (dim $133$)
- $E_8$ (dim $248$)

The gauge groups of the Standard Model:
- $U(1) = U(1)$ (circle, not simple but simple up to center)
- $SU(2) = A_1$ (the simplest non-abelian simple Lie group)
- $SU(3) = A_2$
- $SU(5) = A_4$ (Grand Unified Theory gauge group)
- $SO(10) = D_5$ (another GUT candidate)
- $E_8 \times E_8$ (heterotic string theory)

### 13.2 Roots

For a simple Lie algebra $\mathfrak{g}$ with Cartan subalgebra $\mathfrak{h}$ (maximal abelian diagonalizable
subalgebra), the **roots** are the nonzero eigenvalues of the adjoint action of $\mathfrak{h}$:

If $H \in \mathfrak{h}$ and $E_\alpha \in \mathfrak{g}$ satisfy $[H, E_\alpha] = \alpha(H) E_\alpha$ for all $H \in \mathfrak{h}$,
then $\alpha : \mathfrak{h} \to \mathbb{C}$ is a root and $E_\alpha$ is a root vector (raising/lowering operator).

The root system encodes the entire Lie algebra structure:

$$
[H, E_\alpha] = \alpha(H) E_\alpha
$$

$$
[E_\alpha, E_{-\alpha}] = H_\alpha \quad \text{(coroot)}
$$

$$
[E_\alpha, E_\beta] = N_{\alpha\beta} E_{\alpha + \beta} \quad \text{(if } \alpha + \beta \text{ is a root)}
$$

### 13.3 Weights and Representations

For a representation $\rho: \mathfrak{g} \to \mathfrak{gl}(V)$, the **weights** of $\rho$ are the eigenvalues
of the Cartan subalgebra action on $V$:

$$
H v = \lambda(H) v \quad \text{for all } H \in \mathfrak{h} \text{ and some } v \in V
$$

$\lambda$ is a weight and $v$ is a weight vector.

For $SU(2)$:
- The single Cartan generator is $T_3$
- Weights are the eigenvalues $m = -j, \ldots, +j$ for spin-$j$ representation
- The "highest weight" is $j$ (the maximum eigenvalue)

For $SU(3)$:
- Two Cartan generators $T_3, T_8$
- Weights are 2d vectors (isospin, hypercharge)
- The weight diagram of the fundamental $(1, 0)$ rep is the quark triangle ($u, d, s$)

---

## 14. Homomorphisms, Subgroups, and Quotients

### 14.1 Lie Group Homomorphisms

A **Lie group homomorphism** $\varphi: G \to H$ is a smooth group homomorphism.

Key examples:
- $\det : GL(n) \to \mathbb{R}^{\ast}$ (determinant — surjective onto nonzero reals)
- $\det : U(n) \to U(1)$ (determinant of unitary matrix has $|\det| = 1$)
- $\mathrm{Ad} : G \to GL(\mathfrak{g})$ (adjoint representation)
- $\pi : SU(2) \to SO(3)$ (2:1 covering, kernel $= \mathbb{Z}_2$)
- $\exp : \mathfrak{g} \to G$ (not quite a homomorphism, but the exponential map)

### 14.2 Closed Subgroups

By the **Closed Subgroup Theorem** (Cartan): if $H$ is a closed subgroup of a
Lie group $G$, then $H$ is automatically a Lie group (an embedded submanifold).

This means: to specify a subgroup of a matrix Lie group, you just need a
closed subset that is also a subgroup — the smooth structure comes for free.

Examples of closed subgroups:

$$
U(1) \subset SU(2): \quad U(1) = \left\{ \begin{pmatrix} e^{i\theta} & 0 \\ 0 & e^{-i\theta} \end{pmatrix} : \theta \in \mathbb{R} \right\} \text{ (diagonal matrices)}
$$

$$
SU(2) \subset SU(3): \quad \text{Embed as block } \begin{pmatrix} SU(2) & 0 \\ 0 & 1 \end{pmatrix}
$$

$$
U(1) \times U(1) \subset SU(3): \quad \text{the maximal torus (diagonal matrices in } SU(3)\text{)}
$$

### 14.3 Quotient Groups and Homogeneous Spaces

If $H \subset G$ is a closed subgroup, the **quotient space** $G/H$ is a smooth manifold.
If $H$ is a normal subgroup ($g H g^{-1} = H$ for all $g$), then $G/H$ is a Lie group.

Important examples:
- $SO(3) = SU(2)/\mathbb{Z}_2$ ($\mathbb{Z}_2 = \{I, -I\}$ is normal in $SU(2)$)
- $S^2 = SU(2)/U(1) = SO(3)/SO(2)$ (2-sphere as a homogeneous space)
- $S^n = SO(n+1)/SO(n)$ ($n$-sphere as a homogeneous space)
- $\mathbb{R}P^n = SO(n+1)/O(n)$ (real projective space)
- $G_{SM}/H_{\text{unbroken}}$ (symmetry breaking in particle physics)

The last example is crucial: in the Standard Model, after Higgs symmetry breaking:

$$
U(1)_Y \times SU(2)_L \to U(1)_{EM}
$$

The surviving $U(1)_{EM}$ is a specific $U(1)$ subgroup of the original gauge group,
and the broken symmetries become the massive $W^{\pm}, Z$ bosons.

---

## 15. The Topology of Lie Groups

### 15.1 Homotopy Groups

The topological complexity of Lie groups is captured by their homotopy groups:

| $G$ | $\pi_0$ | $\pi_1$ | $\pi_2$ | $\pi_3$ |
|---|---|---|---|---|
| $U(1)$ | $0$ | $\mathbb{Z}$ | $0$ | $0$ |
| $SU(2)$ | $0$ | $0$ | $0$ | $\mathbb{Z}$ |
| $SU(3)$ | $0$ | $0$ | $0$ | $\mathbb{Z}$ |
| $SO(2)$ | $0$ | $\mathbb{Z}$ | $0$ | $0$ |
| $SO(3)$ | $0$ | $\mathbb{Z}_2$ | $0$ | $\mathbb{Z}$ |
| $Sp(1) \cong SU(2)$ | $0$ | $0$ | $0$ | $\mathbb{Z}$ |
| $G_2$ | $0$ | $0$ | $0$ | $\mathbb{Z}$ |

Key facts:
- All compact connected Lie groups are connected: $\pi_0 = 0$
- $\pi_1(G) = \pi_1(G/[G, G])$ — the fundamental group comes from the abelian part
- For simple simply connected compact Lie groups: $\pi_3 = \mathbb{Z}$ (Bott periodicity)
- $\pi_3(G) = \mathbb{Z}$ is the origin of instantons in Yang-Mills theory

### 15.2 Why $\pi_1(U(1)) = \mathbb{Z}$ and Why It Matters

A loop in $U(1) = S^1$ is a continuous map $\gamma: [0, 1] \to S^1$ with $\gamma(0) = \gamma(1) = 1$.
The homotopy class of $\gamma$ is its **winding number**: how many times $\gamma$ wraps around the circle.

- Winding number $+1$: $\gamma(t) = e^{2\pi i t}$ (once counterclockwise)
- Winding number $+2$: $\gamma(t) = e^{4\pi i t}$ (twice counterclockwise)
- Winding number $-1$: $\gamma(t) = e^{-2\pi i t}$ (once clockwise)
- Winding number $0$: any contractible loop

$\pi_1(U(1)) = \mathbb{Z}$ means loops are classified by their winding number.

**Physical consequences:**
1. **Charge quantization:** Representations of $U(1)$ are $\rho_n(e^{i\theta}) = e^{in\theta}$.
   The integer $n$ is the charge, and it is quantized because $n$ must be an integer
   for the representation to be single-valued (well-defined on $S^1$).

2. **Magnetic monopoles:** The Dirac quantization condition $eg = 2\pi$ (in natural units)
   requires the product of electric charge $e$ and magnetic charge $g$ to be quantized.
   This is because the gauge field on $S^2$ around a monopole is classified by
   $\pi_1(U(1)) = \mathbb{Z}$.

3. **Flux quantization:** In superconductors, the magnetic flux through a loop is
   quantized in units of $\Phi_0 = h/2e$. This comes from the $U(1)$ winding number of
   the Cooper pair condensate.

### 15.3 The Bott Periodicity Theorem

For any simple simply connected compact Lie group $G$:

$$
\pi_3(G) = \mathbb{Z}
$$

This is Bott periodicity applied to Lie groups. For $SU(n)$, the map $S^3 \to SU(n)$
generating $\pi_3$ is the basic instanton. The integer generator of $\pi_3$ is what labels
the topological charge (instanton number) of Yang-Mills gauge fields.

### 15.4 Compact vs. Non-Compact Real Forms

Every complex simple Lie algebra has a unique compact real form and multiple
non-compact real forms. In gauge theory:

Complex algebra $A_1 = \mathfrak{sl}(2, \mathbb{C})$:
- Compact real form: $\mathfrak{su}(2) \to SU(2)$ (gauge group of weak interaction)
- Non-compact: $\mathfrak{sl}(2, \mathbb{R}) \to SL(2, \mathbb{R})$ (appears in 2+1 gravity)

Complex algebra $D_1 = \mathfrak{so}(2, \mathbb{C}) \cong A_1$:
- Compact: $\mathfrak{so}(2) \to SO(2) \cong U(1)$ (rotation in 2D)
- Non-compact: $\mathfrak{so}(1, 1) \to SO(1, 1)$ (Lorentz boosts in 1+1D)

The compact forms give **unitary representations** and are used for internal
symmetries (gauge groups). Non-compact forms appear in spacetime symmetries
(Lorentz group, conformal group).

---

## 16. Lie Groups in Gauge Theory: The Full Picture

### 16.1 The Role of Each Structure

Here is how every aspect of a gauge Lie group appears in physics:

| Lie group structure | Physical meaning |
|---|---|
| Manifold dimension | Number of gauge bosons |
| Group multiplication | Composition of gauge transformations |
| Lie algebra basis $\{T_a\}$ | Internal quantum numbers (color, isospin, hypercharge) |
| Structure constants $f^{abc}$ | Self-interaction of gauge bosons (3-point vertex) |
| $f^{abc} f^{abd}$ | 4-point gauge boson vertex |
| Adjoint representation | How gauge bosons carry their own charge |
| Fundamental representation | How matter fields (quarks, leptons) couple |
| Killing form | Kinetic term for gauge field: $\mathrm{tr}(F^2)$ |
| Compact real form | Unitary, bounded representations (probability conserved) |
| $\pi_1(G)$ | Charge quantization, magnetic monopoles |
| $\pi_3(G)$ | Instantons, topological charge |
| Maximal torus | Conserved charges (Cartan generators) |
| Roots | Particle multiplet structure |
| Weights | Quantum numbers of matter particles |
| Casimir operators | Gauge-invariant mass terms |

### 16.2 $U(1)$ in This Picture

$U(1)$ as a gauge group:
- Manifold: $S^1$ (1-dimensional)
- Dimension: $1$ $\to$ $1$ gauge boson (the photon)
- Lie algebra: $i\mathbb{R}$, trivial bracket $\to$ $[A_\mu, A_\nu] = 0$ $\to$ no photon self-coupling
- Adjoint rep: trivial $\to$ photon is neutral (carries no electric charge)
- $\pi_1(U(1)) = \mathbb{Z}$ $\to$ electric charge quantization
- Representations $\rho_n$: $e^{i\theta} \mapsto e^{in\theta}$ $\to$ charge $n \in \mathbb{Z}$ for each particle
- Compact $\to$ unitary representations $\to$ probability conserved
- Killing form: $B(X, Y) = 0$ (1-dim, trivial) $\to$ kinetic term is $\int F_{\mu\nu} F^{\mu\nu} d^4 x$

### 16.3 Why the Gauge Group Determines the Theory

Given a gauge group $G$, the entire Yang-Mills theory is fixed:

1. **Step 1:** Choose $G$
2. **Step 2:** The Lie algebra $\mathfrak{g}$ gives the gauge bosons (one per generator)
3. **Step 3:** The structure constants $f^{abc}$ determine the self-coupling
4. **Step 4:** The Killing form determines the kinetic term $\mathrm{tr}(F^2)$
5. **Step 5:** Choose representations $\rho$ for matter fields
6. **Step 6:** The covariant derivative $D_\mu = \partial_\mu + g \rho(A_\mu)$ determines matter coupling
7. **Step 7:** The full Lagrangian $\mathcal{L} = -\mathrm{tr}(F^2)/4g^2 + \bar\psi(i \not{D} - m)\psi$ is completely specified

For the Standard Model:

$$
G = U(1)_Y \times SU(2)_L \times SU(3)_c
$$

Matter representations chosen to match experimental particle content
$\to$ Every interaction of every known particle is fixed.

### 16.4 The Non-Abelian Obstruction in Transformer Attention

Connecting to the gauge theory tutorial: for multi-head attention, the force on
hidden state $h$ is:

$$
F^{(h)}(h) = \sum_j \mathrm{softmax}\left(\frac{Q_h h \cdot K_h h_j}{\sqrt{d_h}}\right) V_h h_j \qquad \text{(head } h \text{ contribution)}
$$

Total force: $F(h) = \sum_h F^{(h)}(h) = \sum_h A^{(h)}(h) h$ (schematically).

Each head contributes a connection $A^{(h)}$ valued in a $d_h$-dimensional subspace.
After the output projection $W_O$ mixes all heads:

$$
A_{\text{total}}(h) = W_O ( A^{(1)}(h) \oplus \cdots \oplus A^{(H)}(h) ) W_O^{\dagger} \in \mathrm{End}(\mathbb{R}^d)
$$

The curvature:

$$
F_{ij} = (\partial A / \partial h)_{ij} - (\partial A / \partial h)_{ji} + [A_i, A_j]
$$

contains cross-head commutators $[A^{(h)}, A^{(h')}]$ that are non-zero because
the head subspaces, after $W_O$ mixing, are not orthogonal.

This is a **non-abelian gauge obstruction** in exactly the Yang-Mills sense:
just as $SU(N)$ gauge theory has non-zero $[A_\mu, A_\nu]$ that prevents the gauge
potential from being a gradient, the multi-head attention connection has
non-zero cross-head commutators that prevent any scalar potential from
reproducing the full attention dynamics.

### 16.5 Statement Decoder: The Full Logical Chain

The statement

> *"cross-head commutators $[A^{(h)}, A^{(h')}]$ generate non-abelian curvature
> that obstructs any scalar potential on hidden-state space, regardless of capacity"*

can now be unpacked term by term, with every concept traced to its definition
in this tutorial and the Gauge Theory Tutorial.

---

**"cross-head"**

Each attention head $h = 1, \ldots, H$ produces a force on the hidden state.
The word "cross-head" means we are looking at the **interaction between two
different heads $h \neq h'$** — not a single head acting on itself.
Within a single head, $[A^{(h)}, A^{(h)}] = 0$ always (bracket of anything with
itself vanishes, by §4.1a). The interesting structure comes from pairs $h \neq h'$.

---

**"commutators $[A^{(h)}, A^{(h')}]$"**

$A^{(h)}$ is the **connection 1-form** (gauge potential) contributed by head $h$.
After the output projection $W_O$ recombines all heads, $A^{(h)}$ becomes a
**$d \times d$ matrix** acting on the full hidden-state space $\mathbb{R}^d$.

The commutator (§4.1a) is:

$$
[A^{(h)}, A^{(h')}] = A^{(h)} A^{(h')} - A^{(h')} A^{(h)}
$$

This is a matrix: the product of two $d \times d$ matrices minus the product in
reversed order. It measures how much the force contributions of the two
heads fail to commute as linear operators on $\mathbb{R}^d$.

For the commutator to vanish, $A^{(h)}$ and $A^{(h')}$ must commute as matrices:
$A^{(h)} A^{(h')} = A^{(h')} A^{(h)}$. This would require very special alignment
of the weight matrices $W_Q^h, W_K^h, W_V^h, W_O$ — it is not the generic case.

---

**"generate non-abelian curvature"**

The full connection is $A = \sum_h A^{(h)}$. Its curvature (§7.1–7.2 of the
Gauge Theory Tutorial) is:

$$
F = \partial A - \partial A + [A, A]
$$

- The first two terms give the abelian part (would exist even for a single head).
- $[A, A]$ is the non-abelian part:

$$
[A, A] = \sum_{h \neq h'} [A^{(h)}, A^{(h')}] + \sum_h [A^{(h)}, A^{(h)}]
       = \sum_{h \neq h'} [A^{(h)}, A^{(h')}] \quad \text{(the cross-head commutators)}
$$

"Non-abelian curvature" specifically refers to the $[A, A]$ term — the part of
$F$ that exists **because the gauge group is non-abelian** (§1.4 of this tutorial,
§8.2 of the Gauge Theory Tutorial).

When the cross-head commutators $[A^{(h)}, A^{(h')}]$ are non-zero:

$$
[A, A] = \sum_{h \neq h'} [A^{(h)}, A^{(h')}] \neq 0
\quad \Rightarrow \quad F_{\mu\nu} \neq 0 \text{ (non-zero curvature)}
\quad \Rightarrow \quad G \text{ is effectively non-abelian}
$$

---

**"obstructs any scalar potential"**

A **scalar potential** is any smooth function $V : \mathbb{R}^d \to \mathbb{R}$. If a force field
could be written $F = -\nabla V$, it would be **conservative** (path-independent work,
zero net work around any closed loop).

The **obstruction theorem** (§7.4 of the Gauge Theory Tutorial, proved via
Clairaut's theorem) states:

$$
F = -\nabla V \text{ for some scalar } V \iff F_{\mu\nu} = 0 \text{ everywhere}
$$

Contrapositive:

$$
F_{\mu\nu} \neq 0 \implies F \neq -\nabla V \text{ for ANY scalar } V
$$

The proof is an algebraic identity: if $V$ is any smooth function, then

$$
(\partial_i \partial_j V) - (\partial_j \partial_i V) = 0 \quad \text{(Clairaut/Schwarz theorem: mixed partials commute)}
$$

So the curvature of $-\nabla V$ is always exactly zero — not approximately zero,
not zero for large enough $V$ — **exactly zero**, by pure calculus.

Since $F_{\mu\nu} \neq 0$ (from the cross-head commutators), and any scalar $V$ gives
zero curvature, we have:

$$
F \neq -\nabla V \text{ for ANY scalar } V
$$

This is the obstruction.

---

**"on hidden-state space"**

The base manifold is the hidden-state space $\mathbb{R}^d$ — the $d$-dimensional real
vector space in which transformer hidden states live. The connection $A$ and
its curvature $F$ are defined on this space (§16.1 of the Gauge Theory Tutorial).

---

**"regardless of capacity"**

This phrase means: the failure is not about $V_\psi$ being too small, too shallow,
or having too few parameters. Even if you replace $V_\psi$ with:
- A deeper neural network
- An infinitely wide network
- A universal approximator
- The exact analytical solution to any optimization problem

$\ldots$ it still cannot work. The reason:

Every smooth scalar $V$ — simple or complex — satisfies:

$$
\text{Curvature of } (-\nabla V) = \partial_i(-\partial_j V) - \partial_j(-\partial_i V) = 0 \quad \text{(Clairaut, always)}
$$

The attention force satisfies:

$$
\text{Curvature of } F = [A, A] \neq 0 \quad \text{(from cross-head commutators)}
$$

Therefore: $F \neq -\nabla V$ for ANY smooth $V$, no matter how expressive.

This is an **algebraic/topological obstruction** — it is in the same class as
saying "a continuous map from $S^2$ to $\mathbb{R}^2$ cannot be injective" (a topological fact
that holds regardless of how cleverly you choose the map). The capacity of $V_\psi$
is simply irrelevant because the problem is structural, not quantitative.

---

**The Complete Logical Chain**

1. Multi-head attention has $H \geq 2$ heads.
2. Each head $h$ contributes a matrix-valued connection $A^{(h)} \in \mathrm{End}(\mathbb{R}^d)$.
3. After $W_O$ mixing, different heads act in overlapping subspaces.
4. Cross-head commutators $[A^{(h)}, A^{(h')}] = A^{(h)} A^{(h')} - A^{(h')} A^{(h)} \neq 0$ (because matrices acting in overlapping subspaces generically don't commute).
5. $[A^{(h)}, A^{(h')}] \neq 0 \implies [A, A] = \sum_{h \neq h'} [A^{(h)}, A^{(h')}] \neq 0$ (non-abelian).
6. $F_{\mu\nu} = \partial A - \partial A + [A, A] \neq 0$ (non-zero curvature).
7. Obstruction theorem (Clairaut): any scalar $V$ gives zero curvature. $F_{\mu\nu} \neq 0$ but $\text{Curvature}(-\nabla V) = 0$ $\implies F \neq -\nabla V$.
8. No scalar potential of any form or capacity can represent $F$.
9. The shared-$V_\psi$ test fails for GPT-2 middle layers ($R^2 = 0.04$–$0.20$) regardless of $V_\psi$ architecture or parameter count.

**Where each step is defined:**

| Step | Definition location |
|---|---|
| Connection $A^{(h)}$ | Gauge Theory Tutorial §5.2, §16.1 |
| Matrix commutator $[\cdot, \cdot]$ | Lie Groups Tutorial §4.1a |
| Abelian vs. non-abelian | Lie Groups Tutorial §1.4 |
| Curvature $F = dA + [A, A]$ | Gauge Theory Tutorial §7.1 |
| $[A, A] =$ non-abelian part | Gauge Theory Tutorial §7.2 |
| Obstruction theorem | Gauge Theory Tutorial §7.4 |
| "Regardless of capacity" | Clairaut's theorem (calculus) |
| Hidden-state space = base manifold | Gauge Theory Tutorial §16.1 |
| Cross-head commutators in attention | Gauge Theory Tutorial §16.2 |

---

## 17. Reference: Key Facts About $U(1)$, $SU(2)$, $SU(3)$

### Quick Reference Table

| Property | $U(1)$ | $SU(2)$ | $SU(3)$ |
|---|---|---|---|
| Manifold | $S^1$ | $S^3$ | — (8d manifold) |
| Dimension | $1$ | $3$ | $8$ |
| Rank | $1$ | $1$ | $2$ |
| Compact? | Yes | Yes | Yes |
| Simply connected? | No | Yes | Yes |
| $\pi_1$ | $\mathbb{Z}$ | $0$ | $0$ |
| $\pi_3$ | $0$ | $\mathbb{Z}$ | $\mathbb{Z}$ |
| Lie algebra | $\mathfrak{u}(1) = i\mathbb{R}$ | $\mathfrak{su}(2)$ | $\mathfrak{su}(3)$ |
| Generators | $1$ (trivial) | Pauli $\sigma_a/2$ | Gell-Mann $\lambda_a/2$ |
| Structure constants | $0$ | $\epsilon^{abc}$ | $f^{abc}$ |
| Adjoint rep dim | $1$ | $3$ | $8$ |
| Fundamental rep dim | $1$ | $2$ | $3$ |
| Abelian? | Yes | No | No |
| Physical role | Electromagnetism | Weak force | Strong force |
| Gauge bosons | Photon (1) | $W^{\pm}, Z$ (3) | Gluons (8) |

### Key Formulas

**$U(1)$:**

$$
\text{Group:} \quad e^{i\theta} \cdot e^{i\varphi} = e^{i(\theta + \varphi)}
$$

$$
\text{Algebra:} \quad \mathfrak{u}(1) = i\mathbb{R}, \qquad [i\alpha, i\beta] = 0
$$

$$
\text{Exp map:} \quad \exp(i\theta) = e^{i\theta}, \text{ surjective}
$$

$$
\text{Reps:} \quad \rho_n(e^{i\theta}) = e^{in\theta}, \quad n \in \mathbb{Z}
$$

$$
\text{Cover:} \quad \mathbb{R} \to U(1) = \mathbb{R}/2\pi\mathbb{Z}
$$

**$SU(2)$:**

$$
\text{Group:} \quad
\begin{pmatrix} a & -\bar{b} \\ b & \bar{a} \end{pmatrix}
\begin{pmatrix} c & -\bar{d} \\ d & \bar{c} \end{pmatrix}
= \begin{pmatrix} ac - \bar{b} d & -a\bar{d} - \bar{b}\bar{c} \\ bc + \bar{a} d & -b\bar{d} + \bar{a}\bar{c} \end{pmatrix}
$$

$$
\text{Algebra:} \quad [\sigma_a/2, \sigma_b/2] = i \epsilon_{abc} \sigma_c/2
$$

$$
\text{Exp map:} \quad \exp(i\theta \hat{n} \cdot \sigma/2) = \cos(\theta/2) I + i \sin(\theta/2) \hat{n} \cdot \sigma
$$

$$
\text{Reps:} \quad \text{Labeled by } j = 0, 1/2, 1, \ldots, \text{ dimension } 2j + 1
$$

$$
\text{Cover:} \quad SU(2) \to SO(3) = SU(2)/\mathbb{Z}_2 \quad (2{:}1)
$$

**$SU(3)$:**

$$
\text{Generators:} \quad T_a = \lambda_a/2, \quad a = 1, \ldots, 8
$$

$$
\text{Algebra:} \quad [T_a, T_b] = i f^{abc} T_c
$$

$$
\text{Cartan:} \quad T_3 = \lambda_3/2, \quad T_8 = \lambda_8/2 \quad \text{(diagonal generators)}
$$

$$
\text{Reps:} \quad \text{Labeled by } (p, q), \quad \dim = (p+1)(q+1)(p+q+2)/2
$$

$$
\text{Fundamental:} \quad (1, 0) = 3 \text{ (quarks)}, \quad (0, 1) = \bar{3} \text{ (antiquarks)}
$$

$$
\text{Adjoint:} \quad (1, 1) = 8 \text{ (gluons)}
$$

### Recommended References

- **Hall, "Lie Groups, Lie Algebras, and Representations"** — the best mathematical
  treatment for this level; clear proofs, good exercises, covers all the material here.

- **Fulton & Harris, "Representation Theory: A First Course"** — excellent for the
  representation theory chapters (§11–13 above). Very concrete.

- **Bröcker & tom Dieck, "Representations of Compact Lie Groups"** — comprehensive
  reference for compact groups, Peter-Weyl, maximal tori.

- **Adams, "Lectures on Lie Groups"** — short, elegant, focuses on the topology
  (homotopy groups, Bott periodicity).

- **Georgi, "Lie Algebras in Particle Physics"** — physics perspective, very concrete,
  excellent for the gauge theory applications in §16.

- **Nakahara, "Geometry, Topology and Physics"** — Chapters 5 and 10 cover fiber
  bundles and gauge theory at exactly the level needed to connect this tutorial to
  the Gauge Theory Tutorial.

---

*Tutorial version: April 2026.*
*Designed to complement the Gauge Theory Tutorial.*
*Pitched at graduate level with group theory and calculus of variations prerequisites.*
