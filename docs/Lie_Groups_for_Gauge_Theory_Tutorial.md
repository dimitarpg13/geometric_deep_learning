# Lie Groups for Gauge Theory
### A Graduate Tutorial from Matrix Groups to the Structure of U(1), SU(2), SU(3)

**Prerequisites assumed:** Linear algebra (eigenvalues, matrix exponential, determinants),
multivariable calculus, basic topology (open sets, continuity, compactness),
complex numbers and complex analysis basics.

**Goal:** By the end of this tutorial you will understand U(1) as a Lie group from
every angle — as a manifold, as a matrix group, as an abstract group, as a topological
space, and as the gauge group of electromagnetism. You will understand how U(1) relates
to SU(2) and SU(3) and why the structure of these groups completely determines the
form of the fundamental forces.

---

## Table of Contents

1. [What Is a Lie Group? The Three Faces](#1-what-is-a-lie-group)
   - [1.4 Abelian and Non-Abelian Lie Groups](#14-abelian-and-non-abelian-lie-groups)
2. [The Simplest Examples: ℝ, S¹, and U(1)](#2-simplest-examples)
3. [Matrix Lie Groups: GL, SL, O, U, S](#3-matrix-lie-groups)
4. [The Lie Algebra: Tangent Space at the Identity](#4-lie-algebra)
   - [4.0 Unpacking "Tangent Space at the Identity"](#40-unpacking-tangent-space-at-the-identity)
   - [4.1a The Lie Bracket in Detail: Bilinear, Antisymmetric, Jacobi](#41a-the-lie-bracket-in-detail)
5. [The Exponential Map](#5-exponential-map)
6. [The Adjoint Representation](#6-adjoint-representation)
7. [U(1) in Full Detail](#7-u1-full-detail)
8. [SU(2): The Three-Sphere and Rotations](#8-su2)
9. [SO(3) vs SU(2): Covering Groups and Spinors](#9-so3-vs-su2)
10. [SU(3) and the Gell-Mann Matrices](#10-su3)
11. [Representations of Lie Groups](#11-representations)
12. [The Peter-Weyl Theorem and Harmonic Analysis](#12-peter-weyl)
13. [Roots, Weights, and the Classification of Simple Lie Algebras](#13-classification)
14. [Homomorphisms, Subgroups, and Quotients](#14-homomorphisms)
15. [The Topology of Lie Groups](#15-topology)
16. [Lie Groups in Gauge Theory: The Full Picture](#16-gauge-theory-connection)
17. [Reference: Key Facts About U(1), SU(2), SU(3)](#17-reference)

---

## 1. What Is a Lie Group? The Three Faces

### 1.1 The Three Structures

A **Lie group** is an object with three simultaneous structures that are mutually compatible:

```
Face 1: A GROUP
  - Set G with binary operation · : G × G → G
  - Associativity: (ab)c = a(bc)
  - Identity element: e ∈ G with ge = eg = g
  - Inverses: for each g ∈ G, ∃ g⁻¹ with gg⁻¹ = e

Face 2: A SMOOTH MANIFOLD
  - G is a smooth manifold of dimension n
  - Has an atlas of coordinate charts
  - Transition maps are smooth (C^∞)

Face 3: COMPATIBILITY
  - The group operations are smooth maps:
    μ: G × G → G,   (g, h) ↦ gh   (smooth)
    ι: G → G,       g ↦ g⁻¹        (smooth)
```

The compatibility condition is what makes Lie groups so powerful: the algebraic
structure (group) and the geometric structure (manifold) are woven together.
Every theorem about one structure applies to the other.

### 1.2 Why This Structure Matters for Physics

In physics, we use Lie groups as symmetry groups. The manifold structure means:
- We can talk about **continuous** symmetries (infinitesimal transformations)
- We can differentiate transformations → Lie algebra
- We can integrate → action of the group on fields

The group structure means:
- Compositions of symmetries are symmetries
- Every symmetry has an inverse
- The identity transformation is included

Without the smooth structure, we would only have discrete symmetries (like crystal
symmetries). The smoothness is what gives us conservation laws via Noether's theorem.

### 1.3 The Lie Group–Lie Algebra Correspondence

The most important theorem about Lie groups:

```
Every Lie group G has an associated Lie algebra 𝔤 = T_e G
(the tangent space at the identity).

This correspondence is:
  - Local: 𝔤 determines G near the identity
  - Global (for simply connected G): the exponential map exp: 𝔤 → G
    gives a local diffeomorphism near 0 ∈ 𝔤

Simply connected Lie groups ←1:1→ Lie algebras (Lie's Third Theorem)
```

The Lie algebra is a vector space with a bilinear antisymmetric bracket [·,·]
satisfying the Jacobi identity — each of these terms is defined precisely in
§4.1a. It is often much easier to work with (it is linear!) and then
"exponentiate" back to the group.

### 1.4 Abelian and Non-Abelian Lie Groups

These two terms appear throughout gauge theory and are worth defining precisely
at the outset, because the entire difference between electromagnetism (simple,
linear) and the strong/weak forces (complex, nonlinear) reduces to them.

**Definition (Abelian Lie group).** A Lie group G is **abelian** (or commutative)
if group multiplication commutes:

```
gh = hg    for ALL g, h ∈ G
```

**Definition (Non-abelian Lie group).** A Lie group G is **non-abelian**
(or non-commutative) if there EXIST g, h ∈ G such that:

```
gh ≠ hg    (at least one pair of elements that do not commute)
```

The word "non-abelian" means the group does not have the abelian property —
it is not required to fail commutativity everywhere, only somewhere.

**The Lie algebra test.** Because the Lie bracket [X,Y] measures the
infinitesimal failure of group elements to commute (see §4.1a), a connected
Lie group is abelian if and only if its Lie algebra satisfies:

```
G abelian   ⟺   [X, Y] = 0 for ALL X, Y ∈ 𝔤
```

This means: checking commutativity of the group reduces to checking whether
all brackets vanish in the Lie algebra — a purely linear algebra problem.

**Standard examples:**

| Group | Abelian? | Why |
|---|---|---|
| U(1) = {e^{iθ}} | **Yes** | e^{iα}e^{iβ} = e^{i(α+β)} = e^{iβ}e^{iα} |
| (ℝ, +) | **Yes** | x + y = y + x |
| SU(2) | **No** | σ₁σ₂ = iσ₃ ≠ -iσ₃ = σ₂σ₁ (Pauli matrices don't commute) |
| SU(3) | **No** | Most generators don't commute |
| GL(n), n≥2 | **No** | Matrix multiplication is not commutative |

**Why this matters for gauge theory:**

```
G = U(1)  (abelian):
  [A_μ, A_ν] = 0   for all μ, ν
  F_μν = ∂_μA_ν - ∂_νA_μ + [A_μ, A_ν] = ∂_μA_ν - ∂_νA_μ   (linear in A)
  → Photon does not self-interact
  → Maxwell equations are linear
  → Scalar potential CAN exist (and does: the electromagnetic potential)

G = SU(2) or SU(3)  (non-abelian):
  [A_μ, A_ν] ≠ 0   generically
  F_μν = ∂_μA_ν - ∂_νA_μ + [A_μ, A_ν]   (nonlinear in A)
  → Gauge bosons self-interact (W bosons carry weak charge, gluons carry color)
  → Yang-Mills equations are nonlinear
  → No scalar potential can exist that reproduces F_μν
```

The last bullet is the **non-abelian gauge obstruction**: as soon as [A_μ, A_ν] ≠ 0,
the curvature F_μν ≠ 0, and by the obstruction theorem (§7.4 in the Gauge Theory
Tutorial), no scalar potential V — of any form, any complexity — can satisfy
F = -∇V. This is an algebraic identity, not an approximation limit.

---

## 2. The Simplest Examples: ℝ, S¹, and U(1)

### 2.1 The Real Line (ℝ, +)

The real line with addition is the simplest Lie group:

```
Manifold:  ℝ¹ (a 1-dimensional smooth manifold, the real line)
Group op:  (x, y) ↦ x + y
Identity:  0
Inverse:   x ↦ -x
Lie algebra: ℝ (with trivial bracket [x,y] = 0)
```

The exponential map exp: ℝ → ℝ sends the Lie algebra to itself.
(For additive groups, exp is literally the identity map on ℝ when viewed abstractly,
though numerically exp(t) = eᵗ takes ℝ to ℝ₊ as a map of manifolds.)

### 2.2 The Circle S¹ as a Lie Group

The unit circle in ℝ²:

```
S¹ = {(x,y) ∈ ℝ² : x² + y² = 1}
```

with multiplication defined by angle addition:

```
(x₁,y₁) · (x₂,y₂) = (x₁x₂ - y₁y₂, x₁y₂ + y₁x₂)
```

(this is just complex multiplication of unit complex numbers.)

Written in terms of angles θ ∈ [0, 2π):

```
θ₁ · θ₂ = θ₁ + θ₂ (mod 2π)
```

So S¹ with this multiplication is the group ℝ/2πℤ (the real line modulo
the integers scaled by 2π).

### 2.3 U(1): The Complex Definition

The **unitary group U(1)** is the group of 1×1 unitary matrices — complex numbers
of modulus 1:

```
U(1) = {z ∈ ℂ : |z| = 1} = {e^{iθ} : θ ∈ ℝ}
```

with multiplication = complex multiplication.

**Claim: U(1) ≅ S¹ as Lie groups.**

The isomorphism is: e^{iθ} ↔ (cos θ, sin θ). Under this identification:
- The group operation (complex multiplication) corresponds to angle addition
- The manifold structure (circle) is the same
- The smooth structure is the same

So U(1), S¹, and ℝ/2πℤ are three names for the same Lie group.

### 2.4 Why U(1) Is Called a "Circle Group"

Topologically, U(1) ≅ S¹. This has profound consequences:

```
π₀(U(1)) = 0     (connected: one piece)
π₁(U(1)) = ℤ     (fundamental group: loops wind around the circle by integer amounts)
π₂(U(1)) = 0     (no non-trivial 2-spheres)
```

The fact that π₁(U(1)) = ℤ means there are topologically distinct ways to map a
loop into U(1) — classified by the winding number. This is the origin of
**magnetic monopoles** in physics (magnetic charge is quantized because it is
measured by this ℤ).

---

## 3. Matrix Lie Groups: GL, SL, O, U, S

### 3.1 The General Linear Group GL(n, 𝔽)

The group of invertible n×n matrices over a field 𝔽 (ℝ or ℂ):

```
GL(n, ℝ) = {A ∈ Mat(n×n, ℝ) : det A ≠ 0}
GL(n, ℂ) = {A ∈ Mat(n×n, ℂ) : det A ≠ 0}
```

- Dimension (as a manifold): n² over ℝ, 2n² over ℝ (= n² over ℂ)
- GL(n) is an open subset of Mat(n×n) — the condition det A ≠ 0 is open
- Not compact (matrices can have arbitrarily large or small entries)
- Not connected: GL(n, ℝ) has two connected components (det > 0 and det < 0)

### 3.2 The Special Linear Group SL(n)

```
SL(n, ℝ) = {A ∈ GL(n, ℝ) : det A = 1}
SL(n, ℂ) = {A ∈ GL(n, ℂ) : det A = 1}
```

- Dimension: n²-1 (one constraint det = 1 reduces dimension by 1)
- SL(n) is a closed subgroup of GL(n) (det = 1 is a closed condition)
- SL(n) is connected

### 3.3 The Orthogonal and Special Orthogonal Groups

```
O(n) = {A ∈ GL(n, ℝ) : A^T A = I}   (orthogonal matrices)
SO(n) = {A ∈ O(n) : det A = 1}        (special orthogonal = rotation matrices)
```

- O(n) has two connected components: det = +1 (rotations) and det = -1 (reflections)
- SO(n) = O(n) component containing identity = proper rotations
- Dimension: n(n-1)/2
- Compact: entries bounded by |A_{ij}| ≤ 1 from A^T A = I

Physical examples:
```
SO(2) ≅ U(1) = rotations in 2D  (dim = 1)
SO(3)         = rotations in 3D  (dim = 3)
SO(3,1)       = Lorentz group    (dim = 6, but non-compact due to boosts)
```

### 3.4 The Unitary and Special Unitary Groups

```
U(n) = {A ∈ GL(n, ℂ) : A† A = I}    (unitary matrices, A† = Ā^T)
SU(n) = {A ∈ U(n) : det A = 1}        (special unitary)
```

- Dimension: n² (as a real manifold; n² real constraints from U†U = I reduce 2n² by n²)
  More precisely: dim U(n) = n², dim SU(n) = n²-1
- Both compact: unitary matrices have |eigenvalues| = 1, so entries bounded
- Both connected

**Explicit matrix forms:**

```
U(1) = {e^{iθ} : θ ∈ [0,2π)} — 1×1 unitary matrices, dim=1

SU(2) = {(a   -b̄) : a,b ∈ ℂ, |a|²+|b|² = 1} — dim=3
         (b    ā)

U(2) = {A ∈ GL(2,ℂ) : A†A = I} — dim=4

SU(3) = {A ∈ GL(3,ℂ) : A†A=I, det A=1} — dim=8
```

### 3.5 The Relationship Between These Groups

```
SU(n) ⊂ U(n) ⊂ GL(n, ℂ)
  ↕
U(n) ≅ SU(n) × U(1) / ℤₙ   (almost a direct product)
U(n) ≅ (SU(n) × U(1)) / {(e^{2πik/n} I, e^{-2πik/n}) : k = 0,...,n-1}
```

The exact sequence:

```
1 → SU(n) → U(n) →^det U(1) → 1
```

The determinant map U(n) → U(1) is surjective with kernel SU(n), so U(n)/SU(n) ≅ U(1).
This means U(n) is "SU(n) times U(1)" with a discrete identification.

---

## 4. The Lie Algebra: Tangent Space at the Identity

### 4.0 Unpacking "Tangent Space at the Identity"

Before giving the formal definition, it is worth building up what this phrase means
from scratch — because it contains three separate ideas packed into four words.

---

#### Step 1: What Is a Tangent Vector?

Start with a smooth manifold M (think of a surface in ℝ³, or a more abstract space
like a Lie group). A **tangent vector** at a point p ∈ M is the velocity vector of
a smooth curve passing through p.

Concretely: take any smooth curve γ: (-ε, ε) → M with γ(0) = p.
The tangent vector to γ at p is its velocity:

```
γ̇(0) = dγ/dt|_{t=0}   ∈  T_p M
```

Two curves that pass through p with the same velocity define the same tangent vector.
So a tangent vector is really an **equivalence class of curves** through p,
where two curves are equivalent if they have the same first-order behavior at p.

**In ℝⁿ:** Every tangent vector at every point is literally a vector in ℝⁿ.
The tangent space T_p ℝⁿ = ℝⁿ for all p.

**On a sphere S²:** At a point p ∈ S², the tangent space T_p S² is the plane
tangent to the sphere at p — a 2-dimensional flat vector space sitting in ℝ³.
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

#### Step 2: The Tangent Space T_p M

The **tangent space** at p is the collection of all tangent vectors at p:

```
T_p M = {γ̇(0) : γ smooth curve in M with γ(0) = p}
```

Key facts:
- T_p M is a **vector space** of dimension equal to dim(M)
- Each point has its own tangent space — T_p M and T_q M are different spaces for p ≠ q
- Tangent vectors at p "live at p" — you cannot directly compare tangent vectors
  at different points without a connection (which is exactly what gauge theory provides)

**In coordinates:** If (x¹,...,xⁿ) are local coordinates near p, then
T_p M has basis {∂/∂x¹|_p, ..., ∂/∂xⁿ|_p}.
A general tangent vector is v = vⁱ ∂/∂xⁱ|_p for some real numbers v¹,...,vⁿ.

**For a matrix Lie group G ⊂ GL(n):**
G is a subset of the space Mat(n×n) ≅ ℝⁿ². The tangent space at any point
g ∈ G consists of the velocity vectors of smooth curves in G passing through g:

```
T_g G = {γ̇(0) : γ: (-ε,ε) → G smooth, γ(0) = g}
       ⊂ Mat(n×n)
```

These are matrices — specifically the matrices that are "tangent to G at g."

---

#### Step 3: Why the Identity Is Special

A Lie group G has a distinguished point: the **identity element** e (the matrix I
for matrix groups). The tangent space at the identity:

```
T_e G = {γ̇(0) : γ smooth curve in G, γ(0) = e}
```

is special for one reason: the group structure lets you **move any tangent space
back to T_e G** via left or right multiplication.

For any g ∈ G, left multiplication L_g: G → G defined by L_g(h) = gh is a
diffeomorphism. Its derivative maps tangent spaces:

```
(dL_g)_e : T_e G → T_g G
```

This map is a linear isomorphism — it identifies T_e G with T_g G for every g.
So all tangent spaces are isomorphic to T_e G, and T_e G is the canonical
"home base" for all tangent vectors.

**Physical interpretation:** Every infinitesimal transformation in the group can
be "brought home" to the identity. The Lie algebra T_e G collects all infinitesimal
group transformations in one place.

---

#### Step 4: Computing T_e G for Matrix Groups

For a matrix Lie group G ⊂ GL(n), a curve γ through the identity means:

```
γ: (-ε, ε) → G ⊂ GL(n)    with γ(0) = I
```

The tangent vector at the identity is:

```
γ̇(0) = dγ/dt|_{t=0}   ∈ Mat(n×n)
```

So T_e G consists of all matrices X that are the velocity of some smooth curve in G
starting at I.

**For U(1):** Curves through the identity in U(1) = {e^{iθ}} look like γ(t) = e^{iα(t)}
with α(0) = 0. The velocity is:

```
γ̇(0) = i α̇(0)   ∈ iℝ
```

Any purely imaginary number is achievable (choose α(t) = ct for any c ∈ ℝ).
So T_e U(1) = iℝ — the purely imaginary numbers. ✓

**For SU(2):** Curves through I in SU(2) ⊂ Mat(2×2,ℂ). Write γ(t) = I + tX + O(t²).
For γ(t) ∈ SU(2) we need:
- γ(t)†γ(t) = I:   differentiating → X† + X = 0   (skew-Hermitian)
- det γ(t) = 1:    differentiating → tr X = 0       (traceless)

So T_e SU(2) = {X ∈ Mat(2×2,ℂ) : X† = -X, tr X = 0} = 𝔰𝔲(2). ✓

This computation — differentiate the defining equations of G at the identity — is the
general recipe for finding the Lie algebra of any matrix Lie group.

---

#### Step 5: The Lie Algebra IS T_e G, as a Vector Space

The tangent space T_e G is a vector space of dimension n = dim(G):

```
Dimension:     dim(T_e G) = dim(G) as a manifold

Examples:
  T_e U(1) = iℝ          dim = 1
  T_e SU(2) = 𝔰𝔲(2)     dim = 3
  T_e SU(3) = 𝔰𝔲(3)     dim = 8
  T_e GL(n) = Mat(n×n)   dim = n²
```

The Lie algebra 𝔤 is T_e G **as a vector space**, plus the additional structure of
the Lie bracket [·,·] coming from the group multiplication.

The bracket does not come from the tangent space structure alone — it comes from
differentiating the group commutator ghg⁻¹h⁻¹ twice at the identity. Concretely:

```
If γ₁(t) = e^{tX} and γ₂(t) = e^{tY} are curves through I, then:

e^{tX} e^{tY} e^{-tX} e^{-tY} = I + t²[X,Y] + O(t³)
```

The bracket [X,Y] = XY - YX measures the second-order failure to commute.
It is invisible at first order (the linear term vanishes), and it first appears at
second order — which is why we differentiate the group commutator **twice**
to extract the Lie algebra structure.

---

#### The Complete Picture: Three Things at Once

```
T_e G  as a set         =  equivalence classes of smooth curves through e
T_e G  as a vector space =  infinitesimal directions you can move away from e in G
T_e G  as a Lie algebra  =  infinitesimal generators of G, with bracket [X,Y]
                             encoding the non-commutativity of group multiplication
```

The phrase "tangent space at the identity" refers to the first two.
The full Lie algebra 𝔤 = T_e G adds the third — the bracket — on top.

---

#### A Visual Summary for U(1)

```
The group U(1) = S¹:

          1 ∈ U(1)
          |
    ─────●─────────  ←  T_e U(1) = iℝ (vertical line through 1)
         |↑
         | γ̇(0) = iα̇(0)   (tangent vector = purely imaginary number)
         |
  γ(t) = e^{iα(t)}   (a curve in U(1) starting at e = 1)

The circle U(1) is 1-dimensional.
Its tangent space at e = 1 is a 1-dimensional real vector space: iℝ.
The "vertical line" is the Lie algebra 𝔲(1) = iℝ.
The exponential map wraps this line around the circle: e^{iθ} ↦ U(1).
```

---

### 4.1 Definition

The **Lie algebra** 𝔤 of a Lie group G is:

```
𝔤 = T_e G   (tangent space at the identity element e)
```

As a vector space, 𝔤 has dimension equal to the dimension of G (as a manifold).

But 𝔤 has more structure than just a vector space: it has a **Lie bracket** [·,·]:

```
[·,·]: 𝔤 × 𝔤 → 𝔤
```

satisfying:
```
1. Bilinearity:    [αX + βY, Z] = α[X,Z] + β[Y,Z]
2. Antisymmetry:   [X,Y] = -[Y,X]
3. Jacobi identity: [X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0
```

### 4.1a The Lie Bracket in Detail: Bilinear, Antisymmetric, Jacobi

The three axioms listed above — bilinearity, antisymmetry, and the Jacobi identity —
deserve a careful unpacking. They are not arbitrary: each encodes a precise geometric
or physical requirement.

#### What the bracket IS

The bracket [·,·] is a rule that takes **two vectors in 𝔤 and produces a third**:

```
[·,·] : 𝔤 × 𝔤 → 𝔤

[X, Y] = Z    where X, Y, Z ∈ 𝔤
```

For matrix Lie algebras — the case you encounter in all gauge theory — it is the
**matrix commutator**:

```
[X, Y] = XY - YX
```

The three axioms are then properties of this commutator.

#### Bilinear

"Bilinear" means **linear in each slot separately**, with the other held fixed.

**Linear in the first slot:**
```
[αX + βY, Z] = α[X, Z] + β[Y, Z]
```

**Linear in the second slot:**
```
[X, αY + βZ] = α[X, Y] + β[X, Z]
```

Together: scaling or summing generators works consistently in either argument.
This is the same structure as the determinant (bilinear in rows), the inner product,
or matrix multiplication — "bi-linear" simply means linear in both inputs.

**Verification for the matrix commutator:**
```
[αX + βY, Z] = (αX + βY)Z - Z(αX + βY)
             = αXZ + βYZ - αZX - βZY
             = α(XZ - ZX) + β(YZ - ZY)
             = α[X,Z] + β[Y,Z]   ✓
```

**Why bilinearity is required:** The bracket must respect the vector space structure
of 𝔤. If X generates a rotation by angle θ, then 2X should generate a rotation by
2θ, and the bracket should reflect this scaling consistently.

#### Antisymmetric

"Antisymmetric" means **swapping the two arguments negates the result**:

```
[X, Y] = -[Y, X]    for all X, Y ∈ 𝔤
```

**Immediate consequences:**

*The bracket of anything with itself is zero:*
```
[X, X] = -[X, X]    →    2[X, X] = 0    →    [X, X] = 0
```

*Over ℝ or ℂ (characteristic ≠ 2): "[X,X] = 0 for all X" and "[X,Y] = -[Y,X]"
are equivalent conditions.*

**Verification for the matrix commutator:**
```
[Y, X] = YX - XY = -(XY - YX) = -[X, Y]   ✓
```

**Why antisymmetry is required — the geometric meaning:**

Antisymmetry encodes the fact that the bracket measures **non-commutativity** —
how much two group transformations fail to commute. The Baker-Campbell-Hausdorff
formula makes this precise:

```
e^{tX} e^{tY} = e^{tY} e^{tX} · e^{t²[X,Y] + O(t³)}
```

The bracket [X,Y] is the leading correction when you reverse the order of two
transformations. Swapping X ↔ Y reverses which ordering is "first" and which is
"second", so [Y,X] = -[X,Y]. When [X,Y] = 0, the two transformations commute
exactly — there is no correction.

This is the mathematical explanation of a crucial physical fact:

```
U(1):    [iα, iβ] = 0    →  phase rotations always commute
                           →  photons do not interact with each other
                           →  Maxwell equations are linear

SU(2):   [Tₐ, Tᵦ] ≠ 0  →  weak isospin rotations do not commute
                           →  W bosons carry weak charge and self-interact
                           →  Yang-Mills equations are nonlinear
```

#### The Jacobi Identity

Together with bilinearity and antisymmetry, the bracket must satisfy:

```
[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0
```

This is the condition that makes [·,·] a **Lie bracket** rather than an arbitrary
bilinear antisymmetric map.

**Intuition:** Define ad_X(Y) = [X, Y] (the operator "bracket with X on the left").
The Jacobi identity says:

```
ad_X([Y, Z]) = [ad_X(Y), Z] + [Y, ad_X(Z)]
```

This is a **Leibniz rule** (product rule) for ad_X acting on the bracket.
It says: "differentiating the bracket is like differentiating a product."
This is necessary for the exponential map to be well-behaved and for the Lie algebra
to faithfully represent the group multiplication near the identity.

**Verification for matrix commutators:**
```
[X,[Y,Z]] = X(YZ-ZY) - (YZ-ZY)X = XYZ - XZY - YZX + ZYX
[Y,[Z,X]] =                       = YZX - YXZ - ZXY + XZY
[Z,[X,Y]] =                       = ZXY - ZYX - XYZ + YXZ

Sum: each of the 12 terms appears exactly twice with opposite signs → 0   ✓
```

#### The Complete Definition

A **Lie algebra** is a vector space 𝔤 over a field 𝔽 with a bracket [·,·]: 𝔤×𝔤→𝔤
satisfying:

```
L1. Bilinearity:      [αX + βY, Z] = α[X,Z] + β[Y,Z]
                      [X, αY + βZ] = α[X,Y] + β[X,Z]

L2. Antisymmetry:     [X, Y] = -[Y, X]

L3. Jacobi identity:  [X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0
```

#### A Concrete 3D Example: 𝔰𝔲(2)

The Lie algebra 𝔰𝔲(2) has basis {T₁, T₂, T₃} with brackets:

```
[T₁, T₂] = iT₃
[T₂, T₃] = iT₁
[T₃, T₁] = iT₂
```

Compactly: [Tₐ, Tᵦ] = i ε_{abc} Tc  where ε is the Levi-Civita symbol.

Check **bilinearity**:
```
[2T₁ + 3T₂, T₃] = 2[T₁,T₃] + 3[T₂,T₃]
                 = 2(-iT₂) + 3(iT₁)
                 = 3iT₁ - 2iT₂   ✓  (a linear combination of basis vectors)
```

Check **antisymmetry**:
```
[T₂, T₁] = -[T₁, T₂] = -iT₃   ✓
```

Check **Jacobi**:
```
[T₁,[T₂,T₃]] + [T₂,[T₃,T₁]] + [T₃,[T₁,T₂]]
= [T₁, iT₁]  + [T₂, iT₂]  + [T₃, iT₃]
= i[T₁,T₁]  + i[T₂,T₂]  + i[T₃,T₃]
= i·0 + i·0 + i·0 = 0   ✓    (bracket of anything with itself is zero)
```

#### Summary Table

| Word | Meaning | Formula |
|---|---|---|
| **Bracket** | Binary operation producing a third element | [X,Y] ∈ 𝔤 |
| **Bilinear** | Linear in each argument separately | [αX+βY, Z] = α[X,Z]+β[Y,Z] |
| **Antisymmetric** | Swap arguments → negate result | [X,Y] = -[Y,X] |
| **[X,X] = 0** | Consequence of antisymmetry | [X,X] = 0 for all X |
| **Jacobi identity** | Cyclic sum of nested brackets = 0 | [X,[Y,Z]]+[Y,[Z,X]]+[Z,[X,Y]]=0 |
| **For matrices** | All three hold for the commutator | [X,Y] = XY - YX |
| **For U(1)** | All brackets vanish (abelian) | [X,Y] = 0 for all X,Y |

---

### 4.2 Computing the Lie Algebra for Matrix Groups

For matrix Lie groups, the Lie algebra has a concrete description:

**Theorem.** For a matrix Lie group G ⊂ GL(n), the Lie algebra 𝔤 consists of all
matrices X such that e^{tX} ∈ G for all t ∈ ℝ.

The Lie bracket is simply the matrix commutator:

```
[X, Y] = XY - YX
```

This is automatic from the matrix structure and is the reason matrix Lie groups
are particularly tractable.

### 4.3 The Lie Algebras of the Classical Groups

**𝔤𝔩(n) = T_I GL(n):**
All n×n matrices (the identity is open in GL(n), so the tangent space is the full
matrix space):
```
𝔤𝔩(n, ℝ) = Mat(n×n, ℝ)   with bracket [X,Y] = XY - YX
𝔤𝔩(n, ℂ) = Mat(n×n, ℂ)
```

**𝔰𝔩(n) = T_I SL(n):**
The condition det(e^{tX}) = 1 differentiates to tr(X) = 0:
```
𝔰𝔩(n) = {X ∈ Mat(n×n) : tr X = 0}   (traceless matrices)
```
Proof: d/dt|_{t=0} det(e^{tX}) = tr(X) · det(I) = tr(X).

**𝔬(n) = T_I O(n):**
The condition (e^{tX})^T e^{tX} = I differentiates to X^T + X = 0:
```
𝔬(n) = {X ∈ Mat(n×n, ℝ) : X^T = -X}   (antisymmetric matrices)
```
Dimension: n(n-1)/2 ✓

**𝔲(n) = T_I U(n):**
The condition (e^{tX})† e^{tX} = I differentiates to X† + X = 0:
```
𝔲(n) = {X ∈ Mat(n×n, ℂ) : X† = -X}   (skew-Hermitian matrices)
```

**𝔰𝔲(n) = T_I SU(n):**
Both traceless AND skew-Hermitian:
```
𝔰𝔲(n) = {X ∈ Mat(n×n, ℂ) : X† = -X, tr X = 0}
```
Dimension: n²-1 ✓

### 4.4 The Lie Algebra of U(1)

```
U(1) = {e^{iθ} : θ ∈ ℝ}

Lie algebra 𝔲(1):
  d/dt|_{t=0} e^{itX} = iX must be in T_e U(1)
  e^{tX} ∈ U(1) for all t means |e^{tX}| = 1
  → X must be purely imaginary: X = iθ for some θ ∈ ℝ

𝔲(1) = iℝ = {iθ : θ ∈ ℝ}
```

As a vector space, 𝔲(1) ≅ ℝ (one-dimensional).
The Lie bracket: [iθ₁, iθ₂] = i²θ₁θ₂ - i²θ₂θ₁ = 0.

**U(1) has an abelian Lie algebra — all brackets vanish.**

This is the fundamental reason electromagnetism is simpler than the weak and strong
forces: the gauge group U(1) is abelian, so [A_μ, A_ν] = 0, and the photon
self-interaction term vanishes from the field strength tensor.

---

## 5. The Exponential Map

### 5.1 Definition and Convergence

The **exponential map** exp: 𝔤 → G is defined by:

```
exp(X) = e^X = Σ_{n=0}^∞ X^n / n! = I + X + X²/2! + X³/3! + ...
```

For matrix Lie groups, this is literally the matrix exponential, which converges
absolutely for all matrices X.

### 5.2 Key Properties

```
1. exp(0) = I                           (identity at origin)
2. exp((s+t)X) = exp(sX)exp(tX)        (one-parameter subgroup)
3. d/dt|_{t=0} exp(tX) = X             (X is the velocity at the identity)
4. exp(-X) = (exp X)⁻¹                 (inverse via negation)
5. det(exp X) = e^{tr X}               (Liouville formula)
6. If XY = YX: exp(X+Y) = exp(X)exp(Y) (Baker-Campbell-Hausdorff reduces here)
7. In general: exp(X)exp(Y) ≠ exp(X+Y) for non-commuting X, Y
```

Property 2 means t ↦ exp(tX) is a **one-parameter subgroup** of G — a smooth
homomorphism ℝ → G. Every one-parameter subgroup arises this way.

### 5.3 The Baker-Campbell-Hausdorff Formula

For non-commuting X, Y:

```
exp(X) exp(Y) = exp(X + Y + 1/2[X,Y] + 1/12([X,[X,Y]] - [Y,[X,Y]]) + ...)
```

The full series involves only nested Lie brackets. This shows that the group
multiplication near the identity is completely determined by the Lie bracket.

**Consequence for gauge theory:** The non-abelian structure [A_μ, A_ν] ≠ 0 is
precisely the BCH correction when composing gauge transformations. It is the
reason why Yang-Mills theory has self-interacting gauge bosons.

### 5.4 When Is exp Surjective?

- **For compact, connected Lie groups** (e.g., U(1), SU(n), SO(n)): exp is surjective.
  Every group element can be written as e^X for some X ∈ 𝔤.

- **For non-compact groups** (e.g., GL(n), SL(n)): exp is not surjective in general.
  Example: in SL(2,ℝ), the matrix [-1, 1; 0, -1] is not in the image of exp.

- **For U(1) specifically:**
  exp(iθ) = e^{iθ}, and every z = e^{iφ} ∈ U(1) is hit by iφ ∈ 𝔲(1).
  The map exp: 𝔲(1) = iℝ → U(1) is surjective but not injective:
  exp(iθ) = exp(i(θ + 2πn)) for any integer n.
  The kernel is 2πi ℤ ⊂ iℝ.

### 5.5 The Local Diffeomorphism Property

Near the identity, exp is a local diffeomorphism: there exists a neighborhood
U of 0 in 𝔤 and a neighborhood V of e in G such that exp: U → V is a
diffeomorphism with smooth inverse log: V → U.

This is proved by the inverse function theorem: (d exp)_0 = id.

This is why **working in the Lie algebra is the same as working near the identity
in the Lie group** — they are locally diffeomorphic.

---

## 6. The Adjoint Representation

### 6.1 The Adjoint Action of G on 𝔤

For g ∈ G, the conjugation map C_g: G → G defined by C_g(h) = ghg⁻¹ is a Lie
group automorphism (it preserves the group structure). Its derivative at the
identity gives a linear map on 𝔤:

```
Ad_g: 𝔤 → 𝔤
Ad_g(X) = (dC_g)_e(X)
```

For matrix groups: Ad_g(X) = gXg⁻¹.

**Key property:** Ad_g is a Lie algebra automorphism:
```
Ad_g([X,Y]) = [Ad_g(X), Ad_g(Y)]
```

The map Ad: G → GL(𝔤) defined by g ↦ Ad_g is a **representation** of G
(a Lie group homomorphism from G to GL of some vector space) — specifically
the **adjoint representation**.

### 6.2 The Adjoint Representation of the Lie Algebra

Differentiating Ad: G → GL(𝔤) at the identity gives the **adjoint representation
of the Lie algebra**:

```
ad: 𝔤 → 𝔤𝔩(𝔤)   (linear maps on 𝔤)
ad_X(Y) = [X, Y]
```

This is just the Lie bracket! The Jacobi identity becomes:
```
ad_X ∘ ad_Y - ad_Y ∘ ad_X = ad_{[X,Y]}   (Jacobi = representation property)
```

### 6.3 Why This Matters for Gauge Theory

In Yang-Mills theory, the field strength F_μν transforms in the adjoint representation:

```
F_μν → gF_μνg⁻¹ = Ad_g(F_μν)
```

This is exactly the adjoint action. The gauge bosons (which are sections of the
adjoint bundle) transform under Ad_g, not under the fundamental representation.

For U(1): Ad_g = identity (since U(1) is abelian: gXg⁻¹ = X for all g, X).
This is why the photon is neutral — the U(1) adjoint representation is trivial.

For SU(N): Ad_g is non-trivial. Gluons and W/Z bosons carry their own charge.

### 6.4 The Killing Form

The **Killing form** is the symmetric bilinear form on 𝔤:

```
B(X,Y) = tr(ad_X ∘ ad_Y)   (trace in the adjoint representation)
```

For matrix Lie algebras, this simplifies to:

```
B(X,Y) = 2N · tr(XY)   (for 𝔰𝔲(N) in the fundamental representation)
```

The Killing form is:
- **Negative definite** for compact semisimple Lie algebras (like 𝔰𝔲(N))
- **Non-degenerate** for semisimple Lie algebras (Cartan's criterion)
- Used to raise/lower Lie algebra indices (like the metric tensor)
- Appears in the Yang-Mills action: S_YM = -1/(4g²) ∫ tr(F_μν F^μν) d⁴x

---

## 7. U(1) in Full Detail

### 7.1 Five Equivalent Descriptions

```
Description 1 (Topological):   S¹ = unit circle in ℝ² or ℂ
Description 2 (Algebraic):     ℝ/ℤ or ℝ/2πℤ (quotient of ℝ by integer lattice)
Description 3 (Matrix):        {e^{iθ} : θ ∈ ℝ} ⊂ ℂ ≅ GL(1,ℂ)
Description 4 (Representation):smallest compact connected group whose representations
                                 classify U(1) charges (integers ℤ)
Description 5 (Gauge):         gauge group of electromagnetism
```

All five are isomorphic as Lie groups.

### 7.2 The Manifold Structure of U(1)

U(1) is a 1-dimensional smooth manifold — a circle. The single coordinate chart
(with overlap) is:

```
Chart 1: U₁ = U(1) \ {-1}
  φ₁: U₁ → (-π, π)
  φ₁(e^{iθ}) = θ   for θ ∈ (-π, π)

Chart 2: U₂ = U(1) \ {+1}
  φ₂: U₂ → (0, 2π)
  φ₂(e^{iθ}) = θ   for θ ∈ (0, 2π)

Transition map on U₁ ∩ U₂ = U(1) \ {±1}:
  Two components: upper semicircle θ ∈ (0,π) and lower semicircle θ ∈ (-π,0)
  On each: φ₂ ∘ φ₁⁻¹ is the identity (just the same coordinate)
```

U(1) requires two charts because S¹ is not diffeomorphic to an open interval
(it is compact, while open intervals are not).

### 7.3 The Lie Algebra 𝔲(1)

```
𝔲(1) = T_{1} U(1) = iℝ = {iθ : θ ∈ ℝ}

As a vector space: 1-dimensional over ℝ
Generator:  T = i  (or sometimes written as T = -i by convention)
Bracket: [iα, iβ] = (iα)(iβ) - (iβ)(iα) = 0  (abelian!)
```

The Lie algebra is abelian — the only 1-dimensional Lie algebra over ℝ.

**The standard physics convention** uses Hermitian generators:

```
Instead of working with 𝔲(1) = iℝ, physicists write the generator as:
T = 1  (real, Hermitian)
And write group elements as: g = e^{iθT} = e^{iθ}
So the gauge potential A_μ = A_μ · T is real-valued
```

This is why the gauge potential in electromagnetism is a real-valued 1-form
(not imaginary-valued), and why the field strength F_μν = ∂_μA_ν - ∂_νA_μ is real.

### 7.4 Representations of U(1)

The irreducible representations (irreps) of U(1) are all 1-dimensional
(since U(1) is abelian, by Schur's lemma). They are classified by integers n ∈ ℤ:

```
ρ_n: U(1) → GL(1,ℂ) = ℂ*
ρ_n(e^{iθ}) = e^{inθ}
```

The integer n is the **charge** of the representation. For electromagnetism:
- Electron: charge n = -1 (or +1 depending on sign convention)
- Positron: charge n = +1
- Photon: charge n = 0 (neutral, in the adjoint)
- W boson: charge n = ±1 (but under SU(2)×U(1), not just U(1))

**Why charges are integers:** The map n ↦ ρ_n is a group homomorphism from ℤ
to the set of 1-d representations of U(1). For ρ_n to be a well-defined
representation (single-valued), we need:
```
ρ_n(e^{iθ}) = e^{inθ} with e^{in(θ+2π)} = e^{inθ}
```
This requires n ∈ ℤ. The quantization of electric charge is a direct consequence
of the compactness of U(1) (π₁(U(1)) = ℤ implies charge quantization).

### 7.5 U(1) as a Quotient

The most illuminating way to understand U(1) is as a quotient:

```
U(1) ≅ ℝ/2πℤ

The quotient map q: ℝ → U(1) is:
q(θ) = e^{iθ}

This is a covering map:
  - q is a local homeomorphism (local diffeomorphism)
  - Each point e^{iθ} ∈ U(1) has preimage q⁻¹(e^{iθ}) = {θ + 2πn : n ∈ ℤ}
  - The fiber is countably infinite: ℤ
```

This is the **universal cover** of U(1): ℝ is simply connected, and U(1) = ℝ/2πℤ.

The covering structure explains:
- Why π₁(U(1)) = ℤ (loops in U(1) lift to paths in ℝ; the winding number is the
  endpoint displacement divided by 2π)
- Why charges are quantized (representations of U(1) come from representations of
  ℝ that are periodic with period 2π — i.e., e^{inθ} for n ∈ ℤ)
- Why monopole charges are quantized (Dirac quantization condition)

### 7.6 The Hopf Fibration: U(1) Inside S³

The **Hopf fibration** is a fiber bundle:

```
U(1) → S³ → S²
```

S³ is the 3-sphere {(z₁,z₂) ∈ ℂ² : |z₁|² + |z₂|² = 1}.
The bundle map π: S³ → S² sends (z₁,z₂) ↦ the point on S² given by the
ratio z₁/z₂ ∈ ℂ ∪ {∞} ≅ S² (Riemann sphere).
The fiber over each point is a copy of U(1) = {e^{iθ}(z₁,z₂)}.

Physical significance: the Hopf fibration is the geometric structure of a magnetic
monopole of strength 1/2. The non-trivial topology of the bundle (it is not
S² × S¹) is the reason the monopole exists and its charge is quantized.

---

## 8. SU(2): The Three-Sphere and Rotations

### 8.1 The Matrix Realization

```
SU(2) = {(a  -b̄) : a, b ∈ ℂ, |a|² + |b|² = 1}
         (b   ā)
```

Parametrize with a = x₀ + ix₃, b = x₂ + ix₁, where x₀² + x₁² + x₂² + x₃² = 1:

```
SU(2) = {x₀I + i(x₁σ₁ + x₂σ₂ + x₃σ₃) : x₀² + |**x**|² = 1}

where σ₁, σ₂, σ₃ are the Pauli matrices:
σ₁ = (0 1)    σ₂ = (0 -i)    σ₃ = (1  0)
     (1 0)         (i  0)          (0 -1)
```

**As a manifold:** SU(2) ≅ S³ (the 3-sphere in ℝ⁴ = ℂ²).
This is the key geometric fact about SU(2):

```
SU(2) is the unit 3-sphere:   x₀² + x₁² + x₂² + x₃² = 1
```

### 8.2 The Lie Algebra 𝔰𝔲(2)

```
𝔰𝔲(2) = {X ∈ Mat(2×2, ℂ) : X† = -X, tr X = 0}
        = {(  ia    -b̄ + ic) : a, b, c ∈ ℝ}
           (-b̄ - ic  -ia  )
```

Standard basis (using i/2 times Pauli matrices):

```
T₁ = iσ₁/2 = (0  i/2)    T₂ = iσ₂/2 = (0  1/2)    T₃ = iσ₃/2 = (i/2   0)
              (i/2  0)                  (-1/2  0)                   (0  -i/2)
```

Or in physics convention with Hermitian generators (dropping the i):

```
Lₐ = σₐ/2   (a = 1,2,3)

Lie bracket: [Lₐ, Lᵦ] = i ε_{abc} Lc
             (ε_{abc} is the Levi-Civita symbol)
```

This is the **angular momentum algebra** from quantum mechanics! SU(2) is the
quantum mechanical rotation group.

### 8.3 The Structure Constants of SU(2)

The Lie bracket [Lₐ, Lᵦ] = i ε_{abc} Lc means the **structure constants** are:

```
f^{abc} = ε^{abc}   (Levi-Civita symbol)

f^{123} = f^{231} = f^{312} = 1
f^{132} = f^{213} = f^{321} = -1
all others = 0
```

In gauge theory, the structure constants determine:
- The 3-gluon vertex (in QCD with SU(3), analogous structure)
- The self-interaction of the W bosons (in SU(2) weak theory)
- The non-abelian part of F_μν: F_μν^c = ∂_μA_ν^c - ∂_νA_μ^c - g f^{abc}A_μ^aA_ν^b

### 8.4 The Exponential Map for SU(2)

For a general element X = **n** · **L** = nₐLₐ (unit vector **n**, angle θ):

```
exp(θ **n** · **L**) = cos(θ/2) I + 2i sin(θ/2) **n** · **L**
                     = cos(θ/2) I + i sin(θ/2)(n₁σ₁ + n₂σ₂ + n₃σ₃)
```

This is a rotation by angle θ about axis **n** in the spin-1/2 representation.

**Key property:** exp is surjective for SU(2) (since SU(2) ≅ S³ is simply connected
and compact). Every element of SU(2) is e^X for some X ∈ 𝔰𝔲(2).

### 8.5 SU(2) Is Simply Connected

```
π₀(SU(2)) = 0   (connected)
π₁(SU(2)) = 0   (simply connected — no non-trivial loops!)
π₂(SU(2)) = 0
π₃(SU(2)) = ℤ   (non-trivial 3-spheres — instantons!)
```

Contrast with:
```
π₁(U(1)) = ℤ   (U(1) is NOT simply connected)
π₁(SO(3)) = ℤ₂  (SO(3) is NOT simply connected)
```

The simple connectivity of SU(2) is why:
- SU(2) is the universal cover of SO(3)
- The exponential map exp: 𝔰𝔲(2) → SU(2) determines SU(2) completely
- There are no theta-term complications from π₁ (but there are instantons from π₃)

---

## 9. SO(3) vs SU(2): Covering Groups and Spinors

### 9.1 The 2:1 Covering Map

There is a surjective Lie group homomorphism:

```
φ: SU(2) → SO(3)
φ(U) ψ_R = U **v** U†   for **v** = v^a σ_a ∈ 𝔰𝔲(2) ≅ ℝ³
```

This map is 2:1: both U and -U map to the same rotation R ∈ SO(3).

```
ker(φ) = {I, -I} = ℤ₂
SO(3) = SU(2)/ℤ₂
```

### 9.2 Why ±I Map to the Same Rotation

```
φ(U)**v** = U**v**U†
φ(-U)**v** = (-U)**v**(-U)† = U**v**U†   (the two minus signs cancel!)
```

So φ(U) = φ(-U) — they give the same rotation.

### 9.3 The Physical Consequence: Spinors

A spin-1/2 particle (electron) transforms under SU(2), not SO(3):

```
Under rotation by 2π: U = exp(2π · i σ₃/2) = exp(iπσ₃) = -I
```

So a 2π rotation sends a spin-1/2 state to its negative:

```
|ψ⟩ → -|ψ⟩   under 2π rotation
```

This is observable! In neutron interferometry, the beam acquires a phase of -1
under 360° rotation and returns to its original phase after 720°. This is the
spinor property, and it comes directly from the 2:1 covering SU(2) → SO(3).

### 9.4 The General Pattern: Covering Groups

For every Lie group G, there is a unique **universal cover** G̃:
- G̃ is simply connected (π₁(G̃) = 0)
- G = G̃/Γ for a discrete normal subgroup Γ ≅ π₁(G)

```
Group G          Universal cover G̃    Quotient
U(1) = S¹        ℝ                    ℝ/ℤ
SO(2)            ℝ                    ℝ/ℤ
SO(3)            SU(2)                SU(2)/ℤ₂
SO(n), n≥3       Spin(n)              Spin(n)/ℤ₂
SO(n,1), Lorentz SL(2,ℂ)             SL(2,ℂ)/ℤ₂
```

In gauge theory, the relevant group is always the simply connected cover,
because the path integral and representations are better behaved.

---

## 10. SU(3) and the Gell-Mann Matrices

### 10.1 The Structure of SU(3)

```
SU(3) = {A ∈ GL(3,ℂ) : A†A = I, det A = 1}
Dimension: 8   (as a real manifold)
Rank: 2         (dimension of maximal torus = dimension of maximal abelian subgroup)
```

SU(3) is the gauge group of QCD (quantum chromodynamics). The 8 gauge bosons
are the gluons, corresponding to the 8 generators of 𝔰𝔲(3).

### 10.2 The Gell-Mann Matrices

The standard basis for 𝔰𝔲(3) uses the **Gell-Mann matrices** λ₁,...,λ₈
(analogous to Pauli matrices for SU(2)):

```
λ₁ = (0 1 0)    λ₂ = (0 -i 0)   λ₃ = (1  0  0)
     (1 0 0)         (i  0 0)        (0 -1  0)
     (0 0 0)         (0  0 0)        (0  0  0)

λ₄ = (0 0 1)    λ₅ = (0 0 -i)   λ₆ = (0 0 0)
     (0 0 0)         (0 0  0)        (0 0 1)
     (1 0 0)         (i 0  0)        (0 1 0)

λ₇ = (0  0  0)  λ₈ = 1/√3 (1  0  0)
     (0  0 -i)             (0  1  0)
     (0  i  0)             (0  0 -2)
```

Properties:
```
tr(λₐ) = 0           (traceless)
tr(λₐλᵦ) = 2δₐᵦ     (orthonormality)
[λₐ, λᵦ] = 2i f^{abc} λc   (Lie bracket, f^{abc} = SU(3) structure constants)
```

Generators in physics convention: Tₐ = λₐ/2.

### 10.3 The Structure of SU(3): Cartan Subalgebra

The **Cartan subalgebra** of 𝔰𝔲(3) is the maximal abelian subalgebra:

```
Cartan subalgebra = span{T₃, T₈} = span{λ₃/2, λ₈/2}
```

These are the simultaneously diagonalizable generators — they commute:
[T₃, T₈] = 0.

In QCD: T₃ and T₈ correspond to the two independently conserved color charges.
In the quark model: T₃ is the isospin component, T₈ is related to hypercharge.

### 10.4 Roots and the Root System

The 6 remaining generators (λ₁,...,λ₇ minus λ₃, λ₈) can be organized into
**raising and lowering operators** — ladder operators that shift between
eigenstates of the Cartan subalgebra.

For SU(3), the 6 roots form the A₂ root system:

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
of the Ω⁻ baryon (the missing corner of the decuplet) before it was discovered.

---

## 11. Representations of Lie Groups

### 11.1 Definition

A **representation** of a Lie group G is a smooth group homomorphism:

```
ρ: G → GL(V)
```

where V is a finite-dimensional vector space (the **representation space**) and
GL(V) is the group of invertible linear maps on V.

The **dimension** of the representation is dim V.
The representation is **faithful** if ker(ρ) = {e} (injective).

### 11.2 The Induced Representation of the Lie Algebra

Given a representation ρ: G → GL(V), differentiating at the identity gives:

```
dρ: 𝔤 → 𝔤𝔩(V) = End(V)
dρ(X) = d/dt|_{t=0} ρ(exp(tX))
```

This is a Lie algebra homomorphism: dρ([X,Y]) = [dρ(X), dρ(Y)].

For matrix groups: ρ(exp X) = exp(dρ(X)) (exponential intertwines the two).

### 11.3 Irreducible Representations

A representation ρ: G → GL(V) is **irreducible** (an irrep) if there is no
proper non-zero subspace W ⊂ V with ρ(g)W ⊆ W for all g ∈ G.

**Schur's Lemma:** If ρ₁, ρ₂ are irreps, then any intertwining map T: V₁ → V₂
(with ρ₂(g)T = Tρ₁(g)) is either 0 or an isomorphism.

Consequence for abelian groups (like U(1)):
All irreps are 1-dimensional. (Because any 1-d subspace is invariant for
abelian G, so the only irreducible representation on V is when dim V = 1.)

### 11.4 The Irreps of SU(2)

The irreps of SU(2) are labeled by a half-integer j = 0, 1/2, 1, 3/2, 2, ...:

```
ρⱼ: SU(2) → GL(ℂ^{2j+1})

Basis: |j, m⟩ for m = -j, -j+1, ..., j-1, j

Action: T₃|j,m⟩ = m|j,m⟩
        T±|j,m⟩ = √(j(j+1) - m(m±1)) |j,m±1⟩   (T± = T₁ ± iT₂)
```

In physics:
- j=0: singlet (spin-0, scalar particle)
- j=1/2: doublet (spin-1/2, electron, quarks, neutrinos)
- j=1: triplet (spin-1, vector particles like W bosons)
- j=3/2: quadruplet (spin-3/2, delta baryons)

**Dimension of irrep j:** 2j+1.

**Tensor product:** ρⱼ₁ ⊗ ρⱼ₂ = ⊕_{j=|j₁-j₂|}^{j₁+j₂} ρⱼ (Clebsch-Gordan decomposition)

### 11.5 The Irreps of SU(3)

The irreps of SU(3) are labeled by two non-negative integers (p, q):

```
Dimension = (p+1)(q+1)(p+q+2)/2

Examples:
(0,0): dimension 1   — singlet (colorless)
(1,0): dimension 3   — fundamental representation (quarks)
(0,1): dimension 3̄   — anti-fundamental (antiquarks)
(1,1): dimension 8   — adjoint representation (gluons)
(3,0): dimension 10  — decuplet (baryon resonances)
(2,1): dimension 15
```

The famous quark model prediction: hadrons are color singlets (p=q=0),
built from quarks in the 3 and 3̄ and their bound states.

---

## 12. The Peter-Weyl Theorem and Harmonic Analysis

### 12.1 Functions on a Compact Lie Group

For a compact Lie group G with Haar measure dg (the unique translation-invariant
measure with ∫_G dg = 1), the space L²(G) of square-integrable functions has
a beautiful decomposition.

**Peter-Weyl Theorem:** The matrix coefficients of irreducible representations:

```
ρⱼᵢⱼ: G → ℂ,  g ↦ ⟨eᵢ, ρⱼ(g)eⱼ⟩  (i,j = 1,...,dim ρ)
```

form an orthonormal basis for L²(G):

```
L²(G) = ⊕_ρ irreps V_ρ ⊗ V_ρ*
```

where V_ρ is the representation space of irrep ρ.

### 12.2 Application to U(1): Fourier Series

For G = U(1):

- Irreps: ρ_n(e^{iθ}) = e^{inθ} for n ∈ ℤ
- Matrix coefficients: e^{inθ}
- Peter-Weyl: L²(U(1)) = L²(S¹) = ⊕_{n∈ℤ} ℂe^{inθ}

This is **Fourier series**: every square-integrable function on the circle
has a Fourier expansion f(θ) = Σ_{n∈ℤ} cₙ e^{inθ}.

The Peter-Weyl theorem is the group-theoretic generalization of Fourier analysis
to any compact Lie group.

### 12.3 Application to SU(2): Spherical Harmonics

For G = SU(2) ≅ S³:

- Irreps: ρⱼ for j = 0, 1/2, 1, 3/2, ...
- Peter-Weyl gives an orthonormal basis for L²(S³) labeled by (j, m, m')
- Restricting to functions invariant under U(1) ⊂ SU(2) gives **spherical harmonics**
  on S² = SU(2)/U(1)

The spherical harmonics Y_l^m(θ,φ) are precisely the matrix coefficients of
SU(2) representations restricted to the appropriate quotient. This is why
spherical harmonics have the angular momentum quantum numbers they do.

---

## 13. Roots, Weights, and the Classification of Simple Lie Algebras

### 13.1 The Cartan Classification

Every simple Lie algebra over ℂ belongs to one of:

```
Classical algebras:
  Aₙ (n≥1):  𝔰𝔩(n+1,ℂ)  → SU(n+1) (compact real form)
  Bₙ (n≥2):  𝔰𝔬(2n+1,ℂ) → SO(2n+1)
  Cₙ (n≥3):  𝔰𝔭(2n,ℂ)   → Sp(n) (symplectic)
  Dₙ (n≥4):  𝔰𝔬(2n,ℂ)   → SO(2n)

Exceptional algebras:
  G₂ (dim 14)
  F₄ (dim 52)
  E₆ (dim 78)
  E₇ (dim 133)
  E₈ (dim 248)
```

The gauge groups of the Standard Model:
```
U(1)  = U(1)    (circle, not simple but simple up to center)
SU(2) = A₁      (the simplest non-abelian simple Lie group)
SU(3) = A₂
SU(5) = A₄      (Grand Unified Theory gauge group)
SO(10) = D₅     (another GUT candidate)
E₈ × E₈         (heterotic string theory)
```

### 13.2 Roots

For a simple Lie algebra 𝔤 with Cartan subalgebra 𝔥 (maximal abelian diagonalizable
subalgebra), the **roots** are the nonzero eigenvalues of the adjoint action of 𝔥:

```
If H ∈ 𝔥 and E_α ∈ 𝔤 satisfy [H, E_α] = α(H) E_α for all H ∈ 𝔥
then α: 𝔥 → ℂ is a root and E_α is a root vector (raising/lowering operator)
```

The root system encodes the entire Lie algebra structure:
```
[H, E_α] = α(H) E_α
[E_α, E_{-α}] = H_α   (coroot)
[E_α, E_β] = N_{αβ} E_{α+β}   (if α+β is a root)
```

### 13.3 Weights and Representations

For a representation ρ: 𝔤 → 𝔤𝔩(V), the **weights** of ρ are the eigenvalues
of the Cartan subalgebra action on V:

```
Hv = λ(H)v for all H ∈ 𝔥 and some v ∈ V
```

λ is a weight and v is a weight vector.

For SU(2):
- The single Cartan generator is T₃
- Weights are the eigenvalues m = -j, ..., +j for spin-j representation
- The "highest weight" is j (the maximum eigenvalue)

For SU(3):
- Two Cartan generators T₃, T₈
- Weights are 2d vectors (isospin, hypercharge)
- The weight diagram of the fundamental (1,0) rep is the quark triangle (u, d, s)

---

## 14. Homomorphisms, Subgroups, and Quotients

### 14.1 Lie Group Homomorphisms

A **Lie group homomorphism** φ: G → H is a smooth group homomorphism.

Key examples:
```
det: GL(n) → ℝ*         (determinant — surjective onto nonzero reals)
det: U(n) → U(1)        (determinant of unitary matrix has |det| = 1)
Ad: G → GL(𝔤)           (adjoint representation)
π: SU(2) → SO(3)        (2:1 covering, kernel = ℤ₂)
exp: 𝔤 → G              (not quite a homomorphism, but the exponential map)
```

### 14.2 Closed Subgroups

By the **Closed Subgroup Theorem** (Cartan): if H is a closed subgroup of a
Lie group G, then H is automatically a Lie group (an embedded submanifold).

This means: to specify a subgroup of a matrix Lie group, you just need a
closed subset that is also a subgroup — the smooth structure comes for free.

Examples of closed subgroups:
```
U(1) ⊂ SU(2):   U(1) = {(e^{iθ}  0  ) : θ ∈ ℝ} (diagonal matrices)
                          (0  e^{-iθ})

SU(2) ⊂ SU(3):  Embed as block (SU(2)  0)
                                (0     1)

U(1)×U(1) ⊂ SU(3): The maximal torus (diagonal matrices in SU(3))
```

### 14.3 Quotient Groups and Homogeneous Spaces

If H ⊂ G is a closed subgroup, the **quotient space** G/H is a smooth manifold.
If H is a normal subgroup (gHg⁻¹ = H for all g), then G/H is a Lie group.

Important examples:
```
SO(3) = SU(2)/ℤ₂          (ℤ₂ = {I, -I} is normal in SU(2))
S² = SU(2)/U(1) = SO(3)/SO(2)   (2-sphere as a homogeneous space)
S^n = SO(n+1)/SO(n)        (n-sphere as a homogeneous space)
ℝP^n = SO(n+1)/O(n)        (real projective space)
G_SM/H_unbroken            (symmetry breaking in particle physics)
```

The last example is crucial: in the Standard Model, after Higgs symmetry breaking:
```
U(1)_Y × SU(2)_L → U(1)_EM
```
The surviving U(1)_EM is a specific U(1) subgroup of the original gauge group,
and the broken symmetries become the massive W±, Z bosons.

---

## 15. The Topology of Lie Groups

### 15.1 Homotopy Groups

The topological complexity of Lie groups is captured by their homotopy groups:

```
G          π₀   π₁      π₂    π₃
U(1)        0    ℤ       0     0
SU(2)       0    0       0     ℤ
SU(3)       0    0       0     ℤ
SO(2)       0    ℤ       0     0
SO(3)       0    ℤ₂      0     ℤ
Sp(1)≅SU(2) 0    0       0     ℤ
G₂          0    0       0     ℤ
```

Key facts:
- All compact connected Lie groups are connected: π₀ = 0
- π₁(G) = π₁(G/[G,G]) — the fundamental group comes from the abelian part
- For simple simply connected compact Lie groups: π₃ = ℤ (Bott periodicity)
- π₃(G) = ℤ is the origin of instantons in Yang-Mills theory

### 15.2 Why π₁(U(1)) = ℤ and Why It Matters

A loop in U(1) = S¹ is a continuous map γ: [0,1] → S¹ with γ(0) = γ(1) = 1.
The homotopy class of γ is its **winding number**: how many times γ wraps around the circle.

```
Winding number +1: γ(t) = e^{2πit}   (once counterclockwise)
Winding number +2: γ(t) = e^{4πit}   (twice counterclockwise)
Winding number -1: γ(t) = e^{-2πit}  (once clockwise)
Winding number 0:  any contractible loop
```

π₁(U(1)) = ℤ means loops are classified by their winding number.

**Physical consequences:**
1. **Charge quantization:** Representations of U(1) are ρ_n(e^{iθ}) = e^{inθ}.
   The integer n is the charge, and it is quantized because n must be an integer
   for the representation to be single-valued (well-defined on S¹).

2. **Magnetic monopoles:** The Dirac quantization condition eg = 2π (in natural units)
   requires the product of electric charge e and magnetic charge g to be quantized.
   This is because the gauge field on S² around a monopole is classified by
   π₁(U(1)) = ℤ.

3. **Flux quantization:** In superconductors, the magnetic flux through a loop is
   quantized in units of Φ₀ = h/2e. This comes from the U(1) winding number of
   the Cooper pair condensate.

### 15.3 The Bott Periodicity Theorem

For any simple simply connected compact Lie group G:

```
π₃(G) = ℤ
```

This is Bott periodicity applied to Lie groups. For SU(n), the map S³ → SU(n)
generating π₃ is the basic instanton. The integer generator of π₃ is what labels
the topological charge (instanton number) of Yang-Mills gauge fields.

### 15.4 Compact vs. Non-Compact Real Forms

Every complex simple Lie algebra has a unique compact real form and multiple
non-compact real forms. In gauge theory:

```
Complex algebra A₁ = 𝔰𝔩(2,ℂ):
  Compact real form: 𝔰𝔲(2)  → SU(2) (gauge group of weak interaction)
  Non-compact:       𝔰𝔩(2,ℝ) → SL(2,ℝ) (appears in 2+1 gravity)

Complex algebra D₁ = 𝔰𝔬(2,ℂ) ≅ A₁:
  Compact: 𝔰𝔬(2) → SO(2) ≅ U(1) (rotation in 2D)
  Non-compact: 𝔰𝔬(1,1) → SO(1,1) (Lorentz boosts in 1+1D)
```

The compact forms give **unitary representations** and are used for internal
symmetries (gauge groups). Non-compact forms appear in spacetime symmetries
(Lorentz group, conformal group).

---

## 16. Lie Groups in Gauge Theory: The Full Picture

### 16.1 The Role of Each Structure

Here is how every aspect of a gauge Lie group appears in physics:

```
Lie group structure        →  Physical meaning
─────────────────────────────────────────────────────────────────
Manifold dimension          →  Number of gauge bosons
Group multiplication        →  Composition of gauge transformations
Lie algebra basis {Tₐ}     →  Internal quantum numbers (color, isospin, hypercharge)
Structure constants f^{abc} →  Self-interaction of gauge bosons (3-point vertex)
f^{abc} f^{abd}             →  4-point gauge boson vertex
Adjoint representation      →  How gauge bosons carry their own charge
Fundamental representation  →  How matter fields (quarks, leptons) couple
Killing form                →  Kinetic term for gauge field: tr(F²)
Compact real form           →  Unitary, bounded representations (probability conserved)
π₁(G)                      →  Charge quantization, magnetic monopoles
π₃(G)                      →  Instantons, topological charge
Maximal torus               →  Conserved charges (Cartan generators)
Roots                       →  Particle multiplet structure
Weights                     →  Quantum numbers of matter particles
Casimir operators           →  Gauge-invariant mass terms
```

### 16.2 U(1) in This Picture

```
U(1) as a gauge group:
  Manifold:    S¹ (1-dimensional)
  Dimension:   1  → 1 gauge boson (the photon)
  Lie algebra: iℝ, trivial bracket → [A_μ, A_ν] = 0 → no photon self-coupling
  Adjoint rep: trivial → photon is neutral (carries no electric charge)
  π₁(U(1)) = ℤ → electric charge quantization
  Representations ρ_n: e^{iθ} ↦ e^{inθ} → charge n ∈ ℤ for each particle
  Compact → unitary representations → probability conserved
  Killing form: B(X,Y) = 0 (1-dim, trivial) → kinetic term is ∫ F_μν F^μν d⁴x
```

### 16.3 Why the Gauge Group Determines the Theory

Given a gauge group G, the entire Yang-Mills theory is fixed:

```
Step 1: Choose G
   ↓
Step 2: The Lie algebra 𝔤 gives the gauge bosons (one per generator)
   ↓
Step 3: The structure constants f^{abc} determine the self-coupling
   ↓
Step 4: The Killing form determines the kinetic term tr(F²)
   ↓
Step 5: Choose representations ρ for matter fields
   ↓
Step 6: The covariant derivative D_μ = ∂_μ + gρ(A_μ) determines matter coupling
   ↓
Step 7: The full Lagrangian L = -tr(F²)/4g² + ψ̄(iD̸ - m)ψ is completely specified
```

For the Standard Model:
```
G = U(1)_Y × SU(2)_L × SU(3)_c
Matter representations chosen to match experimental particle content
→ Every interaction of every known particle is fixed
```

### 16.4 The Non-Abelian Obstruction in Transformer Attention

Connecting to the gauge theory tutorial: for multi-head attention, the force on
hidden state h is:

```
F^(h)(h) = Σ_j softmax(Q_h h · K_h h_j / √d_h) V_h h_j   (head h contribution)

Total force: F(h) = Σ_h F^(h)(h) = Σ_h A^(h)(h) h   (schematically)
```

Each head contributes a connection A^(h) valued in a dₕ-dimensional subspace.
After the output projection W_O mixes all heads:

```
A_total(h) = W_O (A^(1)(h) ⊕ ... ⊕ A^(H)(h)) W_O^†  ∈  End(ℝᵈ)
```

The curvature:

```
F_ij = (∂A/∂h)_ij - (∂A/∂h)_ji + [A_i, A_j]
```

contains cross-head commutators [A^(h), A^(h')] that are non-zero because
the head subspaces, after W_O mixing, are not orthogonal.

This is a **non-abelian gauge obstruction** in exactly the Yang-Mills sense:
just as SU(N) gauge theory has non-zero [A_μ, A_ν] that prevents the gauge
potential from being a gradient, the multi-head attention connection has
non-zero cross-head commutators that prevent any scalar potential from
reproducing the full attention dynamics.

### 16.5 Statement Decoder: The Full Logical Chain

The statement

> *"cross-head commutators [A^(h), A^(h')] generate non-abelian curvature
> that obstructs any scalar potential on hidden-state space, regardless of capacity"*

can now be unpacked term by term, with every concept traced to its definition
in this tutorial and the Gauge Theory Tutorial.

---

**"cross-head"**

Each attention head h = 1,...,H produces a force on the hidden state.
The word "cross-head" means we are looking at the **interaction between two
different heads h ≠ h'** — not a single head acting on itself.
Within a single head, [A^(h), A^(h)] = 0 always (bracket of anything with
itself vanishes, by §4.1a). The interesting structure comes from pairs h ≠ h'.

---

**"commutators [A^(h), A^(h')]"**

A^(h) is the **connection 1-form** (gauge potential) contributed by head h.
After the output projection W_O recombines all heads, A^(h) becomes a
**d×d matrix** acting on the full hidden-state space ℝᵈ.

The commutator (§4.1a) is:

```
[A^(h), A^(h')] = A^(h) A^(h') - A^(h') A^(h)
```

This is a matrix: the product of two d×d matrices minus the product in
reversed order. It measures how much the force contributions of the two
heads fail to commute as linear operators on ℝᵈ.

For the commutator to vanish, A^(h) and A^(h') must commute as matrices:
A^(h) A^(h') = A^(h') A^(h). This would require very special alignment
of the weight matrices W_Q^h, W_K^h, W_V^h, W_O — it is not the generic case.

---

**"generate non-abelian curvature"**

The full connection is A = Σ_h A^(h). Its curvature (§7.1-7.2 of the
Gauge Theory Tutorial) is:

```
F = ∂A - ∂A + [A, A]
       ↑            ↑
  abelian part    non-abelian part
  (would exist    ([A,A] = Σ_{h≠h'} [A^(h), A^(h')] + Σ_h [A^(h), A^(h)])
   even for a    = Σ_{h≠h'} [A^(h), A^(h')]   ← the cross-head commutators)
   single head)
```

"Non-abelian curvature" specifically refers to the [A,A] term — the part of
F that exists **because the gauge group is non-abelian** (§1.4 of this tutorial,
§8.2 of the Gauge Theory Tutorial).

When the cross-head commutators [A^(h), A^(h')] are non-zero:

```
[A, A] = Σ_{h≠h'} [A^(h), A^(h')] ≠ 0
→  F_μν ≠ 0   (non-zero curvature)
→  G is effectively non-abelian
```

---

**"obstructs any scalar potential"**

A **scalar potential** is any smooth function V: ℝᵈ → ℝ. If a force field
could be written F = -∇V, it would be **conservative** (path-independent work,
zero net work around any closed loop).

The **obstruction theorem** (§7.4 of the Gauge Theory Tutorial, proved via
Clairaut's theorem) states:

```
F = -∇V for some scalar V     ⟺     F_μν = 0 everywhere

Contrapositive:
F_μν ≠ 0                      ⟹     F ≠ -∇V for ANY scalar V
```

The proof is an algebraic identity: if V is any smooth function, then

```
(∂_i ∂_j V) - (∂_j ∂_i V) = 0     (Clairaut/Schwarz theorem: mixed partials commute)
```

So the curvature of -∇V is always exactly zero — not approximately zero,
not zero for large enough V — **exactly zero**, by pure calculus.

Since F_μν ≠ 0 (from the cross-head commutators), and any scalar V gives
zero curvature, we have:

```
F ≠ -∇V for ANY scalar V
```

This is the obstruction.

---

**"on hidden-state space"**

The base manifold is the hidden-state space ℝᵈ — the d-dimensional real
vector space in which transformer hidden states live. The connection A and
its curvature F are defined on this space (§16.1 of the Gauge Theory Tutorial).

---

**"regardless of capacity"**

This phrase means: the failure is not about V_ψ being too small, too shallow,
or having too few parameters. Even if you replace V_ψ with:
- A deeper neural network
- An infinitely wide network
- A universal approximator
- The exact analytical solution to any optimization problem

...it still cannot work. The reason:

```
Every smooth scalar V — simple or complex — satisfies:
  Curvature of (-∇V) = ∂_i(-∂_j V) - ∂_j(-∂_i V) = 0   (Clairaut, always)

The attention force satisfies:
  Curvature of F = [A,A] ≠ 0   (from cross-head commutators)

Therefore: F ≠ -∇V for ANY smooth V, no matter how expressive.
```

This is an **algebraic/topological obstruction** — it is in the same class as
saying "a continuous map from S² to ℝ² cannot be injective" (a topological fact
that holds regardless of how cleverly you choose the map). The capacity of V_ψ
is simply irrelevant because the problem is structural, not quantitative.

---

**The Complete Logical Chain**

```
Multi-head attention has H ≥ 2 heads
         ↓
Each head h contributes a matrix-valued connection A^(h) ∈ End(ℝᵈ)
         ↓
After W_O mixing, different heads act in overlapping subspaces
         ↓
Cross-head commutators [A^(h), A^(h')] = A^(h)A^(h') - A^(h')A^(h) ≠ 0
  (because matrices acting in overlapping subspaces generically don't commute)
         ↓
[A^(h), A^(h')] ≠ 0  →  [A,A] = Σ_{h≠h'} [A^(h), A^(h')] ≠ 0  (non-abelian)
         ↓
F_μν = ∂A - ∂A + [A,A] ≠ 0   (non-zero curvature)
         ↓
Obstruction theorem (Clairaut): any scalar V gives zero curvature
F_μν ≠ 0 but Curvature(-∇V) = 0  →  F ≠ -∇V
         ↓
No scalar potential of any form or capacity can represent F
         ↓
The shared-V_ψ test fails for GPT-2 middle layers (R² = 0.04–0.20)
regardless of V_ψ architecture or parameter count.
```

**Where each step is defined:**

| Step | Definition location |
|---|---|
| Connection A^(h) | Gauge Theory Tutorial §5.2, §16.1 |
| Matrix commutator [·,·] | Lie Groups Tutorial §4.1a |
| Abelian vs. non-abelian | Lie Groups Tutorial §1.4 |
| Curvature F = dA + [A,A] | Gauge Theory Tutorial §7.1 |
| [A,A] = non-abelian part | Gauge Theory Tutorial §7.2 |
| Obstruction theorem | Gauge Theory Tutorial §7.4 |
| "Regardless of capacity" | Clairaut's theorem (calculus) |
| Hidden-state space = base manifold | Gauge Theory Tutorial §16.1 |
| Cross-head commutators in attention | Gauge Theory Tutorial §16.2 |

---

## 17. Reference: Key Facts About U(1), SU(2), SU(3)

### Quick Reference Table

| Property | U(1) | SU(2) | SU(3) |
|---|---|---|---|
| Manifold | S¹ | S³ | — (8d manifold) |
| Dimension | 1 | 3 | 8 |
| Rank | 1 | 1 | 2 |
| Compact? | Yes | Yes | Yes |
| Simply connected? | No | Yes | Yes |
| π₁ | ℤ | 0 | 0 |
| π₃ | 0 | ℤ | ℤ |
| Lie algebra | 𝔲(1) = iℝ | 𝔰𝔲(2) | 𝔰𝔲(3) |
| Generators | 1 (trivial) | Pauli σₐ/2 | Gell-Mann λₐ/2 |
| Structure constants | 0 | ε^{abc} | f^{abc} |
| Adjoint rep dim | 1 | 3 | 8 |
| Fundamental rep dim | 1 | 2 | 3 |
| Abelian? | Yes | No | No |
| Physical role | Electromagnetism | Weak force | Strong force |
| Gauge bosons | Photon (1) | W±, Z (3) | Gluons (8) |

### Key Formulas

**U(1):**
```
Group:     e^{iθ} · e^{iφ} = e^{i(θ+φ)}
Algebra:   𝔲(1) = iℝ,  [iα, iβ] = 0
Exp map:   exp(iθ) = e^{iθ},  surjective
Reps:      ρ_n(e^{iθ}) = e^{inθ},  n ∈ ℤ
Cover:     ℝ → U(1) = ℝ/2πℤ
```

**SU(2):**
```
Group:     (a  -b̄)(c  -d̄) = (ac-b̄d  -ad̄-b̄c)
           (b   ā)(d   c̄)   (bc+ād   -bd̄+āc̄)
Algebra:   [σₐ/2, σᵦ/2] = i ε_{abc} σc/2
Exp map:   exp(iθ n̂·σ/2) = cos(θ/2)I + i sin(θ/2) n̂·σ
Reps:      Labeled by j = 0, 1/2, 1, ..., dimension 2j+1
Cover:     SU(2) → SO(3) = SU(2)/ℤ₂  (2:1)
```

**SU(3):**
```
Generators: Tₐ = λₐ/2,  a = 1,...,8
Algebra:    [Tₐ, Tᵦ] = i f^{abc} Tc
Cartan:     T₃ = λ₃/2,  T₈ = λ₈/2  (diagonal generators)
Reps:       Labeled by (p,q),  dim = (p+1)(q+1)(p+q+2)/2
Fundamental: (1,0) = 3  (quarks),  (0,1) = 3̄  (antiquarks)
Adjoint:    (1,1) = 8  (gluons)
```

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
