# Mathematical verification

## Objective and equilibrium

For unit-disk centers `c_i` with `||c_i-c_j||>=2`, set

```text
M(c) = (1/2) sum_i ||c_i-c_bar||^2
     = (1/(2n)) sum_{i<j} ||c_i-c_j||^2,
g_ij(c) = (||c_i-c_j||^2-4)/2.
```

The exact contact graph has edges where `g_ij=0` and strict inequalities
on every nonedge. Let `A` contain the Cartesian gradients of its contacts.
The gradient of `M` has components `c_i-c_bar`.
For any velocity `v`,

```text
g_ij(c+t v) = g_ij(c)+t Dg_ij(c)[v]+(t^2/2)||v_i-v_j||^2.
```

Every `v` with `A v>=0` therefore gives a feasible ray for small positive
`t`. At a local minimum, `grad M dot v>=0` on that cone. Linear-cone duality
gives `grad M=A^T lambda` for some `lambda>=0`, without a contact-rank
assumption. Pairing with the centered configuration gives `M=2 sum lambda`.

## Global reduction

After centering, each sublevel `M<=U` is compact: the nonoverlap inequalities
are closed and `sum ||c_i||^2=2M`. Thus a global minimum exists. The included
exact configurations give upper bounds `U_5=34/5` and `U_6=29/3`.

At a disconnected minimum, small independent translations of each contact
component are feasible. Their first variations force every component
centroid to be zero. Rotating a component about zero preserves `M` until a
new contact occurs. To prove that a new contact occurs, choose a center `p`
of smallest radius `r`, in component `C`. If a full rotation of `C` creates
none, every outside center of radius `s>=r` satisfies `s-r>2`. There are at
least two outside centers: a sole outside center would be a singleton at
zero and force coinciding centers. Every other center has
`s>=max(r,2-r)>=|2-r|`. Hence

```text
M > (1/2)[r^2+(n-3)(2-r)^2+2(r+2)^2]
  = 8+(5/2)r^2                  for n=5,
  = 29/3+3(r-1/3)^2             for n=6,
```

contradicting `M<=U_n`. At the first new contact, components merge with all
old contacts preserved and objective unchanged. Repetition produces a
connected global minimum of the same value.

Straight contact segments form a planar embedding. A third center cannot
lie on a segment of length 2. If independent contact segments `ab,cd` crossed,
strict triangle inequalities at their crossing would give
`||a-c||+||b-d||<||a-b||+||c-d||=4`, contradicting separation.
If a vertex `x` were inside a simple enclosing contact cycle, each cycle
edge subtends an angle at `x` at most `pi/3`: its base is 2 and both radii
are at least 2, so `a^2+b^2-ab>=4`. The signed increments wind through
`2pi`; at least six cycle vertices, plus `x`, are required. Thus every
connected contact graph on at most six disks is outerplanar. The
[included census](census.md) covers these classes.

## Angular coverage and symmetry

Choose a deterministic spanning tree and fix its first oriented edge at
`(0,0),(2,0)`. The other `k=n-2` absolute edge directions place all centers
recursively by `c_child=c_parent+2(cos theta_i,sin theta_i)`. Centering is
a translation. This parametrizes every realization modulo rigid motions.
For a graph with leaves, anchor its first labelled leaf edge. Otherwise,
anchor the first contact edge. Complete the tree breadth first in label order.
The same search loop processes both cases.

Choose directions in `[-pi,pi]`. Reflection negates them all, allowing one
coordinate `r` to be nonnegative. The initial rational box is `[0,22/7]`
in that coordinate and `[-22/7,22/7]` elsewhere. Since `22/7>pi`, it covers
every required representative, including angular boundaries.

For graphs with leaves, add only sibling orders justified by checked graph
automorphisms and these registered tree patterns:

| Pattern | Structure and coverage |
|---|---|
| P0 | Ordered siblings are tree leaves and `r` is outside their groups. Reflect to make `r>=0`, then sort the groups. |
| P1 | `r` and every other constrained group are outside exchanged descendant subtrees. Sorting changes only unconstrained descendant directions. |
| P2 | Siblings `u,v` share parent `p`; `r` is edge `u->x`, with `x` adjacent to both siblings. The four centers form a rhombus, so `theta_r=theta_(p,v)` modulo `2pi`. |
| P3 | Every free edge is an interchangeable tree leaf. Reflect to give the group containing `r` a nonnegative member, send it to `r`, then sort the remaining members. |

For P2 put `alpha=theta_(p,u)` and `beta=theta_(p,v)` in `[-pi,pi]`.
Use identity if `alpha<=beta` and `beta>=0`; swap and reflect if
`alpha<=beta<0`; swap if `beta<alpha` and `alpha>=0`; reflect if
`beta<alpha<0`. Each choice gives the required order and `theta_r>=0`.
Closed endpoints allow either representative of the direction `pi`.
The rhombus identity follows from the two radius-2 circle intersections.

[The symmetry routines](../src/search/leaf_filters.jl) check each
transposition against the graph's edge set and verify its registered tree
pattern. Other patterns use reflection alone. Permutations act by
`p[old]=new`, `c'_(p(i))=c_i`; directions are recovered on the oriented tree.
Graphs without leaves use reflection on their first free angle.

## Necessary conditions in a box

Tree contacts hold identically. A non-tree contact interval must contain
zero. A prescribed nonedge must have a strictly positive squared-distance
gap. Configurations with additional contacts belong to their larger exact
contact graph in the census.

For a leaf `i` adjacent to `j`, equilibrium requires

```text
c_i-c_bar = lambda_ij(c_i-c_j),
lambda_ij = ((c_i-c_bar) dot (c_i-c_j))/4 >= 0.
```

Rotate that leaf about `j`, fixing the other centers. Strict noncontacts
remain feasible locally. Differentiation gives
`M''=4((n-1)/n-lambda_ij)`, so `lambda_ij<=(n-1)/n` is necessary.

For a bridge `pq`, let `S` be its side containing `p`. Summed equilibrium
requires `sum_{i in S}(c_i-c_bar)=lambda_pq(c_p-c_q)`, `lambda_pq>=0`.
If `a,b` are the centroids of `S` and its complement, put `v=a-c_q`,
`w=c_q-b`. Rotating `S` about `c_q` changes only
`kappa ||w+R_t v||^2`, where `kappa=|S||S^c|/(2n)>0`. Its necessary
conditions are `w dot Jv=0` and `-w dot v>=0`, with `J(x,y)=(-y,x)`.
A square subsystem of bridge equations can also exclude a root by Krawczyk.

Every surviving box then enters the contact-balance filter, which forms `[A^T | grad M]` and eliminates all
coordinate equations with interval pivots excluding zero. Row operations
enclose every compatible point system; an eliminated entry is zero by the
identity `a-(a/p)p=0`. A zero-left-side row whose right side excludes zero
proves incompatibility. Back substitution encloses all possible multipliers.
Intersect with required force bounds, then check the complete residual and
`2 sum lambda-M`. An uncertain pivot causes subdivision. Multiplier
enclosures express necessary feasibility; root existence is established
by isolation.

## Isolation, transfer and second variation

All graphs assigned to interval computation use the same loop. First apply
symmetry and geometric tests, then leaf and bridge conditions where applicable,
then the contact-balance filter. From depth 18, a square subsystem of bridge
rotation equations can give an additional exclusion. From depth 12, surviving
boxes attempt root isolation whether or not their graph has leaves. Neither
leaf nor bridge conditions are assumed sufficient by themselves to exclude
an entire graph.

The common limits are 3,000,000 visited boxes per graph, depth 140 and minimum
angular width 2^-60. A failed or inconclusive verification causes subdivision;
an operational limit leaves unresolved boxes and prevents closure.

Let `q` be the number of non-tree contacts, with equations `h(theta)=0`.
Outerplanarity gives `q<=k`. For `q=k`, isolate the square contact system.
For `q<k`, use

```text
F(theta,mu) = (h(theta), grad_theta M-Dh(theta)^T mu) = 0.
```

Pulling Cartesian equilibrium through the tree map shows that every
ambient minimum satisfies this system, with `mu` equal to the non-tree
physical multipliers. Their full-equilibrium enclosure supplies the state
box for isolation.

For a square system, an interval Jacobian and a numerical preconditioner
`C` give the outward-rounded Krawczyk image

```text
K(X)=x0-C F(x0)+(I-C DF(X))(X-x0).
```

A residual excluding zero or `K(X)` disjoint from `X` excludes roots.
Strict inclusion in the interior proves a unique root. Otherwise bisect
the widest rational angular interval, breaking ties by first index, and
visit the lower child first. The two closed children cover their parent.
Uncertain arithmetic and resource limits leave unresolved boxes.

Each isolated root must have strict noncontacts. Refinement rechecks and
contracts that same angular root in its stored oriented tree, maps it to centered Cartesian intervals
and recovers the complete, unclipped contact multipliers. Every isolated
box proceeds to analysis. The search makes no comparison with `U_n`.

At equilibrium, on `A v=0`, the intrinsic quadratic form is

```text
Q(v)=sum_i ||v_i-v_bar||^2-sum_ij lambda_ij ||v_i-v_j||^2.
```

Independent contact rows define a smooth local contact manifold; its
tangents are velocities of curves preserving contacts and locally strict
noncontacts. Twice differentiating `M` along such a curve gives `Q`.
A verified negative value excludes a local minimum. A verified negative
multiplier excludes first-order minimality when the multiplier is unique.
The strict bound `inf(M-U_n)>0` excludes a global minimum.

## Two and three disks

For two disks, `M=||c_1-c_2||^2/4>=1`, with equality exactly for a touching
pair. Its graph `A_` is reserved for exact analysis. Centered positions
are `(-1,0),(1,0)`, contact rank 1 and multiplier `1/2`.

For three disks, the pair-distance identity gives `M>=2`; equality requires
all three distances to be 2, hence the equilateral triangle. The connected
classes are the path `BW`, excluded by the interval search, and triangle `Bw`,
verified algebraically. Its contact rank is 3 and each multiplier is `1/3`.

## Four disks

A square gives the upper bound 4. A disconnected local minimum would have
every component centroid zero. For partition `3+1`, the three-vertex component
has a vertex touching both others, at distance at most `4/3` from its centroid,
overlapping the singleton at zero. For `2+2`, the pairs are `+/-u,+/-v`,
with unit `u,v`; the four squared cross distances sum to 8, whereas separation
requires at least 16. The other partitions contain coinciding singletons.
Thus every minimum is connected.

Enumeration of all 64 labeled graphs gives six connected classes.
K4 (`C~`) is impossible: three displacements from one center would have
Gram matrix with diagonal 4 and off-diagonal 2, of rank 3 in the plane.
The interval search excludes `CF,CN,CR`. The remaining classes are the four-cycle
`Cr` and its five-contact endpoint `C^`.

For `Cr`, its two opposite centers can be fixed at `(-a,0),(a,0)`, `a>0`.
The other two are the distinct intersections of their radius-2 circles:

```text
c_1=(-a,0), c_2=(0,b), c_3=(0,-b), c_4=(a,0),
a>0, b>0, a^2+b^2=4.
```

Strict noncontacts are `1<a^2<3`; the closed packing domain is
`1<=a^2<=3`. Throughout, `M=a^2+b^2=4`. The exact verifier works in
`Q[a,b]/(b^2-(4-a^2))`, checks the gaps `4(a^2-1),4(3-a^2)`, cycle
multipliers `1/2`, and a contact minor `-2a^3 b`. The rank is 4 and the
one-dimensional nonrigid tangent has zero intrinsic second variation.

The endpoints are congruent and have graph `C^`. The representative
`(-sqrt(3),0),(sqrt(3),0),(0,1),(0,-1)` has contact rank 5, side
multipliers `1/2` and additional-contact multiplier zero. Exact pair-distance
checks identify it with both family boundaries under relabeling. Coverage
of the cycle family and these endpoints gives the complete equality set.

## Rigid target graphs

Fix an equilateral seed at `(0,0),(2,0),(1,sqrt(3))`. Each attachment to
both endpoints of a placed contact edge has the two circle-intersection
choices. Enumerating these choices covers every realization modulo isometry:

| n | Graph | Seed | Attachments `(vertex,parent,parent)` | Branches |
|---:|---|---|---|---:|
| 5 | `DR{` | `(1,3,5)` | `(4,3,5),(2,5,4)` | 4 |
| 6 | `EDZw` | `(1,4,6)` | `(3,4,6),(5,6,3),(2,6,5)` | 8 |

Exact arithmetic in `Q(sqrt(3))` checks every distance and leaves one labeled
packing branch for each graph. Distance and automorphism checks identify
its input configuration. Exact analysis verifies objectives `34/5,29/3`,
strictly positive multipliers and contact ranks `7,9`, respectively.

These rank and sign checks imply strict local minimality modulo rigid
motions. Fix the centroid and one edge orientation. A sequence of distinct
feasible representatives approaching the target with no larger objective
would have normalized displacement limit `v` with `A v>=0` and
`grad M dot v<=0`. Positive equilibrium forces `A v=0`. Rank `2n-3`
leaves only rigid velocities, removed by the section, contradicting unit
length of `v`.

## The continuous critical family

Graph ``E`dg`` has contacts `(12,15,26,34,36,45,56)`. Fix
`c_5=(-1,0),c_6=(1,0)`. Its two equal-side four-cycles are rhombi, so every
realization has

```text
c_1=(-1,0)+2u, c_2=(1,0)+2u,
c_4=(-1,0)+2w, c_3=(1,0)+2w,  ||u||=||w||=1.
```

Writing `d=u dot w`, the strict noncontact gaps are
`8(1-d-u_x+w_x)`, `4(1-2d)` (twice), `4(1-2u_x)`,
`8(1-d+u_x-w_x)`, `4(1+2u_x)`, `4(1+2w_x)`, `4(1-2w_x)`.
Thus `-1/2<u_x,w_x<1/2`, `d<1/2` and `1-d>|u_x-w_x|` describe the
entire open packing domain. The seven contact rows are independent:
the horizontal and `u` directions at vertex 1 eliminate rows 12 and 15;
vertex 2 then eliminates 26. Vertices 4 and 3 similarly eliminate 34,45,36,
leaving 56.

Here `M=11-(8/3)(1+u dot w)`. Its two angular derivatives vanish exactly
when `det(u,w)=0`. The choice `w=u` overlaps centers; `w=-u` gives the
whole critical set. Set `a=-2u_x`, `s=2u_y`:

```text
-1<a<1, s=+/-sqrt(4-a^2),
c_1=(-1-a,s), c_2=(1-a,s), c_3=(1+a,-s),
c_4=(-1+a,-s), c_5=(-1,0), c_6=(1,0).
```

Both labeled signs are included. Exact polynomial/radical identities check
this domain, contact rank, stationarity and `M=11`. Rotating `u` while
holding `w=-u` fixed preserves contacts and locally strict noncontacts.
Along that curve `M(t)=11-(8/3)(1-cos t)`, so `M''(0)=-8/3`.
Every critical family member is excluded by this feasible negative variation.

## EPVW additional-contact root

Anchor `c_1=(0,0),c_3=(2,0)`. On free tree edges `1->6,3->4,6->2,6->5`,
the registered angles are `(2/3,2/3,2/3,1/3)pi`; reduced multipliers on
`25,45,46` are `(1/2,1,1/2)`. Exact verification gives `M=11`, reduced
KKT determinant `-4608` and additional contact `14`. The old-to-new
permutation `(4,1,2,6,3,5)` maps the full contact graph to EQlw and `14`
to `46`. When a Krawczyk uniqueness box contains this exact witness state,
its unique root belongs to EQlw, whose graph has its own analysis route.

## Closure

For five and six disks the census partitions as follows:

| Route | n=5 | n=6 |
|---|---:|---:|
| Rigid target | 1 | 1 |
| Continuous critical family | 0 | 1 |
| Interval search (including graphs with leaves) | 12 | 44 |
| Total | 13 | 46 |

Complete angular coverage, zero unresolved boxes and terminal analysis
exclusions exhaust every possible connected minimum outside the targets.
The component-rotation reduction then gives the global values `U_n`.
A disconnected equality case would, in its last merging rotation, give
disconnected packings of the same value approaching the connected target.
They are not congruent to it, contradicting its strict local minimality.
The verified targets therefore give every equality case for n=5,6.
