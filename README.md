# Multiwfn-SLF — Stress Localization Function and related stress-trace descriptors

Modified build of [Multiwfn](http://sobereva.com/multiwfn/) adding the Stress
Localization Function (SLF) and a set of related real-space functions derived
from the trace of the quantum-mechanical stress tensor.

This code accompanies:

> M\. Husain and T. Korona, *The Stress Localization Function — a
> Stress-Derived Real-Space Descriptor for Chemical Bonding Analysis*
> (Accepted - Journal of Chemical Theory and Computation).

---

## What is added

All additions are exposed through Multiwfn's user-defined real-space function
mechanism (function index `100`), selected by the `iuserfunc` index.

| `iuserfunc` | Symbol | Quantity |
|---|---|---|
| `201` | `SLF(r)` | Stress Localization Function, `1/(1 + chi_sigma^2)` |
| `207` | `t^s_{-1/8}(r)` | Scaled stress trace, `(G + lambda grad^2 rho)/rho^(5/3)` at `lambda = -1/8` |
| `208` | `t^s_0(r)` | Scaled stress trace at `lambda = 0` |
| `209` | `chi_sigma(r)` | Normalized isotropic stress, `|Tr sigma| / sigma_homog` |
| `210` | `chi_G(r)` | Normalized kinetic contribution |
| `211` | `chi_L(r)` | Normalized Laplacian contribution |
| `212` | `chi_sigma(r)` | Same as 209, evaluated from the decomposition `|chi_G - chi_L|` |

### Definitions

```
Tr sigma(r)    = (1/4) grad^2 rho(r) - 2 G(r)
sigma_homog(r) = 2 C_F rho(r)^(5/3),    C_F = (3/10)(3 pi^2)^(2/3) = 2.871234
chi_sigma(r)   = |Tr sigma(r)| / sigma_homog(r)
SLF(r)         = 1 / (1 + chi_sigma(r)^2)
```

with `G(r)` the positive-definite kinetic energy density. `chi_sigma = 1` and
`SLF = 1/2` for the uniform electron gas.

### Decomposition

```
chi_G(r)      = G(r) / (C_F rho(r)^(5/3))
chi_L(r)      = grad^2 rho(r) / (8 C_F rho(r)^(5/3))
chi_sigma(r)  = |chi_G(r) - chi_L(r)|

SLF(r)        = 1 / (1 + (chi_G(r) - chi_L(r))^2)
t^s_lambda(r) = C_F (chi_G(r) + 8 lambda chi_L(r))
```

Functions `209` and `212` compute `chi_sigma(r)` by independent routes —
directly from the trace, and from the decomposition — and must agree to machine
precision.

---

## Modified files

| File | Change |
|---|---|
| `function.f90` | `SLF_getD`, `Finzel_get_ts`, `SLF_components`; new `case` branches in `userfunc` |
| `topology.f90` | option `8` in the topology menu, *Output SLF values at all CPs* (writes `CPSslf.txt`) |

`settings.ini` is unchanged; `iuserfunc` is set by the user.

---

## Building

The build is unchanged from stock Multiwfn. Intel Fortran and MKL are required.

```bash
cd /path/to/multiwfn-slf
make -j4
```

The binary is written to the Multiwfn root directory. Set `Multiwfnpath` to that
directory so the runtime files are found:

```bash
export Multiwfnpath=/path/to/multiwfn-slf
export PATH=$Multiwfnpath:$PATH
```

---

## Running

### 1. Select the function

Set the index in `settings.ini` before starting Multiwfn:

```ini
iuserfunc= 201
```

Then choose `100` (*User-defined function*) wherever Multiwfn asks which
real-space function to use. The menu echoes the active index:

```
100 User-defined function (iuserfunc=  201), see Section 2.7 of manual
```

### 2. Common analyses

```bash
Multiwfn molecule.wfx
```

| Task | Main-menu path |
|---|---|
| Value at a point | `1` → enter coordinates |
| Profile along a line | `3` → function `100` → define the two endpoints |
| Plane map | `4` → function `100` |
| Grid / cube file | `5` → function `100` → export `.cub` |
| Critical points | `2` → `-11` → `100` → search (below) |

### 3. Topology analysis

```
2            ! Topology analysis       (all lengths in this module are Bohr)
-11          ! Select real space function
100          ! User-defined function  -> picks up iuserfunc = 201
2            ! Search CPs from nuclear positions
3            ! Search CPs from midpoints of atomic pairs
8            ! Output SLF values at all CPs  -> CPSslf.txt
```

Searches are cumulative; run as many as the system needs, then export once.

| Option | Starting points | Use for |
|---|---|---|
| `2` | nuclear positions | nuclear-site CPs |
| `3` | midpoints of atom pairs | bonding attractors |
| `4` | triangle centres of atom triples | three-centre bonds |
| `6` | batch of points within sphere(s) | atomic shell structure, lone pairs |

Switching the real-space function with `-11` resets some CP-search and
path-generation parameters. Once CPs exist the menu entry reads *Delete results
and reselect real space function*, so a further `-11` discards them.

**Search algorithm.** Set under `-1 Set CP searching parameters`. Newton locates
all CP types; steepest ascent finds only maxima and steepest descent only
minima.

**Numerical derivatives.** Multiwfn has analytic derivatives only for a fixed
list of built-in functions, which does not include user-defined ones. Topology
analysis of these functions is therefore numerical: slower, and subject to
Multiwfn's own warning that accuracy cannot be guaranteed. Tighten the
convergence criteria under `-1` if attractor positions look unstable.


---

## Licensing

Multiwfn is the work of Tian Lu. Its license terms permit free distribution of
both the original and modified code, for academic and commercial use alike, and
this repository is distributed on that basis. The terms, reproduced from the
Multiwfn distribution:

> \(a\) Currently Multiwfn is free of charge and open-source for both academic
> and commercial usages, anyone is allowed to freely distribute the original or
> their modified Multiwfn codes to others.
>
> \(b\) Multiwfn can be distributed as a free component of commercial code.
> Selling modified version of Multiwfn may also be granted, however, obtaining
> prior consent from the original author of Multiwfn (Tian Lu) is needed.
>
> \(c\) If Multiwfn is utilized in your work, or your own code incorporated any
> part of Multiwfn code, at least the following original papers of Multiwfn MUST
> BE cited in main text of your paper or code:
> Tian Lu, Feiwu Chen, *J. Comput. Chem.*, **33**, 580–592 (2012);
> Tian Lu, *J. Chem. Phys.*, **161**, 082503 (2024).

The modifications in this repository are made available on the same terms.

---

## Citation

Clause (c) is a condition of the license: **both** Multiwfn papers must be cited
in the main text of any work using this code.

```bibtex
@article{Multiwfn2012,
  author  = {Lu, Tian and Chen, Feiwu},
  title   = {Multiwfn: A multifunctional wavefunction analyzer},
  journal = {J. Comput. Chem.},
  volume  = {33},
  pages   = {580--592},
  year    = {2012},
  doi     = {10.1002/jcc.22885}
}

@article{Multiwfn2024,
  author  = {Lu, Tian},
  title   = {A comprehensive electron wavefunction analysis toolbox for
             chemists, Multiwfn},
  journal = {J. Chem. Phys.},
  volume  = {161},
  number  = {8},
  pages   = {082503},
  year    = {2024},
  doi     = {10.1063/5.0216272}
}

@article{SLF,
  author  = {Husain, Munavvar and Korona, Tatiana},
  title   = {The Stress Localization Function — a Stress-Derived Real-Space
             Descriptor for Chemical Bonding Analysis},
  journal = {accepted - Journal of Chemical Theory and Computation},
  year    = {2026}
}
```

## Contact

Munavvar Husain, Faculty of Chemistry, University of Warsaw.
