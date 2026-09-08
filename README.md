# GTEMPO.jl

GTEMPO is a Julia package that solves quantum impurity problems with
**Grassmann temporal basis representations**, a tensor-network formulation of
the influence-functional (IF) method built on Grassmann-number-valued matrix
product states (GMPS). GTEMPO is a fermionic analogy of the time-evolving
matrix product operator (TEMPO) method, where the latter is aimed for bosonic
impurity problems. It supports imaginary-time (Matsubara), real-time
(Keldysh) and mixed-time (Kadanoff–Baym) contours, multi-orbital impurities,
bosonic (electron–phonon) and superconducting (BCS) baths.

## Features

- **Contours**: imaginary time, real time, mixed (Kadanoff–Baym) time, and
  retarded-interaction orderings; ~25 Grassmann orderings with automatic
  conversion between them.
- **IF construction algorithms**: `PartialIF`, `XTRGIF`, `ExactTTIIF` and
  `TDVPIF`, with exact (zip-up) and boundary-MPS integration of GMPS products.
- **Observables**: Matsubara/real-time Green's functions, greater/lesser GFs,
  occupations, particle and heat currents (also in cached and fast MPO
  variants), with cached evaluation via environment tensors.
- **Impurities**: single-orbital Anderson, multi-orbital Kanamori, IRLM,
  electron–phonon (Fock-lattice) impurities, BCS baths.
- **Reference solvers**: exact-diagonalization references (`Toulouse`,
  `BoundaryDriving`, free-fermion solutions) from
  [ImpurityModelBase.jl](https://github.com/example/ImpurityModelBase.jl).

## Installation

GTEMPO depends on several in-development companion packages that are not yet
registered (`Z2Tensors`, `ImpurityModelBase`, `QuAPI`, `ExpExp`). Clone them
alongside GTEMPO and add GTEMPO with a path-based dependency:

```julia
using Pkg
Pkg.develop(path = "path/to/GTEMPO")
```

A prepared `Manifest.toml` resolving the path dependencies ships with the
repository.

## Quick start

Solve the single-impurity Anderson model on the Keldysh contour and measure
the real-time Green's function:

```julia
using GTEMPO

# bath and impurity
bath  = fermionicbath(semicircular(5.0), β = 10.0, μ = 0.0)
model = AndersonIM(U = 2.0, μ = -1.0)

# Grassmann lattice on the Keldysh contour
lattice = GrassmannLattice(N = 10, δt = 0.05, contour = :real)

# bath correlation function and influence functional
corr = correlationfunction(bath, lattice)
I    = hybriddynamics(lattice, corr, ExactTTIIF(algmult = SVDCompression(truncdimcutoff(D = 100, ϵ = 1e-9))))

# impurity dynamics and boundary condition
K = sysdynamics(lattice, model, trunc = truncdimcutoff(D = 100, ϵ = 1e-9))
K = boundarycondition!(K, lattice)

# observables
gt = [greater(lattice, t, K, I) for t in 1:lattice.k]     # G^>(t, 0)
lt = [lesser(lattice, t, K, I) for t in 1:lattice.k]      # G^<(t, 0)
```

See `docs/tutorials/` for complete single-orbital and multi-orbital examples
(Matsubara, Keldysh and Kadanoff–Baym setups).

## Package structure

```
src/
├── grassmanntensor/     # Grassmann-valued tensors (Z2 symmetric)
├── grassmannmps/        # Grassmann MPS: algebra, canonicalization, multiplication
├── lattices/            # Grassmann lattices and orderings (imag/real/mixed contours)
├── integration/         # exact, BMP (zip-up) and cached integration of GMPS products
├── influencefunctional/ # IF construction: PartialIF / XTRGIF / ExactTTIIF / TDVPIF
├── sysdynamics/         # impurity dynamics (Anderson, Kanamori, IRLM, ...)
├── gvconnections/       # boundary / bulk Grassmann-number connections
├── observables/         # GFs, occupations, currents (plain, cached, fast variants)
├── mpo/                 # MPO Hamiltonians and Schur-type MPO tensors
├── electronphonon/      # Fock-lattice electron–phonon path
└── bcsinfluencefunctional/ # BCS-bath influence functionals
```

## Documentation

Build the documentation with [Documenter.jl](https://documenter.juliadocs.org):

```julia
julia --project=docs docs/make.jl
```

The documentation source lives in `docs/src/`; the working tutorials are in
`docs/tutorials/`.

## Testing

```julia
julia --project=. test/runtests.jl
```

The suite is organized as (i) API unit tests, (ii) functional tests against
exact-diagonalization and analytic references (few-mode and continuous baths,
BCS and electron–phonon paths, multi-orbital impurities), including stepwise
evolution and transport benchmarks.

## License

GTEMPO is distributed under the terms of the license in [LICENSE](LICENSE).

## Citation

If you use GTEMPO in your research, please cite one or some of the following papers:

```bibtex
@article{chen2024gtempo,
  title   = {Grassmann time-evolving matrix product operators for quantum impurity models},
  author  = {Chen, Ruofan and Xu, Xiansong and Guo, Chu},
  journal = {Phys. Rev. B},
  volume  = {109},
  pages   = {045140},
  year    = {2024},
  doi     = {10.1103/PhysRevB.109.045140}
}
@article{chen2024equilibrium,
  title   = {Grassmann time-evolving matrix product operators for equilibrium quantum impurity problems},
  author  = {Chen, Ruofan and Xu, Xiansong and Guo, Chu},
  journal = {New J. Phys.},
  volume  = {26},
  pages   = {013019},
  year    = {2024},
  doi     = {10.1088/1367-2630/ad19fa}
}
@article{chen2024realtime,
  title   = {Real-time impurity solver using Grassmann time-evolving matrix product operators},
  author  = {Chen, Ruofan and Xu, Xiansong and Guo, Chu},
  journal = {Phys. Rev. B},
  volume  = {109},
  pages   = {165113},
  year    = {2024},
  doi     = {10.1103/PhysRevB.109.165113}
}
@article{guo2024influencefunctional,
  title   = {Efficient construction of the Feynman-Vernon influence functional as matrix product states},
  author  = {Guo, Chu and Chen, Ruofan},
  journal = {SciPost Phys. Core},
  volume  = {7},
  pages   = {063},
  year    = {2024},
  doi     = {10.21468/SciPostPhysCore.7.3.063}
}
@article{chen2024kadanoffbaym,
  title   = {Solving quantum impurity problems on the L-shaped Kadanoff-Baym contour},
  author  = {Chen, Ruofan and Guo, Chu},
  journal = {Phys. Rev. B},
  volume  = {110},
  pages   = {165114},
  year    = {2024},
  doi     = {10.1103/PhysRevB.110.165114}
}
@article{xu2024pathintegral,
  title   = {Grassmann time-evolving matrix product operators: An efficient numerical approach for fermionic path integral simulations},
  author  = {Xu, Xiansong and Guo, Chu and Chen, Ruofan},
  journal = {J. Chem. Phys.},
  volume  = {161},
  pages   = {151001},
  year    = {2024},
  doi     = {10.1063/5.0226167}
}
@article{chen2025polaron,
  title   = {Tensor network algorithm to solve polaron impurity problems},
  author  = {Chen, Ruofan and Gu, Lei and Guo, Chu},
  journal = {Chin. Phys. Lett.},
  volume  = {42},
  pages   = {120701},
  year    = {2025},
  doi     = {10.1088/0256-307X/42/12/120701}
}
@article{sun2025scalable,
  title   = {Scalable tensor network algorithm for quantum impurity problems},
  author  = {Sun, Zhijie and Chen, Ruofan and Li, Zhenyu and Guo, Chu},
  journal = {Phys. Rev. B},
  volume  = {112},
  pages   = {155115},
  year    = {2025},
  doi     = {10.1103/PhysRevB.112.155115}
}
@misc{guo2026superconducting,
  title         = {Grassmann time-evolving matrix product operators for fermionic impurities coupled to a superconducting bath},
  author        = {Guo, Chu and Wu, Wei and Xu, Xiansong and Chen, Ping-Xing and Yue, Changming and Jiang, Tian and Chen, Ruofan},
  year          = {2026},
  eprint        = {2604.23301},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.str-el}
}
```
