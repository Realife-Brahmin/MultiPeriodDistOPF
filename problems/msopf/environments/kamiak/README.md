# Recorded sf6937 Kamiak environments

These are the original, byte-preserved environments used for the sf6937 reference campaign, with Julia **1.12.7**. They are separate because the jobs selected different environments:

- `sdp_hsl`: AC edge searches and idle-certificate phases (A, AD, C and V).
- `sdp_pardiso`: whole-day certificate and audit phases (CD and VD), using Clarabel's `mkl` backend.

Both environments pin JuMP 1.30.0, MathOptInterface 1.50.0, Ipopt 1.14.1, Clarabel 0.11.1 and OpenDSSDirect 0.9.9. The whole-day environment adds Pardiso.jl 1.1.2, MKL_jll 2025.2.0+0 and IntelOpenMP_jll 2025.2.0+0. The manifests contain no machine-local package paths. Exact file hashes and the full package inventory are in the [canonical configuration](../../results/kamiak/farkas_day148b_nov03/reference_config/canonical_configuration.json).

## Install the recorded package environment

From the repository root, in Bash or Git Bash:

```bash
julia --startup-file=no --project=problems/msopf/environments/kamiak/sdp_hsl -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=problems/msopf/environments/kamiak/sdp_pardiso -e 'using Pkg; Pkg.instantiate()'
```

Use Julia 1.12.7 and preserve the manifests. Instantiation may download package artifacts. Solver libraries and their applicable licence/setup requirements remain separate from the package snapshot.

## Runtime details outside the Manifest

The original Ipopt jobs selected a private MA57 **3.11.3** build through `hsllib`, with AMD ordering and a preload library for METIS symbols. Although HSL.jl 0.5.6 / HSL_jll 4.0.8+0 are installed, that JLL was not the MA57 backend loaded by the recorded jobs. No HSL binaries are redistributed here. Configure your own compatible, appropriately licensed library rather than copying the author's home-directory paths.

The OpenDSS engine recorded in the reference audit is DSS C-API **0.14.3**. OpenDSSDirect's version alone does not identify the loaded engine. Gurobi 1.9.2 / Gurobi_jll 13.0.0 are installed but were not used in these reference phases.

CD selected `FB_ENV=sdp_pardiso FB_LDL=mkl NTH=8`. Julia, OpenBLAS and OpenMP thread environment variables were 1; MKL_NUM_THREADS was 8. Effective Pardiso thread count was not logged, so eight allocated CPUs must not be described as eight busy solver threads.

These files preserve the recorded software state. Release preparation checked their hashes and dependency UUID consistency; it did not reinstall the environments or rerun the original optimization campaign.
