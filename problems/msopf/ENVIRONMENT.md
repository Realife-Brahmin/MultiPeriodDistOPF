# MSOPF environments

The canonical reference is the configuration used for **sf6937 on Kamiak**. Its original Project/Manifest pairs are tracked under [`environments/kamiak/`](environments/kamiak/README.md): `sdp_hsl` for AC and idle-certificate phases, and `sdp_pardiso` for whole-day CD/VD phases. Julia is 1.12.7. Follow that README for package versions, installation and external-library qualifications.

The root `problems/msopf/Project.toml` is retained for the general/historical drivers. It is not the complete environment used for the November certificate campaign. The older branch-flow project remains under `archive/legacy-environment/`. Do not mix these environments when claiming exact reproduction.

Earlier work used the environment copied from `problems/mpopf/tadmm` in October 2026. Its historical Manifest pins Julia 1.12.5, JuMP 1.30.0, Ipopt 1.14.1 and OpenDSSDirect 0.9.9, but does not provide the November Clarabel/PARDISO environment.

The supplied 309 snapshot matches Kamiak's `sdp_hsl` byte-for-byte. The ten idle-edge directions run on that Windows host retain their original provenance; they do not redefine the canonical Kamiak setup. The reference audit records a Windows MA57 build for those tasks.

Do not update or resolve a recorded campaign environment as part of reproducing it. Use a separate environment for experiments with newer packages. Hardware, solver-library choices, actual runtime backends and the scope of exact checks are recorded in [RELEASE.md](RELEASE.md) and the archived canonical configuration.
