# IEEE 9500-node test system, original files (unmodified)

Downloaded 2026-10-20 from the GridAPPS-D CIMHub repository (PNNL),
https://github.com/GRIDAPPSD/CIMHub, branch `feature/SETO`, folder `ieee9500/`:
`README.rst` and `base/` as they are there. Copyright Battelle Memorial Institute;
redistribution with or without modification is permitted under the BSD-style terms in
`LICENSE.md` (CIMHub's `license.md`), which must stay with these files.

Three radial 12.47 kV feeders, each fed by its own 69/12.47 kV substation (S1, S2, S3),
joined by normally-open switches; 2,702 primary buses at 12.47 kV, 1,275 customer loads
(13.7 MW balanced), 9 single-phase regulator banks, 10 capacitor banks.

`Generators.dss` uses `PVSystem ... model=7`, which the OpenDSS in OpenDSSDirect 0.9.9 rejects;
the converter leaves the built-in DERs out anyway (we place our own).

The single-phase deck made from these is `../ieee9500_1ph`, written by
`envs/multi_poi/root_level/convert_to_1ph.jl ieee9500` and checked by `check_1ph.jl ieee9500`.
