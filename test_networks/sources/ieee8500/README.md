# IEEE 8500-node test feeder (OpenDSS)

- Source: https://github.com/GRIDAPPSD/Powergrid-Models, branch `cimhub_refactor`, folder
  `models/feeders/OpenDSS/IEEE/IEEE8500` (the default branch keeps only CIM XML), downloaded
  2026-10-20, unmodified.
- Masters used: `Master.dss` (balanced loads) and `Master-unbal.dss`.
- File fixes: none. The solve needs more than OpenDSS's default 15 iterations (67 here);
  `Run_8500Node.dss` itself sets `Maxiterations=20`.
- Not downloaded: IEEE8500_3subs, published only as CIM XML (GRIDAPPSD/CIMHub
  `model_output_tests/IEEE8500_3subs.xml`); its successor, the 9500-node system, is in
  `rawData/ieee9500_source`.
- Citation: R. F. Arritt and R. C. Dugan, "The IEEE 8500-node test feeder," IEEE PES T&D 2010.
- Snapshot results: `summary.csv`.
