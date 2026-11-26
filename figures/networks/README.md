# Network diagrams: editable SVGs and PDF exports

These illustrations accompany the feeder models in ../../test_networks/. SVGs are
editable sources (Inkscape); PDFs are vector exports for manuscripts.

- `ieee123C_1ph`: ieee123, 128 modeled buses.
- `ieee2522C_1ph`: med2522, 2,522 buses. The raw dataset filename is retained.
- `large10kC_1ph_top_level`: large10k overview.
- `large10kC_1ph_low_level`: detail of a constituent feeder.

The diagrams illustrate topology and study device placement; they
are not optimized dispatches. In particular, the medium diagram retains
root devices that are excluded from the current exported optimization.
Illustrations and original feeder datasets have distinct provenance.

## Reuse and citation

Please cite this repository when using the editable Inkscape SVGs or their
exported figures:

> Aryan Ritwajeet Jha and contributors, *MultiPeriodDistOPF*,
> https://github.com/Realife-Brahmin/MultiPeriodDistOPF.

Retain applicable notices; see [LICENSE](../../LICENSE). The medium and large
network source publications are listed below.

Exploratory illustration utilities are under ../../scratch/. Editing a
committed SVG and exporting its drawing area is separate from regenerating
a layout with those utilities.

## Underlying network publications

- med2522: R. R. Jha, A. Inaolaji, B. D. Biswas, A. Suresh, A. Dubey,
  S. Paudyal and S. Kamalasadan, "Distribution Grid Optimal Power Flow
  (D-OPF): Modeling, Analysis, and Benchmarking," IEEE Transactions on Power
  Systems, 38(4), 3654–3668, 2023.
  https://doi.org/10.1109/TPWRS.2022.3204227
- large10k: R. Sadnan and A. Dubey, "Distributed Computing for Scalable
  Optimal Power Flow in Large Radial Electric Power Distribution Systems
  with Distributed Energy Resources," HICSS-56, 2023.
  https://aisel.aisnet.org/hicss-56/es/monitoring/8

The ieee123 diagram depicts this repository's balanced single-phase adaptation
of the IEEE 123-node feeder. The filename is not a claim that its modeled
bus count or added DER placements equal the original feeder specification.
Figure credit and test-system citation serve different purposes.
