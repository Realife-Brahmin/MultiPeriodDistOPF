# OpenDSS test networks

You can use these networks independently of the optimization code. They are
balanced single-phase research adaptations with loads, solar PV and batteries;
they are **not unmodified copies** of the original IEEE/EPRI feeders.

## Pick a network

- **Small starting case:** [ieee123 — 128 modeled buses](../rawData/ieee123C_1ph/).
  Balanced single-phase adaptation of the IEEE 123-node feeder, with added DERs.
- **Medium case:** [med2522 — 2,522 buses](../rawData/ieee2522C_1ph/).
  The stable filename is `ieee2522C_1ph`; this is not an official IEEE 2,522-node
  test feeder. Its provenance is a single-phase IEEE 8500 primary-network
  adaptation used in the D-OPF study cited below.
- **Large case:** [large10k — 10,321 buses](../rawData/large10kC_1ph/).
  Synthetic aggregate single-phase system. Read its
  [electrical corrections](../rawData/large10kC_1ph/CORRECTION.md) before comparing
  it with another version of the source system.
- **Multiple sources:** [small two-source example](../rawData/small2poi_1ph/)
  or [five-source IEEE123 adaptation](../rawData/ieee123_5poi_1ph/).
- **More networks:** [all adapted decks](../rawData/), including converted
  EPRI, IEEE 8500/9500 and SMART-DS cases. Their
  [published source decks and provenance](../test_systems/README.md) are separate.

The bus counts above describe these repository models, not the official feeder
names. Device placement, voltage bases, source settings and one-phase equivalents
are study choices. Check each deck and the relevant study before reusing it.

## Try one in OpenDSS

[Download the repository ZIP](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/archive/refs/heads/master.zip)
and extract it, or clone the repository. OpenDSS alone is sufficient; Julia,
Python and optimization-solver licenses are not needed to run the deck.

In OpenDSS, replace the path below with your extracted directory:

```text
Compile [C:/path/to/MultiPeriodDistOPF/rawData/ieee123C_1ph/Master.dss]
Solve
```

This runs the deck's configured power flow, not a battery-scheduling optimization.
Most decks redirect shared profiles from their parent directory: **do not
download only `Master.dss` or move its network folder on its own**.

To package just one adapted network with its profiles and citation, Python 3.10+
can create a self-contained ZIP (no extra packages):

```bash
python tools/package_test_network.py ieee123C_1ph --output ieee123-network.zip
```

Extract the ZIP, then compile `test_networks/ieee123C_1ph/Master.dss` inside it.
The packaged DSS/profile files are byte-identical to the repository files.
Existing output files are never overwritten.

### Why are the files still in rawData?

`test_networks/` is the human-facing catalogue; `rawData/` remains the physical
storage path used by existing parsers and recorded studies. This avoids breaking
working code or historical reproduction merely to improve discoverability.
The standalone ZIP uses the clearer `test_networks/` layout. Neither this guide
nor the packager changes the electrical data.

## Cite the networks you use

Please cite this repository:

> Aryan Ritwajeet Jha and contributors, *MultiPeriodDistOPF*,
> https://github.com/Realife-Brahmin/MultiPeriodDistOPF.

Also cite the underlying feeder source where applicable:

- IEEE123: the original IEEE 123-node test feeder and its applicable source notices.
- med2522: R. R. Jha et al., “Distribution Grid Optimal Power Flow (D-OPF):
  Modeling, Analysis, and Benchmarking,” *IEEE Transactions on Power Systems*,
  38(4), 3654–3668, 2023. [DOI](https://doi.org/10.1109/TPWRS.2022.3204227).
- large10k: R. Sadnan and A. Dubey, “Distributed Computing for Scalable Optimal
  Power Flow in Large Radial Electric Power Distribution Systems with Distributed
  Energy Resources,” HICSS-56, 2023.
  [Source publication](https://aisel.aisnet.org/hicss-56/es/monitoring/8/).
- Converted EPRI/IEEE/SMART-DS decks: see the
  [source provenance and citations](../test_systems/README.md).

[CITATION.cff](../CITATION.cff) supplies GitHub's repository citation metadata.
Please retain applicable upstream notices. The
citation request does not replace or add restrictions to the repository's
[MIT license](../LICENSE) or to third-party source terms.
