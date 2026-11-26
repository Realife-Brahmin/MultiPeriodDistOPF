# Supporting network tools

These optional helpers retain their original paths and inputs:

- `fetch_test_systems.sh`: downloads published upstream test-system sources into `test_networks/sources/`; requires Bash, curl, and an authenticated GitHub CLI for repository listings.
- `fetch_smartds_day.py`: retrieves one day of public SMART-DS P1U load profiles; Python standard library only. Default: day 148.
- `package_test_network.py`: exports an adapted network with shared profiles and attribution into a new ZIP without changing input files; Python standard library only.

Run from the repository root. Downloads depend on upstream availability. Review the [network provenance](../test_networks/README.md) and each script before use. These helpers are not required to browse the included simulation results.
