# How to Use Zonal Plotting Functions

## Quick Start

Add these lines to your `multi_poi_mpopf.jl` after solving all slack configurations:

```julia
# Include the zonal plotter
include(joinpath(@__DIR__, "envs", "multi_poi", "ZonalPlotter.jl"))
using .ZonalPlotter

# Print zone statistics
print_zone_statistics(data)

# Generate zonal plots for each slack configuration
for slack_sub in slack_substations
    result = results_by_slack[slack_sub]
    
    # PV dispatch by zone
    plot_zonal_pv_dispatch(
        result, data, slack_sub,
        showPlots=showPlots,
        savePlots=savePlots,
        filename=joinpath(plots_dir, "zonal_pv_slack$(slack_sub).png")
    )
    
    # Battery actions by zone
    plot_zonal_battery_actions(
        result, data, slack_sub,
        showPlots=showPlots,
        savePlots=savePlots,
        filename=joinpath(plots_dir, "zonal_battery_slack$(slack_sub).png")
    )
end
```

## What You Get

### 1. Zone Statistics Table
Printed to console showing:
- Number of PV units per zone
- Number of battery units per zone
- Load pattern characteristics

### 2. Zonal PV Dispatch Plot
- Single plot with 5 curves (one per zone)
- Y-axis: PV Real Power [kW]
- X-axis: Time Period (t)
- Colors match your substation color scheme
- Shows which zones have PV resources

### 3. Zonal Battery Actions Plot
- Two-panel vertical layout (like your voltage plots)
- **Top panel**: Battery charging/discharging power by zone
  - Positive = discharging
  - Negative = charging
- **Bottom panel**: Aggregate state of charge (SOC) by zone
- Colors match your substation color scheme

## Understanding the Fluctuations

### Why Substation 1s has huge swings:
- **Zone 1**: 5 PV + 8 batteries + residential sinusoidal load
- Day: High PV export → large negative power
- Night: Pure load import → large positive power
- Battery optimization amplifies these swings

### Why Substations 2s and 3s are moderate:
- **Zone 2**: 3 PV + 4 batteries (moderate resources)
- **Zone 3**: 4 PV + 6 batteries (moderate resources)
- Load patterns (commercial sinusoid + industrial shifts) create moderate variability

### Why Substations 4s and 5s are stable:
- **Zone 4**: 5 PV + 8 batteries BUT office load perfectly matches PV timing (hours 8-18)
  → Local generation-load balance → minimal net exchange
- **Zone 5**: 0 PV + 0 batteries
  → Pure load-driven → very stable import

## Color Scheme

All plots use your consistent 5-zone color palette:
- Zone 1 (1s): Deep Blue #0072B2
- Zone 2 (2s): Amber #E69F00
- Zone 3 (3s): Teal #009E73
- Zone 4 (4s): Magenta #CC79A7
- Zone 5 (5s): Rust Red #D55E00

## File Locations

- **Zone mapping**: `envs/multi_poi/zone_analysis.md`
- **Plotting functions**: `envs/multi_poi/ZonalPlotter.jl`
- **Generated plots**: `envs/multi_poi/processedData/ieee123_5poi_1ph_T24/plots/zonal_*.png`
