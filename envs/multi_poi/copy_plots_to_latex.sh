#!/bin/bash
# Script to copy plots from multi_poi simulation to LaTeX figures folder

# Define source and destination directories
PLOTS_SOURCE="/c/Users/aryan/Documents/documents_general/MultiPeriodDistOPF/envs/multi_poi/processedData/ieee123_5poi_1ph_T24/plots"
LATEX_DEST="/c/Users/aryan/Documents/documents_general/documentsCreated/PESGM2026_Multi-Source-Multi-Period-OPF/figures"

echo "=================================="
echo "Copying plots to LaTeX figures folder"
echo "=================================="
echo ""
echo "Source: $PLOTS_SOURCE"
echo "Destination: $LATEX_DEST"
echo ""

# Create destination if it doesn't exist
mkdir -p "$LATEX_DEST"

# Copy angle-only plots
echo "Copying angle-only plots..."
cp "$PLOTS_SOURCE"/angle_only_slack*.png "$LATEX_DEST/" && echo "  ✓ Copied 5 angle-only plots"

# Copy combined angle+voltage plots
echo "Copying angle+voltage plots..."
cp "$PLOTS_SOURCE"/angle_voltage_slack*.png "$LATEX_DEST/" && echo "  ✓ Copied 5 angle+voltage plots"

echo ""
echo "=================================="
echo "✓ All plots copied successfully!"
echo "=================================="
echo ""
echo "Latest plot timestamps:"
ls -lh "$LATEX_DEST"/angle_only*.png | awk '{print "  " $6, $7, $8, $9}' | head -5
