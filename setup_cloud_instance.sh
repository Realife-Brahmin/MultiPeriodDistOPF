#!/bin/bash
# Cloud Instance Setup Script for tADMM
# Run this on a fresh AWS/GCP/Azure instance after SSH login

set -e  # Exit on error

echo "============================================"
echo "Setting up Julia + Gurobi for tADMM"
echo "============================================"

# 1. Install Julia
echo "[1/4] Installing Julia 1.10.5..."
cd ~
wget -q https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.5-linux-x86_64.tar.gz
tar -xzf julia-1.10.5-linux-x86_64.tar.gz
sudo ln -sf ~/julia-1.10.5/bin/julia /usr/local/bin/julia
rm julia-1.10.5-linux-x86_64.tar.gz
echo "✓ Julia installed: $(julia --version)"

# 2. Install Gurobi
echo "[2/4] Installing Gurobi 11.0..."
wget -q https://packages.gurobi.com/11.0/gurobi11.0.0_linux64.tar.gz
tar -xzf gurobi11.0.0_linux64.tar.gz
rm gurobi11.0.0_linux64.tar.gz

# Add to bashrc for permanent
echo 'export GUROBI_HOME="$HOME/gurobi1100/linux64"' >> ~/.bashrc
echo 'export PATH="${PATH}:${GUROBI_HOME}/bin"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"' >> ~/.bashrc
source ~/.bashrc

echo "✓ Gurobi installed"
echo "⚠  IMPORTANT: Get your academic license key from https://www.gurobi.com/academia/"
echo "   Then run: grbgetkey YOUR-LICENSE-KEY"

# 3. Clone repository (optional, uncomment if needed)
# echo "[3/4] Cloning repository..."
# git clone https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
# cd MultiPeriodDistOPF

# 4. Install Julia packages
echo "[3/4] Installing Julia packages..."
julia -e 'using Pkg; Pkg.add(["JuMP", "Gurobi", "Ipopt", "PowerModelsDistribution", "Plots", "Printf", "LinearAlgebra"])'
echo "✓ Julia packages installed"

# 5. Configure threading
NPROC=$(nproc)
echo "[4/4] Configuring for $NPROC cores..."
echo "export JULIA_NUM_THREADS=$NPROC" >> ~/.bashrc
source ~/.bashrc

echo ""
echo "============================================"
echo "Setup complete! ✓"
echo "============================================"
echo ""
echo "Next steps:"
echo "1. Get Gurobi license: grbgetkey YOUR-KEY"
echo "2. Upload your code: scp -r MultiPeriodDistOPF user@instance:~/"
echo "3. In tadmm_socp.jl, set: USE_THREADING = true"
echo "4. Run: julia tadmm_socp.jl"
echo ""
echo "Current threads: $JULIA_NUM_THREADS"
echo "Available cores: $NPROC"
echo ""
