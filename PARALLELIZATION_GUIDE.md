# Parallelization Guide for tADMM Solver

## Quick Start (Local Machine)

### 1. Enable Threading
```bash
# Set number of threads (use your CPU core count)
export JULIA_NUM_THREADS=16

# Run Julia
julia tadmm_socp.jl
```

### 2. Enable Threading in Code
In `tadmm_socp.jl`, line ~10:
```julia
const USE_THREADING = true  # Change from false to true
```

### Expected Speedup
- 8-16× faster on typical workstation
- Linear scaling up to your CPU core count

---

## HPC Cluster Setup

### Check if Julia is Available
```bash
ssh your_hpc_cluster
module avail julia       # Check if Julia module exists
which julia              # Check system-wide installation
```

### Option A: HPC Has Julia Module
```bash
#!/bin/bash
#SBATCH --job-name=tadmm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=02:00:00
#SBATCH --mem=128GB

module load julia/1.10   # Load Julia module
export JULIA_NUM_THREADS=48

julia tadmm_socp.jl
```

### Option B: Install Julia Yourself (No Admin Needed!)
```bash
# 1. Download Julia
cd ~
wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.5-linux-x86_64.tar.gz

# 2. Extract
tar -xzf julia-1.10.5-linux-x86_64.tar.gz

# 3. Add to PATH (add to ~/.bashrc for permanent)
export PATH="$HOME/julia-1.10.5/bin:$PATH"

# 4. Verify
julia --version

# 5. Install packages (one time)
julia -e 'using Pkg; Pkg.add(["JuMP", "Gurobi", "PowerModelsDistribution"])'
```

### SLURM Script Example
```bash
#!/bin/bash
#SBATCH --job-name=tadmm_parallel
#SBATCH --nodes=1
#SBATCH --ntasks=48          # One task per time period
#SBATCH --time=01:00:00
#SBATCH --mem=256GB
#SBATCH --partition=compute  # Your cluster's partition name

# Load Julia (method depends on your HPC)
export PATH="$HOME/julia-1.10.5/bin:$PATH"
# OR: module load julia/1.10

# Set threading
export JULIA_NUM_THREADS=48

# Run
cd $SLURM_SUBMIT_DIR
julia tadmm_socp.jl

# Save output
echo "Job completed: $(date)"
```

Submit: `sbatch run_tadmm.sh`

---

## AWS Cloud Setup

### 1. Launch EC2 Instance
```bash
# Choose instance type based on cores needed
# c7i.24xlarge: 96 vCPUs ($3.40/hour, spot: $1.02/hour)
# c7i.48xlarge: 192 vCPUs ($6.80/hour, spot: $2.04/hour)

# Use AWS CLI or console to launch
aws ec2 run-instances \
  --image-id ami-0c55b159cbfafe1f0 \
  --instance-type c7i.24xlarge \
  --key-name your-key \
  --spot-instance-type "one-time"  # For cost savings
```

### 2. Setup Julia on EC2
```bash
# SSH to instance
ssh -i your-key.pem ubuntu@ec2-instance-ip

# Install Julia
wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.5-linux-x86_64.tar.gz
tar -xzf julia-1.10.5-linux-x86_64.tar.gz
sudo ln -s ~/julia-1.10.5/bin/julia /usr/local/bin/julia

# Install Gurobi (download from gurobi.com with academic license)
wget https://packages.gurobi.com/11.0/gurobi11.0.0_linux64.tar.gz
tar -xzf gurobi11.0.0_linux64.tar.gz
export GUROBI_HOME="$HOME/gurobi1100/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

# Get your academic license key (from gurobi.com)
grbgetkey YOUR-LICENSE-KEY

# Install Julia packages
julia -e 'using Pkg; Pkg.add(["JuMP", "Gurobi"])'
```

### 3. Transfer Your Code
```bash
# From your local machine
scp -i your-key.pem -r MultiPeriodDistOPF ubuntu@ec2-ip:~/

# OR use git
ssh -i your-key.pem ubuntu@ec2-ip
git clone https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
```

### 4. Run with Threading
```bash
export JULIA_NUM_THREADS=96  # Or whatever your instance has
cd MultiPeriodDistOPF
julia tadmm_socp.jl
```

### 5. Cost Optimization
- **Spot Instances**: 60-80% cheaper, can be interrupted
- **Auto-shutdown**: Set up automatic termination after job completes
  ```bash
  julia tadmm_socp.jl && sudo shutdown -h now
  ```
- **S3 for results**: Store outputs in S3, terminate instance immediately
  ```bash
  aws s3 cp results/ s3://your-bucket/results/ --recursive
  ```

---

## JuliaHub (Easiest for Julia)

### 1. Sign Up
Visit https://juliahub.com and create account

### 2. Upload Code
- Zip your project folder
- Upload via web interface

### 3. Configure Job
```julia
using JuliaHub

job = JuliaHub.submit_job(
    "tadmm_socp.jl",
    ncpu=48,
    memory=128,  # GB
    time=2  # hours
)

# Monitor
JuliaHub.job_status(job)

# Download results
JuliaHub.download_job_results(job, "results/")
```

**Advantages:**
- Pre-configured Julia environment
- No setup needed
- Pay per use
- Great for occasional large runs

---

## Performance Expectations

### Current (Sequential)
- T=48: ~10 seconds wall time (311s total subproblem time)
- T=96: ~20 seconds wall time (600s+ total subproblem time)

### With 48 Threads
- T=48: ~7-8 seconds (6.5s subproblem time, near-linear speedup)
- T=96: ~14-15 seconds (13s subproblem time, near-linear speedup)

### Why Near-Linear?
- Subproblems are **embarrassingly parallel** (completely independent)
- Only synchronization point is consensus update (negligible time)
- Gurobi solve is most expensive operation → parallelizes perfectly

---

## Troubleshooting

### "Thread count is 1"
```bash
# Check current threads
julia -e 'println(Threads.nthreads())'

# Set before starting Julia
export JULIA_NUM_THREADS=16
```

### "Gurobi license error"
- Academic license is machine-locked
- Need to re-register on HPC/cloud instance
- Get your license key: https://www.gurobi.com/academia/academic-program-and-licenses/
- Run: `grbgetkey YOUR-KEY`

### "Out of memory"
- Each subproblem needs ~1-2 GB
- For T=48 with 48 threads: Need 96+ GB RAM
- Reduce thread count or use higher-memory instance

### "Slower with threading"
- Gurobi might be using multiple threads per subproblem
- Set single-threaded Gurobi: In your code, add:
  ```julia
  set_optimizer_attribute(model, "Threads", 1)
  ```

---

## Next Steps

1. **Test locally first**: `USE_THREADING=true` with your laptop
2. **Contact HPC support**: Ask about Julia, most will help you install
3. **Start with small runs**: Test with T=24 before full T=96
4. **Benchmark**: Compare sequential vs parallel times
5. **Scale up**: Once working, use more cores for larger problems

---

## Questions?

Common issues are usually:
1. Thread count not set → No speedup
2. Gurobi using multiple threads → Oversubscription
3. Memory limits → Reduce parallel count

Good luck with your parallelization! 🚀
