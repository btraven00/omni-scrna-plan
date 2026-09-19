#!/bin/bash
# ==============================================================================
# SLURM JOB CONFIGURATION
# Edit the #SBATCH lines below to control resources.
# ==============================================================================

#SBATCH --job-name=mark-testing      # Job name (shows in squeue)
#SBATCH --output=logs/%x_%j.out      # stdout log (%x=job name, %j=job ID)
#SBATCH --error=logs/%x_%j.err       # stderr log

# --- Time & partition ---
#SBATCH --time=04:00:00              # Max wall time (HH:MM:SS)
#SBATCH --partition=main             # Partition/queue name

# --- CPUs & memory ---
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=1                   # Number of tasks (MPI ranks)
#SBATCH --cpus-per-task=50           # CPU cores per task
#SBATCH --mem=100G                   # Total memory per node (e.g. 8G, 32G)
# --mem-per-cpu=4G                   # Alternative: memory per CPU core

# --- GPUs (comment out if not needed) ---
#SBATCH --gpus 1                     # Number of GPUs (e.g. gpu:1, gpu:a100:2)
# #SBATCH --constraint=L4              # Specific GPU type (cluster-dependent)

# --- Notifications (optional) ---
# #SBATCH --mail-type=END,FAIL       # Email on job end or failure
# #SBATCH --mail-user=mark.robinson@mls.uzh.ch  # Your email address

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

# Exit immediately if any command fails
set -euo pipefail

echo "=============================================="
echo "Job:       $SLURM_JOB_NAME ($SLURM_JOB_ID)"
echo "Node:      $SLURMD_NODENAME"
echo "Started:   $(date)"
echo "Directory: $(pwd)"
echo "=============================================="

# --- Conda activation ---
# Slurm jobs often start with a bare environment, so we source conda explicitly.
# Adjust the path below to match your conda/mamba installation.

CONDA_ENV="omni-gather"                   # <-- set your environment name here

# Try common conda installation paths; adjust if yours differs
CONDA_BASE="$HOME/miniforge3"

if [ -z "$CONDA_BASE" ]; then
    echo "ERROR: Could not find a conda installation. Set CONDA_BASE manually."
    exit 1
fi

# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

echo "Python:    $(which python)"
echo "Conda env: $CONDA_DEFAULT_ENV"
echo "=============================================="

# ==============================================================================
# CREATE LOG DIRECTORY (if it doesn't exist)
# ==============================================================================

mkdir -p logs

# ==============================================================================
# YOUR COMMANDS GO HERE
# ==============================================================================

# Example: run a script
ob run benchmark_minimal.yaml -c 48 --with-capability gpu --unpinned

# ==============================================================================
echo "Finished:  $(date)"
# ==============================================================================

