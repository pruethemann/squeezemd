#!/bin/bash
# Run a command pinned to one free GPU (lock-file pool).
# Usage: gpu_run.sh <command> [args...]
# Number of GPUs: $SQUEEZE_NGPU, else auto-detected via nvidia-smi.
NGPU=${SQUEEZE_NGPU:-$(nvidia-smi -L 2>/dev/null | wc -l)}
if [ "$NGPU" -lt 1 ]; then
    exec "$@"   # no GPU found: run unpinned (OpenMM falls back to CPU)
fi
while true; do
    for ((i = 0; i < NGPU; i++)); do
        exec 9>"/tmp/squeezemd_gpu_${USER:-$(id -un)}_$i.lock"
        if flock -n 9; then
            export CUDA_VISIBLE_DEVICES=$i
            echo "[gpu_run] using GPU $i" >&2
            exec "$@"   # fd 9 is inherited, so the lock is held until the job exits
        fi
    done
    sleep 10
done
