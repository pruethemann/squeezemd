#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import tempfile

import pandas as pd

# import your parser helpers
from .transform_contact_data import parse_lines  # adjust import path if needed



def run_one_frame(extract_bin, posco_bin, topo, traj, frame, sequence_path, complex, mutation, seed):

    sequence = pd.read_parquet(sequence_path)
    if not sequence.index.is_unique:
        sequence = sequence[~sequence.index.duplicated(keep='first')]


    tmpdir = tempfile.mkdtemp(prefix="posco_stream_")
    lig_fifo = os.path.join(tmpdir, "lig.pdb")
    rec_fifo = os.path.join(tmpdir, "rec.pdb")
    os.mkfifo(lig_fifo)
    os.mkfifo(rec_fifo)

    try:
        # Start po-sco in background (it will read FIFOs)
        psco = subprocess.Popen(
            [posco_bin, rec_fifo, lig_fifo, "-b"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        # Write the frame’s receptor/ligand PDB into the FIFOs (no files persisted)
        extract = subprocess.run(
            [
                extract_bin,
                "--topo", topo,
                "--traj", traj,
                "--frame", str(frame),
                "--lig_frame", lig_fifo,
                "--rec_frame", rec_fifo,
                "--sequence", sequence_path,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        if extract.returncode != 0:
            stderr = extract.stderr.strip()
            raise RuntimeError(f"extract-contact-frames failed for frame {frame}:\n{stderr}")

        out, err = psco.communicate()
        if psco.returncode != 0:
            err = (err or "").strip()
            raise RuntimeError(f"po-sco failed for frame {frame}:\n{err}")

        lines = out.splitlines(True)  # keep line endings
        df = parse_lines(lines, sequence, complex, mutation, seed, frame)
        return df

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

def parse_arugments():
    ap = argparse.ArgumentParser()
    ap.add_argument("--topo", required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--n-frames", type=int, required=True)
    ap.add_argument("--complex", dest="complex_name", required=True)
    ap.add_argument("--mutation", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--extract-bin", default="extract-contact-frames")
    ap.add_argument("--posco-bin", default="po-sco")
    return ap.parse_args()


def main():

    args = parse_arugments()

    args.sequence = "sequence.parquet"

    dfs = []
    for frame in range(args.n_frames):
        dfs.append(
            run_one_frame(
                args.extract_bin,
                args.posco_bin,
                args.topo,
                args.traj,
                frame,
                args.sequence,
                args.complex_name,
                args.mutation,
                args.seed,
            )
        )

    res = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    res.to_parquet(args.out, index=False)


if __name__ == "__main__":
    main()