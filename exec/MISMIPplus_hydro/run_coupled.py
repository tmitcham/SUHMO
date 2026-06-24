#!/usr/bin/env python3
"""
Coupling script for the MISMIP+ BISICLES–SUHMO coupled hydrology example.

Workflow
--------
This script orchestrates an iterative coupling between BISICLES (ice dynamics)
and SUHMO (subglacial hydrology).  Each coupling iteration proceeds as follows:

  1. BISICLES advances the ice geometry for `coupling_steps` timesteps
     (restarting from its own checkpoint on iterations > 0).
     BISICLES writes plot files containing ice thickness and velocity.

  2. The latest BISICLES plot file is symlinked to a fixed filename
     so SUHMO can always find the current ice geometry at a known path.

  3. SUHMO integrates the subglacial hydrology equations for ~1 yr of
     simulated time (set by maxStep in input.hydro), converging to a
     near-steady drainage state for the current ice geometry.
     SUHMO writes the effective pressure field N to a fixed output file.

  4. On the next BISICLES iteration the coupling script injects N from
     SUHMO via the LevelData effective-pressure option.

Data exchange files
-------------------
  BISICLES -> SUHMO  :  bisicles_for_suhmo.2d.hdf5
                        (symlink to the latest BISICLES plot file;
                         contains: thickness, xVel, yVel)

  SUHMO -> BISICLES  :  suhmo_for_bisicles.2d.hdf5
                        (written by SUHMO each iteration;
                         contains: effectivePressure)

Spin-up
-------
Before running this script, BISICLES should be spun up to a near-steady
grounded-ice geometry using hydrostatic effective pressure:

    mpirun -n <N> <bisicles_exe> inputs.bisicles

Run for ~10 000 yr (adjust main.maxTime in inputs.bisicles) until the
grounding-line position converges.  Then use the final checkpoint as the
--bisicles-restart argument to this script.

Usage
-----
    # Quick test (4 MPI ranks, 5 coupling iters of 5 BISICLES steps each):
    python run_coupled.py --nprocs 4 --coupling-steps 5 --num-coupling-iters 5

    # Production run with custom executables and spin-up restart:
    python run_coupled.py \\
        --nprocs 16 \\
        --bisicles-exe /path/to/driver2d.ex \\
        --suhmo-exe    /path/to/Suhmo2d.ex  \\
        --coupling-steps 20 \\
        --num-coupling-iters 50 \\
        --bisicles-restart chk.bisicles.000100.2d.hdf5
"""

import argparse
import glob
import os
import shutil
import sys
import subprocess


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_latest_file(prefix, ext="2d.hdf5"):
    """Return the lexicographically last file matching <prefix>*.<ext>."""
    files = sorted(glob.glob(f"{prefix}*.{ext}"))
    return files[-1] if files else None


def symlink_force(src, dst):
    """Create or replace symlink dst -> src."""
    if os.path.islink(dst) or os.path.exists(dst):
        os.remove(dst)
    os.symlink(src, dst)
    print(f"[COUPLED] linked {dst} -> {src}")


def run_cmd(cmd, label):
    """Print and execute a shell command; abort on non-zero exit."""
    print(f"\n{'='*60}")
    print(f"[COUPLED] {label}")
    print(f"  cmd: {cmd}")
    print(f"{'='*60}\n")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"[COUPLED] ERROR: '{label}' exited with code {result.returncode}")
        sys.exit(result.returncode)


# ---------------------------------------------------------------------------
# BISICLES
# ---------------------------------------------------------------------------

def write_bisicles_input(args, coupling_iter, restart_chk):
    """
    Write a per-iteration BISICLES input file by copying the base file and
    appending coupling overrides.  Returns the temporary file path.
    """
    tmp = f"inputs.bisicles.iter{coupling_iter:04d}"
    shutil.copy2(args.bisicles_input, tmp)

    with open(tmp, "a") as f:
        f.write(f"\n# --- coupling iteration {coupling_iter} overrides ---\n")

        # Advance to the end of this coupling window
        max_step = args.coupling_steps * (coupling_iter + 1)
        f.write(f"main.maxStep = {max_step}\n")

        # Write a checkpoint at the end of each window so SUHMO can restart
        f.write(f"amr.check_interval = {args.coupling_steps}\n")

        # Restart from the previous BISICLES checkpoint (not on first iter)
        if restart_chk is not None:
            f.write(f"amr.restart_file = {restart_chk}\n")

        # Use effective pressure from SUHMO once it has produced output.
        # On iteration 0 the base input file defaults to hydrostatic N.
        if os.path.exists(args.s2b_file):
            f.write(f"main.effectivePressure = LevelData\n")
            f.write(f"LevelDataEffectivePressure.file = {args.s2b_file}\n")
            f.write(f"LevelDataEffectivePressure.variable = effectivePressure\n")

        # Per-iteration pout so logs are not overwritten between restarts
        f.write(f"main.poutBaseName = pout.bisicles.iter{coupling_iter:04d}\n")

    return tmp


def run_bisicles(args, coupling_iter, restart_chk=None):
    """Run BISICLES for one coupling window; return path of latest plot file."""
    tmp = write_bisicles_input(args, coupling_iter, restart_chk)
    cmd = f"mpirun -n {args.nprocs} {args.bisicles_exe} {tmp}"
    run_cmd(cmd, f"BISICLES iteration {coupling_iter}")

    latest = find_latest_file(args.bisicles_plot_prefix)
    if latest is None:
        print(f"[COUPLED] ERROR: no BISICLES plot file found "
              f"(prefix='{args.bisicles_plot_prefix}')")
        sys.exit(1)

    symlink_force(latest, args.b2s_file)
    return latest


# ---------------------------------------------------------------------------
# SUHMO
# ---------------------------------------------------------------------------

def write_suhmo_input(args, coupling_iter):
    """
    Write a per-iteration SUHMO input file by copying the base file and
    appending coupling overrides.  Returns the temporary file path.
    """
    tmp = f"input.hydro.iter{coupling_iter:04d}"
    shutil.copy2(args.suhmo_input, tmp)

    with open(tmp, "a") as f:
        f.write(f"\n# --- coupling iteration {coupling_iter} overrides ---\n")
        f.write(f"suhmo.coupled_to_bisicles = true\n")
        f.write(f"suhmo.bisicles_input_file = {args.b2s_file}\n")
        f.write(f"suhmo.output_N_file       = {args.s2b_file}\n")
        # Per-iteration output prefixes so files are not overwritten
        f.write(f"main.poutBaseName        = pout.suhmo.iter{coupling_iter:04d}\n")
        f.write(f"AmrHydro.plot_prefix     = plot.suhmo.iter{coupling_iter:04d}.\n")
        f.write(f"AmrHydro.check_prefix    = chk.suhmo.iter{coupling_iter:04d}.\n")

    return tmp


def run_suhmo(args, coupling_iter):
    """Run SUHMO to near-steady state for the current BISICLES geometry."""
    tmp = write_suhmo_input(args, coupling_iter)
    cmd = f"mpirun -n {args.nprocs} {args.suhmo_exe} {tmp}"
    run_cmd(cmd, f"SUHMO steady-state solve, iteration {coupling_iter}")

    if not os.path.exists(args.s2b_file):
        print(f"[COUPLED] ERROR: SUHMO did not produce '{args.s2b_file}'")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Orchestrate iteratively coupled BISICLES–SUHMO MISMIP+ run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--nprocs", type=int, default=4,
                        help="MPI process count for both models")
    parser.add_argument("--coupling-steps", type=int, default=5,
                        help="BISICLES timesteps per coupling window")
    parser.add_argument("--num-coupling-iters", type=int, default=10,
                        help="Total number of coupling iterations to run")

    # Executables
    parser.add_argument(
        "--bisicles-exe", type=str,
        default=(
            "/home/tm17544/BISICLES/bisicles-uob-SUHMO-sep-exec/code/exec2D/"
            "driver2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex"
        ),
        help="Path to the BISICLES driver executable",
    )
    parser.add_argument(
        "--suhmo-exe", type=str,
        default=(
            "/home/tm17544/BISICLES/SUHMO-BISICLES-sep-exec/exec/BisiclesCoupled/"
            "Suhmo2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex"
        ),
        help="Path to the SUHMO executable (BisiclesCoupled build)",
    )

    # Input files
    parser.add_argument("--bisicles-input", type=str,
                        default="inputs.bisicles",
                        help="BISICLES ParmParse input file")
    parser.add_argument("--suhmo-input", type=str,
                        default="input.hydro",
                        help="SUHMO ParmParse input file")

    # Optional spin-up restart checkpoint
    parser.add_argument("--bisicles-restart", type=str, default=None,
                        help="BISICLES checkpoint file from the spin-up phase "
                             "(if omitted, BISICLES starts from scratch)")

    # Plot prefix used to locate the most recent BISICLES output
    parser.add_argument("--bisicles-plot-prefix", type=str,
                        default="plot.bisicles.",
                        help="Must match amr.plot_prefix in inputs.bisicles")

    # Coupling exchange filenames
    parser.add_argument("--b2s-file", type=str,
                        default="bisicles_for_suhmo.2d.hdf5",
                        help="BISICLES -> SUHMO coupling file (symlink)")
    parser.add_argument("--s2b-file", type=str,
                        default="suhmo_for_bisicles.2d.hdf5",
                        help="SUHMO -> BISICLES effective-pressure file")

    args = parser.parse_args()

    bisicles_chk = args.bisicles_restart   # None on first run; set by spin-up

    for iteration in range(args.num_coupling_iters):
        print(f"\n{'#'*60}")
        print(f"# COUPLING ITERATION {iteration} / {args.num_coupling_iters - 1}")
        print(f"{'#'*60}")

        # ── BISICLES ─────────────────────────────────────────────────────────
        run_bisicles(args, iteration, restart_chk=bisicles_chk)

        # Locate checkpoint written at end of this window for the next restart
        bisicles_chk = find_latest_file("chk.bisicles.")
        if bisicles_chk is None and iteration < args.num_coupling_iters - 1:
            print("[COUPLED] WARNING: no BISICLES checkpoint found for restart")

        # ── SUHMO ────────────────────────────────────────────────────────────
        run_suhmo(args, iteration)

    print(f"\n{'#'*60}")
    print(f"# COUPLED SIMULATION COMPLETE")
    print(f"# {args.num_coupling_iters} iterations, "
          f"{args.num_coupling_iters * args.coupling_steps} BISICLES steps total")
    print(f"{'#'*60}")


if __name__ == "__main__":
    main()
