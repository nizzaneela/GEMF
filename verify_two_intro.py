#!/usr/bin/env python3

"""
Takes an integer 0..109999 from the command line,
runs one_intro twice and then two_intro in a unique directory for verification.
"""

import sys
import os
import shutil
import numpy as np

# Import from your main "replicate" script:
from replicate import one_intro, two_intro, NUM_SIMS

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <parameters_file_name> <index (0..109999)>")
        sys.exit(1)

    parameters_file_name = sys.argv[1]
    i = int(sys.argv[2])

    # Validate the index
    if i < 0 or i > 109999:
        print("Error: index must be between 0 and 109999 inclusive.")
        sys.exit(1)

    # Reproduce the main script’s seed logic
    seed = 42
    # spawn(2) yields two child SeedSequence objects
    one_intro_seedseq, two_intro_seedseq = np.random.SeedSequence(seed).spawn(2)
    # Then spawn 110,000 seeds for one_intro:
    one_intro_seedseq_list = one_intro_seedseq.spawn(NUM_SIMS)
    # And 110,000 seeds for two_intro:
    two_intro_seedseq_list = two_intro_seedseq.spawn(NUM_SIMS)

    # Create the RNG for the i-th two_intro run
    rng = np.random.default_rng(two_intro_seedseq_list[i])

    # Create or switch to a verification directory
    verification_dir  = f"{parameters_file_name.split('_')[0]}_verifications"
    if os.path.basename(os.getcwd()) != verification_dir:
        os.makedirs(verification_dir, exist_ok=True)
        # Optionally copy the parameters file for reference
        shutil.copy2(parameters_file_name, os.path.join(verification_dir, parameters_file_name))
        os.chdir(verification_dir)

    # 1) Pick two single-intro runs at random and run them first
    one_intro_runs = []
    for _ in range(2):
        # Randomly pick an index from 0..(NUM_SIMS-1)
        chosen_idx = rng.integers(NUM_SIMS)
        one_intro_runs.append(chosen_idx)

        # Output directory for the single-intro run
        single_intro_outdir = os.path.join("one_intro", str(chosen_idx))

        # Run the single-intro simulation
        one_intro(
            single_intro_outdir,
            parameters_file_name,
            one_intro_seedseq_list[chosen_idx]
        )

    # 2) Now combine those two single-intro simulations using two_intro
    two_intro_outdir = os.path.join("two_intro", str(i))
    os.makedirs(two_intro_outdir, exist_ok=True)

    # Call two_intro with the same rng
    two_intro(
        two_intro_outdir,
        parameters_file_name,
        rng,
        one_intro_runs[0],
        one_intro_runs[1]
    )

    print(f"Verification run {i} completed.")
    print(f"Results are in: {os.path.join(verification_dir, 'two_intro', str(i))}")

if __name__ == "__main__":
    main()

