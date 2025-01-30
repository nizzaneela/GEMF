#!/usr/bin/env python3

"""
Takes an integer 0..109999 from the command line,
and runs one_intro to verify that run in a unique directory.
"""

import sys
import os
import shutil
import numpy as np
from replicate import one_intro, NUM_SIMS

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <parameters_file_name> <index (0..109999)>")
        sys.exit(1)
    parameters_file_name = sys.argv[1]
    i = int(sys.argv[2])
    if i < 0 or i > 109999:
        print("Error: index must be between 0 and 109999 inclusive.")
        sys.exit(1)

    # Reproduce replicate's SeedSequence list
    seed = 42
    one_intro_seedseq, two_intro_seedseq = np.random.SeedSequence(seed).spawn(2)
    seedseq_list = one_intro_seedseq.spawn(NUM_SIMS)

    # Get the i-th seed from the list
    selected_seedseq = seedseq_list[i]

    # ensure verification directory is cwd
    verification_dir  = f"{parameters_file_name.split('.')[0]}_verifications"
    if os.path.basename(os.getcwd()) != verification_dir:
        os.makedirs(verification_dir, exist_ok=True)
        shutil.copy2(parameters_file_name, os.path.join(verification_dir, parameters_file_name))
        os.chdir(verification_dir)

    # set output directory
    output_dir = os.path.join("one_intro", str(i))

    # Run the single-introduction simulation using the i-th SeedSequence
    one_intro(output_dir, parameters_file_name, selected_seedseq)

    print(f"Verification run {i} completed. Results are in: {os.path.join(verification_dir, "one_intro", str(i))}")


if __name__ == "__main__":
    main()

