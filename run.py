import numpy as np
import concurrent.futures
import sys
import os
import shutil
import subprocess
import treeswift
import hashlib

def filecontents_to_seed(file_path):
    # generate a seed from a hash of the GEMF input paramters file

    # Read the contents of the file
    with open(file_path, 'rb') as f:
        file_contents = f.read()
    # Hash the contents using SHA-256
    hash_object = hashlib.sha256(file_contents)
    # Convert the hash to an integer
    hash_integer = int(hash_object.hexdigest(), base=16)
    # Use a part of the hash to generate a seed
    seed_value = hash_integer % (2**32)  # Ensure the seed fits into a 32-bit integer
    return seed_value

def setup_GEMF_param(simulation_dir, params_file_name, new_seed):
    # copy the GEMF input paramters file into the simulation directory,
    # with a new seed

    new_params_file_name = os.path.join(simulation_dir, params_file_name)
    with open(params_file_name, 'r') as file:
        lines = file.readlines()
    for i, line in enumerate(lines):
        if line.strip() == "[RANDOM_SEED]":
            lines[i + 1] = str(new_seed) + "\n"
            break
    with open(new_params_file_name, 'w') as file:
        file.writelines(lines)


def get_primary_case(simulation_dir):
    # get the primary case
    file_path = os.path.join(simulation_dir, 'transmission_network.txt')
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            primary_case = parts[1]
            break
    return primary_case

def generate_samples(simulation_dir, params_file_name, rng):
    # Simulate sampling of the virus from individuals that were infected
    # in the simulated epidemic.
    # Also return a set "active", which contains each infection_id of cases
    # that had not recovered by the end of the simulation

    # get the id of the primary case from transmission network
    primary_case = get_primary_case(simulation_dir)

    # setup dict for disease info of ascertained cases and primary case
    # and  time of first hospitalisation
    infections = {}
    # add the primary case
    infections[primary_case] = {'start_time': 0, 'end_time': 0}
    first_hospitalisation_time = 0
    # read disease info from output.txt
    file_path = os.path.join(simulation_dir, 'output.txt')
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            infection_id = parts[1]
            event_time = float(parts[0])
            event_type = parts[3]
            if infection_id == primary_case and event_type == '5': # primary case start times
                infections[infection_id] = {'start_time': event_time, 'end_time': 0}
            elif event_type == '3':  # ascertainement event
                infections[infection_id] = {'start_time': event_time, 'end_time': 0}
            elif event_type == '7':  # Recovery event
                if infection_id in infections:
                    infections[infection_id]['end_time'] = event_time
            elif event_type == '6':  # Hospitalisation event
                # record first hospitalisation time
                if first_hospitalisation_time == 0:
                    first_hospitalisation_time = event_time
                # ensure all hospitalized cases are sampled by making the sure their start time doesn't precede the first hsopitlazation time'
                elif infections[infection_id]['start_time'] < first_hospitalisation_time:
                    infections[infection_id]['start_time'] = first_hospitalisation_time

    # get max time from the params file
    max_time = 0
    file_path = os.path.join(simulation_dir, params_file_name)
    with open(file_path, 'r') as file:
        lines = file.readlines()
    for i, line in enumerate(lines):
        if line.strip() == "[MAX_TIME]":
            max_time = float(lines[i + 1].strip())
            break

    # generate sample times and write them
    # store set of infections that were not recovered at end of simulation
    active = set()
    file_path = os.path.join(simulation_dir, 'sample_times.txt')
    with open(file_path, 'w') as file:
        for infection_id in infections:
            start_time = infections[infection_id]['start_time']
            # set end time to sim end time if not recovered, and add to active
            if infections[infection_id]['end_time'] == 0:
                infections[infection_id]['end_time'] = max_time
                active.add(infection_id)
            end_time = infections[infection_id]['end_time']
            sample_time = rng.uniform(start_time, end_time)
            # write sample time if it is the primary case or after first hospitalisation
            if sample_time > first_hospitalisation_time or infection_id == primary_case:
                file.write(f"{infection_id}\t{sample_time}\n")
    return active, rng

def get_clade_info(parent, rng):
    lineages = 0
    one_mutation_subclade_roots = []
    two_mutation_subclade_roots = []
    for child in parent.child_nodes():
        if child.is_leaf():
            lineages += 1
        else:
            mutations = rng.poisson(29903*0.00092*child.get_edge_length())
            if mutations > 0:
                lineages += 1
                one_mutation_subclade_roots.append(child)
                if mutations > 1:
                    two_mutation_subclade_roots.append(child)
            else:
                derived_lineages, derived_one_mutation_subclade_roots, derived_two_mutation_subclade_roots, rng = get_clade_info(child, rng)
                lineages += derived_lineages
                one_mutation_subclade_roots.extend(derived_one_mutation_subclade_roots)
                two_mutation_subclade_roots.extend(derived_two_mutation_subclade_roots)
    return lineages, one_mutation_subclade_roots, two_mutation_subclade_roots, rng

def simulate(simulation_no, seed, params_file_name):

    rng = np.random.default_rng(seed)

    simulation_dir = params_file_name.split('.')[0] + f"/{simulation_no}"
    os.makedirs(simulation_dir, exist_ok=True)

    # run GEMF to generate transmission network
    GEMF_rng_seed = rng.integers(low=0, high=2_147_483_648)
    setup_GEMF_param(simulation_dir, params_file_name, GEMF_rng_seed)
    command = f"GEMF {params_file_name}"
    logfile_path = os.path.join(simulation_dir, f"GEMF_log.txt")
    with open(logfile_path, 'w') as logfile:
        subprocess.run(command, cwd=simulation_dir, shell=True, stdout=logfile, stderr=subprocess.STDOUT)

    # generate samples and get ids of infections that are active at sim end
    active, rng = generate_samples(simulation_dir, params_file_name, rng)

    # run CoaTran to generate tree
    coatran_rng_seed = rng.integers(low=0, high=2_147_483_648)
    command = "coatran_constant transmission_network.txt sample_times.txt 1"
    env = os.environ.copy()
    env["COATRAN_RNG_SEED"] = str(coatran_rng_seed)
    with open(os.path.join(simulation_dir, "tree.time.nwk"), 'w') as output_file:
        subprocess.run(command, cwd=simulation_dir, shell=True, stdout=output_file, env=env)

    # get tree without primary case sample
    tree = treeswift.read_tree_newick(os.path.join(simulation_dir, "tree.time.nwk"))
    primary_case= get_primary_case(simulation_dir)
    excluded_labels = set()
    for label in tree.labels(leaves=True, internal=False):
        if label.split('|')[1] == primary_case:
            excluded_labels.add(label)
    tree = tree.extract_tree_without(excluded_labels)
    tree.suppress_unifurcations()

    # get stable_coalescence using labels of infections that are active at sim end
    stable_coalescence_labels = set()
    for label in tree.labels(leaves=True, internal=False):
        if label.split('|')[1] in active:
            active.remove(label.split('|')[1])
            stable_coalescence_labels.add(label)
    # start with MRCA of those active infections as stable_coalescence
    stable_coalescence = tree.mrca(stable_coalescence_labels)
    # shift to parent if there is no time difference (and there is a parent)
    while stable_coalescence.get_edge_length() == 0 and not stable_coalescence.is_root():
        stable_coalescence = stable_coalescence.get_parent()

    # simulate mutations along branches from stable coalescent and classify clade topology
    lineages, one_mutation_subclade_roots, two_mutation_subclade_roots, rng = get_clade_info(stable_coalescence, rng)
    polytomy= False
    AB = False
    CC = False
    # check for basal polytomy
    if lineages >= 100:
        polytomy = True
        size = sum([1 for _ in stable_coalescence.traverse_leaves()])
        # check for two_mutation_subclade polytomy
        for subclade_root in two_mutation_subclade_roots:
            subclade_lineages, _, _, rng = get_clade_info(subclade_root, rng)
            if subclade_lineages >= 100:
                # check for relative sizes
                subclade_size = sum([1 for _ in subclade_root.traverse_leaves()])
                if 0.3 <= (subclade_size/size) <= 0.7:
                    AB = True
    # check if two one_mutation clades
    elif lineages == 2 and len(one_mutation_subclade_roots) == 2:
        # check for one_mutation_subclade polytomies
        first_subclade_lineages, _, _, rng = get_clade_info(one_mutation_subclade_roots[0], rng)
        if first_subclade_lineages >= 100:
            second_subclade_lineages, _, _, rng = get_clade_info(one_mutation_subclade_roots[1], rng)
            if second_subclade_lineages >= 100:
                # check relative sizes
                size = sum([1 for _ in stable_coalescence.traverse_leaves()])
                subclade_size = sum([1 for _ in one_mutation_subclade_roots[0].traverse_leaves()])
                if 0.3 <= (subclade_size/size) <= 0.7:
                    CC = True

    # write topologies
    file_path = os.path.join(simulation_dir, 'poly_AB_CC.txt')
    with open(file_path, 'w') as file:
        file.write(str(polytomy) + ' ' + str(AB) + ' ' + str(CC))
    # read failures
    file_path = os.path.join(simulation_dir, 'failures.txt')
    with open(file_path, 'r') as file:
        failures = int(file.readline().split()[0])

    return failures, polytomy, AB, CC

def main(parameters_file_name, num_processors, num_sims):
    # seed based on file contents
    seed_value = filecontents_to_seed(parameters_file_name)
    seeds = np.random.SeedSequence(seed_value).spawn(num_sims)
    file_names = [parameters_file_name] * num_sims
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_processors) as executor:
        results = list(executor.map(simulate, range(num_sims), seeds, file_names))

    # Initialize counters
    failure_count, success_count, polytomy_count, AB_count, CC_count = 0, 0, 0, 0, 0

    # Process results
    for failures, polytomy, AB, CC in results:
        failure_count += failures
        success_count += 1
        if polytomy:
            polytomy_count += 1
        if AB:
            AB_count += 1
        if CC:
            CC_count += 1

    # Bayes Factor calculations
    uncon_posteriors = np.array([1.68, 80.85, 10.32, 0.92])
    recCA_posteriors = np.array([77.28, 8.18, 10.49, 3.71])

    uncon_BF = (0.25 * uncon_posteriors.sum() * (polytomy_count/success_count)**2) / \
               (0.5 * uncon_posteriors[:2].sum() * (AB_count/success_count) + 0.5 * uncon_posteriors[2:].sum() * (CC_count/success_count))
    recCA_BF = (0.25 * recCA_posteriors.sum() * (polytomy_count/success_count)**2) / \
               (0.5 * recCA_posteriors[:2].sum() * (AB_count/success_count) + 0.5 * recCA_posteriors[2:].sum() * (CC_count/success_count))

    # write results summary
    with open(file_name.split('.')[0] + "_results.txt", 'w') as file:
        file.write('Parameter file: ' + file_name + '\n')
        file.write(f'Seed: {seed_value}\n')
        file.write(f'Successful epidemics: {success_count}\n')
        file.write(f'Failure rate: {failure_count/(failure_count+success_count)}\n')
        file.write(f'Basal polytomy rate: {polytomy_count/success_count}\n')
        file.write(f'AB topology rate: {AB_count/success_count}\n')
        file.write(f'CC topology rate: {CC_count/success_count}\n')
        file.write(f'unconstrained Bayes factor: {uncon_BF}\n')
        file.write(f'recCA Bayes factor: {recCA_BF}\n')

    # display results summary
    print('Parameter file: ' + file_name)
    print(f'Seed: {seed_value}')
    print(f'Successful epidemics: {success_count}')
    print(f'Failure rate: {failure_count/(failure_count+success_count)}')
    print(f'Basal polytomy rate: {polytomy_count/success_count}')
    print(f'AB topology rate: {AB_count/success_count}')
    print(f'CC topology rate: {CC_count/success_count}')
    print(f'unconstrained Bayes factor: {uncon_BF}')
    print(f'recCA Bayes factor: {recCA_BF}')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 script.py <params_file> <num_processors> <num_sims>")
        sys.exit(1)
    file_name = sys.argv[1]
    num_processors = int(sys.argv[2])
    num_sims = int(sys.argv[3])
    main(file_name, num_processors, num_sims)
