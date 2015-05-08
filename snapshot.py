"""
A module containing functions for storing a population
to a file, and loading it again.
"""
from __future__ import print_function
import datetime
import csv, tarfile
from collections import deque
try:
    import cPickle as pickle
except ImportError:
    import pickle
from mutation import Mutation
from subpopulation import Subpopulation
from population import Population
from utils import delete_local_file


def save_population_to_file(popn, fpath):
    """Save population snapshot to file.

    Take a snapshot of the current simulation
    population and save it to timestamped CSV files,
    then compress the files into a .tar archive.
    """
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M")
    param_fname = "params_{}.pkl".format(timestamp)
    mut_fname = "mutations_{}.csv".format(timestamp)
    clone_fname = "clones_{}.csv".format(timestamp)
    snapshot_fnames = [param_fname, mut_fname, clone_fname]

    root_clone = popn.subpop
    all_muts = popn.all_mutations

    # create snapshot files
    save_parameters_to_file(popn, param_fname)
    save_muts_to_file(all_muts, mut_fname)
    save_clones_to_file(root_clone, clone_fname)

    # store files in tar archive
    archive_name = "{}/population_{}.tar.gz".format(fpath, timestamp)
    popn_archive = tarfile.open(archive_name, "w:gz")
    for fname in snapshot_fnames:
        popn_archive.add(fname)
    popn_archive.close()

    # now delete the individual files
    for fname in snapshot_fnames:
        delete_local_file(fname)


def save_parameters_to_file(popn, param_fname):
    """
    Save population parameters to file.

    Use Python's Pickle module to simplify reloading from file.
    Assumes population is being stored *before* treatment introduced,
    so we don't record, for instance, popn.selective_pressure.
    """
    # get specific population parameters we want to record,
    # mainly so we don't have to waste time recalculating
    # when we load the population
    popn_params = {'tumoursize': popn.tumoursize,
                   'clonecount': popn.clonecount,
                   'avg_pro_rate': popn.avg_pro_rate,
                   'avg_mut_rate': popn.avg_mut_rate,}

    # pickle the global parameter set along with these
    # population parameters
    pickled_params = pickle.dumps((popn.opt, popn_params))

    # write to file
    with open(param_fname, 'w') as param_file:
        param_file.write(pickled_params)


def save_muts_to_file(all_muts, mut_fname):
    """
    Save the entire mutation dictionary to a CSV file.

    NOTE: there is an implicit, currently unenforced contract
    between this function and Mutation.init_from_file(). That
    method relies on this function storing valid mutation
    attributes, with attribute names as column headers.

    The best way to enforce this contract, down the track,
    might be to set up a test suite. Trying to enforce it
    in code would likely be ugly, time-consuming and inflexible.
    """
    with open(mut_fname, 'w') as mut_file:
        writer = csv.writer(mut_file)
        header = ['mut_id', 'mut_type',
                  'prolif_rate_effect', 'mut_rate_effect',
                  'resist_strength', 'original_clone_id', 's_time']
        writer.writerow(header)

        for mut_type in all_muts:
            # don't write dead mutations to file
            if mut_type == 'dead':
                continue
            for mut in all_muts[mut_type]:
                mut_data = [mut.mut_id, mut.mut_type,
                            mut.prolif_rate_effect, mut.mut_rate_effect,
                            mut.resist_strength, mut.original_clone.clone_id,
                            mut.s_time]
                # convert all data to strings, to ensure, in
                # particular, that None values are written in such a way
                # that they can later be parsed by literal_eval
                writer.writerow([repr(item) for item in mut_data])


def save_clones_to_file(root_clone, clone_fname):
    """
    Save the entire clone tree to a CSV file.

    This function stores the clones to file in
    breadth-first search order, to facilitate reconstruction
    of the clone tree when loading the population from file.

    NOTE: there is an implicit, currently unenforced contract
    between this function and Subpopulation.init_from_file(). That
    method relies on this function storing valid clone
    attributes, with attribute names as column headers.
    See note in save_muts_to_file above.
    """
    with open(clone_fname, 'w') as clone_file:
        writer = csv.writer(clone_file)
        header = ['clone_id', 'parent_id', 'num_children',
                  'prolif_rate', 'mut_rate', 'size', 'precrash_size',
                  'depth', 's_time', 'd_time', 'branch_length',
                  'is_resistant', 'resist_strength', 'mutations',]
        writer.writerow(header)

        # we use a LIFO queue to store the clones to file in BFS order
        queue = deque([root_clone])

        # a dictionary to map clones to their parent IDs
        parents = {}

        parents[root_clone.clone_id] = None

        while queue:
            curr_clone = queue.popleft()
            # convert dictionary of Mutations
            # to dictionary of mutation IDs
            mut_id_dict = {}
            for mut_type in curr_clone.mutations:
                mut_id_dict[mut_type] = [mut.mut_id
                                         for mut
                                         in curr_clone.mutations[mut_type]]
            clone_data = [curr_clone.clone_id,
                          # parent_id - BFS write order ensures this is set
                          # by the time we write each clone to file
                          parents[curr_clone.clone_id],
                          # num_children
                          len(curr_clone.nodes),
                          curr_clone.prolif_rate, curr_clone.mut_rate,
                          curr_clone.size, curr_clone.precrash_size,
                          curr_clone.depth,
                          curr_clone.s_time, curr_clone.d_time,
                          curr_clone.branch_length,
                          curr_clone.is_resistant, curr_clone.resist_strength,
                          # mutations
                          mut_id_dict]
            # convert all data to string rep, to ensure, in
            # particular, that None vals are written in such a way
            # that it can later be parsed by literal_eval
            writer.writerow([repr(item) for item in clone_data])
            # now add the children of this clone to the queue,
            # if there are any
            for child in curr_clone.nodes:
                parents[child.clone_id] = curr_clone.clone_id
                queue.append(child)


def load_population_from_file(archive_path):
    """Load a population from a snapshot.

    Load a population from a snapshot,
    which should be stored in a .tar archive.
    """
    popn_archive = tarfile.open(archive_path)

    # get filenames for extraction
    # NB this may break if filenames in save_population_to_file() are changed
    clone_fname = mut_fname = param_fname = None

    for filename in popn_archive.getnames():
        if filename.startswith('clones_'):
            clone_fname = filename
        elif filename.startswith('mutations_'):
            mut_fname = filename
        elif filename.startswith('params_'):
            param_fname = filename
        else:
            raise Exception("population archive contains invalid files")

    try:
        # extract snapshot files to current working directory
        popn_archive.extract(clone_fname)
        popn_archive.extract(mut_fname)
        popn_archive.extract(param_fname)
    except:
        # one of our snapshot files hasn't been found
        raise
    finally:
        popn_archive.close()

    # load parameters, mutations and clones
    opt, popn_params = load_parameters_from_file(param_fname)
    all_muts, mutation_map = load_muts_from_file(opt, mut_fname)
    root_clone = load_clones_from_file(opt, mutation_map, clone_fname)

    # construct population from parameter set,
    # clone tree and mutation dictionary
    new_popn = Population.init_from_file(opt, popn_params, root_clone, all_muts)

    # now delete the individual snapshot files,
    # as we will always load popn from an archive
    for fname in [clone_fname, mut_fname, param_fname]:
        delete_local_file(fname)

    return new_popn


def load_parameters_from_file(param_fname):
    """
    Load global parameters from population snapshot.

    Assumes parameters have been pickled. There
    is a contract between this function and
    save_parameters_to_file() above.
    """
    with open(param_fname) as param_file:
        data = param_file.read()

    opt, popn_params = pickle.loads(data)

    return opt, popn_params


def load_muts_from_file(opt, mut_fname):
    """
    Load all mutations from a population snapshot.

    Load the entire population of mutations from
    a population snapshot file.

    Inputs
    ------
    opt: global parameter set, used in initialising each mutation
    mut_fname: path to CSV file containing mutation snapshot (which
               must have been written using save_muts_to_file above)

    Returns
    -------
    2-tuple of dictionaries: (all_muts, mutation_map)

    all_muts: a dictionary of mutations of the form
              {mut_type: list of muts of that type}, passed
              back to simulation as master list of mutations.

    mutation_map: a dictionary of mutations of the form
                  {mut_id: Mutation}, used to reconstruct
                  the relationships between clones and mutations when
                  loading population snapshot from file.
    """
    all_muts = {'b':[], 'n':[], 'd':[], 'r':[]}
    mutation_map = {}

    with open(mut_fname) as mut_file:
        mut_reader = csv.DictReader(mut_file)
        for row in mut_reader:
            # initialise new Mutation object
            new_mut = Mutation.init_from_file(opt, all_muts, row)

            # add it to ID map
            mutation_map[new_mut.mut_id] = new_mut

            # add it to main mutation dictionary
            try:
                all_muts[new_mut.mut_type].append(new_mut)
            except KeyError:
                raise KeyError("invalid mutation type (from file): {}".format(new_mut.mut_type))

    return all_muts, mutation_map


def load_clones_from_file(opt, mutation_map, clone_fname):
    """
    Load all clones from a population snapshot.

    Load the entire population of clones from
    a population snapshot file. This function
    loads the clones in breadth-first search
    order in order to efficiently reconstruct
    the clone tree; consequently, it assumes
    the population was *stored* in BFS order.

    Inputs
    ------
    opt:          global parameter set, used
                  in clone initialisation
    mutation_map: dictionary mapping mutation IDs
                  to Mutation objects
    clone_fname:  path to CSV file containing clone
                  population snapshot

    Returns
    -------
    The root of the clone tree (i.e. the initial subpopulation)
    """
    # as this is a BFS algorithm, we use a LIFO queue
    parent_queue = deque([])

    with open(clone_fname) as clone_file:
        reader = csv.DictReader(clone_file)

        # root clone should always be first clone in file
        root_row = next(reader)
        root_clone = Subpopulation.init_from_file(opt, mutation_map, root_row)
        parent_queue.append(root_clone)

        for row in reader:
            # initialise clone (including associating it with its mutations)
            new_clone = Subpopulation.init_from_file(opt, mutation_map, row)

            # this clone's parent should be first parent in queue
            curr_parent = parent_queue[0]
            if curr_parent.clone_id != new_clone.parent_id:
                # something has gone wrong
                raise Exception("clones stored in incorrect order; abort loading population.")

            # connect new clone to its parent
            curr_parent.nodes.append(new_clone)

            if len(curr_parent.nodes) == curr_parent.num_children:
                # this clone's parent is now connected to all
                # its children, so remove it from the queue
                del curr_parent.num_children
                parent_queue.popleft()

            if new_clone.num_children > 0:
                parent_queue.append(new_clone)

    # entire tree has been restored
    return root_clone
