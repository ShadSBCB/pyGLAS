# -*- coding: utf-8 -*-
import os
import sys
import textwrap
import itertools
from argparse import ArgumentParser as ArgP

chicos_residues = {
    1: [1.57, 2.44],
    2: [2.42, 3.46],
    3: [2.47, 1.53],
    4: [2.47, 1.50],
    5: [1.50, 2.50],
    6: [2.50, 7.46],
    7: [1.50, 7.46],
    8: [3.34, 4.57],
    9: [4.53, 3.34],
    10: [4.53, 3.38],
    11: [3.38, 4.50],
    12: [3.44, 5.54],
    13: [5.54, 6.41],
    14: [3.47, 5.57],
    15: [5.57, 3.51],
    16: [3.51, 5.60],
    17: [2.43, 7.53],
    18: [3.36, 6.48],
    19: [3.40, 6.44],
    20: [1.46, 7.47],
    21: [1.49, 7.50],
    22: [6.51, 7.39],
    23: [6.51, 7.38],
    24: [6.47, 7.45]
}

# --------------------------------------------------------------------------- #
parser = ArgP(
        description="pyGLAS takes a GPCR structure and determines it's \
        GPCR-Likeness Assessment Score (GLAS). Given a trajectory, it will \
        calculate the GLAS timeseries.")

parser.add_argument(
        '-r', '--receptor', dest='recept', 
        default=None, help='Receptor file in PDB or GRO format.')

parser.add_argument(
        '-x', '--trajectory', dest='trj', default=None,
        help='Single trajectory file (compatible with MDAnalysis).')

# =============================================================================
# parser.add_argument(
#         '--chicos', dest='chicos_file', default=None,
#         help='List of residue numbers for CHICOs (see \
#         examples/Single_Topology for an example).')
# =============================================================================

parser.add_argument(
        '--cons', dest='cons_file', default=None,
        help='File of residues in position 50 in each helix.')

parser.add_argument(
        '-l', '--log', dest='log_file', default='pyGLAS.log',
        help='Name for output log file. Defaults to "pyGLAS.log".')

parser.add_argument(
        '--format', dest='run_type', default="short",
        help='Specifies run type. Defaults to short. \
        Please see documentation for more.')

parser.add_argument(
        '--cutoff', dest='cutoff', default=4.0, 
        help='Defines cutoff for distance calculations. Default is 4A.')

parser.add_argument(
        '--debug', dest='debug', default=False,
        help='Triggers debugging checks.')

args = parser.parse_args()

# --------------------------------------------------------------------------- #

#####################      Pre-processing functions      ######################

def load_from_position50_file(filepath):
    """
    Loads data into a dictionary for later use.
    
    :param filepath: path to position50.txt file.
    """
    with open(filepath) as f:
        position50 = {numpy.around(numpy.float(line.split()[0]), 2): 
            numpy.int(line.split()[1]) for line in f}
    
    return position50


def gen_bw_numbers(position50):
    """
    Generates a dictionary of resids that match the residues from position 
    X.35 to X.65, where X is the helix number
    """
    
    bw_numbers = position50.copy()
    for resid in sorted(bw_numbers):
        for i in numpy.linspace(start=-0.16, stop=0.16, num=33):
            bw_numbers[numpy.float('{:.2f}'.format(float(resid) + i))] = int(bw_numbers[resid] + i*100)

    return bw_numbers


def convert_pairs_to_fulllength(chicos_residues, bw_numbers):
    """
    Returns a dictionary of chicos pairs with resids.
    
    :param chicos_residues: hardcoded chicos pairs (with BW numbers).
    :param bw_numbers: dictionary of BW numbers (key: BW number, value: resid)
    
    :return chicos: same as chicos_residues, with resids.
    """
    
    chicos = chicos_residues.copy()
    for k,v in chicos_residues.items():
        chicos[k] = [bw_numbers[v[0]], bw_numbers[v[1]]]
    
    return chicos


######################    Functions to output results    ######################

def short_static_GLAS(system, GLAS):
    """
    Simple. 
    Writes GLAS to a file that starts with input filename and ends in .short.
    
    :param system: MDAnalysis Universe with a structure.
    :param GLAS: GLAS.
    """
    # Use in not long form.
    filename = "{}_GLAS.short".format(
            os.path.splitext(os.path.basename(system.filename))[0])

    backup_file(filename)

    with open(filename, 'w') as f:
        f.write("##################################################\n")
        f.write("#           Thank you for using pyGLAS           #\n")
        f.write("#         The GLAS of your receptor is {}        #\n".format(
                                                                        GLAS))
        f.write("#  Please ensure you cite all relevant sources.  #\n")
        f.write("##################################################\n")


def long_static_GLAS(system, chicos_dict, array, GLAS):
    """
    Writes long form GLAS result. 
    Writes output in the form of a tab separated value file with 
    pair number and minimum distance and ends in .long.
    
    :param system: MDAnalysis Universe with a structure.
    :param dct:
    :param array:
    :param GLAS: GLAS.
    """
    filename = "{}_GLAS.long".format(
            os.path.splitext(os.path.basename(system.filename))[0])

    with open(filename, 'w') as f:
        f.write("##################################################\n")
        f.write("#           Thank you for using pyGLAS           #\n")
        f.write("#         The GLAS of your receptor is {}        #\n".format(
                                                                        GLAS))
        f.write("#  Please ensure you cite all relevant sources.  #\n")
        f.write("##################################################\n")
        f.write("\n")
        f.write("# Since you have requested a long form output, \n")
        f.write("# below are the 24 pairs with the minimum distance \n")
        f.write("# between them.\n")
        f.write("\n")
        f.write("# Pair\tMinimum Distance\n")

        for key, dist in itertools.izip(chicos_dict, array):
            f.write("{0}\t{1:.3f}\n".format(key, dist))


def short_dynamic_GLAS(system, dGLAS):
    """
    Plots time dependent GLAS in a single plot.
    Also outputs time varied GLAS in file.
    
    :param system: MDAnalysis Universe with a structure and a trajectory.
    :param dGLAS: array of GLAS timeseries.
    """
    
    filename = "{}_GLAStimeseries.csv".format(
            os.path.splitext(
                    os.path.basename(system.filename))[0])
    
    with open(filename, "w") as f:
        f.write("Frame,GLAS\n")
        for entry in dGLAS.transpose():
            f.write("{},{}\n".format(entry[0], entry[1]))
    
    plotname = "{}_GLAStimeseries.svg".format(
            os.path.splitext(
                    os.path.basename(system.filename))[0])
    
    fig = pyplot.figure(figsize=(4,6))
    pyplot.plot(dGLAS[0], dGLAS[1], 'k--')
    pyplot.title('GLAS timeseries')
    pyplot.xlabel('GLAS')
    pyplot.xlim(numpy.min(dGLAS[0]), numpy.max(dGLAS[0]))
    pyplot.ylabel('Frame Number')
    pyplot.ylim(0, 24)
    fig.savefig(plotname)


def long_dynamic_GLAS(system, dGLAS, distsDict):
    """
    Plots time dependent GLAS in a single plot, and the 24 individual pair 
    timeseries in a 4x6 grid of plots.
    
    :param system: MDAnalysis Universe with a structure and a trajectory.
    :param dGLAS: array of GLAS timeseries.
    :param distsDict: dictionary with format key (pair number) and 
    minimum distance (value).
    """
    
    short_dynamic_GLAS(system, dGLAS)
    
    for k,v in distsDict.items():
        with open("Pair_{}_timeseries.csv".format(k), "w") as f:
            for i in v:
                f.write("{}\n".format(i))
    
    plotsname = "{}_GLAS_pairstimeseries.svg".format(
            os.path.splitext(
                    os.path.basename(system.filename))[0])
    
    rows, cols = 4, 6
    fig, axes = pyplot.subplots(ncols=cols, nrows=rows, 
                                sharex=True, sharey=True,
                                figsize=(30, 20))
    
    for idx,pair in numpy.ndenumerate(
            numpy.reshape(sorted(distsDict.keys()), (rows,cols))):
        axes[idx[0]][idx[1]].plot(distsDict[pair])
        axes[idx[0]][idx[1]].set_title("Pair {}".format(pair))
    
    axes[-1][-1].set_xlim(numpy.min(dGLAS[0]), numpy.max(dGLAS[0]))
    
    fig.savefig(plotsname)


###################    Functions to calculate results    ######################


def solo_protocol(system, chicos_dict, cutoff = 4.0, debug = False):
    """
    Uses system without trajectory to return GLAS.
    
    :param system: MDAnalysis Universe with a structure.
    :param chicos_dict: dictionary with pair and resids.
    :param cutoff: cutoff for distance calculation 
    (default is 4.0; can be changed with args.cutoff)
    :param debug: Used for debugging 
    (default is False; can be changed with args.debug)
    """
    
    minDistArray = []

    for key in chicos_dict:

        # Get resids for MDAnalysis; excludes all hydrogen atoms.
        resid1 = "resid {} and not name H*".format(chicos_dict[key][0])
        resid2 = "resid {} and not name H*".format(chicos_dict[key][1])

        # Debug statement
        if debug:
            # Check if code is assigning strings correctly
            print resid1, resid2

        # Grab coordinates from resid variables.
        coords1 = system.select_atoms(resid1).positions
        coords2 = system.select_atoms(resid2).positions

        # Append minimum to minDistArray
        minDistArray.append(numpy.amin(cdists.distance_array(coords1, coords2)))

    GLAS = numpy.sum(numpy.array(minDistArray) < cutoff)

    return numpy.array(minDistArray), GLAS


def dynamic_GLAS(system, chicos_dict, cutoff = 4.0, debug = False):
    """
    Uses the system with trajectory to calculate the time resolved GLAS.
    
    :param system: MDAnalysis Universe with a structure and a trajectory.
    :param chicos_dict: dictionary with pair and resids.
    :param cutoff: cutoff for distance calculation 
    (default is 4.0; can be changed with args.cutoff)
    :param debug: Used for debugging 
    (default is False; can be changed with args.debug)
    """

    dGLAS = numpy.zeros((2, system.trajectory.n_frames))

    # To store minimum distances
    minDistArray = []

    distsDict = {key: [] for key in chicos_dict}

    for ts in system.trajectory:

        # Assign frame number to 'x'-axis
        dGLAS[0][ts.frame] = ts.frame

        for key in chicos_dict:
            # Assign row values to variables and convert to strings
            resid1 = "resid {} and not name H*".format(chicos_dict[key][0])
            resid2 = "resid {} and not name H*".format(chicos_dict[key][1])

            # Debug statement
            # Check if code is assigning strings correctly
            if debug:
                print resid1, resid2

            # Grab coordinates from resid variables.
            coords1 = system.select_atoms(resid1).positions
            coords2 = system.select_atoms(resid2).positions

            # Pass coords variables to dists.distance_array function
            ContactsMatrix = cdists.distance_array(coords1, coords2)

            distsDict[key].append(numpy.amin(ContactsMatrix))

            # Append minimum to minDistArray
            minDistArray.append(numpy.amin(ContactsMatrix))

        GLAS = numpy.sum(numpy.array(minDistArray) < cutoff)

        # Reset minDistArray
        minDistArray = []

        # Assign GLAS to 'y'-axis
        dGLAS[1][ts.frame] = GLAS

    return dGLAS, distsDict


######################      Miscellaneous functions      ######################

def backup_file(filename):
    """
    Creates a backup file if filename exists.
    
    :param filename: name for the file created.
    """
    if os.path.exists(filename):
        target_name = filename + ".bkup"
        failure = True
        if not os.path.exists(target_name):
            os.rename(filename, target_name)
            failure = False
        else:
            for i in range(20):
                alt_target_name = "{}.{}".format(target_name, str(i))
                if os.path.exists(alt_target_name):
                    continue
                else:
                    os.rename(filename, alt_target_name)
                    failure = False
                    break
        if failure:
            raise IOError("Too many backups. Clean up and try again")


def _log_and_print(log, message):
    log.write(textwrap.fill(message, 80))
    log.write("\n")
    print textwrap.fill(message, 80)


# Log file
if os.path.isfile('{}'.format(args.log_file)):
    backup_file(args.log_file)
    log_file = open("{}".format(args.log_file), 'w')
else:
    log_file = open("{}".format(args.log_file), 'w')

# import external libraries and packages
# MDAnalysis for Universe and distance calculation
try:
    import MDAnalysis
    _log_and_print(
            log=log_file,
            message='Found MDAnalysis version {}.\n'.format(
                    MDAnalysis.__version__))
except ImportError:
    _log_and_print(
            log=log_file,
            message='Cannot find MDAnalysis in current Python environment. Please ensure you have a recent (14 or above) version of MDAnalysis, then start again. It can usually be installed using conda or from github.')
    sys.exit()

from MDAnalysis import Universe
from MDAnalysis.lib import distances as cdists

# pandas for data handling and plotting.
try:
    import pandas
    _log_and_print(log=log_file,
                   message='Found pandas version {}.\n'.format(
            pandas.__version__))
except ImportError:
    _log_and_print(log=log_file,
                   message='Cannot find Pandas in current Python environment. Please ensure you have a recent version of Pandas, then start again. It can usually be installed using conda or from github.')

# numpy for data handling
try:
    import numpy
    _log_and_print(log=log_file,
                   message='Found numpy version {}.\n'.format(
            pandas.__version__))
except ImportError:
    _log_and_print(log=log_file,
                   message='Cannot find numpy in current Python environment. Please ensure you have a recent version of numpy, then start again. It can usually be installed using conda or from github.')

# numpy for data handling
try:
    import matplotlib
    _log_and_print(log=log_file,
                   message='Found matplotlib version {}.\n'.format(
            matplotlib.__version__))
    from matplotlib import pyplot
except ImportError:
    _log_and_print(log=log_file,
                   message='Cannot find matplotlib in current Python environment. Please ensure you have a recent version of matplotlib, then start again. It can usually be installed using conda or from github.')
    if not args.trj:
        _log_and_print(log=log_file,
                       message='Since there is no trajectory input, matplotlib is not required to run the full script.')
    else:
        sys.exit()

if __name__ == '__main__':

    # Sanity check
    log_file.write('Sanity check...\n')
    
    # Check if:
    
    # Structure exists:
    if args.recept is None:
        _log_and_print(log=log_file,
                       message='Error. Is the input file correct?')
        sys.exit()
    
    # If more than one output type is selected
    if args.run_type not in ["short", "long"]:
        _log_and_print(log=log_file, 
                       message="Please pick only one output type. Consult the documentation for more.")
        sys.exit()
    
    try:
        struct = Universe(args.recept)
    except Exception:
        _log_and_print(log=log_file,
                       message='There is a problem with your input structure. Please check your file and try again.')
        sys.exit()
    
    position50_data = load_from_position50_file(args.cons_file)
    bw_numbers = gen_bw_numbers(position50_data)
    chicos = convert_pairs_to_fulllength(chicos_residues=chicos_residues,
                                         bw_numbers=bw_numbers)
    
    if args.trj is None:
        minDistArray, GLAS = solo_protocol(system=struct,
                                           chicos_dict=chicos,
                                           cutoff=args.cutoff,
                                           debug=args.debug)
        
        if args.run_type == "short":
            short_static_GLAS(system=struct, GLAS=GLAS)
            _log_and_print(log=log_file,
                           message="Script completed. Thank you for using pyGLAS!")
            sys.exit()
        
        elif args.run_type == "long":
            long_static_GLAS(system=struct,
                             chicos_dict=chicos, 
                             array=minDistArray,
                             GLAS=GLAS)
            _log_and_print(log=log_file,
                           message="Script completed. Thank you for using pyGLAS!")
            sys.exit()
        
    elif args.trj is not None:
        try:
            struct.load_new(args.trj)
        except Exception:
            _log_and_print(log=log_file,
                           message='There is a problem with your input trajectory. Please check your file and try again.')
            sys.exit()

        dGLAS, distsDict = dynamic_GLAS(system=struct,
                                        chicos_dict=chicos,
                                        cutoff=args.cutoff,
                                        debug=args.debug)
        
        if args.run_type == "short":
            short_dynamic_GLAS(system=struct,
                               dGLAS=dGLAS)
            _log_and_print(log=log_file,
                           message="Script completed. Thank you for using  pyGLAS!")
            sys.exit()
        
        elif args.run_type == "long":
            long_dynamic_GLAS(system=struct,
                              dGLAS=dGLAS,
                              distsDict=distsDict)
