from argparse import ArgumentParser as ArgP

import os

import sys

import MDAnalysis
from MDAnalysis import Universe
from MDAnalysis.lib import distances as cdists

import numpy as np

import csv

import itertools

from ggplot import *

import pandas as pd
from pandas import Series

from Bio.PDB import *

import math

import pickle

###############################################################################

parser = ArgP(
    description="pyGLAS takes a GPCR structure and determines it's GPCR-Likeness Assessment Score (reference).\n \
    Additionally, given a trajectory, it will plot the state of the NACHOs during your trajectory.")
parser.add_argument('-r', '--receptor', dest='recept', default=None, help='Receptor file in PDB format.')
parser.add_argument('-x', '--trajectory', dest='trj', nargs="*", default=None, )
parser.add_argument('--chicos', dest='chicosfile', default=None, help='Add list of residue numbers '
                                                                      '(see examples/Single_Topology for an example.')
parser.add_argument('--nachos', dest='nachosfile', default=None, help='Add list of residue numbers '
                                                                      '(see examples/Single_Topology for an example.')
parser.add_argument('--no_nachos', dest='no_nachos', action='store_true', help='NACHOs will not be calculated.')
parser.add_argument('--no_chicos', dest='no_chicos', action='store_true', help='CHICOs will not be calculated.')
parser.add_argument('--format', dest='run_type', default="short", help='Specifies run type. Please see documentation'
                                                                       'for more information.')
parser.add_argument('--cutoff', dest='cutoff', default=4.0, help='Defines cutoff. Default is 4A.')

args = parser.parse_args()

###############################################################################

# Location check
print "Output type check...\n"

# If more than one output type is selected
if 0 > len(args.run_type) > 1:
    print "Please pick only one output type. Consult the documentation for more."
    sys.exit()

# Check for compatibility of input/output
if args.run_type == "plain" and args.no_chicos == True:
    print "The 'plain' option is meant to calculate GLAS and output a simple text file. Very useful when used " \
          "alongside modelling software such as Modeller (reference). Please ensure you do not flag '--no_chicos' " \
          "when using the plain option."
    sys.exit()

if args.nachosfile is True and args.no_nachos is True:
    print "You are receiving this error message because you have requested I simultaneously calculate the NACHOs " \
    "and don't calculate the NACHOs. I will now return to the dark library from which I was called and disturb " \
    "everybody with the crunching sound of these delicious triangles."
    sys.exit()

if args.chicosfile is True and args.no_chicos is True:
    print "You are receiving this error message because you have requested I simultaneously calculate the CHICOs " \
    "and don't calculate the CHICOs. Facing this existential crisis, I will now return to my extravagant life as a " \
    "singing sensation. It's pyGLAS time!"

# Otherwise, proceed
print "Checks out!\n"

# Location check
print "Checking input...\n"

# Warn user that we will attempt to get NACHOs and CHICOs from GPCRDb.
if args.chicosfile is None and args.nachosfile is None:
    print "Since you didn't input any CHICOs or NACHOs files, I will attempt to calculate them from the input " \
          "structure. I will pass the input receptor to the GPCRDb and retrieve the results (coming soon...). " \
          "Please ensure you cite the GPCRDb (see references directory) if this shows up in the screen."

###############################################################################

###############################################################################
#                                                                             #
#                 Functions used in the script listed here                    #
#                                                                             #
###############################################################################

###############################################################################

chicosTemplate = {
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

nachosTemplate = {
    1: [3.43, 6.40],
    2: [3.43, 6.41],
    3: [3.46, 6.37],
    4: [3.50, 6.37],
    5: [6.36, 7.53],
    6: [6.40, 7.49],
    7: [1.53, 7.53],
    8: [3.50, 6.40],
    9: [3.43, 7.49],
    10: [3.43, 7.53],
    11: [3.46, 7.53],
    12: [5.55, 6.41],
    13: [5.58, 6.40],
    14: [5.62, 6.37],
    15: [2.43, 7.54]
}


#######################    Functions to load input    #########################

def system_load(topology, trajectory):
    """Load topology into MDAnalysis Universe with or without trajectory info"""

    # Define system as a global variable; I don't particularly like this idea.
    global system

    # Load into system topologies only
    if args.trj is None:
        system = MDAnalysis.Universe(args.recept)
        # For flow purposes
        print "Successfully loaded topology.\n"
        return system
    else:
        system = MDAnalysis.Universe(args.recept, args.trj)
        # For flow purposes
        print "Successfully loaded topology and trajectories.\n"
        return system


def chicos_from_structure(chicosTemplate, receptor):
    """Parses file and extracts BW numbers from bfactor column"""

    print "Attempting to get CHICOs from input file...\n"

    # Start a parser instance in BioPython
    parser = PDBParser()

    # Load structure into parser
    recept = parser.get_structure('User_Receptor', receptor)

    # Define receptor_chicos as a global variable; I don't particularly like this idea.
    global receptor_chicos

    receptor_chicos = {}

    bw_residues_dict = {}

    # This parses the input receptor and outputs a list of residues with all required information
    residue_list = Selection.unfold_entities(recept, 'R')

    # And now with a simple for loop...
    for residue in residue_list:

        if residue.get_id()[0] == ' ':

            # We get the bfactor from the file (which is already the BW number)
            bwnumber = residue['N'].bfactor

            if 9 > bwnumber > 0:
                # And add it to a dictionary
                bw_residues_dict[bwnumber] = residue.get_id()[1]

    for pair in chicosTemplate:
        receptor_chicos[pair] = [bw_residues_dict[chicosTemplate[pair][0]], bw_residues_dict[chicosTemplate[pair][1]]]

    print "Since you asked me to find the CHICOs, I will write them to a file so you can verify and confirm them.\n "
    "If any are wrong, please change them and re-run pyGLAS using the flag --chicos followed by the CHICOs file.\n"

    with open("{}_CHICOs.txt".format(os.path.splitext(os.path.basename(system.filename))[0]), 'w') as chicos_file:
        for key in receptor_chicos:
            chicos_file.write("{}\t{}\t{}\n".format(key, receptor_chicos[key][0], receptor_chicos[key][1]))

    return receptor_chicos


def nachos_from_structure(nachosTemplate, receptor):
    """Parses file and extracts BW numbers from bfactor column"""
    print "Attempting to get NACHOs from input file...\n"

    # Start a parser instance in BioPython
    parser = PDBParser()

    # Load structure into parser
    recept = parser.get_structure('User_Receptor', receptor)

    # Define receptor_nachos as a global variable; I don't particularly like this idea.
    global receptor_nachos
    receptor_nachos = {}

    bw_residues_dict = {}

    # This parses the input receptor and outputs a list of residues with all required information
    residue_list = Selection.unfold_entities(recept, 'R')

    for residue in residue_list:
        if residue.get_id()[0] == ' ':
            bwnumber = residue['N'].bfactor
            if 9 > bwnumber > 0:
                bw_residues_dict[bwnumber] = residue.get_id()[1]

    for pair in nachosTemplate:
        receptor_nachos[pair] = [bw_residues_dict[nachosTemplate[pair][0]], bw_residues_dict[nachosTemplate[pair][1]]]

    print "Since you asked me to find the NACHOs, I will write them to a file so you can verify and confirm them.\n"
    "If any are wrong, please change them and re-run pyGLAS using the flag --nachos followed by the NACHOs file.\n"

    with open("{}_NACHOs.txt".format(os.path.splitext(os.path.basename(system.filename))[0]), 'w') as nachos_file:
        for key in receptor_nachos:
            nachos_file.write("{}\t{}\t{}\n".format(key, receptor_nachos[key][0], receptor_nachos[key][1]))

    return receptor_nachos


def conf_parse_dictionary(file, type):
    """Parse user input (works for both CHICOs and NACHOs)"""

    # Must ensure that data is in a tab separated file; simpler to parse.
    data = csv.reader(open(file), delimiter="\t")

    pairs = {}

    for line in data:
        pairs[int(line[0])] = [line[1], line[2]]

    if type == "chicos":
        global chicosDict
        chicosDict = pairs
        return chicosDict
    elif type == "nachos":
        global nachosDict
        nachosDict = pairs
        return nachosDict


######################    Functions to output results    ######################

def plain_output(system, GLAS):
    """Trial output format. To be used in conjunction with Modeller."""

    # Simply put, it will output a very plain file with a number (GLAS).
    filename = str("%s_GLAS.txt" % os.path.splitext(os.path.basename(system.filename))[0])
    file = open(filename, 'w')
    file.write("{}".format(GLAS))


def short_static_GLAS(system, GLAS):
    """Simple. Writes GLAS to a file that starts with input filename and ends in .short. """
    # Use in not long form.
    filename = str("%s_GLAS.short" % os.path.splitext(os.path.basename(system.filename))[0])

    file = open(filename, 'w')
    file.write("##################################################\n")
    file.write("#           Thank you for using pyGLAS           #\n")
    file.write("#         The GLAS of your receptor is {}        #\n".format(GLAS))
    file.write("#  Please ensure you cite all relevant sources.  #\n")
    file.write("##################################################\n")


def short_static_SAS(system, SAS):
    """Simple. Writes SAS to a file that starts with input filename and ends in .short. """
    # Use in not long form.
    filename = str("%s_SAS.short" % os.path.splitext \
        (os.path.basename(system.filename))[0])

    file = open(filename, 'w')
    file.write("##################################################\n")
    file.write("#           Thank you for using pyGLAS           #\n")
    file.write("#         The SAS of your receptor is {}         #\n".format(SAS))
    file.write("#  Please ensure you cite all relevant sources.  #\n")
    file.write("##################################################\n")


def long_static_GLAS(system, dict, array, GLAS):
    """Writes long form GLAS result. Writes output in the form of a tab separated value file with pair number and
    minimum distance."""
    filename = str("%s_GLAS.long" % os.path.splitext(os.path.basename(system.filename))[0])

    file = open(filename, 'w')

    file.write("##################################################\n")
    file.write("#           Thank you for using pyGLAS           #\n")
    file.write("#         The GLAS of your receptor is {}        #\n".format(GLAS))
    file.write("#  Please ensure you cite all relevant sources.  #\n")
    file.write("##################################################\n")
    file.write("\n")
    file.write("Since you have requested a long form output, below\n")
    file.write("are the 24 pairs with the minimum distance between\n")
    file.write("them.\n")
    file.write("\n")
    file.write("\tPair #\tMinimum Distance\n")

    for key, dist in itertools.izip(dict, array):
        file.write("\t{0}\t{1:.3f}\n".format(key, dist))


def long_static_SAS(system, dict, SAS):
    """Writes long form SAS result. Writes output in the form of a tab separated value file with pair number and
    minimum distance."""
    filename = str("%s_SAS.long" % os.path.splitext(os.path.basename \
                                                        (system.filename))[0])

    file = open(filename, 'w')

    file.write("##################################################\n")
    file.write("#           Thank you for using pyGLAS           #\n")
    file.write("#         The GLAS of your receptor is {}        #\n".format(SAS))
    file.write("#  Please ensure you cite all relevant sources.  #\n")
    file.write("##################################################\n")
    file.write("\n")
    file.write("Since you have requested a long form output, below\n")
    file.write("are the 15 pairs with the minimum distance between\n")
    file.write("them.\n")

    file.write("\tPair #\tMinimum Distance\n")

    for key in dict:
        file.write("\t{0}\t{1:.3f}\n".format(key, dict[key]))


def short_dynamic_GLAS(system, dGLAS):
    """Plots time dependent GLAS in a single plot. Also outputs time varied GLAS as a pickle file."""

    # Pickle data
    pickle.dump(dGLAS, open('GLAS_change_{}.p'.format(os.path.splitext(os.path.basename(system.filename))[0]), 'wb'))

    # Prepare dGLAS for ggplot
    dGLAS_prep = pd.DataFrame(dGLAS)
    # This creates a pandas DataFrame for ggplot. It needs to be transposed.

    dGLAS_prep_trans = dGLAS_prep.transpose()

    dGLAS_prep_trans.columns = ['Frame Number', 'GLAS']

    timed_GLAS = ggplot(aes(x = 'Frame Number', y = 'GLAS'), data=dGLAS_prep_trans) + \
                            ggtitle(element_text(text="GLAS for {}".format(system.filename), size=20)) + \
                            xlab(element_text(text="Frame Number", size=16)) + \
                            ylab(element_text(text="GLAS", size=16)) + \
                            xlim(0, int(max(dGLAS[0]))) + \
                            ylim(int(min(dGLAS[1])) - 1, 24) + \
                            geom_point()

    timed_GLAS.save("%s_GLAS_overTime.svg" % os.path.splitext(os.path.basename(system.filename))[0])


def long_dynamic_GLAS(system, dGLAS, dict):
    """Plots time dependent GLAS in a single plot and all minDists for all pairs for all times. Pickles data."""

    # Pickle data
    pickle.dump(distsDict,
                open('CHICOs_change_{}.p'.format(os.path.splitext(os.path.basename(system.filename))[0]), 'wb'))

    # Prepare dictionary containing contact info for plotting.
    # Make pandas DataFrame
    chicos_prep = pd.DataFrame.from_dict(distsDict)

    # Add index to DataFrame
    chicos_prep['index'] = Series(range(len(chicos_prep)))

    # Change variable to reflect pair number
    chicos_prep.columns = ['%s-%s' %(x,y) for x,y in chicosTemplate.values()] + ['index']

    # Melt DataFrame by index
    chicos_prep_molten = pd.melt(chicos_prep, id_vars=['index'])

    # Change columns to match plots
    chicos_prep_molten.columns = ['Frame Number', 'Pair Number', 'Min Distance (A)']

    # Plot
    fachicos = ggplot(aes(x='Frame Number', y='Min Distance (A)'), data=chicos_prep_molten) + \
                    theme_gray() + \
                    ggtitle(element_text(text="Minimum Distances for All Pairs", size=20)) + \
                    xlab(element_text(text="Frame Number", size=16)) + \
                    ylab(element_text(text="Min Distance (A)", size=16)) + \
                    ylim(0, math.ceil(chicos_prep_molten.max()[2])) + \
                    xlim(0, int(max(dGLAS[0]))) + \
                    facet_wrap(x='Pair Number', nrow=4) + \
                    geom_point()

    fachicos.save('CHICOs_MinimumDistances_for_{}.svg'.format(os.path.splitext(os.path.basename(system.filename))[0]))


def long_dynamic_nachos(system, dict):
    """Plots time dependent NACHOs i.e. all minDists for all NACHO pairs for all times. Pickles data."""

    # Pickle data
    pickle.dump(dict, open('NACHOs_change_{}.p'.format(os.path.splitext(os.path.basename(system.filename))[0]), 'wb'))

    # Prep data for plotting
    nachos_prep = pd.DataFrame.from_dict(dict)

    # Add index to DataFrame
    nachos_prep['index'] = Series(range(len(nachos_prep)))

    # Change variable to reflect pair number
    nachos_prep.columns = ['%s-%s' % (x, y) for x, y in nachosTemplate.values()] + ['index']

    # Melt DataFrame
    nachos_prep_molten = pd.melt(nachos_prep, id_vars=['index'])

    # Change column names
    nachos_prep_molten.columns = ['Frame Number', 'Pair Number', 'Min Distance (A)']

    facnachos = ggplot(aes(x='Frame Number', y='Min Distance (A)'), data=nachos_prep_molten) + \
                    theme_gray() + \
                    ggtitle(element_text(text="Minimum Distances for All NACHO Pairs", size=20)) + \
                    xlab(element_text(text="Frame Number", size=16)) + \
                    ylab(element_text(text="Min Distance (A)", size=16)) + \
                    xlim(0, math.ceil(nachos_prep_molten.max()[0])) + \
                    ylim(0, math.ceil(nachos_prep_molten.max()[2])) + \
                    facet_wrap(x='Pair Number', nrow=5) + \
                    geom_point()

    facnachos.save("{}_NACHOs_overTime.svg".format(os.path.splitext(os.path.basename(system.filename))[0]))


pass


######################    Functions to calculate GLAS    ######################

def solo_protocol(system, dict):
    """Uses system without trajectory to return GLAS. No NACHOs input."""

    # Location check
    print "Pre-initialising relevant variables and arrays...\n"

    # Create variable GLAS
    global GLAS
    GLAS = 0

    # To store minimum distances
    global minDistArray
    minDistArray = []

    for key in dict:

        # Assign row values to variables and convert to strings; excludes all hydrogen atoms,
        # such that distances are between heavy atoms only
        resid1 = str("resid %s and not name H*") % dict[key][0]
        resid2 = str("resid %s and not name H*") % dict[key][1]

        # Debug statement
        # Check if code is assigning strings correctly
        # print resid1, resid2

        # Grab coordinates from resid variables. WARNING: get_positions() will change in MDAnalysis 16.0!
        coords1 = system.select_atoms(resid1).get_positions()
        coords2 = system.select_atoms(resid2).get_positions()

        # Pass coords variables to dists.distance_array function
        ContactsMatrix = cdists.distance_array(coords1, coords2)

        # Append minimum to minDistArray
        minDistArray.append(np.amin(ContactsMatrix))

        # Restore ContactsMatrix to 0, or it will get confused.
        ContactsMatrix = 0

    # Take output from chicosMinDistArray
    if minDistArray is not None:
        for dist in minDistArray:
            if dist < args.cutoff:
                GLAS += 1

    return minDistArray, GLAS


def dynamic_GLAS(system, dict):
    """Uses the system with trajectory to calculate the time resolved GLAS"""

    # Location check
    print "Pre-initialising relevant variables and arrays...\n"

    # Create variable GLAS
    global GLAS
    GLAS = 0

    # Create array dGLAS, of size (2, n of frames)
    global dGLAS
    dGLAS = np.zeros((2, system.trajectory.n_frames))

    # To store minimum distances
    global minDistArray
    minDistArray = []

    global distsDict
    distsDict = {}

    for key in dict:
        distsDict[key] = []

    for ts in system.trajectory:

        # Assign frame number to 'x'-axis
        dGLAS[0][ts.frame] = ts.frame

        # Debug statement
        # Check where we are
        # print "Calculating GLAS for frame {}\n".format(ts)

        for key in dict:
            # Assign row values to variables and convert to strings
            resid1 = str("resid {} and not name H*").format(dict[key][0])
            resid2 = str("resid {} and not name H*").format(dict[key][1])

            # Debug statement
            # Check if code is assigning strings correctly
            # print resid1, resid2

            # Grab coordinates from resid variables. WARNING: get_positions()
            # will change in MDAnalysis 16.0!
            coords1 = system.select_atoms(resid1).get_positions()
            coords2 = system.select_atoms(resid2).get_positions()

            # Debug statement
            # Check if code is coordinates correctly
            # print coords1, coords2

            # Pass coords variables to dists.distance_array function
            ContactsMatrix = cdists.distance_array(coords1, coords2)

            distsDict[key].append(np.amin(ContactsMatrix))

            # Append minimum to minDistArray
            minDistArray.append(np.amin(ContactsMatrix))

            # Restore ContactsMatrix to 0, or it will get confused.
            ContactsMatrix = 0

        # Take output from minDistArray
        if minDistArray is not None:
            for dist in minDistArray:
                if dist < args.cutoff:
                    GLAS += 1

        # Reset minDistArray
        minDistArray = []

        # Assign GLAS to 'y'-axis
        dGLAS[1][ts.frame] = GLAS

        # Restart GLAS
        GLAS = 0

    return dGLAS, distsDict


def static_SAS(system, dict):
    """Uses system without trajectory to return SAS. No CHICOs input."""

    # Location check
    print "Pre-initialising relevant variables and arrays...\n"
    print "Warning: SAS is a score based on NACHOs. Please make sure you understand the concept before using " \
          "any potential results.\n"

    # Create variable SAS
    global SAS
    SAS = 0

    global minDistDict
    minDistDict = {}

    for key in dict:

        # Assign row values to variables and convert to strings
        resid1 = str("resid %s and not name H*") % dict[key][0]
        resid2 = str("resid %s and not name H*") % dict[key][1]

        # Debug statement
        # Check if code is assigning strings correctly
        # print resid1, resid2

        # Grab coordinates from resid variables. WARNING: get_positions() will change in MDAnalysis 16.0!
        coords1 = system.select_atoms(resid1).get_positions()
        coords2 = system.select_atoms(resid2).get_positions()

        # Pass coords variables to cdists.distance_array function
        ContactsMatrix = cdists.distance_array(coords1, coords2)

        # Append minimum to minDistArray
        minDistDict[key] = np.amin(ContactsMatrix)

        # Restore ContactsMatrix to 0, or it will get confused.
        ContactsMatrix = 0

    # Take output from minDistArray
    if minDistDict is not None:
        for key in minDistDict:
            # Enter the arbitrary score: SAS (System Activation Score)
            if key < 8:
                if minDistDict[key] < args.cutoff:
                    SAS -= 1
            elif key > 7:
                if minDistDict[key] < args.cutoff:
                    SAS += 1

    return minDistDict, SAS


def calc_NACHOs(system, dict):
    """Uses the system with trajectory to calculate the state of the NACHOs."""

    # Array to store minimum distances
    global minDistArray
    minDistArray = []

    global distsDict
    distsDict = {}

    for key in dict:
        distsDict[key] = []

    for ts in system.trajectory:

        for key in dict:
            # Assign row values to variables and convert to strings
            resid1 = str("resid {} and not name H*").format(dict[key][0]);
            resid2 = str("resid {} and not name H*").format(dict[key][1]);

            # Debug statement
            # Check if code is assigning strings correctly
            # print resid1, resid2

            # Grab coordinates from resid variables. WARNING: get_positions()
            # will change in MDAnalysis 16.0!
            coords1 = system.select_atoms(resid1).get_positions();
            coords2 = system.select_atoms(resid2).get_positions();

            # Debug statement
            # Check if code is coordinates correctly
            # print coords1, coords2

            # Pass coords variables to dists.distance_array function
            ContactsMatrix = cdists.distance_array(coords1, coords2)

            distsDict[key].append(np.amin(ContactsMatrix))

            # Append minimum to minDistArray
            minDistArray.append(np.amin(ContactsMatrix))

            # Restore ContactsMatrix to 0, or it will get confused.
            ContactsMatrix = 0

    return distsDict


pass


# Running GLAS

# Load system
system = system_load(args.recept, args.trj)

# Load CHICOs
if args.chicosfile is not None:

    # Load CHICOs using
    chicosDict = conf_parse_dictionary(args.chicosfile, "chicos")
    print "CHICOs loaded from user input.\n"

elif args.chicosfile is None:

    # Else look for chicos from input
    chicosDict = chicos_from_structure(chicosTemplate, args.recept)
    print "CHICOs loaded from input structure.\n"

# Load NACHOs
if args.nachosfile is not None:

    # Load from user input
    nachosDict = conf_parse_dictionary(args.nachosfile, "nachos")
    print "NACHOs loaded from user input.\n"

elif args.nachosfile is None and args.no_nachos is not True:
    # Load from GPCRDb output
    nachosDict = nachos_from_structure(nachosTemplate, args.recept)
    print "NACHOs loaded from input structure.\n"

# Solo GLAS run

# Run trajectory-less GLAS
if args.trj is None:

    if args.no_chicos is not True:
        # Activate solo protocol with loaded system and chicosDict
        solo_protocol(system, chicosDict)

    elif args.no_nachos is not True:
        # Calculate SAS with loaded system and nachosDict
        static_SAS(system, nachosDict)

    else:

        static_SAS(system, nachosDict)

        solo_protocol(system, chicosDict)

    # Output data; activates different protocols for different run_type
    if args.run_type == "long":

        if args.no_chicos is not True:
            # Outputs long version (GLAS + min distances)
            long_static_GLAS(system, chicosDict, minDistArray, GLAS)

        elif args.no_nachos is not True:
            # Outputs long version (SAS + min distances)
            long_static_SAS(system, minDistDict, SAS)

        else:

            long_static_GLAS(system, chicosDict, minDistArray, GLAS)

            long_static_SAS(system, minDistDict, SAS)

    elif args.run_type == "plain":

        # Outputs plain version (GLAS in a text file)
        plain_output(system, GLAS)

    else:

        if args.no_chicos is not True:
            # Default output: short.
            short_static_GLAS(system, GLAS)

        elif args.no_nachos is not True:

            short_static_SAS(system, SAS)

        else:

            short_static_GLAS(system, GLAS)

            short_static_SAS(system, SAS)

# Dynamic GLAS run

# Run GLAS with trajectory
if args.trj is not None:

    if args.no_chicos is not True:
        # Runs dynamic_GLAS and outputs dGLAS, a time resolved GLAS for additional analysis of protein stability.
        dynamic_GLAS(system, chicosDict)

    # If NACHOs are given, this returns time_nachos; useful in simulations that bias receptors to certain states.
    if args.no_nachos is not True:
        time_nachos = calc_NACHOs(system, nachosDict)

    if args.run_type == "short":

        if args.no_chicos is not True:
            # Returns 1 plot with time resolved GLAS
            short_dynamic_GLAS(system, dGLAS)

        if args.no_nachos is not True:
            # Returns 1 plot with time resolved SAS
            short_static_SAS(system, time_nachos)

    elif args.run_type == "long":

        if args.no_chicos is not True:
            # Returns 24 plots with time resolved contacts
            long_dynamic_GLAS(system, dGLAS, distsDict)

        elif args.no_nachos is not True:
            # Returns 15 plots with time resolved contacts
            long_dynamic_nachos(system, time_nachos)