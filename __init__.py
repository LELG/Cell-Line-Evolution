"""
POPLN: a tumour evolution and genetic instability simulation
============================================================

This python program simulates the growth and mutation of a
(generic) cancer tumour. The smallest simulated unit is a
genetically homogeneous cancer clone; the tumour is comprised
of a number of such clones.

The program simulates the introduction, and in some cases
maintenance, of medical treatment, in the form of a global
penalty on proliferation rates.

Further documentation is available in module and function docstrings.

Contents
--------

:: Python modules

    analytics     --- Track / analyse data
    compact       --- ???
    dropdata      --- Export data to do w/ heterogeneous populations
    main          --- Parse parameters and run simulation
    plotdata      --- Plot results
    population    --- Class and functions for entire tumour
    simulator     --- High-level simulation control and logic
    subpopulation --- Class, functions for individual clones
    treatment     --- Class, functions for treatment
    tree_to_xml   --- Export phylogenetic tree as XML file
    snapshot      --- Store and load population 'snapshots'
    utils         --- Various utility functions
    constants     --- Various constants used throughout the package

:: Scripts etc.

    run_simulation.sh --- Parse config file, create directories, run sim
    default.conf      --- Sample config file
"""
