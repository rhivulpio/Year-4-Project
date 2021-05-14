# Introduction to the project

This is a Year 4 MSCi Particle Physics Project in which a statistical analysis of Monte Carlo simulated data of the LHeC electron-proton collider was conducted. The goal is to ultimately determine the precision with which the mass of the Higgs boson can be measured.

It was determined that the focus for this project would be the H->ZZ*->4l (l=e or mu) decay mode; these are the events that you will find in the signal input file.

In this repository you will find my ROOT/C++ analysis which was used to verify the validity of the Monte Carlo simulated data, and then plot the reconstructed invariant 4 lepton mass. Following this, a statistical analysis was conducted in order to obtain the final numerical result. It was found that, based on this simulated data of the LHeC, and at the integrated luminosity goal of 1 ab^{-1}, the Higgs boson mass can be measured with a precision of 1.69 ± 0.32 GeV.

# Explanation of the contents of this repository

The code "Process.cxx" takes as input a file from the Delphes simulation framework (https://cp3.irmp.ucl.ac.be/projects/delphes) and was used to implement a physics analysis event selection and fill some histograms. This file was used to create histograms for the signal and background events, which could then be compared by another .cxx file.

The "SignalBkgd.cxx" file was used to compare the signal and background outputs from the "Process.cxx" file. The main aim of this was to calculate the significance and signal to background ratio for the data, and create mass reconstructions that contained both the signal and the background data on one plot. This file was used to obtain the final stacked mass reconstruction.

The "Utilities.cxx" file was used to store global constants and function declarations that needed to be used across both "Process.cxx" and "SignalBkgd.cxx".

The "Process.h" file was used to declare all of the histograms that were plotted, and to store function declarations that were used in the "Process.cxx" file.

The "output.log" file was created as a way to view the particle data for a selected number of events in order to help with debugging the code.

# How to compile the code

Set up the necessary environment:

cc7

source /disk/moose/general/asc/Build_cc7_Delphes342_Py8243/SetupEnv.sh

Compile the code:

make

# How to run the code

The program takes two arguments, the first is the input simulation file and the second is the name of a new output ROOT file which will contain the histograms produced by the code. To run do:

./Process [path/to/inputfile] [output.root]

(where you will need to substitute these dummy arguments with the appropriate path to the input file and your chosen name for the output file)

Input file for signal events: /disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_H4l.Delphes.root
Input file for ZZ*->4l background events: /disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_Zll.100k.v2.Delphes.root
Input file for Z->4l background events: /disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_Z4l.10k.Delphes.root

# Making changes to the code

Whenever you update the C++ code, make sure you re-compile the code before re-running it by issuing the following command again:

make

# Whenever you open a new terminal

Repeat the steps to setup the environment whenever you open a new terminal window, before you attempt to compile or run the code:

cc7

source /disk/moose/general/asc/Build_cc7_Delphes342_Py8243/SetupEnv.sh

