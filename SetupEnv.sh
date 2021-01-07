#!/bin/bash


export HCDIR=$(readlink -f $(dirname ${BASH_SOURCE}))
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase

source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

lsetup "root 6.10.04-x86_64-slc6-gcc62-opt"
lsetup "boost boost-1.62.0-python2.7-x86_64-slc6-gcc62"
lsetup "cmake"
export BOOST_ROOT=$ALRB_BOOST_ROOT

#setup consistant root version
BASEDIR=/disk/moose/atlas/hsg2/PowhegPythiaDelphes/sl6

export FASTJET=$BASEDIR/fastjet-install
export OPENMPI=$BASEDIR/openmpi-install
export CPPFLAGS=" -I$OPENMPI/include"
export MPICXX=$OPENMPI/bin/mpic++
export PYTHIA8=$BASEDIR/PYTHIA8-install
export DELPHES=$BASEDIR/Delphes-3.4.1
export LHAPDF=$BASEDIR/LHAPDF-install
export POWHEG=$BASEDIR/POWHEG-BOX-V2
export OPENLOOPS=$BASEDIR/OpenLoops-1.3.1
export HEPMC=$BASEDIR/HepMC-install
export SHERPA=$BASEDIR/Sherpa-install-MPI

# ASC Hack!
export PYTHIA8=/disk/moose/general/asc/PYTHIA8-install/
export DELPHES=/disk/moose/general/asc/Delphes-3.4.1/

export PATH=$HCDIR/Software/Build/bin:$SHERPA/bin:$LHAPDF/bin:$PYTHIA8/bin:$FASTJET/bin:$OPENMPI/bin:$DELPHES:$PATH
export PYTHONPATH=$DELPHES/python:$PYTHONPATH
export LD_LIBRARY_PATH=$SERPA/lib/SHERPA-MC:$HEPMC/lib:$PYTHIA8/lib:$OPENMPI/lib:$DELPHES:$LD_LIBRARY_PATH
