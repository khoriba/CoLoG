#!/bin/bash

if [ -z $CONDA_PREFIX ]; then
  echo 'Please activate conda environment and try again.'
  exit 1
fi

set -eu

echo Starting setup.

echo Downloading picard.
mkdir -p $CONDA_PREFIX/share/picard-2.10.10/
wget -O $CONDA_PREFIX/share/picard-2.10.10/picard.jar https://github.com/broadinstitute/picard/releases/download/2.10.10/picard.jar
echo Done.

echo Cloning SimSeq.
git clone https://github.com/jstjohn/SimSeq $CONDA_PREFIX/share/SimSeq
git clone https://github.com/nasasaki/Alt_nCov2019_primers $CONDA_PREFIX/share/Alt_nCov2019_primers
chmod +x $CONDA_PREFIX/share/Alt_nCov2019_primers/tools/plot_depth.py
echo Done.

echo Copying scripts.
cp COLOG/scripts/* $CONDA_PREFIX/bin/
echo Done.

echo Preparing reference.
mkdir -p $CONDA_PREFIX/opt/reference
cp COLOG/reference/* $CONDA_PREFIX/opt/reference/
echo Done.

echo Setup completed successfull.
