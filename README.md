## Loop-extruder mediated stiffening can globally order bacterial chromosomes

This repository contains code for running simulations of loop-extrusion in an *E. coli* like system for the manuscript [Loop-extruder mediated rigidity can globally order bacterial chromosomes](https://www.biorxiv.org/content/10.1101/2024.10.10.617531v1), and for analyzing and plotting the resulting data.

The repository is structured as follows:

### looplib bacterial 
1D loop extrusion code forked from Anton Goloborodko's [github.com/open2c/looplib](https://github.com/open2c/looplib)

Edits include adding periodic boundary conditions for circular bacterial chromosomes, as well as implementing loop-extruder bypassing upon collisions

Running a loop-extruder simulation with the default length typically takes approximately a minute

### main simulations
Python code forked from Hugo Brandao's [github.com/hbbrandao/bacterialSMCtrajectories](https://github.com/hbbrandao/bacterialSMCtrajectories)

To run simulations, you need to install [polychrom](https://github.com/open2c/polychrom).

The code first runs 1D loop extrusion simulations, and then a 3D polymer simulation where the saved loop-extruder trajectories are used to constrain the polymer.

Edits include a container class for storing information about bacteria to simulate, as well as the ability to load 1D loop extrusion data from files

Depending on the number of loop-extruders, running a set of 5 simulation trajectories can take 2-6 days on an average GPU

Simulation files can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.13908509)

### analysis
Based on my code for analysis in a previous paper [github.com/PLSysGitHub/loop-extrusion_with_replication_analysis](https://github.com/PLSysGitHub/loop-extrusion_with_replication_analysis)

Analyze simulation data and then make plots

The directory also includes output figures. Processed data files can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.13908509)

Segments of the analysis pipeline take 1-10 minutes
