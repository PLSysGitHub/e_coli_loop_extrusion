"""
This file runs analysis on simulations in an infinite tube.
It's in a separate script because some statistics don't make sense without a finite confinement. We're mainly interested in analysing what the rest length of the simulated chromosome is.
"""

include("src/Unreplicating_analysis.jl")
include("src/Unreplicating_plots.jl")

sample_from=1000 #which block number sampling should start from
skip_done=false #whether or not exisiting stat files should be overwritten
monomer_size=33/1000 #μm
figure_type="pdf" #can also save as PDF; just change to "pdf"

#where to save data and plots
config_dir="../Results_3D/Combined_infinite_tube/" #Directory where simulation output .h5 are stored
out_dir= "Stats/Infinite_tube_from_$sample_from/" #statistics are saved here
plot_dir= "Plots/Infinite_tube_from_$sample_from/" #plots are saved here
convergence_dir="Stats/Infinite_tube_convergence/" #stats as a function of simulation time

#first do mean z-positions, segregated fractions
calc_infinite_tube_stats(config_dir,out_dir, skip_done=false, sample_from=sample_from)
compare_infinite_tube_stats_to_predictions(out_dir, plot_dir, figure_type, monomer_size)

check_unreplicating_convergence(config_dir,convergence_dir)
make_convergence_plots(convergence_dir, plot_dir, monomer_size,filetype=figure_type)

#Save a file with the number of simulations for each condition
check_num_samples_steady_state(config_dir, convergence_dir, sample_from)

