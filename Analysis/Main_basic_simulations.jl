"""
This is the pipeline that was used to analyse simulations with varying confinement dimensions, M, λ.

It does the following:
- Calculate mean and std of monomer long axis positions etc
- Calculate the orientation angles corresponding to different trajectories
- Calculate the rates at which the L3 and R3 loci switch positions
- Make plots both for individual simulations and comparing simulations

The script was used to generate plots for Figures 2 and 3.
"""

include("src/Unreplicating_analysis.jl")
include("src/Unreplicating_plots.jl")

sample_from=1000 #which block number sampling should start from
skip_done=false #whether or not exisiting stat files should be overwritten
monomer_size=33/1000 #μm
figure_type="pdf" #can also save as png, svg etc 

#where to save data and plots
config_dir="../Results_3D/Combined_steady_state/" #Directory where simulation output .h5 are stored
out_dir= "Stats/Steady_state_from_$sample_from/" #statistics are saved here
plot_dir= "Plots/Steady_state_from_$sample_from/" #plots are saved here
convergence_dir="Stats/Steady_state_convergence/" #stats as a function of simulation time

#first do mean z-positions
calc_fast_unreplicating_stats(config_dir,out_dir, skip_done=skip_done, sample_from=sample_from)
plot_fast_unreplicating_stats(out_dir, plot_dir, monomer_size, figure_type)

#Do cosine fits to z-axis data to get orientation angles
calc_angle_trajectories(config_dir, out_dir, skip_done=skip_done)
calc_angle_trajectories(config_dir, out_dir, skip_done=skip_done, second_order=true)

plot_examples_angle_fits(config_dir, plot_dir, figure_type; sampled_ts=[2000])
calc_mean_angle_traj_stats(out_dir)
plot_all_angle_trajectories(out_dir, plot_dir, figure_type)
plot_all_angle_trajectories(out_dir, plot_dir, figure_type, true)
plot_all_angle_histograms(out_dir, plot_dir, figure_type)

#Like Makela et al 2021, calculate number of times theres a swtich LR->RL
calc_lr_flips(config_dir, out_dir)

#And then some plots comparing different conditions
compare_unreplicating_stats(out_dir, plot_dir, figure_type, nu=0.5)
compare_unreplicating_stats(out_dir, plot_dir, figure_type, only_original_dims=true, nu=0.5)

for params in [(50,200), (100, 125), (100, 100), (30, 250), (30, 333)]
    M=params[1]
    λ=params[2]
    compare_varying_size_stats(out_dir, plot_dir, figure_type, num_smcs=M, loop_size=λ)
end

compare_varying_ter_stats(out_dir, plot_dir, figure_type)

compare_smc_mutants(out_dir, plot_dir, figure_type)
compare_smc_mutants(out_dir, plot_dir, figure_type, λ=150)

compare_LE_stats_to_predictions(out_dir, plot_dir, figure_type)

#Check convergence
check_unreplicating_convergence(config_dir,convergence_dir)
make_convergence_plots(convergence_dir, plot_dir, monomer_size,filetype=figure_type)

#Save a file with the number of simulations for each condition
check_num_samples_steady_state(config_dir, convergence_dir, sample_from)
