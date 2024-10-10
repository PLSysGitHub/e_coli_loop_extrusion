"""
This file runs the same pipeline for simulations without bypassing, with a different starting position, with loop-extruders fixed in space, or with a monomer tethered to the cell pole.

It calculates:
- Statistics like the mean and std of the monomer long axis positions
- The orientation angles for all configurations
- Convergence plots for several statistics
- The number of simulation trajectories for each set of simulation parameters
- Makes plots for individual simlations

"""

include("src/Unreplicating_analysis.jl")
include("src/Unreplicating_plots.jl")
figure_type="pdf" #can also save as png, svg etc 

sample_from=1000 #which block number sampling should start from
skip_done=true #whether or not exisiting stat files should be overwritten
monomer_size=33/1000 #μm

for config_dir in ["../Results_3D/Combined_start_top_monomer_0", "../Results_3D/Steady_state_tether_monomer_1150","../Results_3D/Combined_no_bypassing", "../Results_3D/Steady_state_fixed_smcs"]
    #where to save results
    out_dir=replace(config_dir, "../Results_3D/"=>"Stats/")*"_from_$sample_from/"
    plot_dir=replace(out_dir, "Stats/"=>"Plots/")
    convergence_dir=replace(config_dir, "../Results_3D/"=>"Stats/")*"_convergence/"

    #first do mean z-positions etc
    calc_fast_unreplicating_stats(config_dir,out_dir, skip_done=skip_done, sample_from=sample_from)
    plot_fast_unreplicating_stats(out_dir, plot_dir, monomer_size, figure_type)
    
    #Do cosine fits to z-axis data to get orientation angles
    calc_angle_trajectories(config_dir, out_dir, skip_done=skip_done)
    plot_examples_angle_fits(config_dir, plot_dir, figure_type; sampled_ts=[2000])
    calc_mean_angle_traj_stats(out_dir)
    plot_all_angle_trajectories(out_dir, plot_dir, figure_type)
    plot_all_angle_histograms(out_dir, plot_dir, figure_type)

    #We need some additional plots for simulations with tethers or no bypassing
    if contains(config_dir, "tether")
        compare_z_stds(out_dir, plot_dir, figure_type)
    elseif contains(config_dir,"no_bypassing") 
        calc_lr_flips(config_dir, out_dir)
        compare_unreplicating_stats(out_dir, plot_dir, figure_type, only_original_dims=true, nu=1/2)
        compare_LE_stats_to_predictions(out_dir, plot_dir, figure_type, nu_loop=1/2)
    end

    #Check convergence
    check_unreplicating_convergence(config_dir,convergence_dir, orient_arms=contains(config_dir, "top_monomer_0"))
    make_convergence_plots(convergence_dir, plot_dir, monomer_size,filetype=figure_type)

    #Save a file with the number of simulations for each condition
    check_num_samples_steady_state(config_dir, convergence_dir, sample_from)
end

