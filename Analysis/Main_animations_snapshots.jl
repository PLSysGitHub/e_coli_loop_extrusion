"""
Make animations and snapshots of simulations
"""

include("src/Trajectory_animations.jl")

#Some parameters
sample_from=1000 #which block number sampling should start from
skip_done=true #whether or not exisiting stat files should be overwritten
monomer_size=33/1000 #μm
figure_type="png" #can also save as PDF; just change to "pdf"

#where to save data and plots
config_dir="../Results_3D/Combined_steady_state/" #Directory where simulation output .h5 are stored
plot_dir= "Plots/Steady_state_from_$sample_from/" #plots are saved here

#where to save data and plots
snapshot_dir= joinpath(plot_dir, "Snapshots/") #plots are saved here
animation_dir=joinpath(plot_dir,"Animations/")

#First some snapshots
snapshots_per_sim(config_dir, snapshot_dir; skip_done=skip_done, file_type=figure_type, monomer_size=monomer_size, which_files=[1, 5])
#Now make the animations
one_animation_per_sim(config_dir, animation_dir; skip_done=skip_done, which_files=[1, 5])
