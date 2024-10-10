"""
This file was used to check whether the diffusion of the chromosome segments in the simulations
occurrs at the expected rate.

The diffusion constant in the presence of loop-extruders was of the same order of magnitude
as that found by Weber et.al. 2012.
"""

include("src/Unreplicating_plots.jl")
include("src/Unreplicating_analysis.jl")


fig_type=".pdf" #can also save png or other formats
nt=500 #total 500 s, sample every 1 s

#Constants
monomer_size=33/1000 #μm
N=4600
data_dir="../Results_3D/Short_times/"
out_dir="Stats/Short_times/"
base_plot_dir="Plots/Check_time_scales/"
skip_done=true #whether or not exisiting stat files should be overwritten

mkpath(out_dir)

#For each directory, calculate the mean squared distance of the origin and all monomers
for subdir in subdirs(data_dir)
    initial_stage_files=readdir(subdir, join=true)
    initial_stage_files=filter(contains("blocks"), initial_stage_files)
    plot_directory=replace(subdir,data_dir=>base_plot_dir)*"/"
    L_0=parse_after(subdir, "_L_")*monomer_size
    stats_file=replace(subdir, data_dir=>out_dir)*".h5"

    if !skip_done || !isfile(stats_file)
        mean_rsq, error_rsq=fetch_msd_trajectories(initial_stage_files)
        h5open(stats_file, "w") do f
            @write f mean_rsq
            @write f error_rsq
        end
    else
        mean_rsq=h5read(stats_file, "mean_rsq")
        error_rsq=h5read(stats_file, "error_rsq")
    end
    mkpath(plot_directory)

    ts=findall(mean_rsq[1,:].>0)
    mean_rsq=mean_rsq[:,ts]

    plot_mean_rsq(ts, mean_rsq)
    savefig(plot_directory*"diffusion_early_all_loci"*fig_type)

    plot_mean_rsq(ts.*3, mean_rsq)
    savefig(plot_directory*"diffusion_early_triple_time_all_loci"*fig_type)

    plot_mean_rsq(ts.*2, mean_rsq)
    savefig(plot_directory*"diffusion_early_double_time_all_loci"*fig_type)
end

