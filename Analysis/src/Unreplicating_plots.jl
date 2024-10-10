"""
This file contains functions that loop through saved results and save plots.

Helper functions that actually make individual plots are stored in Plot_functions.jl

"""

using ProgressMeter
include("Plot_functions.jl")

"""
Fetch the flip frequency from the specified directory.

# Arguments
- `dir::String`: The directory containing the data.
- `out_dir::String`: The output directory where the flip frequencies file is located.

# Returns
- `flip_freq::Any`: The flip frequency read from the file.
"""
function fetch_flip_frequency(dir, out_dir)
    f_name=joinpath(out_dir, "flip_frequencies.h5")
    if dir[end]=='/'
        dir=dir[1:end-1]
    end
    h5open(f_name, "r") do f
        ks=keys(f)
        wanted_key=findfirst(contains(basename(dir)), ks)
        wanted_key=ks[wanted_key]
        flip_freq=read(f[wanted_key], "mean_freq")
        return flip_freq
    end
end

"""
Plot segregated fraction and long axis position of oris over time for steady state simulations.

Plots can be used to assess whether simulations have converged.
"""
function make_convergence_plots(data_dir, base_plot_dir, monomer_size; spacer=1, filetype="png")
    if filetype[1]!='.'
        filetype="."*filetype
    end
    plot_dir=joinpath(base_plot_dir, "Convergence/")
    mkpath(plot_dir)

    files=readdir(data_dir, join=true)
    files=files[contains.(files, "_convergence.h5")]
    for f in files
        out_plot=replace(f, data_dir=>plot_dir)
        out_plot=replace(out_plot, ".h5"=>"")
        mkpath(out_plot)

        zs=h5read(f, "zs")*monomer_size
        region_distances=h5read(f, "region_extensions")*monomer_size

        plot(zs[1,:], xlabel="Simulation time", ylabel="Mean long axis position [μm]", label=L"ori")
        plot!(zs[1150, :], label="90°")
        plot!(zs[3450,:], label="270°")
        savefig(out_plot*"/mean_zs_convergence"*filetype)

        plot(region_distances[1,:], xlabel="Simulation time", ylabel="Region extent [μm]", label=L"ori")
        plot!(region_distances[3, :], label=L"ter")
        savefig(out_plot*"/region_extent_convergence"*filetype)
    end
end

"""
Plot histogram of the sum squared residuals and the orientation angle for all simulations.
"""
function plot_all_angle_histograms(data_dir, base_plot_dir, filetype="png", second_order=false, first_time=1000, max_time=2000, max_res=500, nbins=30)
    if filetype[1]!='.'
        filetype="."*filetype
    end
    plot_dir=joinpath(base_plot_dir, "Fast_stats_all/")
    for dir in subdirs(data_dir)
        mkpath(replace(dir, data_dir=>plot_dir))
        if second_order
            filename=dir*"/angle_trajectories_second_order.h5"
        else
            filename=dir*"/angle_trajectories.h5"
        end
        out_plot=replace(dir, data_dir=>plot_dir)
        ang_trajs_all=[]
        residuals_all=[]
        h5open(filename, "r") do file
            for k in keys(file)
                angs=read(file[k], "angles")[first_time:end]
                res=read(file[k], "residuals")[first_time:end]
                angs=angs[res.<max_res]
                push!(residuals_all, res...)
                push!(ang_trajs_all,angs...)
            end
        end
        #residuals to get a feel for how good the fits are
        lbl="Mean $(round(mean(residuals_all), sigdigits=2))\nStd $(round(std(residuals_all), sigdigits=2))"
        histogram(residuals_all, xlabel="Sum squared residuals [μm^2]", norm=true, ylabel="Probability", color=8, label=lbl, legend=:topright)
        if second_order
            savefig(out_plot*"/ssr_histogram_second_order"*filetype)
        else
            savefig(out_plot*"/ssr_histogram"*filetype)
        end

        #then the angles
        if length(ang_trajs_all)>0
            histogram(ang_trajs_all, xlabel="Orientation angle", norm=:probability, ylabel="Probability", xlims=(-3.15, 3.15), ylims=(0,0.2), color=8, bins=nbins)
            annotate!(2, 0.17, "Std $(round(std(ang_trajs_all), sigdigits=2))")
            if second_order
                savefig(out_plot*"/angle_histogram_second_order"*filetype)
            else
                savefig(out_plot*"/angle_histogram"*filetype)
            end

            #Now modulo pi/2
            ang_trajs_all[ang_trajs_all.>pi/2].-=pi
            ang_trajs_all[ang_trajs_all.<-pi/2].+=pi
            histogram(ang_trajs_all, xlabel="Orientation angle mod π", norm=:probability, ylabel="Probability", xlims=(-3.15/2, 3.15/2), color=8, ylims=(0,0.2), bins=nbins)
            annotate!(1, 0.17, "Std $(round(std(ang_trajs_all), sigdigits=2))")
            if second_order
                savefig(out_plot*"/angle_mod_histogram_second_order"*filetype)
            else
                savefig(out_plot*"/angle_mod_histogram"*filetype)
            end
        end
    end
end

"""
Plot the angle trajectories and residuals for all simulations.
"""
function plot_all_angle_trajectories(data_dir, base_plot_dir, filetype="png", second_order=false, max_num_trajectories=5, max_time=1999)
    if filetype[1]!='.'
        filetype="."*filetype
    end
    plot_dir=joinpath(base_plot_dir, "Fast_stats_all/")
    for dir in subdirs(data_dir)
        mkpath(replace(dir, data_dir=>plot_dir))
        if second_order
            filename=dir*"/angle_trajectories_second_order.h5"
        else
            filename=dir*"/angle_trajectories.h5"
        end
        out_plot=replace(dir, data_dir=>plot_dir)
        ang_trajs_all=[]
        residuals_all=[]
        if second_order
            ang_trajs_ratio_amps=[]
        end
        h5open(filename, "r") do file
            for k in keys(file)
                angs=read(file[k], "angles")
                res=read(file[k], "residuals")
                push!(residuals_all, res)
                if second_order
                    amps=read(file[k], "amplitudes")
                    push!(ang_trajs_all, angs[:,1])
                    push!(ang_trajs_ratio_amps, abs.(amps[:,2])./abs.(amps[:,1]))
                else
                    push!(ang_trajs_all,angs)
                end
            end
        end
        long_enough=filter(x->length(x)>=max_time, ang_trajs_all)
        long_enough=[a[1:max_time] for a in long_enough]
        num_trajs=length(long_enough)
        ang_trajs=hcat(long_enough...)[1:max_time,:]

        if num_trajs<max_num_trajectories
            sampled=1:num_trajs
        else
            sampled=sample(1:num_trajs, max_num_trajectories, replace=false)
        end

        plot_angle_trajectories(ang_trajs[:,sampled])
        if second_order
            savefig(out_plot*"/angle_trajectories_second_order"*filetype)
            plot(ang_trajs_ratio_amps, ylabel=L"|A_2|/|A_1|", xlabel="Simulation time")
            savefig(out_plot*"/angle_trajectories_second_order_compare_amps"*filetype)
        else
            savefig(out_plot*"/angle_trajectories"*filetype)
        end
        plot_residual_trajectories(residuals_all)
        if second_order
            savefig(out_plot*"/ssr_second_order"*filetype)
        else
            savefig(out_plot*"/ssr"*filetype)
        end
    end
end

"""
Plot statistics like mean long axis positions for all simulations
"""
function plot_fast_unreplicating_stats(data_dir, base_plot_dir,monomer_size=0.033, filetype="png")
    if filetype[1]!='.'
        filetype="."*filetype
    end
    plot_dir=joinpath(base_plot_dir, "Fast_stats_all/")
    for dir in subdirs(data_dir)
        mkpath(replace(dir, data_dir=>plot_dir))
        filename=dir*"/steady_state_stats.h5"
        out_plot=replace(dir, data_dir=>plot_dir)
        mean_zs=h5read(filename, "mean_zs")*monomer_size
        N=length(mean_zs)
        L=time_to_height(0)
        ter_length=parse_after(dir, "ter_size_")
        ter_width=Int(ter_length/2)

        plot_unreplicating_av_z(mean_zs, ter_width; size=(550,400))
        savefig(out_plot*"/mean_zs"*filetype)
        
        std_zs=h5read(filename, "zs_std")*monomer_size
        plot_unreplicating_std_z(std_zs, ter_width)
        savefig(out_plot*"/std_zs"*filetype)

        region_extents=h5read(filename, "mean_region_extents").*monomer_size
        std_region_extents=h5read(filename, "std_region_extents").*monomer_size

        scatter(region_extents, yerr=std_region_extents, xlabel="Genomic position", ylabel="Mean region extent, [μm]", xticks=([1,2,3,4], ["0°", "90°", "180°", "270°"]))
        savefig(out_plot*"/region_extents"*filetype)

        #read the distributions of distances and plot histograms

        filename=dir*"/steady_state_distance_distributions.h5"
        ori_region=h5read(filename, "ori_extents").*monomer_size
        ter_region=h5read(filename, "ter_extents").*monomer_size
        l_region=h5read(filename, "left_extents").*monomer_size
        r_region=h5read(filename, "right_extents").*monomer_size

        density_plot(l_region, 50, left_color, xlabel="Region extension [μm]", ylabel="Probability", label="Left", xlims=(0,1.8), ylims=(-0.005, 0.08))
        density_plot!(r_region, 50, right_color, label="Right")
        density_plot!(ori_region, 50, ori_color, label=L"ori")
        density_plot!(ter_region, 50, ter_color, label=L"ter")
        savefig(out_plot*"/distance_distributions"*filetype)

        filename=dir*"/steady_state_position_distributions.h5"
        ori_region=h5read(filename, "ori_pos").*monomer_size
        ori_plus_region=h5read(filename, "ori_plus_pos").*monomer_size
        ter_region=h5read(filename, "ter_pos").*monomer_size
        ter_plus_region=h5read(filename, "ter_plus_pos").*monomer_size
        L0=time_to_height(0)

        density_plot(ori_region./L0, 40, ori_color, label=L"ori", xlabel="Relative position along nucleoid", ylabel="Probability",size=(800, 400))
        density_plot!(ori_plus_region./L0, 40, left_color, label=L"ori+500\mathrm{kb}")
        density_plot!(ter_region./L0, 40, ter_color, label=L"ter")
        density_plot!(ter_plus_region./L0, 40, right_color, label=L"ter+200\mathrm{kb}")
        savefig(out_plot*"/position_distributions"*filetype)
    end
end

"""
Compare simulations where the confinement size varies
"""
function compare_varying_size_stats(data_dir, base_plot_dir, filetype="png"; monomer_size=0.033, num_smcs=50, loop_size=150, L_0=54.42240234227263, r0=15.151515151515152)
    if filetype[1]!='.'
        filetype="."*filetype
    end
    filetype="_M_$(num_smcs)_lambda_$(loop_size)"*filetype
    plot_dir=joinpath(base_plot_dir, "Compare_varying_size/")
    mkpath(plot_dir)

    radii=Vector{Float64}()
    lengths=Vector{Float64}()
    ds_ter=Vector{Float64}() #800 kb; bounds at 190 and 270
    std_ter=Vector{Float64}()
    sem_lifetimes=Vector{Float64}()
    std_angle=Vector{Float64}()
    flip_freqs=Vector{Float64}()

    relevant_dirs=filter(contains("num_smcs_$(num_smcs)_loop_size_$(loop_size)"), subdirs(data_dir))
    relevant_dirs=filter(contains("ter_size_800_ter_strength_100"), relevant_dirs)

    for dir in relevant_dirs
        filename=dir*"/steady_state_stats.h5"
        mean_extents=h5read(filename, "mean_region_extents").*monomer_size
        std_extents=h5read(filename, "std_region_extents").*monomer_size
        r=parse_after(dir, "_r_")*monomer_size
        L=parse_after(dir, "_L_")*monomer_size

        push!(radii,r)
        push!(lengths, L)

        push!(ds_ter, mean_extents[3])
        push!(std_ter, std_extents[3])

        filename=dir*"/mean_angle_traj_stats.h5"
        push!(std_angle, h5read(filename, "std_angle_mod"))
        
        flip_f=fetch_flip_frequency(dir,data_dir)
        push!(flip_freqs, flip_f)
    end

    colors=:black
    ms=:circle
    x=1:length(radii)
    Ds=round.(radii.*2, sigdigits=2)
    Ls=round.(lengths, sigdigits=2)
    
    sortp=sortperm(std_angle, rev=true)
    Ds=Ds[sortp]
    Ls=Ls[sortp]
    flip_freqs=flip_freqs[sortp]
    std_angle=std_angle[sortp]
    ds_ter=ds_ter[sortp]
    std_ter=std_ter[sortp]

    xlab= ["($(L),$(D))" for (L,D) in Pair.(Ls, Ds)]

    suffix="varying_dims"*filetype

    scatter(std_angle, 1:length(std_angle), markershape=ms, yticks=(x,xlab), xlabel="Std angle mod π [rad]", ylabel=L"(L,D)\ [\mathrm{\mu m}]", color=colors, size=(550,400))
    savefig(plot_dir*"/std_angle_mod_pi"*suffix)

    scatter(ds_ter, 1:length(std_angle), xerr=std_ter, markershape=ms, yticks=(x,xlab), xlabel="Extension ter [μm]", ylabel=L"(L,D)\ [\mathrm{\mu m}]", color=colors, size=(550,400))
    savefig(plot_dir*"/extension_ter"*suffix)

    scatter(flip_freqs, xticks=(x,xlab), markershape=ms, color=colors,ylabel="Flipping frequency [1/sim time]", size=(550,400))
    savefig(plot_dir*"/flip_frequency"*suffix)
end

"""
Compare simulations where the ter region size varies
"""
function compare_varying_ter_stats(data_dir, base_plot_dir, filetype="png", monomer_size=0.033, num_smcs=50, loop_size=200.0, L_0=54.42240234227263, min_ters=[0,200], r0=15.151515151515152)
    if filetype[1]!='.'
        filetype="."*filetype
    end
    plot_dir=joinpath(base_plot_dir, "Compare_varying_ter_size/")
    mkpath(plot_dir)

    ter_lengths=Vector{Float64}()
    ds_ter=Vector{Float64}() #800 kb; bounds at 190 and 270
    std_ter=Vector{Float64}()
    std_angle=Vector{Float64}()
    flip_freqs=Vector{Float64}()

    relevant_dirs=filter(contains("e_coli_r_$(r0)_L_$(L_0)_num_smcs_$(num_smcs)_loop_size_$(loop_size)"), subdirs(data_dir))
    relevant_dirs=filter(contains("ter_strength_100"), relevant_dirs)

    for dir in relevant_dirs
        filename=dir*"/steady_state_stats.h5"
        mean_extents=h5read(filename, "mean_region_extents").*monomer_size
        std_extents=h5read(filename, "std_region_extents").*monomer_size
        loop_l=h5read(filename, "mean_loop_l")
        ter_length=parse_after(dir, "ter_size_")

        push!(ter_lengths, ter_length)

        if ter_length>0
            push!(ds_ter, mean_extents[3])
            push!(std_ter, std_extents[3])
        else
            push!(ds_ter,0)
            push!(std_ter, 0)
        end

        filename=dir*"/mean_angle_traj_stats.h5"        
        push!(std_angle, h5read(filename, "std_angle_mod"))
        
        flip_f=fetch_flip_frequency(dir,data_dir)
        push!(flip_freqs, flip_f)
    end

    colors=:black
    ms=:circle
    L=L_0*monomer_size
    x=ter_lengths
    xlab="ter length [kb]"
    suffix="_vs_ter"*filetype

    hline([L], ribbon=(L,0),fillcolor=cell_color, linealpha=0, fillalpha=0.3, xlims= (minimum(x)*0.9, maximum(x)*1.1)) 
    scatter!(x, ds_ter, yerr=std_ter, color=colors, xlabel=xlab, ylabel="Extension ter region [μm]", markershape=ms)
    savefig(plot_dir*"/ter_dist"*suffix)

    for min_ter in min_ters
        mask=x.>=min_ter
        suffix="_vs_ter_leq_$min_ter"*filetype

        scatter(x[mask], std_angle[mask], markershape=ms, xlabel=xlab, ylabel="Std angle mod π [rad]", color=colors)
        savefig(plot_dir*"/std_angle_mod_pi"*suffix)

        scatter(x[mask], flip_freqs[mask], xlabel=xlab, markershape=ms, color=colors,ylabel="Flipping frequency [1/sim time]")
        savefig(plot_dir*"/flip_frequency"*suffix)
    end
end

"""
Compare the standard deviations and means of the long axis positions for different simulations
"""
function compare_z_stds(data_dir, base_plot_dir, figtype="pdf"; num_smcs=50, λ=200, monomer_size=0.033, N=4600, ter_size=800, wt_dir="Stats/Steady_state_from_1000/e_coli_r_15.151515151515152_L_54.42240234227263_num_smcs_50_loop_size_200.0_colrate_0.03_trunc_5.0_ter_size_800_ter_strength_100")
    if figtype[1]!='.'
        figtype="."*figtype
    end
    figtype="_M_$(num_smcs)_lambda_$(λ)"*figtype

    plot_dir=joinpath(base_plot_dir, "Compare_z_std/")
    mkpath(plot_dir)

    relevant_dirs=subdirs(data_dir)

    #filter the MukBEF mutant
    no_smcs_dir=filter(contains("No_smcs_"), relevant_dirs)[1]
    #and the "MatP" mutant
    relevant_dirs=filter(contains("e_coli_r_15.151515151515152_L_54.42240234227263_num_smcs_$(num_smcs)_loop_size_$λ"), relevant_dirs)

    no_ter_dir=filter(contains("ter_size_$(ter_size)_ter_strength_1.0"), relevant_dirs)[1]

    std_in_pos=Vector{Vector{Float64}}()
    mean_pos=Vector{Vector{Float64}}()

    for dir in [wt_dir, no_smcs_dir, no_ter_dir]
        filename=dir*"/steady_state_stats.h5"
        push!(std_in_pos,h5read(filename, "zs_std").*monomer_size)
        push!(mean_pos,h5read(filename, "mean_zs").*monomer_size)
    end

    labels=["WT", "No SMCs", "No ter"]

    plot_unreplicating_std_z(std_in_pos[1];label=labels[1])
    plot_unreplicating_std_z!(std_in_pos[2];label=labels[2], alpha=0.3)
    plot_unreplicating_std_z!(std_in_pos[3]; label=labels[3], alpha=0.6, linestyle=:dash)
    savefig(plot_dir*"/std_long_axis_positions"*figtype)

    plot_unreplicating_av_z(mean_pos[1]; label=labels[1], size=(600,400))
    plot_unreplicating_av_z!(mean_pos[2]; label=labels[2], alpha=0.3)
    plot_unreplicating_av_z!(mean_pos[3]; linestyle=:dash, label=labels[3], alpha=0.6)
    savefig(plot_dir*"/mean_long_axis_positions"*figtype)
end

"""
Compare simulations without SMCs, without a ter region, and with both
"""
function compare_smc_mutants(data_dir, base_plot_dir, figtype="pdf"; num_smcs=50, λ=200, monomer_size=0.033, N=4600, ter_size=800)
    if figtype[1]!='.'
        figtype="."*figtype
    end
    figtype="_M_$(num_smcs)_lambda_$(λ)"*figtype

    plot_dir=joinpath(base_plot_dir, "Compare_SMC_mutants/")
    mkpath(plot_dir)

    relevant_dirs=subdirs(data_dir)

    #filter the MukBEF mutant
    no_smcs_dir=filter(contains("No_smcs_"), relevant_dirs)[1]
    #and the "MatP" mutant
    relevant_dirs=filter(contains("e_coli_r_15.151515151515152_L_54.42240234227263_num_smcs_$(num_smcs)_loop_size_$λ"), relevant_dirs)

    no_ter_dir=filter(contains("ter_size_$(ter_size)_ter_strength_1.0"), relevant_dirs)[1]

    #filter out the smcs file
    smcs_dir=filter(contains("ter_size_$(ter_size)_ter_strength_100"), relevant_dirs)[1]

    L_R_separate_halves=Vector{Float64}()
    L_R_separate_halves_error=Vector{Float64}()
    ds_ter=Vector{Float64}() #800 kb; bounds at 190 and 270
    std_ter=Vector{Float64}()
    flip_freqs=Vector{Float64}()

    std_in_pos=Vector{Vector{Float64}}()

    for dir in [no_smcs_dir, no_ter_dir, smcs_dir]
        filename=dir*"/steady_state_stats.h5"
        push!(L_R_separate_halves,h5read(filename, "L_R_separate_halves"))
        push!(std_in_pos,h5read(filename, "zs_std").*monomer_size)
        push!(L_R_separate_halves_error,h5read(filename, "L_R_separate_halves_error"))
        mean_extents=h5read(filename, "L_R_distance").*monomer_size
        std_extents=h5read(filename, "L_R_distance_std").*monomer_size
        num_samp=h5read(filename, "num_samp")
        push!(ds_ter, mean_extents)
        push!(std_ter, std_extents/sqrt(num_samp))

        flip_f=fetch_flip_frequency(dir,data_dir)
        push!(flip_freqs, flip_f)
    end

    labels=["No SMCs", "No ter", "WT"]

    plot_unreplicating_std_z(std_in_pos[1];label=labels[1], alpha=0.3)
    plot_unreplicating_std_z!(std_in_pos[2]; label=labels[2], alpha=0.6, linestyle=:dash)
    plot_unreplicating_std_z!(std_in_pos[3]; label=labels[3])
    savefig(plot_dir*"/std_long_axis_positions"*figtype)

    makela_data=[56.6 ,65.7, 97.8]
    makela_errors=[3.2, 0.8, 0.6]
    sim_std=sqrt.(L_R_separate_halves.*(1 .-L_R_separate_halves))
    sim_color=:black
    makela_1=:white
    makela_2=:lightgray

    scatter((1:3).-0.1, L_R_separate_halves*100, yerr=sim_std,yticks=40:10:105, ylabel="L3-R3 separate halves [%]", xticks=(1:3,labels), size=(450,400), ylims=(45,103), label="Simulations", legend=:topleft, color=sim_color)
    scatter!((1:3).+0.1, makela_data, yerr=makela_errors, label="Mäkelä et al. 2021", color=makela_1)
    hline!([50], color=:black, linestyle=:dash)
    savefig(plot_dir*"/L3_R3_sep"*figtype)

    flip_freqs./=flip_freqs[3]
    makela_data=[NaN, 0.78,0.081]
    makela_errors=[NaN, 0.023,0.013]
    #calculate the relative frequency of flipping, and do error propagation
    rel_errors=makela_errors./makela_data
    makela_data./=makela_data[3]
    #for f=a/b, relative errors squared add
    makela_errors=map(i->sqrt(rel_errors[i]^2+rel_errors[1]^2), 1:3).*makela_data
    makela_errors[1]=0

    max_y_val=max(maximum(flip_freqs),maximum(makela_data.+makela_errors))
    scatter((1:3).-0.1, flip_freqs, ylabel="Relative flip freq", xticks=(1:3,labels), size=(450,400), label="Simulations", legend=:bottomleft, color=sim_color)
    scatter!((1:3).+0.1, makela_data, yerr=makela_errors, label="Mäkelä et al. 2021", color=makela_1, ylims=[0, max_y_val+1])
    savefig(plot_dir*"/flip_freq"*figtype)

    #d L-R
    scatter((1:3).-0.15, ds_ter, yerr=std_ter, ylabel="Mean L3-R3 distance [μm]", xticks=(1:3,labels), size=(450,400), label="Simulations", legend=:topleft, ylims=(0,2.1), color=sim_color)
    #data from 2021 paper
    makela_data=[1.18, 0.82, 1.57]
    makela_errors=[0.07, 0.03, 0.12]
    scatter!((1:3), makela_data, yerr=makela_errors, label="Mäkelä et al. 2021", color=makela_1, legend=:topleft)
    #data from 2020 paper
    makela_data=[NaN, 0.69, 1.01]
    makela_errors=[NaN, 0.01, 0.03]
    scatter!((1:3).+0.15, makela_data, yerr=makela_errors, label="Mäkelä & Sherratt 2020", color=makela_2)
    savefig(plot_dir*"/L3_R3_extension"*figtype)
end

"""
Compare measures of LoR order
"""
function compare_unreplicating_stats(data_dir, base_plot_dir, filetype="png"; monomer_size=0.033, N=4600, ter_size=800, only_original_dims=false, nu=0.42)
    nu_loop=2/5
    if filetype[1]!='.'
        filetype="."*filetype
    end
    filetype="nu_$nu"*filetype
    if only_original_dims
        plot_dir=joinpath(base_plot_dir, "Compare_varying_lambda_d/")
    else
        plot_dir=joinpath(base_plot_dir, "Compare_varying_lambda_d_and_dims/")
    end
    mkpath(plot_dir)
    Ms=Vector{Int64}()
    λs=Vector{Float64}()
    ls=Vector{Float64}()
    Ls=Vector{Float64}()
    Ds=Vector{Float64}()
    N_bb=Vector{Float64}()
    std_N_bb=Vector{Float64}()
    R_g_loop=Vector{Float64}()
    L_R_separate_halves=Vector{Float64}()
    ds_ter=Vector{Float64}() #800 kb; bounds at 190 and 270
    std_ter=Vector{Float64}()
    std_theta=Vector{Float64}()
    flip_freqs=Vector{Float64}()
    std_z_ori=Vector{Float64}()

    relevant_dirs=filter(contains("_ter_size_$(ter_size)_ter_strength_100"), subdirs(data_dir))

    #no_smcs_folder=filter(contains("No_smcs"), relevant_dirs)

    for dir in relevant_dirs
        if !contains(dir, "No_smcs")
            L_0=parse_after(dir, "_L_")*monomer_size
            radius=parse_after(dir, "_r_")*monomer_size
            if (radius==0.5 && L_0==54.42240234227263*monomer_size) || !only_original_dims
                filename=dir*"/steady_state_stats.h5"
                lr_sep=h5read(filename, "L_R_separate_halves")
                loop_size=parse_after(dir, "loop_size_")
                num_smcs=parse_after(dir, "num_smcs_")
                bb_l=h5read(filename, "mean_bb_l")
                std_bb_l=h5read(filename, "std_bb_l")
                mean_extents=h5read(filename, "mean_region_extents").*monomer_size
                std_extents=h5read(filename, "std_region_extents").*monomer_size
                loop_l=h5read(filename, "mean_loop_l")
                mean_loop_extents=h5read(filename, "mean_loop_ext").*monomer_size
                std_z=h5read(filename, "/zs_std")[1]*monomer_size

                push!(R_g_loop, mean_loop_extents)
                push!(std_z_ori, std_z./L_0)
                push!(L_R_separate_halves,lr_sep)
                push!(λs, loop_size)
                push!(Ls, L_0)
                push!(Ds, 2*radius)
                push!(ls, loop_l)
                push!(ds_ter, mean_extents[3])
                push!(std_ter, std_extents[3])
                push!(Ms, num_smcs)
                push!(N_bb, bb_l)
                push!(std_N_bb, std_bb_l)

                filename=dir*"/mean_angle_traj_stats.h5"
                push!(std_theta, h5read(filename, "std_angle_mod"))

                flip_f=fetch_flip_frequency(dir,data_dir)
                push!(flip_freqs, flip_f)
            end
        end
    end

    #rescale things to dimensionless units
    C=N_bb.*monomer_size.^(1/nu).*Ds.^(1-1/nu).*(N_bb./N).^((1-1/nu)/2)./Ls
    R_g_loop.*=4 ./Ds
   
    #Make "phase diagrams"
    xs=R_g_loop
    ys=C
    plot_phase_diagram(xs,ys,ds_ter./Ls, Ms, "Ter extension/cell length")
    savefig(plot_dir*"/phase_diagram_ter_extension"*filetype)

    plot_phase_diagram(xs,ys,L_R_separate_halves, Ms, "Prob. L3 R3 separate halves")
    savefig(plot_dir*"/phase_diagram_lr_separate"*filetype)

    plot_phase_diagram(xs,ys,std_theta, Ms, "Std. θ mod π")
    savefig(plot_dir*"/phase_diagram_std_theta_mod_pi"*filetype)

    plot_phase_diagram(xs,ys,flip_freqs, Ms, "Flip frequency")
    savefig(plot_dir*"/phase_diagram_flip_frequency"*filetype)

    scatter_plot_with_colors(std_theta, std_z_ori, Ms, "Std. θ mod π","Std. rel. ori position", scale=:identity)
    savefig(plot_dir*"/std_theta_vs_std_z_ori"*filetype)

    scatter_plot_with_colors(L_R_separate_halves, std_z_ori, Ms, "Prob. L3 R3 separate halves","Std. rel. ori position", scale=:identity)
    savefig(plot_dir*"/prob_lr_separate_vs_std_z_ori"*filetype)

    #First, plot against backbone length
    for x in [N_bb, ls, R_g_loop]
        if x==N_bb
            xlab="Backbone length [kb]"
            suffix="vs_bb_length"*filetype
        elseif x==C
            xlab="Compaction factor C"
            suffix="vs_compaction_factor"*filetype
        elseif x==ls
            xlab="Mean loop size"
            suffix="_vs_loop_size"*filetype
        elseif x==R_g_loop
            xlab=L"4 R_{\mathrm{g}}^{\mathrm{loop}}/D"
            suffix="_vs_loop_extension"*filetype
        else
            error("Unknown x-vector for plots")
        end

        scatter_plot_with_colors(x,ds_ter./Ls, Ms, xlab, "Extension ter/L")
        savefig(plot_dir*"/ter_dist"*suffix)

        scatter_plot_with_colors(x,std_theta, Ms, xlab, "Std angle mod π")
        savefig(plot_dir*"/std_angle_mod_pi"*suffix)

        scatter_plot_with_colors(x,L_R_separate_halves, Ms, xlab, "Prob. L3-R3 separate")
        savefig(plot_dir*"/L3_R3_separate"*suffix)

        scatter_plot_with_colors(x,flip_freqs, Ms, xlab, "Flip frequency")
        savefig(plot_dir*"/flip_frequency"*suffix)
    end    

end

"""
Compare backbone length etc to analytical predictions
"""
function compare_LE_stats_to_predictions(data_dir, base_plot_dir, filetype="png", monomer_size=0.033, N=4600, ter_size=800, only_original_dims=true; nu_loop=2/5)
    if filetype[1]!='.'
        filetype="."*filetype
    end
    plot_dir=joinpath(base_plot_dir, "Compare_LE_stats_to_predictions/")
    mkpath(plot_dir)
    loop_sizes=Vector{Float64}()
    actual_loop_sizes=Vector{Float64}()
    nums_smcs=Vector{Int64}()
    bb_ls=Vector{Float64}()
    std_bb_ls=Vector{Float64}()
    lengths=Vector{Float64}()
    radii=Vector{Float64}()
    loop_extensions=Vector{Float64}()

    relevant_dirs=filter(contains("_ter_size_$ter_size"), subdirs(data_dir))

    for dir in relevant_dirs
        if !contains(dir, "No_smcs") && !contains(dir, "ter_strength_1.0")
            filename=dir*"/steady_state_stats.h5"
            L_0=parse_after(dir, "_L_")*monomer_size
            if contains(dir, "_r_")
                radius=parse_after(dir, "_r_")*monomer_size
            else
                radius=0.5 #in microns
            end
            if (radius==0.5 && L_0==54.42240234227263*monomer_size) || !only_original_dims
                loop_size=parse_after(dir, "loop_size_")
                num_smcs=parse_after(dir, "num_smcs_")
                bb_l=h5read(filename, "mean_bb_l")
                std_bb_l=h5read(filename, "std_bb_l")
                loop_l=h5read(filename, "mean_loop_l")
                num_samp=h5read(filename, "num_samp")

                mean_loop_extents=h5read(filename, "mean_loop_ext").*monomer_size
                push!(loop_extensions, mean_loop_extents)
                push!(loop_sizes, loop_size)
                push!(lengths, L_0)
                push!(radii, radius)
                push!(actual_loop_sizes, loop_l)
                push!(nums_smcs, num_smcs)
                push!(bb_ls, bb_l)
                push!(std_bb_ls, std_bb_l/sqrt(num_samp))
            end
        end
    end
    ds=(N-ter_size)./nums_smcs

    #show that backbone length scales as expected
    scatter_plot_with_colors(loop_sizes./ds, bb_ls, nums_smcs, L"\lambda/d", L"N_{\mathrm{bb}}", scale=:identity)
    plot!(x->(N-ter_size)*exp(-x), color=:black, linestyle=:dash, xlims=(0,3.5), label=L"N\times e^{-\lambda/d}", legend=:topright)
    savefig(plot_dir*"/backbone_lengths_vs_loop_size"*filetype)

    scatter_plot_with_colors(actual_loop_sizes./ds, bb_ls, nums_smcs, L"\ell/d", L"N_{\mathrm{bb}}", scale=:identity)
    plot!(x->(N-ter_size)*exp(-x), color=:black, linestyle=:dash, xlims=(0,3.5), label=L"N\times e^{-\ell/d}", legend=:topright)
    savefig(plot_dir*"/backbone_lengths_vs_actual_loop_size"*filetype)

    #Loop length
    scatter_plot_with_colors(loop_sizes, actual_loop_sizes,nums_smcs, L"\lambda", L"\ell", scale=:identity)
    plot!(x->x, color=:black, linestyle=:dash, label=L"\ell=\lambda", legend=:topleft, xlims=(0, 430))
    savefig(plot_dir*"/actual_loop_sizes_vs_lambda"*filetype)

    scatter_plot_with_colors(actual_loop_sizes, loop_extensions, nums_smcs, L"\ell",L"R_\mathrm{g}^{\mathrm{loop}}", fitted_power=true)
    savefig(plot_dir*"/loop_r_g_vs_loop_size_fit"*filetype)
end

"""
Compare simulations of chromosome in infinite tube to predictions
"""
function compare_infinite_tube_stats_to_predictions(data_dir, base_plot_dir, filetype="png", monomer_size=0.033, N=4600, ter_size=800, nu=0.42)
    if filetype[1]!='.'
        filetype="."*filetype
    end
    filetype="nu_$(nu)_"*filetype
    plot_dir=joinpath(base_plot_dir, "Compare_LE_stats_to_predictions/")
    mkpath(plot_dir)

    fracs_bb=Vector{Float64}()
    seg_lengths=Vector{Float64}()
    gs=Vector{Float64}()
    bb_ls=Vector{Float64}()
    N_bb_w_smcs=Vector{Float64}()
    radii=Vector{Float64}()
    extents=Vector{Float64}()
    std_extents=Vector{Float64}()
    ter_extents=Vector{Float64}()
    nucleoid_extents=Vector{Float64}()
    LR_extents=Vector{Float64}()
    lambdas=Vector{Float64}()
    nums_smcs=Vector{Float64}()
    loop_extensions=Vector{Float64}()
    loop_lengths=Vector{Float64}()

    relevant_dirs=filter(contains("_ter_size_$ter_size"), subdirs(data_dir))

    for dir in relevant_dirs
        if !contains(dir, "No_smcs")
            filename=dir*"/steady_state_stats.h5"
            if contains(dir, "_r_")
                radius=parse_after(dir, "_r_")*monomer_size
            else
                radius=0.5 #in microns
            end
            ds=h5read(filename, "ds")
            nucleoid_d=h5read(filename, "nucleoid_extension")*monomer_size
            extent_per_d=h5read(filename,"mean_extensions").*monomer_size
            std_per_d=h5read(filename,"std_extensions").*monomer_size
            num_samp=h5read(filename, "num_samp")
            std_per_d./=sqrt(num_samp)
            LR_d=h5read(filename, "LR_extension")*monomer_size

            bb_l=h5read(filename, "mean_bb_l")
            bb_l_w_smcs=h5read(filename, "N_bb_w_smcs")
            g=bb_l/(N-ter_size)
            frac_bb=bb_l/N
            mean_loop_l=h5read(filename, "mean_loop_l")
            lambda=parse_after(filename, "loop_size_")
            num_smcs=parse_after(filename, "num_smcs_")
            ter_d=h5read(filename, "ter_extension")*monomer_size
            loop_extent=h5read(filename, "mean_loop_ext")*monomer_size
            for (index, d) in enumerate(ds)
                push!(N_bb_w_smcs, bb_l_w_smcs)
                push!(extents, extent_per_d[index])
                push!(loop_extensions, loop_extent)
                push!(std_extents, std_per_d[index])
                push!(lambdas, lambda)
                push!(loop_lengths, mean_loop_l)
                push!(ter_extents, ter_d)
                push!(LR_extents, LR_d)
                push!(nucleoid_extents, nucleoid_d)
                push!(radii, radius)
                push!(fracs_bb, frac_bb)
                push!(seg_lengths, d)
                push!(gs, g)
                push!(nums_smcs, num_smcs)
                push!(bb_ls, bb_l)
            end
        end
    end
    Ds=2 .*radii
    ds=(N-ter_size)./nums_smcs
    nu_loop=2/5
    seg_lengths.*=gs

    scatter_plot_with_colors(loop_lengths, loop_extensions, nums_smcs,L"\ell\ \mathrm{[kb]}" ,L"R_\mathrm{g}^{\mathrm{loop}}", fitted_power=true)
    savefig(plot_dir*"lambda_vs_loop_rg"*filetype)

    scatter_plot_with_colors(loop_extensions, LR_extents./bb_ls, nums_smcs,L"R_\mathrm{g}^{\mathrm{loop}}", L"z_{\mathrm{bb}}/N_{\mathrm{bb}}")
    savefig(plot_dir*"h_vs_LR_extension_N_bb"*filetype)

    scatter_plot_with_colors(bb_ls, LR_extents./bb_ls, nums_smcs,L"N_{\mathrm{bb}}", L"z_{\mathrm{bb}}/N_{\mathrm{bb}}", fitted_power=true)
    savefig(plot_dir*"N_bb_vs_LR_extension_N_bb"*filetype)
    
    scatter_plot_with_colors(N_bb_w_smcs, LR_extents./N_bb_w_smcs, nums_smcs,L"N_{\mathrm{bb}}b+M_{\mathrm{base}}\ell_{\mathrm{loop}}", L"z_{\mathrm{bb}}/\ell_{\mathrm{bb}}", fitted_power=true)
    savefig(plot_dir*"N_bb_plus_smcs_vs_LR_extension_N_bb"*filetype)


    scatter(N_bb_w_smcs, nucleoid_extents, label="Chromosome", color=:black,
            xlabel="Backbone length [μm]", ylabel="Extension [μm]", ylims=(0, 8))
    scatter!(N_bb_w_smcs, LR_extents, label="Looped region", shape=:rect, color=:gray)
    scatter!(N_bb_w_smcs, ter_extents, label="End-to-end", shape=:xcross, color=:white)
    savefig(plot_dir*"compare_different_extents"*filetype)
end

"""
For all sims, plot examples of fitted cosine curves
"""
function plot_examples_angle_fits(config_dir, base_plot_dir, filetype="png"; monomer_size=0.033, sampled_ts=[2000])
    if filetype[1]!='.'
        filetype="."*filetype
    end
    plot_dir=joinpath(base_plot_dir, "Fast_stats_all/")
    for dir in subdirs(config_dir)
        ter_width=floor(Int, parse_after(dir, "ter_size_"))
        out_path=replace(dir, config_dir=>plot_dir)
        mkpath(out_path)
        files=filter(contains("blocks"), readdir(dir, join=true))
        for file in files
            h5open(file) do f
                for key in keys(f)
                    i=0
                    try
                        i=read_attribute(f[key], "block")+1
                    catch
                        println("bad index! file $f, index $key")
                        continue
                    end
                    if i in sampled_ts
                        sim_run=read_attribute(f[key], "simulation_run")+1
                        pos=read(f, "$key/pos")
                        plot_z_vs_fitted_cosine(pos, monomer_size, ter_width, false)
                        savefig(joinpath(out_path,"fitted_angle_first_order_sim_$(sim_run)_time_point_$i"*filetype))

                        plot_z_vs_fitted_cosine(pos, monomer_size, ter_width, true)
                        savefig(joinpath(out_path,"fitted_angle_second_order_sim_$(sim_run)_time_point_$i"*filetype))
                    end
                end
            end
        end
    end
end

"""
    plot_mean_rsq(ts, mean_rsq)

Plots the mean squared displacement (MSD) over time for different labels.

# Arguments
- `ts::Array`: The simulation time points.
- `mean_rsq::Array`: The mean squared displacement data.
"""
function plot_mean_rsq(ts, mean_rsq)
    scatter(ts, mean_rsq[2,:], label=L"ter", color=ter_color, markerstrokewidth=0, alpha=0.5, xlabel="Simulation time [s]", ylabel="MSD [μm²]", size=(800, 400), scale=:log10)
    scatter!(ts, mean_rsq[1,:], label=L"ori", color=ori_color, markerstrokewidth=0, alpha=0.5)
    scatter!(ts, mean_rsq[3,:], label=L"90^{\circ}", color=right_color, markerstrokewidth=0, alpha=0.5)
    scatter!(ts, mean_rsq[4,:], label=L"-90^{\circ}", color=left_color, markerstrokewidth=0, alpha=0.5)
    plot!(ts, 10^(-1.5)*ts.^0.4, color=:gray, linestyle=:dash, alpha=0.5, label="Weber et al., PRL, 2010")
    plot!(ts, 10^(-2)*ts.^0.4, color=:black, linestyle=:dot, alpha=0.5, label="Experiment mean")
    plot!(ts, 10^(-3)*ts.^0.4, color=:gray, linestyle=:dash, alpha=0.5)
end


