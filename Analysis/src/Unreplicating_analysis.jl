true
include("Statistics_functions.jl")
using DelimitedFiles, LsqFit, ImageFiltering, ProgressMeter

"""
Calculate frequency of LR flips in the trajectories. (L3 and R3 locus switching positions)
"""
function calc_lr_flips(config_dir, out_dir; sample_from=1000, L=-128, R=122)
    p=Progress(length(subdirs(config_dir)), "Calculating number LR flips trajectories")
    out_file_name=joinpath(out_dir, "flip_frequencies.h5")
    if isfile(out_file_name)
        print("Found previous file $out_file_name. overwriting.")
        rm(out_file_name)
    end
    
    for dir in subdirs(config_dir)
        config_files = filter(x -> contains(x, "blocks"), readdir(dir, join=true))
        N=4600-(800-parse_after(dir, "ter_size_"))
        L_ind=floor(Int, N*((360+L)%360)/360)
        R_ind=floor(Int, N*((360+R)%360)/360)

        trajectory_data = Dict{Int, Vector{Vector{Int}}}() #sim id to vector of LR orientations

        for f in config_files
            if contains(f, "sim_")
                sim_run = parse_after(f, "sim_")
            else
                sim_run=0
            end
            h5open(f) do file
                for i in keys(file)
                    time_p=0
                    try
                        pos = read(file, "$i/pos")
                        time_p = read_attribute(file[i], "block")
                    catch
                        @warn "Corrupt index! file $f, index $i"
                        continue
                    end
                    
                    if time_p > sample_from
                        pos = read(file, "$i/pos")[3,[R_ind, L_ind]]
                        L_before_R=pos[2]<pos[1]
                        L_before_R=-1+2*L_before_R #-1 vs +1
                        time = read_attribute(file[i], "block")
                        sim = read_attribute(file[i], "simulation_run") + 5 * sim_run
                        if !haskey(trajectory_data, sim)
                            trajectory_data[sim] = []
                        end
                        push!(trajectory_data[sim],vcat(time,L_before_R) )
                    end
                end
            end
        end

        if !isempty(trajectory_data)
            nums_flips=Vector{Int}()
            lengths_traj=Vector{Int}()
            for (sim, data) in trajectory_data
                sorted_data = sort(data, by = x -> x[1])  # Sort by time
                LR_traj=getindex.(sorted_data, 2)
                num_flips=0
                prev_lr=LR_traj[1]
                for i in 2:length(LR_traj)
                    new_lr=LR_traj[i]
                    if new_lr!=prev_lr
                        num_flips+=1
                    end
                    prev_lr=new_lr
                end
                push!(nums_flips, num_flips)
                push!(lengths_traj, length(LR_traj))
            end

            mean_freq=sum(nums_flips)/sum(lengths_traj)
            h5open(out_file_name, "cw") do f
                g=create_group(f,basename(dir))
                g["mean_freq"]=mean_freq
                g["nums_flips"]=nums_flips
                g["lengths_traj"]=lengths_traj
            end
        end
    next!(p)
    end
    finish!(p)
end
"""
Calculate the angles as a function of time for all simulations
"""
function calc_angle_trajectories(config_dir, out_dir; skip_done=false, sample_from=0, second_order=false, monomer_size=0.033)
    Threads.@threads for dir in subdirs(config_dir)
        out_path = replace(dir, config_dir => out_dir)
        mkpath(out_path)

        out_file_name = joinpath(out_path, second_order ? "angle_trajectories_second_order.h5" : "angle_trajectories.h5")

        config_files = filter(x -> contains(x, "blocks"), readdir(dir, join=true))

        if !isfile(out_file_name) || !skip_done
            trajectory_data = Dict{Int, Vector{Vector{Float64}}}() #sim id to vector of t,angle, res etc

            for f in config_files
                if contains(f, "sim_")
                    sim_run = parse_after(f, "sim_")
                else
                    sim_run=0
                end
                h5open(f) do file
                    for i in keys(file)
                        time_p=0
                        try
                            pos = read(file, "$i/pos")
                            time_p = read_attribute(file[i], "block")
                        catch
                            @warn "Corrupt index! file $f, index $i"
                            continue
                        end
                        
                        if time_p > sample_from
                            zs = read(file, "$i/pos")[3,:].*monomer_size
                            fit_for_angle = fit_angle_to_traj(zs, second_order)
                            params=fit_for_angle.param
                            residual = sum(fit_for_angle.resid.^2)
                            time = read_attribute(file[i], "block")
                            sim = read_attribute(file[i], "simulation_run") + 10 * sim_run

                            if !haskey(trajectory_data, sim)
                                trajectory_data[sim] = []
                            end
                            push!(trajectory_data[sim],vcat(time, params, residual) )
                        end
                    end
                end
            end

            if !isempty(trajectory_data)
                h5open(out_file_name, "w") do file
                    for (sim, data) in trajectory_data
                        sorted_data = sort(data, by = x -> x[1])  # Sort by time
                        times = getindex.(sorted_data, 1)
                        if second_order
                            amplitudes=getindex.(sorted_data, 2)
                            angles=getindex.(sorted_data, 3)
                            amps_2=getindex.(sorted_data, 4)
                            angles_2=getindex.(sorted_data, 5)
                            residuals=getindex.(sorted_data, 6)

                            amplitudes=hcat(amplitudes, amps_2)
                            angles=hcat(angles, angles_2)
                            θ_0 = angles[1,1]
                        else
                            amplitudes=getindex.(sorted_data, 2)
                            angles=getindex.(sorted_data, 3)
                            residuals = getindex.(sorted_data, 4)
                            θ_0 = angles[1]
                        end
                        if abs(θ_0) > pi/3
                            angles .+= pi
                            angles[angles .> pi] .-= 2 * pi
                        end
                        
                        write_with_compression(file, "angles", angles, sim)
                        write_with_compression(file, "residuals", residuals, sim)
                        write_with_compression(file, "amplitudes", amplitudes, sim)
                    end
                end
            end
        end
    end
end

"""
Discretize the angle trajectory into three states; lor, rol and in between
"""
function three_state_angle_trajectory(trajectory, smoothening_width=20)
    smoothened_traj=smoothened_curve(abs.(trajectory), smoothening_width) #between 0 and pi
    smoothened_traj[smoothened_traj.<pi/3].=-1 #lor; angle close to 0
    smoothened_traj[smoothened_traj.>2*pi/3].=0 #rol angle close to pi
    smoothened_traj[smoothened_traj.>0].=1 #in between
    return smoothened_traj
end

"""
Given an angle trajectory, calculate stats like fraction of time in each state and lifetime of initial state
"""
function angle_traj_stats(trajectory, smoothening_width=20, sample_from=1000)
    three_state_traj=three_state_angle_trajectory(trajectory, smoothening_width)
    three_state_traj_converged=three_state_traj[sample_from:end]
    lor=sum(three_state_traj_converged.==-1)
    rol=sum(three_state_traj_converged.==0)
    rotating=sum(three_state_traj_converged.==1)
    fractions=[lor,rol,rotating]./length(three_state_traj_converged)
    first_lor=findfirst(three_state_traj.==-1) #probably 1
    if isnothing(first_lor)
        return fractions, 0
    end
    first_rot=findfirst(three_state_traj[first_lor:end].>-1)
    lifetime=length(trajectory)
    if !isnothing(first_rot)
        lifetime=first_rot
    end
    return fractions, lifetime
end

"""
For all directories, calculate the mean angle trajectory stats
"""
function calc_mean_angle_traj_stats(stats_folder, smoothening_width=20, sample_from=500, max_res=500)
    for f in subdirs(stats_folder)
        traj_file=f*"/angle_trajectories.h5"
        mean_stat_file=f*"/mean_angle_traj_stats.h5"
        mean_fracs=zeros(3)
        mean_lifetime=0
        std_lifetime=0
        std_fracs=zeros(3)
        mean_angle=0
        std_angle=0
        mean_angle_mod=0
        std_angle_mod=0
        mean_ssr=0
        std_ssr=0

        num_traj=0
        num_samp=0
        h5open(traj_file) do infile
            num_traj=length(keys(infile))
            for k in keys(infile)
                traj=read(infile[k], "angles")
                if length(traj)>sample_from
                    fracs, lt=angle_traj_stats(traj, smoothening_width, sample_from)
                    mean_lifetime+=lt
                    mean_fracs.+=fracs
                    std_lifetime+=lt^2
                    std_fracs.+=fracs.^2
                    # for the stats, only beyond sample_from and with good residual
                    rss=read(infile[k], "residuals")
                    rss=rss[sample_from:end]
                    mean_ssr+=mean(rss)
                    std_ssr+=mean(rss)^2

                    traj=traj[sample_from:end]
                    traj=traj[rss.<max_res]
                    mean_angle+=sum(traj)
                    std_angle+=sum(traj.^2)
                    #now modulo pi. better std measure
                    traj[traj.>pi/2].-=pi
                    traj[traj.<-pi/2].+=pi
                    mean_angle_mod+=sum(traj)
                    std_angle_mod+=sum(traj.^2)
                    num_samp+=length(traj)
                end
            end
        end
        
        if num_samp>1
            mean_fracs./=num_traj
            std_fracs./=num_traj
            std_fracs=sqrt.(std_fracs.-mean_fracs.^2)
            sem_fracs=std_fracs./sqrt(num_traj)

            mean_lifetime/=num_traj
            std_lifetime/=num_traj
            std_lifetime=sqrt(std_lifetime-mean_lifetime^2)
            sem_lifetime=std_lifetime/sqrt(num_traj)

            mean_angle/=num_samp
            std_angle/=num_samp
            std_angle=sqrt(std_angle-mean_angle^2)

            mean_angle_mod/=num_samp
            std_angle_mod/=num_samp
            std_angle_mod=sqrt(std_angle_mod-mean_angle_mod^2)
                

            mean_ssr/=num_traj
            std_ssr/=num_traj
            std_ssr=sqrt(std_ssr-mean_ssr^2)

            h5open(mean_stat_file, "w") do outfile
                @write outfile mean_lifetime
                @write outfile std_lifetime
                @write outfile mean_fracs
                @write outfile std_fracs
                @write outfile sem_lifetime
                @write outfile sem_fracs
                @write outfile mean_angle
                @write outfile std_angle
                @write outfile mean_angle_mod
                @write outfile std_angle_mod
                @write outfile mean_ssr
                @write outfile std_ssr
            end
        else
            print("Got $num_samp samples for $mean_stat_file")
        end
    end
end

"""
Extract SMC loop lengths
"""
function extract_loop_lengths(smcs, ter_size=800, N_non_ter=3800)
    N=N_non_ter+ter_size
    in_bb=trues(N) #also considers ter region
    ter=floor(Int,N/2)
    ter_width=Int(ter_size/2)
    loop_lengths=[]
    for smc in eachcol(smcs)
        inds=copy(smc)
        sort!(inds)
        if norm(ter-inds[1])<ter_width && norm(ter-inds[2])<ter_width
            #both in ter region.
            continue
        elseif inds[1]<ter && inds[2]>ter #on different arms; loop across ori
            loop_l=inds[1]+N-inds[2]
            if loop_l<N_non_ter
                push!(loop_lengths,loop_l)
            end
        else
            push!(loop_lengths, inds[2]-inds[1])
        end
    end
    return loop_lengths 
end

"""
Calculate the radius of gyration of a given set of monomers
"""
function radius_of_gyration(pos)
    com=mean(pos, dims=2)
    norm_pos=pos.-com
    return sqrt(sum(norm_pos.^2)/size(pos)[2])
end

"""
Calculate mean radius of gyration of loops
"""
function loop_rgs(smcs, pos, ter_size=800, N_non_ter=3800)
    N=N_non_ter+ter_size
    ter=floor(Int,N/2)
    ter_width=Int(ter_size/2)
    loop_sizes=[]
    for smc in eachcol(smcs)
        inds=copy(smc)
        sort!(inds)
        if norm(ter-inds[1])<ter_width && norm(ter-inds[2])<ter_width
            #both in ter region.
            continue
        elseif inds[1]<ter && inds[2]>ter #on different arms; loop across ori
            inside_loop=vcat(1:inds[1]-1,inds[2]+1:N)
        else
            inside_loop=inds[1]+1:inds[2]-1
        end
        if length(inside_loop)>5
            loop_monomers=pos[:,inside_loop]
            loop_rg=radius_of_gyration(loop_monomers)
            push!(loop_sizes, loop_rg)
        end
    end
    return loop_sizes
end

"""
Extract the backbone, ie part of looped region not in loops
"""
function extract_backbone(smcs, ter_size=800, N_non_ter=3800)
    N=N_non_ter+ter_size
    in_bb=trues(N) #also considers ter region
    ter=floor(Int,N/2)
    ter_width=floor(Int, ter_size/2)
    for smc in eachcol(smcs)
        inds=copy(smc)
        sort!(inds)
        if norm(ter-inds[1])<ter_width && norm(ter-inds[2])<ter_width
            #both in ter region.
            continue
        elseif inds[1]<ter && inds[2]>ter #on different arms; loop across ori
            in_bb[1:inds[1]-1].=false
            in_bb[inds[2]+1:N].=false
        else
            in_bb[inds[1]+1:inds[2]-1].=false
        end
    end
    circshift!(in_bb, -ter-ter_width)
    return in_bb[1:N-ter_size]
end

"""
Function that calculates contour length of backbone, including contributions from loop-extruders.
"""
function true_backbone_length(bb, monomer_size=0.033, condensin_size=0.05)
    if all(bb)
        return length(bb)
    end

    spacings, loopsizes=backbone_spacings_and_loopsizes(bb)
    smc_contribution=length(loopsizes)*condensin_size

    btwn_smcs=bb[findfirst(.!bb):findlast(.!bb)]
    if any(btwn_smcs)
        btwn_smcs=btwn_smcs[findfirst(btwn_smcs):findlast(btwn_smcs)]
        return sum(btwn_smcs)*monomer_size+smc_contribution
    else
        return smc_contribution
    end
end

"""
Calculate the spacings between loops and loop sizes
"""
function backbone_spacings_and_loopsizes(bb)
    spacings=[]
    loopsizes=[]
    L=length(bb)
    prev_type=bb[1]
    current_l=1
    for i in 2:L
        new_type=bb[i]
        if new_type!=prev_type || i==L
            if prev_type
                push!(spacings, current_l)
            else
                push!(loopsizes, current_l)
            end
            prev_type=new_type
            current_l=1
        else
            current_l+=1
        end
    end

    return spacings,loopsizes
end

"""
Calculate end-to-end distances of regions in the polymer
"""
function calc_region_distances(pos, ter_size=800)
    N=size(pos)[2]
    #if there's no ter region, we consider the WT value
    region_size = ter_size==0 ? 800 : ter_size
    #we take four regions; ori, ter, left, right.
    centers=floor.(Int, [1, 0.25, 0.5, 0.75].*N)
    locus_pairs=map(center-> periodic_ind.([center-region_size/2, center+region_size/2], N), centers)

    ds=map(p-> abs(pos[3,p[1]]-pos[3,p[2]]), locus_pairs)
    return ds
end

"""
For simulations in infinite tube, calculate statistics and save to .h5 file
"""
function calc_infinite_tube_stats(config_dir,out_dir; N_non_ter=3800, skip_done=true, sample_from=1000, ds=[1000, 1500, 2000, 2500])
    println("Calculating fast stats for unconfined polymers")
    for dir in subdirs(config_dir)
        out_path=replace(dir, config_dir=>out_dir)
        if !isfile(out_path*"/steady_state_stats.h5") || !skip_done
            ter_size=parse_after(dir, "ter_size_")
            N=N_non_ter+ter_size
            looped_region=vcat(1:Int((N-ter_size)/2), Int((N+ter_size)/2):N)
            has_smcs=!contains(dir, "No_smcs_")

            #arrays for stats
            num_samp=0
            mean_bb_l=0
            std_bb_l=0
            N_bb_w_smcs=0
            mean_loop_l=0
            mean_loop_ext=0
            std_loop_l=0
            num_loops=0
            mean_extensions=zeros(length(ds))
            std_extensions=zeros(length(ds))
            ter_extension=0
            LR_extension=0
            nucleoid_extension=0
            nucleoid_std=0
            ter_std=0
            LR_std=0

            for f in filter(contains("blocks"), readdir(dir, join=true))
                h5open(f) do file
                    for i in keys(file)
                        time_p=0
                        try
                            pos=read(file, "$i/pos")
                            if has_smcs
                                smcs=read(file, "$i/SMCs").+1
                            end
                            time_p=read_attribute(file[i], "block")
                        catch
                            @warn "Corrupt index! file $f, index $i"
                            continue
                        end
                        if time_p>sample_from
                            pos=read(file, "$i/pos")
                            zs=pos[3,:]
                            com_z=mean(zs)
                            zs.-=com_z
                            if has_smcs
                                smcs=read(file, "$i/SMCs").+1
                                bb=extract_backbone(smcs, ter_size, N_non_ter)
                                loop_lengths=extract_loop_lengths(smcs,ter_size, N_non_ter)
                                loop_ext=loop_rgs(smcs, pos, ter_size, N_non_ter)
                            else
                                bb=trues(N_non_ter)
                                loop_lengths=[0]
                                loop_ext=0
                            end
                            spacings, lls=backbone_spacings_and_loopsizes(bb)
                            N_bb=sum(bb)
                            N_bb_w_smcs+=true_backbone_length(bb)
                            mean_bb_l+=N_bb
                            std_bb_l+=N_bb^2
                            mean_loop_l+=sum(loop_lengths)
                            std_loop_l+=sum(loop_lengths.^2)
                            num_loops+=length(loop_lengths)
                            mean_loop_ext+=mean(loop_ext)
                            ter_d=abs(zs[Int(N/2-ter_size/2)]-zs[Int(N/2+ter_size/2)])
                            ter_d=abs(zs[Int(N/2-ter_size/2)]-zs[Int(N/2+ter_size/2)])
                            nucleoid_d=maximum(zs)-minimum(zs)
                            LR_d=maximum(zs[looped_region])-minimum(zs[looped_region])

                            ter_extension+=ter_d
                            ter_std+=ter_d^2
                            LR_extension+=LR_d
                            LR_std+=LR_d^2
                            nucleoid_extension+=nucleoid_d
                            nucleoid_std+=nucleoid_d^2

                            for (index, d) in enumerate(ds)
                                d_half=floor(Int, d/2)
                                extension=abs(zs[end-d_half]-zs[d_half])
                                mean_extensions[index]+=extension
                                std_extensions[index]+=extension^2
                            end

                            num_samp+=1
                        end
                    end
                end
            end

            if num_samp>0
                mean_loop_ext/=num_samp
                mean_bb_l/=num_samp
                N_bb_w_smcs/=num_samp
                std_bb_l/=num_samp
                mean_loop_l/=num_loops
                std_loop_l/=num_loops
                num_loops/=num_samp
                std_bb_l=sqrt(std_bb_l-mean_bb_l^2)
                std_loop_l=sqrt(std_loop_l-mean_loop_l^2)
                mean_extensions./=num_samp
                ter_extension/=num_samp
                LR_extension/=num_samp
                nucleoid_extension/=num_samp
                std_extensions./=num_samp
                std_extensions=sqrt.(std_extensions.-mean_extensions.^2)
                #save to file
                if !ispath(out_path)
                    mkpath(out_path)
                end
                h5open(out_path*"/steady_state_stats.h5", "w") do file
                    #single 
                    @write file num_samp
                    @write file N_bb_w_smcs
                    @write file ds
                    @write file mean_extensions
                    @write file std_extensions
                    @write file mean_bb_l
                    @write file std_bb_l
                    @write file mean_loop_l
                    @write file std_loop_l
                    @write file mean_loop_ext
                    @write file num_loops
                    @write file ter_extension
                    @write file LR_extension
                    @write file nucleoid_extension
                end
            end
        end
    end
end


"""
Given a directory of configurations performed at a fixed replication stage, calculate...

- mean long axis positions

averaged over configurations from simulation step sample_from onwards
"""
function calc_fast_unreplicating_stats(config_dir,out_dir="Stats/Equilibrium_stats/"; N_non_ter=3800, skip_done=true, sample_from=0)
    println("Calculating fast stats")
    mkpath(out_dir)
    Threads.@threads for dir in subdirs(config_dir)
        out_path=replace(dir, config_dir=>out_dir)
        if !isfile(out_path*"/steady_state_stats.h5") || !isfile(out_path*"/steady_state_distance_distributions.h5") || !skip_done
            ter_size=parse_after(dir, "ter_size_")
            N=N_non_ter+ter_size
            has_smcs=!contains(dir, "No_smcs_")

            #arrays for stats
            num_samp=0
            mean_zs=zeros(N)
            zs_std=zeros(N)
            ori_extents=Vector{Float64}()
            ori_pos=Vector{Float64}()
            ori_plus_pos=Vector{Float64}()
            ter_extents=Vector{Float64}()
            ter_pos=Vector{Float64}()
            ter_plus_pos=Vector{Float64}()
            left_extents=Vector{Float64}()
            right_extents=Vector{Float64}()
            L_R_separate_halves=0
            L_R_distance=0
            L_R_distance_std=0
            mean_bb_l=0
            std_bb_l=0
            mean_loop_l=0
            mean_loop_ext=0
            std_loop_l=0
            num_loops=0

            for f in filter(contains("blocks"), readdir(dir, join=true))
                h5open(f) do file
                    for i in keys(file)
                        time_p=0
                        try
                            pos=read(file, "$i/pos")
                            if has_smcs
                                smcs=read(file, "$i/SMCs").+1
                            end
                            time_p=read_attribute(file[i], "block")
                        catch
                            @warn "Corrupt index! file $f, index $i"
                            continue
                        end
                        if time_p>sample_from
                            pos=read(file, "$i/pos")
                            if has_smcs
                                smcs=read(file, "$i/SMCs").+1
                                bb=extract_backbone(smcs, ter_size, N_non_ter)
                                loop_lengths=extract_loop_lengths(smcs,ter_size, N_non_ter)
                                loop_ext=loop_rgs(smcs, pos, ter_size, N_non_ter)
                            else
                                bb=trues(N_non_ter)
                                loop_lengths=[0]
                                loop_ext=0
                            end
                            mean_bb_l+=sum(bb)
                            std_bb_l+=sum(bb)^2
                            mean_loop_l+=sum(loop_lengths)
                            std_loop_l+=sum(loop_lengths.^2)
                            num_loops+=length(loop_lengths)
                            
                            mean_loop_ext+=mean(loop_ext)

                            z=pos[3,:]
                            mean_zs.+=z
                            zs_std.+=z.^2

                            LR_is_sep, LR_sep =fetch_L_R_sep(pos)
                            L_R_separate_halves+=LR_is_sep
                            L_R_distance+=LR_sep
                            L_R_distance_std+=LR_sep^2

                            region_ds=calc_region_distances(pos, ter_size)

                            push!(ori_extents, region_ds[1])
                            push!(left_extents, region_ds[2])
                            push!(ter_extents, region_ds[3])
                            push!(right_extents, region_ds[4])
                            push!(ori_pos, z[1])
                            push!(ori_plus_pos, z[1+500])
                            push!(ter_pos, z[Int(length(z)/2)])
                            push!(ter_plus_pos, z[Int(length(z)/2)+200])

                            num_samp+=1
                        end
                    end
                end
            end

            if num_samp>0
                #divide by number of samples to get average
                for array in [mean_zs, zs_std]
                    array ./= num_samp
                end

                mean_loop_ext/=num_samp
                mean_bb_l/=num_samp
                std_bb_l/=num_samp
                mean_loop_l/=num_loops
                std_loop_l/=num_loops
                L_R_separate_halves/=num_samp
                L_R_distance_std/=num_samp
                L_R_distance/=num_samp
                #for binary variable, the error can be estimated:
                L_R_separate_halves_error=sqrt(L_R_separate_halves*(1-L_R_separate_halves)/num_samp)
                #calculate stds
                zs_std.=sqrt.(zs_std.-mean_zs.^2)
                std_bb_l=sqrt(std_bb_l-mean_bb_l^2)
                mean_region_extents=[mean(ori_extents), mean(left_extents), mean(ter_extents), mean(right_extents)]
                std_region_extents=[std(ori_extents), std(left_extents), std(ter_extents), std(right_extents)]
                std_loop_l=sqrt(std_loop_l-mean_loop_l^2)
                L_R_distance_std=sqrt(L_R_distance_std-L_R_distance^2)
                #save to file
                if !ispath(out_path)
                    mkpath(out_path)
                end
                h5open(out_path*"/steady_state_stats.h5", "w") do file
                    #first the arrays in compressed format
                    write_with_compression(file, "mean_zs", mean_zs)
                    write_with_compression(file, "zs_std", zs_std)
                    #single 
                    @write file num_samp
                    @write file L_R_separate_halves
                    @write file L_R_separate_halves_error
                    @write file L_R_distance
                    @write file L_R_distance_std
                    @write file mean_region_extents
                    @write file std_region_extents
                    @write file mean_bb_l
                    @write file std_bb_l
                    @write file mean_loop_l
                    @write file std_loop_l
                    @write file mean_loop_ext
                    @write file num_loops
                end

                #also the distance distributions for making some histograms
                h5open(out_path*"/steady_state_distance_distributions.h5", "w") do file
                    write_with_compression(file, "ori_extents", ori_extents)
                    write_with_compression(file, "left_extents", left_extents)
                    write_with_compression(file, "ter_extents", ter_extents)
                    write_with_compression(file, "right_extents", right_extents)
                end
                h5open(out_path*"/steady_state_position_distributions.h5", "w") do file
                    write_with_compression(file, "ori_pos", ori_pos)
                    write_with_compression(file, "ter_pos", ter_pos)
                    write_with_compression(file, "ori_plus_pos", ori_plus_pos)
                    write_with_compression(file, "ter_plus_pos", ter_plus_pos)
                end
            end
        end
    end
end

"""
Given a directory of configurations, calculate...

- mean long axis positions
- mean region extensions

as a function of simulation step
"""
function check_unreplicating_convergence(config_dir="Simulation_data/Steady_state/", out_dir="Stats/Convergence/"; max_it=2000, spacer=1, orient_arms=false)
    folders=subdirs(config_dir)
    mkpath(out_dir)
    
    Threads.@threads for folder in folders
        ter_size=parse_after(folder, "ter_size_")
        N=3800+ter_size
        files=filter(contains("blocks"), readdir(folder, join=true))
        no_orient=!(contains(folder, "no_tether")) #if tethered, dont relabel oris
        ter_distinct=!(contains(folder, "No_smcs")) #if smcs, ter distinct from ori distinct even if R=N/2

        zs=zeros(N,max_it)
        counts=zeros(max_it)
        region_extensions=zeros(4, max_it)
        
        #loop over all files and increment arrays
        for file in files
            R=0#parse_after(file)
            h5open(file) do f
                for key in keys(f)
                    i=0
                    try
                        i=read_attribute(f[key], "block")+1
                        pos=read(f, "$key/pos")
                    catch
                        @warn "Corrupt index! file $file, index $key"
                        continue
                    end
                    if i<=max_it
                        pos=read(f, "$key/pos")
                        if orient_arms && mean(pos[3,1:Int(N/2)])<mean(pos[3,Int(N/2):N])
                            #flip so that orientation is the same for averaged simulations
                            pos[3,:].*=-1
                        end

                        zs[:,i].+=pos[3,:]
                        counts[i]+=1
                        region_extensions[:,i].+=calc_region_distances(pos, ter_size)
                    end
                end
            end
        end

        #divide by counts to get averages
        for r in eachrow(zs)
            r./=counts
        end
        for r in eachrow(region_extensions)
            r./=counts
        end
        
        #save the data
        out_path=replace(folder, config_dir=>out_dir)
        if !isdir(dirname(out_path))
            mkpath(dirname(out_path))
        end
        h5open(out_path*"_convergence.h5", "w") do file
            write_with_compression(file, "zs", zs)
            write_with_compression(file, "region_extensions", region_extensions)
        end
    end
end

"""
Given a directory of simulation results, calculate how many simulations were run for each set of conditions
"""
function check_num_samples_steady_state(config_dir="../Results_3D/Steady_state/", out_dir="Stats/", sample_from=1000)
    folders=subdirs(config_dir)
    num_samples=zeros(Int, length(folders))
    traj_lengths=zeros(Int, length(folders))
    for (index,folder) in enumerate(folders)
        files=filter(contains("blocks"), readdir(folder, join=true))
        #loop over all files and increment arrays
        for file in files
            h5open(file) do f
                for key in keys(f)
                    i=0
                    try
                        i=read_attribute(f[key], "block")+1
                    catch
                        @warn "Corrupt index! file $f, index $key"
                        continue
                    end
                    if i==sample_from+1
                        num_samples[index]+=1
                    end
                    traj_lengths[index]=max(traj_lengths[index], i)
                end
            end
        end
    end
    mkpath(out_dir)
    writedlm(out_dir*"num_steady_state_simulations.txt", hcat(replace.(folders, config_dir=>""), num_samples, traj_lengths))
end

"""
Given set of files, calculate mean squared distance over time
"""
function fetch_msd_trajectories(initial_stage_files, tracked_loci=[1, 2300, 1150, 3450])
    #Make arrays and set parameters
    r_squared=Dict{Int, Array{Float64}}() #for each sim, store an array x,y,t
    n_loci=length(tracked_loci)

    for (index, f) in enumerate(initial_stage_files)
        h5open(f) do file
            for i in keys(file)
                block=read_attribute(file[i], "block")+1
                sim = read_attribute(file[i], "simulation_run")
                pos=read(file[i], "pos")[2:3,tracked_loci].*monomer_size #2D diffusion
                if haskey(r_squared, sim)
                    r_squared[sim][:,:,block].=pos
                else
                    r_squared[sim]=zeros(2,n_loci, nt)
                    r_squared[sim][:, :, block].=pos
                end
            end
        end
    end
    #For each simulation, calculate r^2(t)
    for sim in keys(r_squared)
        all_pos=r_squared[sim]
        all_pos.-=all_pos[:,:,1] #subtract first time point
        all_pos.=all_pos.^2
        r_squared[sim]=(all_pos[1,:,:].+all_pos[2,:,:]) #sum y and z
    end

    mean_rsq=zeros(n_loci, nt)
    error_rsq=zeros(n_loci, nt)
    for sim in keys(r_squared)
        mean_rsq.+=r_squared[sim]
        error_rsq.+=r_squared[sim].^2
    end
    mean_rsq./=length(keys(r_squared))
    error_rsq./=length(keys(r_squared))
    error_rsq.=sqrt.(error_rsq.-mean_rsq.^2)
    error_rsq./=sqrt(length(keys(r_squared)))

    return mean_rsq, error_rsq
end
