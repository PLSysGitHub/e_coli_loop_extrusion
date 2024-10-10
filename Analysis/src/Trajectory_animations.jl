using Plots, ImageFiltering, Plots.PlotMeasures, StatsBase
include("helpers.jl")

#gr is much faster than pythonplot for animations...
gr(palette=:seaborn_colorblind, label="")

ori_color=7
ter_color=:black
left_color=3
right_color=2
cell_color=10
cell_alpha=0.2

"""
    plot_strands_unreplicated(pos, ter_size, adjust_size=false)

Plots the strands in an unreplicated state.

# Arguments
- `pos::Array`: A 3xN array containing the positions of the strands.
- `ter_size::Int`: The size of the termination region.
- `adjust_size::Bool`: Whether to adjust the size of the plot based on the z-coordinates. Default is `false`.

# Returns
- The current plot.
"""
function plot_strands_unreplicated(pos, ter_size, adjust_size=false)
    x=pos[1,:]
    y=pos[2,:]
    z=pos[3,:]
    for v in [x,y,z]
        v.-=mean(v)
    end
    if adjust_size
        h=2*maximum(abs.(z))
    else
        h=3.0 #μm
    end

    N=length(x)
    ter=N÷2
    R=ter-ter_size÷2
    L=ter+ter_size÷2

    right=1:R
    ter_region=R:L
    left=L:N

    regions=[right, ter_region, left]
    colors=[right_color, ter_color, left_color]

    plot3d(xlims=(-h/2,h/2), ylims=(-h/2,h/2),zlims=(-h/2,h/2), xlabel="", 
        framestyle=:none, ticks=false, grid=:hide, showaxis=:hide, 
        size=(600,600), margin=0mm, aspect_ratio=1)

    for (region,color) in zip(regions, colors) 
        plot3d!(z[region], x[region], y[region], color=color, line_width=1.5)
    end
    scatter3d!([z[1]], [x[1]], [y[1]], color=ori_color, markersize=8)

    return current()
end

"""
    animate_file(in_file_name::String, out_file_name::String; monomer_size=0.033, adjust_size=false)

Creates an animation from the data in the input file and saves it to the output file.

# Arguments
- `in_file_name::String`: The name of the input file containing the data.
- `out_file_name::String`: The name of the output file to save the animation.
- `monomer_size::Float64`: The size of the monomers. Default is `0.033`.
- `adjust_size::Bool`: Whether to adjust the size of the plot based on the z-coordinates. Default is `false`.
"""
function animate_file(in_file_name::String, out_file_name::String; monomer_size=0.033, adjust_size=false)
    h5open(in_file_name) do file
        ter_size=parse_after(in_file_name, "ter_size_")
        ks=parse.(Int, keys(file))
        if length(ks)>1
            sort!(ks)
            anim = @animate for time_ind in ks
                pos=read(file, "$time_ind/pos").*monomer_size
                plot_strands_unreplicated(pos,ter_size, adjust_size)
            end
            mp4(anim, out_file_name, fps=10)
        end
    end
end

"""
    snapshots_file(in_file_name::String, out_file_name::String; monomer_size=0.033, sample_every=25, file_type=".png", adjust_size=false)

Creates snapshots from the data in the input file and saves them to the output file.

# Arguments
- `in_file_name::String`: The name of the input file containing the data.
- `out_file_name::String`: The name of the output file to save the snapshots.
- `monomer_size::Float64`: The size of the monomers. Default is `0.033`.
- `sample_every::Int`: The interval at which to sample the data. Default is `25`.
- `file_type::String`: The file type for the snapshots. Default is `.png`.
- `adjust_size::Bool`: Whether to adjust the size of the plot based on the z-coordinates. Default is `false`.
"""
function snapshots_file(in_file_name::String, out_file_name::String; monomer_size=0.033, sample_every=25, file_type=".png", adjust_size=false)
    if file_type[1]!='.'
        file_type="."*file_type
    end

    h5open(in_file_name) do file
        ter_size=parse_after(in_file_name, "ter_size_")
        ks=parse.(Int, keys(file))
        sort!(ks)
        for time_ind in ks[1:sample_every:end]
            pos=read(file, "$time_ind/pos").*monomer_size

            plot_strands_unreplicated(pos,ter_size, adjust_size)
            savefig(out_file_name*"_$time_ind"*file_type)
        end
    end
end

"""
    one_animation_per_sim(data_parent_dir::String, out_dir::String; which_files=[1,20], skip_done=true, monomer_size=0.033)

Creates one animation per simulation and saves it to the output directory.

# Arguments
- `data_parent_dir::String`: The parent directory containing the simulation data.
- `out_dir::String`: The output directory to save the animations.
- `which_files::Vector{Int}`: The indices of the files to create animations for. Default is `[1,20]`.
- `skip_done::Bool`: Whether to skip files that have already been processed. Default is `true`.
- `monomer_size::Float64`: The size of the monomers. Default is `0.033`.
"""
function one_animation_per_sim(data_parent_dir::String, out_dir::String; which_files=[1,20], skip_done=true, monomer_size=0.033)
    if !isdir(out_dir)
        mkpath(out_dir)
    end
    adjust_size=contains(data_parent_dir, "infinite_tube") || contains(data_parent_dir, "no_confinement")
    for dir in subdirs(data_parent_dir) 
        for which in which_files
            out_file_name=replace(dir, data_parent_dir=>out_dir)*"_$which.mp4"
            if skip_done && isfile(out_file_name)
                @info "Found file $out_file_name. Skipping animation."
            else
                all_files=filter(contains(".h5"),readdir(dir, join=true))
                if length(all_files) >= which
                    filename=all_files[which]
                    animate_file(filename, out_file_name, monomer_size=monomer_size, adjust_size=adjust_size)
                else
                    @warn "$dir has only $(length(readdir(dir))), no animation for $which..."
                end
            end
        end
    end
end

"""
    snapshots_per_sim(data_parent_dir::String, out_parent_dir::String; skip_done=true, file_type=".png", monomer_size=0.033, which_files=[1,20])

Creates snapshots for each simulation and saves them to the output directory.

# Arguments
- `data_parent_dir::String`: The parent directory containing the simulation data.
- `out_parent_dir::String`: The output directory to save the snapshots.
- `skip_done::Bool`: Whether to skip files that have already been processed. Default is `true`.
- `file_type::String`: The file type for the snapshots. Default is `.png`.
- `monomer_size::Float64`: The size of the monomers. Default is `0.033`.
- `which_files::Vector{Int}`: The indices of the files to create snapshots for. Default is `[1,20]`.
"""
function snapshots_per_sim(data_parent_dir::String, out_parent_dir::String; skip_done=true, file_type=".png", monomer_size=0.033, which_files=[1,20])
    adjust_size=contains(data_parent_dir, "infinite_tube") || contains(data_parent_dir, "no_confinement")
    for dir in subdirs(data_parent_dir)
        out_dir=replace(dir, data_parent_dir=>out_parent_dir)
        if !isdir(out_dir)
            mkpath(out_dir)
        end
        for which in which_files
            out_file_name=joinpath(out_dir,"file_$(which)_snapshot")

            if isfile(out_file_name) && skip_done #already done
                continue
            else
                all_files=filter(contains(".h5"),readdir(dir, join=true))
                if length(all_files) >= which #there are enough files to make snapshots for requested run
                    filename=all_files[which]
                    snapshots_file(filename, out_file_name, monomer_size=monomer_size, file_type=file_type, adjust_size=adjust_size)
                else
                    @warn "$dir has only $(length(readdir(dir))), no snapshots for $which..."
                end
            end
        end
    end
end
