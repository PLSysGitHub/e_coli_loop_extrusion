"""
Given that statistics have been calculated, this file can be used to make the plots
"""

using Plots, DelimitedFiles, HDF5, LaTeXStrings, StatsBase
include("helpers.jl")

#plotting defaults
pythonplot(grid=false,label="",framestyle=:box,colorbar=true, linewidth=2,
    guidefontsize=16, tickfontsize=16,colorbar_tickfontsize=14,legend=:outerright, markersize=7,
    colorbar_titlefontsize=16, legendfontsize=14,palette=:seaborn_colorblind, size=(600,400))

ori_color=7
ter_color=:black
left_color=3
right_color=2
cell_color=10
cell_alpha=0.2


"""
Plot histograms but with lines and points indicating bin heights.
Allows easier plotting of multiple histograms on one plot
"""
function density_plot(data::Vector, nbins::Int, color=:black; kwargs...)
    # Fit the histogram
    bin_edges=range(minimum(data), stop=maximum(data), length=nbins+1)
    h = StatsBase.fit(Histogram, data, bin_edges)

    # Normalize the weights to get the density
    normalized_weights = h.weights / sum(h.weights)

    # Calculate bin midpoints
    bin_midpoints = 0.5 .* vec(h.edges[1][1:end-1] + h.edges[1][2:end])

    # Create the plot (line plot for density)
    plot(bin_midpoints, normalized_weights, color=color, alpha=0.5; kwargs...)

    # Add scatter plot with x-crosses
    scatter!(bin_midpoints, normalized_weights, marker=:x, color=color)
end

"""
Plot histograms but with lines and points indicating bin heights.
Allows easier plotting of multiple histograms on one plot
"""
function density_plot!(data::Vector, nbins::Int, color=:black; kwargs...)
    # Fit the histogram
    bin_edges=range(minimum(data), stop=maximum(data), length=nbins+1)
    h = StatsBase.fit(Histogram, data, bin_edges)

    # Normalize the weights to get the density
    normalized_weights = h.weights / sum(h.weights) 

    # Calculate bin midpoints
    bin_midpoints = 0.5 .* vec(h.edges[1][1:end-1] + h.edges[1][2:end])

    # Create the plot (line plot for density)
    plot!(bin_midpoints, normalized_weights,color=color, alpha=0.5; kwargs...)

    # Add scatter plot with x-crosses
    scatter!(bin_midpoints, normalized_weights, marker=:x, color=color)
end

"""
Helper function that returns xticks in Mb
"""
function xticks_Mb(N, ori_center=true, total_Mb=4.6)
        num_Mb=floor(Int, total_Mb)
        ter=floor(Int, total_Mb/2)
        increment=round(Int, N/num_Mb)
        tick_pos=0:increment:N

        if ori_center
                tick_labels=string.(vcat(ter:num_Mb-1, 0:ter))
        else
                tick_labels=string.(0:num_Mb)
        end
        return (tick_pos, tick_labels)
end


"""
Helper function that returns xticks in Mb
"""
function yticks_Mb(N, ori_center=true, total_Mb=4.6)
        num_Mb=floor(Int, total_Mb)
        ter=floor(Int, total_Mb/2)
        increment=round(Int, N/num_Mb)
        tick_pos=increment:increment:N

        if ori_center
                tick_labels=string.(vcat(ter+1:num_Mb-1, 0:ter))
        else
                tick_labels=string.(1:num_Mb)
        end
        return (tick_pos, tick_labels)
end

"""
For plotting "phase diagrams" using different measures for LoR order
"""
function plot_phase_diagram(xs, ys, zs,nums_smcs, zlabel, xlabel=L"4 R_{\mathrm{g}}^{\mathrm{loop}}/D", ylabel=L"C(N_{\mathrm{bb}},D,L)")
    alphas=0.6
    ms=map(i->shape_per_num(i), nums_smcs)
    scatter(xs, ys, zcolor=zs,colorbartitle=zlabel, xlabel=xlabel, ylabel=ylabel, markershape=ms, xscale=:log10, yscale=:log10, xlims=(10^-0.61, 10^.01), ylims=(10^-0.24, 10^0.38))
    return current()
end

"""
Plot the orientation angles over time
"""
function plot_angle_trajectories(angles)
    scatter(angles, xlabel="Simulation time",markersize=3, markerstrokealpha=0, ylabel="Chromosome orientation, rad", ylims=[-pi*1.1, 1.1*pi])
    hline!([0, pi, -pi], color=:black, linestyle=:dash)
end

"""
Plot the sum squared residuals of the fits to the orientation angles over time
"""
function plot_residual_trajectories(res)
    scatter(res, xlabel="Simulation time",markersize=3, markerstrokealpha=0, ylabel="Sum squared residual")
end

"""
Helper function that returns a specific shape for a given number of SMCs
"""
function shape_per_num(num)
    if num==200
        return :circle
    elseif num==100
        return :cross
    elseif num==50
        return :diamond
    elseif num==30
        return :star4
    elseif num==0
        return :x
    else
        return :diamond
    end
end

"""
For plotting curves of LoR measures against a single variable
"""
function scatter_plot_with_colors(xs, ys, nums, xlabel, ylabel; scale=:log10, yerr=NaN, xerr=NaN, α=NaN, fitted_power=false)
    if any(ys.<=0)
        plot(xscale=scale, xlabel=xlabel, ylabel=ylabel)
    else
        plot(scale=scale, xlabel=xlabel, ylabel=ylabel)
    end
    unique_nums=unique(nums)
    sort!(unique_nums)
    for num in unique_nums
        shape=shape_per_num(num)
        inds=nums.==num

        if any(!isnan, xs[inds]) && any(!isnan, ys[inds])
            clr=color_for_M(num)

            if length(yerr)==length(nums)
                ybars=yerr[inds]
            elseif length(yerr)==1
                ybars=yerr
            end

            if length(xerr)==length(nums)
                xbars=xerr[inds]
            elseif length(xerr)==1
                xbars=xerr
            else
                xbars=NaN
            end

            scatter!(xs[inds], ys[inds], shape=shape, color=clr, label=L"M=%$num", alpha=0.8, yerr=ybars, xerr=xbars)
        end
    end
    if !isnan(α)
        i=findfirst(nums.==50)
        plot!( x->(x/xs[i])^α*ys[i], color=:black, linestyle=:dash, label=L"y\sim x^{%$(round(α, sigdigits=2))}")
    end

    if fitted_power
        p, σs = fit_power_law_with_errors(xs,ys)
        plot!( x->p[1]*x^p[2], color=:gray, linestyle=:dash, label=L"y\sim x^{%$(round(p[2], sigdigits=2))\pm %$(round(σs[2], sigdigits=1))}")
    end

    return current()
end

"""
Given min and max value, return a palette with this range
"""
function scaled_colors(x, pal=:plasma)
    rng=maximum(x)-minimum(x)
    colors=map(i->palette(pal)[round(Int, i/rng*255+1)], x.-minimum(x))
    return colors
end

"""
Given a value M, return a color from a palette
"""
function color_for_M(M, pal=:terrain, M_max=201, M_min=25)
    if M==0 || M==1
        return :black
    end
    rng=M_max-M_min
    return palette(pal)[round(Int, (M-M_min)/rng*101+1)]
end

"""
Plot the mean long axis positions of loci
"""
function plot_unreplicating_av_z(av, ter_width=400; kwargs...)
        N=length(av)
        ter=ceil(Int,N/2)
        L=time_to_height(0)

        ter_r_z=av[ter:ter+ter_width]
        ter_l_z=av[ter-ter_width:ter]

        r_arm=ter+ter_width:N
        l_arm=1:ter-ter_width

        right_arm_z=av[r_arm]
        left_arm_z=av[l_arm]

        #plot the cell height as an outline
        hline([L/2], ribbon=(L, 0), color=cell_color, fillalpha=cell_alpha)
        hline!([-L/2], color=cell_color, xticks=xticks_Mb(N))
        hline!([0], color=:grey, linealpha=0.3)

        #plot means and errors for the first set of simulations
        plot!(ter_width:ter_width+length(right_arm_z)-1, right_arm_z, ylabel="Mean long axis position [μm]",xlabel="Distance to ori [Mb]", color=right_color; kwargs...)
        plot!(ter:ter+length(left_arm_z)-1, left_arm_z, color=left_color; kwargs...)
        plot!(1:ter_width+1,ter_r_z, color=ter_color; kwargs...)
        plot!(N-ter_width:N, ter_l_z, color=ter_color;kwargs...)
        scatter!([ter], [av[1]], color=ori_color)
end

"""
Plot the mean long axis positions of loci
"""
function plot_unreplicating_av_z!(av, ter_width=400; kwargs...)
        N=length(av)
        ter=ceil(Int,N/2)
        L=time_to_height(0)

        ter_r_z=av[ter:ter+ter_width]
        ter_l_z=av[ter-ter_width:ter]

        r_arm=ter+ter_width:N
        l_arm=1:ter-ter_width

        right_arm_z=av[r_arm]
        left_arm_z=av[l_arm]

        plot!(ter_width:ter_width+length(right_arm_z)-1, right_arm_z, ylabel="Mean long axis position [μm]",xlabel="Distance to ori [Mb]", color=right_color; kwargs...)
        plot!(ter:ter+length(left_arm_z)-1, left_arm_z, color=left_color; kwargs...)
        plot!(1:ter_width+1,ter_r_z, color=ter_color; kwargs...)
        plot!(N-ter_width:N, ter_l_z, color=ter_color;kwargs...)
        scatter!([ter], [av[1]], color=ori_color)
end

"""
Plot the std of locus long axis positions
"""
function plot_unreplicating_std_z(std_z, ter_width=400; kwargs...)
        N=length(std_z)
        ter=ceil(Int,N/2)
        L=time_to_height(0)

        std_z./=L

        ter_r_z=std_z[ter:ter+ter_width]
        ter_l_z=std_z[ter-ter_width:ter]

        r_arm=ter+ter_width:N
        l_arm=1:ter-ter_width

        right_arm_z=std_z[r_arm]
        left_arm_z=std_z[l_arm]

        #plot means and errors for the first set of simulations
        plot(ter_width:ter_width+length(right_arm_z)-1, right_arm_z, ylabel="Std. rel. position",xlabel="Distance to ori [Mb]", xticks=xticks_Mb(N), color=right_color, ylims=(0, 0.5), size=(600,400); kwargs...)
        plot!(ter:ter+length(left_arm_z)-1, left_arm_z, color=left_color; kwargs...)
        plot!(1:ter_width+1,ter_r_z, color=ter_color; kwargs...)
        plot!(N-ter_width:N, ter_l_z, color=ter_color; kwargs...)
        scatter!([ter], [std_z[1]], color=ori_color; kwargs...)
end


"""
Plot the std of locus long axis positions
"""
function plot_unreplicating_std_z!(std_z, ter_width=400;kwargs...)
        N=length(std_z)
        ter=ceil(Int,N/2)
        L=time_to_height(0)

        std_z./=L

        ter_r_z=std_z[ter:ter+ter_width]
        ter_l_z=std_z[ter-ter_width:ter]

        r_arm=ter+ter_width:N
        l_arm=1:ter-ter_width

        right_arm_z=std_z[r_arm]
        left_arm_z=std_z[l_arm]

        #plot means and errors for the first set of simulations
        plot!(ter_width:ter_width+length(right_arm_z)-1, right_arm_z, color=right_color; kwargs...)
        plot!(ter:ter+length(left_arm_z)-1, left_arm_z, color=left_color; kwargs...)
        plot!(1:ter_width+1,ter_r_z, color=ter_color; kwargs...)
        plot!(N-ter_width:N, ter_l_z, color=ter_color; kwargs...)
        scatter!([ter], [std_z[1]], color=ori_color)
end

"""
Given long vectors of the long axis positions of monomers as a function of simulation time, plot the long axis positions
        of the oris and the ter
"""
function plot_convergence_oris(long_axis_positions, R)
        L=R_to_height(R)
        L_max=R_to_height(404)*1.1

        hline([L/2], ribbon=(L, 0), color=:turquoise, fillalpha=0.2)
        hline!([-L/2], color=:turquoise)
        plot!(long_axis_positions[1,:], xlabel="Simulation time", ylabel="Mean long axis position [μm]", label="Ori 1", color=1)
        plot!(long_axis_positions[N+1,:], label="Ori 2", ylims=(-L_max/2,L_max/2), legend=:outerright, color=2)
        plot!(long_axis_positions[Int(N/2),:], label="Ter", xticks=:auto, color=:black)
end

"""
Plot the long axis positions of a single chromosome configuration, together with the fitted cosine curve used to define the orientation angle
"""
function plot_z_vs_fitted_cosine(pos, monomer_size=0.033, ter_width=400, second_order=false)
    zs=pos[3,:].*monomer_size
    N=length(zs)
    fitted_angles=fit_angle_to_traj(zs, second_order)

    A=fitted_angles.param[1]
    theta=fitted_angles.param[2]+3*pi/2 #merely for plotting
    residual=round(sum(fitted_angles.resid.^2), sigdigits=3)
    
    plot_unreplicating_av_z(zs, ter_width)
    ylabel!("Monomer long axis position")
    if second_order
        B=fitted_angles.param[3]
        theta_2=fitted_angles.param[4]
        plot!(1:N, x-> A*cos(x*2*pi/N+theta)+B*cos(x*4*pi/N+theta_2), color=4, linestyle=:dash,label="SSR $residual", legend=:topright)
    else
        plot!(1:N, x-> A*cos(x*2*pi/N+theta), color=4, linestyle=:dash,label="SSR $residual", legend=:topright)
    end
end
