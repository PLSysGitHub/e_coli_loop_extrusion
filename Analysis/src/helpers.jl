"""
File contains functions for parsing, loading simulation files, and converting between replicated length, time, etc.

"""

using HDF5, LsqFit

periodic_ind(i,N)= i>1 ? Int((i-1)%N+1) : periodic_ind(i+N,N)

model(x,p,N) = p[1].*cos.(2*pi .* x./ N .+ p[2] .+pi/2)

model_second_order(x,p,N) = p[1].*cos.(2*pi .* x./ N .+ p[2] .+pi/2)+p[3].*cos.(4*pi .*x ./N .+p[4] .+pi/2)

"""
We fit a cosine function to the long axis positions of loci.

An angle of 0 corresponds to the initial left-ori-right organization
"""
function fit_angle_to_traj(zs, second_order=false)
    N=length(zs)
    max_zs=maximum(abs.(zs))
    
    if second_order
        model_with_fixed_N=(x,p)->model_second_order(x,p,N)
        fit_angle=LsqFit.curve_fit(model_with_fixed_N, collect(1:N), zs, [max_zs, 0.0, 0, 0], lower=[max_zs/2, -pi, 0, -pi], upper=[1.5*max_zs, pi, 50, pi])
    else
        model_with_fixed_N=(x,p)->model(x,p,N)
        fit_angle=LsqFit.curve_fit(model_with_fixed_N, collect(1:N), zs, [max_zs, 0.0], lower=[max_zs/2, -pi], upper=[1.5*max_zs, pi])
    end
    
    return fit_angle
end

@. power_law_model(x,p)=p[1]*x^p[2]

"""
For studying diffusive behavior; fit power law to a given (time) trajectory
"""
function fit_power_law_with_errors(x,y)
    p0=[1. ,0.]
    lb=[0.,-100.]
    ub=[Inf, 100.]

    fitted_power_law=LsqFit.curve_fit(power_law_model, x,y,p0, lower=lb, upper=ub)

    sigma=margin_error(fitted_power_law, 0.05)

    return coef(fitted_power_law), sigma
end

"""
Write data to a h5 file with compression
"""    
function write_with_compression(file::HDF5.File, dataset_name::String, data, group_name=""; compression_level=4)
    group_name=string(group_name)
    chunk_size=size(data)
    if !isempty(group_name)
        if !haskey(file, group_name)
            create_group(file, group_name)
        end
        # Create the dataset with compression
        write(file, "$group_name/$dataset_name", data, chunk=chunk_size, compress=compression_level)
    else
        write(file, dataset_name, data, chunk=chunk_size, compress=compression_level)
    end
end

"""
Helper function for parsing strings; find first occurency of one of set of characters. If not found, return end of string
"""
function find_first_or_end(s::String, chars)
    # Find the first occurrence of any of the characters
    index = findfirst(chars, s)
    
    # If no characters are found, return the length of the string + 1
    return index !== nothing ? index[1] : length(s) + 1
end

"""
Given a string like "..._L_21_...", return 21 if what="L_"
"""
function parse_after(dir_name, what="_L_")
    from=findfirst(what, dir_name)[end]+1
    last_chars=r"(/|\.h5|_)"

    to=from+find_first_or_end(dir_name[from:end], last_chars)-2
    try
        return parse(Int,dir_name[from:to])
    catch
        try
            return parse(Float64, dir_name[from:to])
        catch
            @error "Could not parse $what from $(dir_name)"
        end
    end
end

# Helper function for getting length at given time
time_to_height(t,L_0=(1600*exp(10*log(2)/60))/1000, growth_exponent=log(2)/60)=L_0*exp(growth_exponent*t)

"""
Given a directory name, return list of subdirectories (full path)
"""
function subdirs(dir::String)
    sdirs=readdir(dir, join=true)
    return sdirs[isdir.(sdirs)]
end

"""
Given a vector, return smoothened curve (rolling average)
"""
function smoothened_curve(x, average_over=5)
    N=length(x)                                                                 
    smoothened=map(i->mean(x[max(1,i-average_over):min(N,i+average_over)]), 1:N) 
    return smoothened[average_over:N-average_over]
end
