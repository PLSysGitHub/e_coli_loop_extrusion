include("helpers.jl")

#helper functions for fetching time point data from time trajectory simulations
using LinearAlgebra, StatsBase, Distances

"""
Return true if left and right index are in opposite cell halves

Based on the measure used by Makela et al. 2021 PNAS

Function returns a boolean (are the loci in different cell halves?)
as well as the 2D projected distance between the loci
"""
function fetch_L_R_sep(pos, L=-128, R=122)
    N=size(pos)[2]
    L_ind=floor(Int, N*((360+L)%360)/360)
    R_ind=floor(Int, N*((360+R)%360)/360)
    separation=norm(pos[:,L_ind].-pos[:,R_ind])
    return sign(pos[3,L_ind])!=sign(pos[3,R_ind]), separation
end

"""
Return MSD curves for a given position trajectory
"""
function fetch_mean_squared_dist_array(pos, linear=true)
    N=size(pos)[2]
    ds=pairwise(Euclidean(), pos, dims=2) #efficient distance matrix calculation

    if linear
        mean_sq=zeros(N-1)
        counts=zeros(N-1)
        for i in 1:N
            for j in i+1:N
                s=j-i
                mean_sq[s]+=ds[j,i]^2
                counts[s]+=1
            end
        end
        mean_sq./=counts
    else
        mean_sq=zeros(floor(Int,N/2))
        counts=zeros(floor(Int,N/2))
        for i in 1:N
            for j in i+1:N
                s=min(j-i, N-(j-i))
                mean_sq[s]+=ds[j,i]^2
                counts[s]+=1
            end
        end
        mean_sq./=counts
    end
    return mean_sq
end
