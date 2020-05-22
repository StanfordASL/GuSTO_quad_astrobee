# -----------------------------------------------------------------------------------------------
# -   Initializing script for the freeflyer example - Thomas Lew and Riccardo Bonalli 12/2019   -
# -----------------------------------------------------------------------------------------------



using LinearAlgebra
using Ipopt
using JuMP
using DifferentialEquations
using NLsolve
using Plots
# using PyPlot
import PyPlot; const plt = PyPlot

include("./Models/astrobee_se3.jl")
include("./SCP/gusto_problem.jl")



# Plotting Functions



# 2D plots

function plot_solutions(scp_problem::GuSTOProblem, model, X_all, U_all)
    N = length(X_all)

    idx = [1,2]
    local fig
    fig = plot(X_all[1][idx[1],:], X_all[1][idx[2],:],
        label="iter = 0",linewidth=2,
        xlims=(-0.5,6.),ylims=(-0.5,6.5),
        xlabel="x Position",ylabel="y Position")
    
    for iter = 2:length(X_all)
        X = X_all[iter]
        plot!(fig, X[idx[1],:], X[idx[2],:],
            label="iter = $(iter - 1)",linewidth=2)
    end

    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plot_circle(p_obs[idx], obs_radius; color=:red, fig=fig)
    end

    return fig
end



# 3D plots

function plot3D_solutions(scp_problem::GuSTOProblem, model, X_all, U_all)
    N = length(X_all)

    pyplot(grid=true)

    idx = [1,2,3]
    local fig
    fig = plot(X_all[1][idx[1],:],X_all[1][idx[2],:],X_all[1][idx[3],:],
        linewidth=2,label="iter = 0",
        xlims=(-0.1,0.5),ylims=(-0.3,0.3),zlims=(-0.2,0.2),
        xlabel="x Position",ylabel="y Position",zlabel="z Position")

    for iter = 2:length(X_all)
        X = X_all[iter]
        plot!(fig,X[idx[1],:],X[idx[2],:],X[idx[3],:],
            linewidth=2,label="iter = $(iter - 1)")
    end

    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plot_sphere(p_obs[idx], obs_radius; color=:red, fig=fig)
    end

    return fig
end



# Additional plotting functions
# Adapted from https://discourse.julialang.org/t/plot-a-circle-with-a-given-radius-with-plots-jl/23295

function plot_circle(position_2d, radius; color=:blue, fig=-1, lab="")
    function circleShape(h, k, r)
        theta = LinRange(0, 2*pi, 500)
        h .+ r*sin.(theta), k .+ r*cos.(theta)
    end

    if fig == -1 # undefined, plot new
        fig = plot(circleShape(position_2d[1],position_2d[2],radius),
                seriestype = [:shape,], lw = 0.5,
                c = color, linecolor = :black,
                fillalpha = 0.5, aspect_ratio = 1, label = lab)
    return fig
    else
        plot!(fig, circleShape(position_2d[1],position_2d[2],radius),
                    seriestype = [:shape,], lw = 0.5,
                    c = color, linecolor = :black,
                    fillalpha = 0.5, aspect_ratio = 1, label = lab)
        return fig
    end
end



function plot_sphere(position_3d, radius; color=:blue, fig=-1, lab="")
    if fig == -1 # undefined, plot new
        u = LinRange(-pi/2.0,pi/2.0,50)
        v = LinRange(0,2*pi,100)
        x, y, z = [], [], []

        for i = 1:length(u)
            for j = 1:length(v)
                push!(x,position_3d[1] + radius*cos(u[i])*cos(v[j]))
                push!(y,position_3d[2] + radius*cos(u[i])*sin(v[j]))
                push!(z,position_3d[3] + radius*sin(u[i]))
            end
        end

        fig = plot(x,y,z,linetype=:surface,label=lab,colorbar=false,
            seriestype = [:shape,],lw = 0.5,
            c = color,fillalpha = 0.5)
    return fig
    else
        u = LinRange(-pi/2.0,pi/2.0,50)
        v = LinRange(0,2*pi,100)
        x, y, z = [], [], []

        for i = 1:length(u)
            for j = 1:length(v)
                push!(x,position_3d[1] + radius*cos(u[i])*cos(v[j]))
                push!(y,position_3d[2] + radius*cos(u[i])*sin(v[j]))
                push!(z,position_3d[3] + radius*sin(u[i]))
            end
        end

        plot!(fig,x,y,z,linetype=:surface,label=lab,colorbar=false,
            seriestype = [:shape,],lw = 0.5,
            c  = color,fillalpha = 0.5)
        return fig
    end
end



