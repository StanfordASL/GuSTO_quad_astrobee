# -----------------------------------------------------------------------------------------------
# -   Initializing script for the freeflyer example - Thomas Lew and Riccardo Bonalli 12/2019   -
# -----------------------------------------------------------------------------------------------



using LinearAlgebra
using Ipopt
using JuMP
using DifferentialEquations
using NLsolve
using Plots

# Python plotting with matplotlib
using PyCall, LaTeXStrings
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





# --------------------
# Plotting with PyPlot
# --------------------
function plt_circle(ax, pos, radius; color="k", alpha=1., label="None")
    # Filled circle
    circle = plt.matplotlib.patches.Circle(pos, radius=radius,
                    color=color, fill=true, alpha=alpha)
    ax.add_patch(circle)
    # Edge of circle, with alpha 1.
    if label=="None"
        c_edge = plt.matplotlib.patches.Circle(pos, radius=radius, 
                        color=color, alpha=1, fill=false)
    else
        c_edge = plt.matplotlib.patches.Circle(pos, radius=radius, 
                        color=color, alpha=1, fill=false, label=label)
    end
    ax.add_patch(c_edge)
    return ax
end

function plt_solutions(scp_problem::GuSTOProblem, model, X_all, U_all;
                        xlims=[-0.5,3.], ylims=[0.0,6.], figsize=(8,6), B_plot_labels=true)
    N = length(X_all)
    idx = [1,2]


    fig = plt.figure(figsize=figsize)
    ax  = plt.gca()

    # Plot SCP solutions
    plt.plot(X_all[1][idx[1],:], X_all[1][idx[2],:],
                    label="Initializer", linewidth=2)
    for iter = 2:length(X_all)
        X = X_all[iter]
        ax.plot(X[idx[1],:], X[idx[2],:],
                    label="Iterate $(iter - 1)", linewidth=2)
    end

    # Plot obstacles
    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plt_circle(ax, p_obs[idx], obs_radius; color="r", alpha=0.3)
    end

    # Settings / Style / Parameters
    # csfont = ['fontname':'Comic Sans MS']
    # hfont = ['fontname':'Helvetica']
    # plt.title('title',**csfont)
    # plt.xlabel('xlabel', **hfont)
    # PyPlot.rc("text", usetex=true)
    # rcParams = PyDict(plt.matplotlib["rcParams"])
    # rcParams["font.size"] = 15
    # rcParams["font.family"] = "Helvetica"
    # plt.xlim(xlims)
    # plt.ylim(ylims)
    if B_plot_labels
        plt.title("Open-Loop Astrobee Trajectories")
        plt.xlabel("x")
        plt.ylabel("y")    
        plt.legend(loc="lower right")
    end
    plt.grid(alpha=0.3)

    plt.draw()

    return fig
end

function plt_final_solution(scp_problem::GuSTOProblem, model, X, U)
    N = length(X_all)
    idx = [1,2]

    fig = plt.figure(figsize=(4.5, 7.5))
    ax  = plt.gca()

    # Plot final solution
    plt.plot(X[idx[1],:], X[idx[2],:], "bo-", 
                linewidth=2, markersize=6)
    plt.plot(Inf*[1,1],Inf*[1,1], "b-", label="Trajectory") # for legend

    # Plot Thrust
    for k = 1:(length(X[1,:])-1)
        xk, uk =  X[:,k], U[:,k]

        uMax, mag = 23.2, 1000.5

        plt.arrow(xk[idx[1]], xk[idx[2]], 
                    mag*(uk[idx[1]]/uMax), mag*(uk[idx[2]]/uMax),
                    color="g")
    end
    plt.plot(Inf*[1,1],Inf*[1,1], "g-", label="Thrust") # for legend

    # Plot obstacles
    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plt_circle(ax, p_obs[idx], obs_radius; color="r", alpha=0.4)
    end
    plt.plot(Inf*[1,1],Inf*[1,1], "r-", label="Obstacle") # for legend

    # Settings / Style / Parameters
    # csfont = ['fontname':'Comic Sans MS']
    # hfont = ['fontname':'Helvetica']
    # plt.title('title',**csfont)
    # plt.xlabel('xlabel', **hfont)
    # PyPlot.rc("text", usetex=true)
    # rcParams = PyDict(plt.matplotlib["rcParams"])
    # rcParams["font.size"] = 15
    # rcParams["font.family"] = "Helvetica"
    plt.title("Final Trajectory")
    # plt.xlim([-0.5,3.])
    # plt.ylim([ 0.0,6.])
    plt.xlabel("x")
    plt.ylabel("y")    
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    
    plt.draw()

    return fig
end

# function plt_final_angle_accel(scp_problem::GuSTOProblem, model, X, U)
#     N = length(X_all)
#     t_max = 2.

#     times = collect(range(0,stop=(SCPproblem.N-1)*SCPproblem.dt,length=SCPproblem.N))[1:SCPproblem.N-1]
#     norms_U = sqrt.(U[1,:].^2+U[2,:].^2+U[3,:].^2)

#     fig = plt.figure(figsize=(4.5, 7.5))

#     # -------------
#     # Tilt Angle
#     plt.subplot(2,1,1)
#     tilt_angles = U[3,:]./norms_U
#     tilt_angles = (180. / pi) * tilt_angles
#     plt.plot(times, tilt_angles, "bo-", 
#                 linewidth=1, markersize=4)
#     # max tilt angle
#     theta_max = (180/pi) * (pi/3.0)
#     plt.plot([0,t_max], theta_max*ones(2), 
#                     color="r", linestyle="dashed", linewidth=2)
#     plt.fill_between([0,t_max], theta_max*ones(2), 90, 
#                         color="r", alpha=0.2)
#     # Params
#     plt.title("Tilt Angle")
#     plt.xlim([0, t_max])
#     plt.ylim([0, 65   ])
#     # plt.xlabel("Time [s]")
#     plt.ylabel(L"$\theta$ [deg]")
#     plt.grid(alpha=0.3)
#     plt.draw()

#     # -------------
#     # Acceleration
#     plt.subplot(2,1,2)
#     fig.tight_layout(pad=2.0)

#     plt.title("Cmd. Acceleration")
#     plt.plot(times, norms_U, "bo-", 
#                 linewidth=1, markersize=4)

#     # max/min acceleration
#     a_min, a_max = 0.6, 23.2
#     plt.plot([0,t_max], a_max*ones(2), 
#                     color="r", linestyle="dashed", linewidth=2)
#     plt.fill_between([0,t_max], a_max*ones(2), 90, 
#                         color="r", alpha=0.2)
#     plt.plot([0,t_max], a_min*ones(2), 
#                     color="r", linestyle="dashed", linewidth=2)
#     plt.fill_between([0,t_max], a_min*ones(2), -5, 
#                         color="r", alpha=0.2)

#     plt.xlim([0, t_max])
#     plt.ylim([-1, 25  ])
#     plt.xlabel("Time [s]")
#     plt.ylabel(L"$u [m/s^2]$")
#     plt.grid(alpha=0.3)
#     plt.draw()

#     return fig
# end