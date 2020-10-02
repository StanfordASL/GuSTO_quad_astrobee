# -------------------------------------------------------------------------------------------
# -   Plotting script for the quadrotor example - Thomas Lew and Riccardo Bonalli 12/2019   -
# -------------------------------------------------------------------------------------------

# Python plotting with matplotlib
using PyCall, LaTeXStrings, Colors
import PyPlot; const plt = PyPlot

include("./Models/quadrotor.jl")
include("./SCP/gusto_problem.jl")



# -------------------------------------------------------
# Plotting misc - Taylor Patrick Reynolds - UW - Sep 2020
# ---
gsp = pyimport("matplotlib.gridspec")
slice(i,j) = pycall(pybuiltin("slice"), PyObject, i,j)
struct CSMPlotCol
    red
    blue
    darkblue
    green
    gray
    purple
    gold
    darkgold
    magenta
    cyan
end
function CSMPlotCol()
    red = [234;61;37]/255
    blue = [0;32;91]/255
    darkblue = [4;28;44]/255
    green = [10;134;61]/255
    gray = [153;153;154]/255
    purple = [51;0;111]/255
    gold = [232;211;162]/255
    darkgold = [145;123;76]/255
    magenta = [1;0;1]
    cyan = [0;1;1]
    return CSMPlotCol(red,blue,darkblue,green,gray,purple,gold,darkgold,magenta,cyan)
end

struct CSMPlotFmt
    col::CSMPlotCol
    circle::Array{Float64,2}
    markersize::Integer
    gridalpha::Float64
    figsize::Tuple{Int64,Int64}
    dblwide::Tuple{Int64,Int64}
    lw::Integer
    labelsize::Integer
    fontsize::Integer
    titlesize::Integer
end
function CSMPlotFmt()
    col     = CSMPlotCol()
    angles  = LinRange(0,2*pi,100)
    circle  = zeros(2,100)
    for k = 1:100
        circle[:,k] = [cos(angles[k]);sin(angles[k])]
    end
    markersize  = 5
    gridalpha   = 0.3
    figsize     = (8,6)
    dblwide     = (12,6)
    linewidth   = 2
    labelsize   = 12
    fontsize    = 14
    titlesize   = 16
    return CSMPlotFmt(col,circle,markersize,gridalpha,
                        figsize,dblwide,linewidth,
                        labelsize,fontsize,titlesize)
end
# ---
# Plotting misc - Taylor Patrick Reynolds - UW - Sep 2020
# -------------------------------------------------------




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
    # Colors of each iteration
    fmt = CSMPlotFmt()
    colors = zeros(3,length(X_all))
    for i = 1:3
        colors[i,:] = LinRange(fmt.col.cyan[i],fmt.col.magenta[i],length(X_all))
    end

    N = length(X_all)
    idx = [1,2]

    fig = plt.figure(figsize=figsize)
    ax  = plt.gca()

    # Plot SCP solutions
    # plt.plot(X_all[1][idx[1],:], X_all[1][idx[2],:],
    #                 label="Initializer", linewidth=2,
    #                 marker="o",color=colors[:,1],linestyle="-")
    # for iter = 2:length(X_all)
    #     X = X_all[iter]
    #     ax.plot(X[idx[1],:], X[idx[2],:],
    #                 label="Iteration $(iter - 1)", linewidth=2,
    #                 marker="o",color=colors[:,iter],linestyle="-")
    # end
    for iter = 1:length(X_all)
        X = X_all[iter]
        ax.plot(X[idx[1],:], X[idx[2],:],
                    label="Iteration $(iter)", linewidth=2,
                    marker="o",color=colors[:,iter],linestyle="-")
    end

    # Plot obstacles
    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plt_circle(ax, p_obs[idx], obs_radius; color="r", alpha=0.3)
    end

    # Settings / Style / Parameters
    PyPlot.rc("text", usetex=true)
    rcParams = PyDict(plt.matplotlib["rcParams"])
    rcParams["font.size"] = 20
    rcParams["font.family"] = "Helvetica"
    plt.xlim(xlims)
    plt.ylim(ylims)
    if B_plot_labels
        plt.title("Open-Loop Quadcopter Trajectories", pad=10)
        plt.xlabel("E [m]")
        plt.ylabel("N [m]")    
        plt.legend(loc="top left", fontsize=18, 
                                      labelspacing=0.1)
    end
    plt.grid(alpha=0.3)

    plt.draw()

    return fig
end

function plt_final_solution(scp_problem::GuSTOProblem, model, X, U)
    # color
    fmt = CSMPlotFmt()

    N = length(X_all)
    idx = [1,2]

    fig = plt.figure(figsize=(4.5, 7.5))
    ax  = plt.gca()

    # Plot final solution
    plt.plot(X[idx[1],:], X[idx[2],:], "bo-", 
                linewidth=2, markersize=6, color=fmt.col.blue)
    plt.plot(Inf*[1,1],Inf*[1,1], "b-", label="Trajectory", color=fmt.col.blue) # for legend


    # Plot obstacles
    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plt_circle(ax, p_obs[idx], obs_radius; color=fmt.col.red, alpha=0.1)
    end
    plt.plot(Inf*[1,1],Inf*[1,1], "r-", label="Obstacle") # for legend

    # Plot Thrust
    for k = 1:(length(X[1,:])-1)
        xk, uk =  X[:,k], U[:,k]

        uMax, mag = 23.2, 1.5

        plt.arrow(xk[idx[1]], xk[idx[2]], 
                    mag*(uk[idx[1]]/uMax), mag*(uk[idx[2]]/uMax),
                    color="g")
    end
    plt.plot(Inf*[1,1],Inf*[1,1], "g-", label="Thrust") # for legend

    # Settings / Style / Parameters
    PyPlot.rc("text", usetex=true)
    rcParams = PyDict(plt.matplotlib["rcParams"])
    rcParams["font.size"] = 20
    rcParams["font.family"] = "Helvetica"
    plt.title("Final Quadcopter Trajectory", pad=15)
    plt.xlim([-0.5,3.])
    plt.ylim([ 0.0,6.])
    plt.xlabel("E [m]")
    plt.ylabel("N [m]")    
    plt.legend(loc="lower right", fontsize=18)
    plt.grid(alpha=0.3)
    
    plt.draw()

    return fig
end

function plt_final_angle_accel(scp_problem::GuSTOProblem, model, X, U)
    # color
    fmt = CSMPlotFmt()


    N = length(X_all)
    t_max = 2.5

    times = collect(range(0,stop=(SCPproblem.N-1)*SCPproblem.dt,length=SCPproblem.N))[1:SCPproblem.N-1]

    NN = length(U[2,:])
    norms_U     = zeros(NN)
    tilt_angles = zeros(NN)
    for k = 1:NN
        norms_U[k]     = norm(U[1:3,k])
        tilt_angles[k] = acosd(U[3,k]/norms_U[k])
    end

    fig = plt.figure(figsize=(5.5, 7.5))

    # -------------
    # Tilt Angle
    plt.subplot(2,1,1)
    plt.plot(times, tilt_angles, "bo-", 
                linewidth=1, markersize=4, color=fmt.col.blue)
    # max tilt angle
    theta_max = (180/pi) * (pi/3.0)
    plt.plot([0,t_max], theta_max*ones(2), 
                    color=fmt.col.red, linestyle="dashed", linewidth=2)
    plt.fill_between([0,t_max], theta_max*ones(2), 90, 
                        color=fmt.col.red, alpha=0.1)
    # Params
    plt.title("Tilt Angle", pad=10)
    plt.xlim([0, t_max])
    plt.ylim([0, 65   ])
    # plt.xlabel("Time [s]")
    plt.ylabel(L"$\theta$ [deg]")
    plt.grid(alpha=0.3)
    plt.draw()

    # -------------
    # Acceleration
    plt.subplot(2,1,2)
    fig.tight_layout(pad=2.0)

    plt.plot(times, norms_U, "bo-", 
                linewidth=1, markersize=4, color=fmt.col.blue)

    # max/min acceleration
    a_min, a_max = 0.6, 23.2
    plt.plot([0,t_max], a_max*ones(2), 
                    color=fmt.col.red, linestyle="dashed", linewidth=2)
    plt.fill_between([0,t_max], a_max*ones(2), 90, 
                        color=fmt.col.red, alpha=0.1)
    plt.plot([0,t_max], a_min*ones(2), 
                    color=fmt.col.red, linestyle="dashed", linewidth=2)
    plt.fill_between([0,t_max], a_min*ones(2), -5, 
                        color=fmt.col.red, alpha=0.1)

    # Parameters / Settings / Labels
    PyPlot.rc("text", usetex=true)
    rcParams = PyDict(plt.matplotlib["rcParams"])
    rcParams["font.size"] = 20
    rcParams["font.family"] = "Helvetica"
    plt.title("Cmd. Acceleration", pad=10)
    plt.xlim([0, t_max])
    plt.ylim([-1, 25  ])
    plt.xlabel("Time [s]")
    plt.ylabel(L"$\|u\|_2\ [m/s^2]$")
    plt.grid(alpha=0.3)
    plt.draw()

    return fig
end