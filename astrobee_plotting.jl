# -----------------------------------------------------------------------------------------------
# -   Initializing script for the freeflyer example - Thomas Lew and Riccardo Bonalli 12/2019   -
# -----------------------------------------------------------------------------------------------


using Plots

# Python plotting with matplotlib
using PyCall, LaTeXStrings, Colors
import PyPlot; const plt = PyPlot

include("./Models/astrobee.jl")
include("./SCP/gusto_problem.jl")

include("Models/polygonal_obstacles.jl")
include("Models/ISS.jl")


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
        # xlims=(-0.1,0.5),ylims=(-0.3,0.3),zlims=(-0.2,0.2),
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

function plot3D_final_solution(scp_problem::GuSTOProblem, model, X, U; B_save_png=false)
    N = length(X_all)

    font = Plots.font("Helvetica", 16)
    titlefont = Plots.font("Helvetica", 28)
    smallfont = Plots.font("Helvetica", 1)
    pyplot(grid=true, guidefont=font, titlefont=titlefont,
        xtickfont=font, ytickfont=font, ztickfont=font, leg=false)
        # xtickfont=smallfont, ytickfont=smallfont, ztickfont=smallfont, leg=false)

    idx = [1,2,3]
    # camera = (40, 40)
    camera = (40, 20)

    xlims = (9.75,12.25)
    ylims = (1.,  4.5)
    # zlims = (4.2, 5.1)
    zlims = (4.2, 5.3)

    X = copy(X)
    idx_in =            X[idx[1],:].<= xlims[2]
    idx_in = idx_in .& (X[idx[1],:].>= xlims[1])
    idx_in = idx_in .& (X[idx[2],:].<= ylims[2])
    idx_in = idx_in .& (X[idx[2],:].>= ylims[1])
    X = X[:, idx_in]

    xticks = collect(10.:1:12.)
    yticks = collect(round(ylims[1],digits=1):1.:round(ylims[2],digits=1))
    zticks = collect(4.2:0.4:round(zlims[2],digits=1))

    local fig
    # Main trajectory plot
    fig = plot(X[idx[1],:],X[idx[2],:],X[idx[3],:],
        linewidth=2,
        # label="iter = 0",
        title="Close View",
        color=:blue,
        xlims=xlims, ylims=ylims,zlims=zlims,
        xticks=xticks,yticks=yticks,zticks=zticks,
        # xlabel="X", ylabel="Y", zlabel="Z",
        camera = camera,
        fontfamily="Helvetica",margin=5Plots.mm,left_margin = 8Plots.mm,
        dpi=300)


    # Spherical obstacles
    # color = :red
    color = RGB(255. / 255,  120. / 255,   120. / 255)
    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plot_sphere(p_obs[idx], obs_radius; color=color, fig=fig)
    end
    for obs_i = 1:length(model.obstacles)
        # sphere projection
        p_obs_2d, obs_rad = model.obstacles[obs_i][1][1:2], model.obstacles[obs_i][2]
        u = LinRange(-pi/2.0,pi/2.0,50)
        v = LinRange(0,2*pi,100)
        x, y, z = [], [], []
        for i = 1:length(u)
            for j = 1:length(v)
                push!(x,p_obs_2d[1] + obs_rad*cos(u[i])*cos(v[j]))
                push!(y,p_obs_2d[2] + obs_rad*cos(u[i])*sin(v[j]))
                push!(z,zlims[1])
            end
        end
        plot!(fig,x,y,z,linetype=:surface,colorbar=false,
            seriestype = [:shape,],
            c  = color,alpha = 0.2, dpi=300)
    end

    # Trajectory
    fig = plot!(X[idx[1],:],X[idx[2],:],zlims[1]*ones(length(X[idx[3],:])),
        linewidth=2,color=:gray,alpha=0.7, dpi=300)
    fig = scatter!(X[idx[1],:],X[idx[2],:],zlims[1]*ones(length(X[idx[3],:])),
        linewidth=4,color=:gray,markerstrokecolor=:darkgray,
        markersize=6, alpha=0.7, dpi=300)
    fig = scatter!(X[idx[1],:],X[idx[2],:],X[idx[3],:],
        linewidth=4,color=:blue,markerstrokecolor=:darkblue,
        markersize=6, alpha=1., dpi=300)

    if B_save_png
        Plots.savefig("figs/astrobee/astro_close3d.png")#,bbox_inches="tight",dpi=300)
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
# ----------------------------------------------------
function plt_circle(ax, pos, radius; 
                        color="k", alpha=1., label="None")
    # Filled circle
    circle = plt.matplotlib.patches.Circle(pos, radius=radius,
                    color=color, fill=true, alpha=alpha)
    ax.add_patch(circle)
    # Edge of circle, with alpha 1.
    if label=="None"
        c_edge = plt.matplotlib.patches.Circle(pos, radius=radius, 
                        color=color, alpha=0.8, fill=false)
    else
        c_edge = plt.matplotlib.patches.Circle(pos, radius=radius, 
                        color=color, alpha=0.8, fill=false, label=label)
    end
    ax.add_patch(c_edge)
    return ax
end
function plt_rectangle(ax, center, widths, additional_w=0; 
                           color="b", alpha=0.1, noFaceColor=false, 
                           label="None")
    """
    Plots a rectangle with a given center and total widths
    arguments:  - center    - (2,) center
                - widths    - (2,) total widths  from and to edges of rectangle
    """
    facecolor = color

    deltas = [widths[1]+additional_w, widths[2]+additional_w]
    bottom_left = (center[1] - deltas[1]/2., center[2] - deltas[2]/2.)
    rect = plt.matplotlib.patches.Rectangle((bottom_left[1],bottom_left[2]),
                                             deltas[1],deltas[2],
                                             # color="g")
                                linewidth=1, edgecolor=color,facecolor=facecolor,alpha=alpha)
    ax.add_patch(rect)
    return ax
end

function plt_obstacles(ax, obstacles, keepin_zones, poly_obs, idx=[1,2])
    for obs in obstacles
        pos, radius = Vector([obs[1][idx[1]],obs[1][idx[2]]]), obs[2]
        ax = plt_circle(ax, pos, radius, color="r", alpha=0.1)
    end
    for obs in keepin_zones
        center, widths = obs.c, 2. * Vector([obs.dx,obs.dy,obs.dz])
        green = [10;134;61]/255
        ax = plt_rectangle(ax, center[idx], widths[idx], color=green, alpha=0.15)
        
    end
    for obs in poly_obs
        center, widths = obs.c, 2. * Vector([obs.dx,obs.dy,obs.dz])
        ax = plt_rectangle(ax, center[idx], widths[idx], color="r", alpha=0.4)
    end
    return ax
end

# ---------------------------------------------

function plt_ISS()
    lims_btm, lims_up = Vector([8.8,-2., 3.5]), Vector([12.2, 8., 6.5])

    obstacles, poly_obs = [], []
    keepin_zones, keepout_zones = get_ISS_zones()


    # *********************************************
    # TRAJECTORY - X-Y
    idx = [1,2]
    fig, ax = plt.subplots(figsize=(6, 10))

    ax = plt_obstacles(ax, obstacles, keepin_zones, poly_obs, idx)

    # Plot start & goal positions
    ax.tick_params("both", labelsize=24) 
    ax.set_xlabel("X", fontsize=24)
    ax.set_ylabel("Y", rotation="horizontal",fontsize=24)
    plt.xlim([5.,15.])
    plt.ylim([-8.,10.])  
    # *********************************************

    plt.draw()
    plt.show()
    # *********************************************
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

    # fig = plt.figure(figsize=figsize)
    fig = plt.figure(figsize=(10,8))
    ax  = plt.gca()

    # # Plot SCP solutions
    # plt.plot(X_all[1][idx[1],:], X_all[1][idx[2],:],
    #                 label="Initializer", linewidth=2)
    # for iter = 2:length(X_all)
    #     X = X_all[iter]
    #     ax.plot(X[idx[1],:], X[idx[2],:],
    #                 label="Iterate $(iter - 1)", linewidth=2)
    # end
    for iter = 1:length(X_all)
        X = X_all[iter]
        ax.plot(X[idx[1],:], X[idx[2],:],
                    label="Iteration $(iter)", linewidth=2,
                    marker="o",color=colors[:,iter],linestyle="-")
    end

    ## ---------------------- ##
    ## ----- ENVIRONMENT ---- ##

    ## ----- obstacles ------ ##
    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plt_circle(ax, p_obs[idx], obs_radius; color="r", alpha=0.1)
    end
    plt.plot(Inf*[1,1],Inf*[1,1], "r-") # for legend

    ## -------- ISS ---------- ##
    lims_btm, lims_up = Vector([5.5,-1.5, 3.5]), Vector([12.2, 7.5, 6.5])
    obstacles, poly_obs = [], []#model.poly_obstacles
    keepin_zones, keepout_zones = get_ISS_zones()
    ax = plt_obstacles(ax, obstacles, keepin_zones, poly_obs, idx)
    # *********************************************

    # Settings / Style / Parameters
    PyPlot.rc("text", usetex=true)
    rcParams = PyDict(plt.matplotlib["rcParams"])
    rcParams["font.size"] = 22
    rcParams["font.family"] = "Helvetica"
    plt.xlim([5.,13])
    plt.ylim([-3,8.])
    if B_plot_labels
        plt.title("All Free-Flyer Trajectories")
        ax.title.set_position([.5, 1.01])
        # plt.xlabel("X", fontsize=26)
        # plt.ylabel("Y", fontsize=26, rotation="horizontal")    
        ax.set_xlabel(L"x_{\mathcal{I}}\ [m]",fontsize=26)
        ax.set_ylabel(L"y_{\mathcal{I}}\ [m]",fontsize=26)
        plt.legend(loc="upper left", labelspacing=0.1)
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
        plt_circle(ax, p_obs[idx], obs_radius; color="r", alpha=0.1)
    end
    plt.plot(Inf*[1,1],Inf*[1,1], "r-", label="Obstacle") # for legend

    # Settings / Style / Parameters
    PyPlot.rc("text", usetex=true)
    rcParams = PyDict(plt.matplotlib["rcParams"])
    rcParams["font.size"] = 20
    rcParams["font.family"] = "Helvetica"
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


function plt_final_solution_ISS(scp_problem::GuSTOProblem, model, X, U,
                                idx=[1,2])
    # color
    fmt = CSMPlotFmt()

    N = length(X_all)
    # idx = [1,2]
    # fig, ax = plt.subplots(figsize=(6, 10))
    fig, ax = plt.subplots(figsize=(8, 8))

    ax  = plt.gca()

    # Plot final solution
    plt.plot(X[idx[1],:], X[idx[2],:], "bo-", 
                linewidth=2, markersize=6, color=fmt.col.blue)
    plt.plot(Inf*[1,1],Inf*[1,1], "b-", label="Trajectory", color=fmt.col.blue) # for legend

    # Plot Thrust
    for k = 1:(length(X[1,:])-1)
        xk, uk =  X[:,k], U[:,k]

        uMax, mag = 23.2, 0.3*1000.5

        plt.arrow(xk[idx[1]], xk[idx[2]], 
                    mag*(uk[idx[1]]/uMax), mag*(uk[idx[2]]/uMax),
                    color="g", linewidth=2.)
    end
    plt.plot(Inf*[1,1],Inf*[1,1], "g-", label="Thrust") # for legend

    ## ---------------------- ##
    ## ----- ENVIRONMENT ---- ##

    ## ----- obstacles ------ ##
    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plt_circle(ax, p_obs[idx], obs_radius; color=fmt.col.red, alpha=0.1)
    end
    plt.plot(Inf*[1,1],Inf*[1,1], "r-", label="Obstacle") # for legend

    ## -------- ISS ---------- ##
    lims_btm, lims_up = Vector([5.5,-2.5, 3.5]), Vector([12.2, 8.0, 6.5])
    obstacles, poly_obs = [], []#model.poly_obstacles
    keepin_zones, keepout_zones = get_ISS_zones()
    ax = plt_obstacles(ax, obstacles, keepin_zones, poly_obs, idx)
    # *********************************************



    # Settings / Style / Parameters
    PyPlot.rc("text", usetex=true)
    rcParams = PyDict(plt.matplotlib["rcParams"])
    rcParams["font.size"] = 26
    rcParams["font.family"] = "Helvetica"
    plt.title("Final Trajectory", fontsize=38)
    ax.title.set_position([.5, 1.01])

    plt.xlabel(L"x_{\mathcal{I}}\ [m]", fontsize=28)
    plt.ylabel(L"y_{\mathcal{I}}\ [m]", fontsize=28)    
    # plt.ylabel(L"y_{\mathcal{I}}\ [m]", fontsize=28, rotation="horizontal")    
    plt.xlim([lims_btm[idx[1]],lims_up[idx[1]]])
    plt.ylim([lims_btm[idx[2]],lims_up[idx[2]]])
    # plt.zlim([lims_btm[3],lims_up[3]])
    plt.legend(loc="center left", bbox_to_anchor=(0.,0.7), handletextpad=0.1)
    plt.grid(alpha=0.3)
    
    plt.draw()

    return fig
end