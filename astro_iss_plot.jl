include("Models/polygonal_obstacles.jl")
include("Models/ISS.jl")

# Python plotting with matplotlib
using PyCall, LaTeXStrings
import PyPlot; const plt = PyPlot

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
                        color=color, alpha=1, fill=false)
    else
        c_edge = plt.matplotlib.patches.Circle(pos, radius=radius, 
                        color=color, alpha=1, fill=false, label=label)
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
        ax = plt_circle(ax, pos, radius, color="r", alpha=0.4)
    end
    for obs in keepin_zones
        center, widths = obs.c, 2. * Vector([obs.dx,obs.dy,obs.dz])
        ax = plt_rectangle(ax, center[idx], widths[idx], color="g")
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
    ax.text(10.45, -0.2, 
            L"$\mathcal{X}_{obs}$", fontsize=24)
    ax.text(9.2, 1.8, 
            L"$\mathcal{X}_{obs}$", fontsize=24)

    # # ----------------
    ax.tick_params("both", labelsize=24) 
    ax.set_xlabel("X", fontsize=24)
    ax.set_ylabel("Y", rotation="horizontal",fontsize=24)
    plt.xlim([lims_btm[idx[1]], lims_up[idx[1]]])
    plt.ylim([lims_btm[idx[2]], lims_up[idx[2]]])  
    # *********************************************

    plt.draw()
    plt.show()
    # *********************************************
end