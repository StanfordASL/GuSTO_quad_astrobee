include("Models/polygonal_obstacles.jl")
include("Models/ISS.jl")

# Python plotting with matplotlib
using PyCall, LaTeXStrings
# using PyPlot
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
    # if center.shape[0] != 2 or widths.shape[0] != 2:
    #     assert False, 'plt_rectangle function can only plot in 2d.'
    facecolor = color
    # if noFaceColor:
    #     facecolor = None

    deltas = [widths[1]+additional_w, widths[2]+additional_w]
    bottom_left = (center[1] - deltas[1]/2., center[2] - deltas[2]/2.)
    rect = plt.matplotlib.patches.Rectangle((bottom_left[1],bottom_left[2]),
                                             deltas[1],deltas[2],
                                             # color="g")
                                linewidth=1, edgecolor=color,facecolor=facecolor,alpha=alpha)
    ax.add_patch(rect)
    # if label=="None"
    #     r_edge = plt.matplotlib.patches.Rectangle((bottom_left[0],bottom_left[1]),
    #                                          deltas[0],deltas[1], 
    #                     color=color, alpha=1, fill=false)
    # else
    #     r_edge = plt.matplotlib.patches.Circle(pos, radius=radius, 
    #                     color=color, alpha=1, fill=false, label=label)
    # end
    # ax.add_patch(r_edge)
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
# def plot(X_in, U_in, V_in, m, Xs_true, Us_true):


    lims_btm, lims_up = Vector([8.8,-2., 3.5]), Vector([12.2, 8., 6.5])

#     rcParams['font.family'] = 'serif'
#     rcParams['font.size'] = 14
#     rc('text', usetex=True)

#     X, U, V = X_in[-1, :, :], U_in[-1, :, :], V_in[-1, :, :, :]
#     V = np.moveaxis(V, -1, 0) # (N, n_x,n_x)

    # obstacles, poly_obs = m.obstacles, m.poly_obstacles[-2:]
    obstacles, poly_obs = [], []
    keepin_zones, keepout_zones = get_ISS_zones()


    # *********************************************
    # TRAJECTORY - X-Y
    idx = [1,2]
    fig, ax = plt.subplots(figsize=(6, 10))

    # plot_mean_traj_with_gaussian_uncertainty_ellipses(X.T, V, 
    #                             additional_radius=m.robot_radius,
    #                             color='blue', probability=m.prob, idx=idx,
    #                             alpha = 0.2)
    ax = plt_obstacles(ax, obstacles, keepin_zones, poly_obs, idx)
    ax.text(10.45, -0.2, 
            L"$\mathcal{X}_{obs}$", fontsize=24)
    ax.text(9.2, 1.8, 
            L"$\mathcal{X}_{obs}$", fontsize=24)

    # Plot start & goal positions
    # ax.plot(m.x_init[idx[0]], m.x_init[idx[1]], '+', 
    #             markersize=16, markeredgewidth=3,color='k')
    # ax.plot(m.x_final[idx[0]], m.x_final[idx[1]], '+', 
    #             markersize=16, markeredgewidth=3,color='k')
    # init_text_pos = [m.x_init[idx[0]]-0.,   m.x_init[idx[1]]-0.4]
    # goal_text_pos = [m.x_final[idx[0]]+0., m.x_final[idx[1]]+0.1]
    # ax.text(init_text_pos[0]-0., init_text_pos[1]-0.3, 
    #         r'$\mathbf{x}_{\textrm{init}}$', fontsize=22)
    # ax.text(goal_text_pos[0]-0.4, goal_text_pos[1]+0.4, 
    #         r'$\mathbf{x}_{\textrm{goal}}$', fontsize=22)
    # # ----------------
    ax.tick_params("both", labelsize=24) 
    ax.set_xlabel("X", fontsize=24)
    ax.set_ylabel("Y", rotation="horizontal",fontsize=24)
    plt.xlim([lims_btm[idx[1]], lims_up[idx[1]]])
    plt.ylim([lims_btm[idx[2]], lims_up[idx[2]]])  
    # *********************************************


    # *********************************************
    # TRAJECTORY - Z-Y
    # idx = [2,1]
    # fig, ax = plt.subplots(figsize=(4, 10))

    # plot_mean_traj_with_gaussian_uncertainty_ellipses(X.T, V, 
    #                             additional_radius=m.robot_radius,
    #                             color='blue', probability=m.prob, idx=idx,
    #                             alpha = 0.2)
    # ax = plot_obstacles(ax, obstacles, keepin_zones, poly_obs, idx)

    # # Plot start & goal positions
    # ax.plot(m.x_init[idx[0]], m.x_init[idx[1]], '+', 
    #             markersize=16, markeredgewidth=3,color='k')
    # ax.plot(m.x_final[idx[0]], m.x_final[idx[1]], '+', 
    #             markersize=16, markeredgewidth=3,color='k')
    # init_text_pos = [m.x_init[idx[0]]-0.,   m.x_init[idx[1]]-0.7]
    # goal_text_pos = [m.x_final[idx[0]]+0.12, m.x_final[idx[1]]+0.1]
    # ax.text(init_text_pos[0]-0.2, init_text_pos[1]-0.1, 
    #         r'$\mathbf{x}_{\textrm{init}}$', fontsize=22)
    # ax.text(goal_text_pos[0]-0.4, goal_text_pos[1]+0.35, 
    #         r'$\mathbf{x}_{\textrm{goal}}$', fontsize=22)
    # # ----------------
    # ax.tick_params("both", labelsize=24) 
    # ax.set_xlabel('Z', fontsize=24)
    # plt.xlim([lims_btm[idx[0]], lims_up[idx[0]]])
    # plt.ylim([lims_btm[idx[1]], lims_up[idx[1]]])
    # # *********************************************


    # # *********************************************
    # # CONTROLS
    # fig, ax = plt.subplots(figsize=(10, 6))
    # ax.plot(list(range(0,U.shape[1])),U[0,:], '--o', color='tab:green')
    # ax.plot(list(range(0,U.shape[1])),U[1,:], '--o', color='tab:orange')
    # ax.plot(list(range(0,U.shape[1])),U[2,:], '--o', color='tab:blue')
    # if m.B_feedback:
    #     Sigmas_u = compute_variances_controls(m, V)
    #     plot_mean_var(U[0,:], Sigmas_u[:,0,0], prob=0.9, color='tab:green')
    #     plot_mean_var(U[1,:], Sigmas_u[:,1,1], prob=0.9, color='tab:orange')
    #     plot_mean_var(U[2,:], Sigmas_u[:,2,2], prob=0.9, color='tab:blue')
    # ax.legend([r'$F_x$',r'$F_y$',r'$F_z$'], fontsize=24)
    # ax.tick_params("both", labelsize=24) 
    # ax.set_xlabel('k', fontsize=24)
    # ax.set_ylabel('F', rotation="horizontal",fontsize=24)

    # fig, ax = plt.subplots(figsize=(10, 6))
    # ax.plot(list(range(0,U.shape[1])),U[3,:], '--o', color='tab:green')
    # ax.plot(list(range(0,U.shape[1])),U[4,:], '--o', color='tab:orange')
    # ax.plot(list(range(0,U.shape[1])),U[5,:], '--o', color='tab:blue')
    # if m.B_feedback:
    #     Sigmas_u = compute_variances_controls(m, V)
    #     plot_mean_var(U[3,:], Sigmas_u[:,3,3], prob=0.9, color='tab:green')
    #     plot_mean_var(U[4,:], Sigmas_u[:,4,4], prob=0.9, color='tab:orange')
    #     plot_mean_var(U[5,:], Sigmas_u[:,5,5], prob=0.9, color='tab:blue')
    # ax.legend([r'$M_x$',r'$M_y$',r'$M_z$'], fontsize=24)
    # ax.tick_params("both", labelsize=24) 
    # ax.set_xlabel('k', fontsize=24)
    # ax.set_ylabel('M', rotation="horizontal",fontsize=24)
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    plt.draw()
    plt.show()
    # *********************************************
end