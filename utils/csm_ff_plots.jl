using PyCall, LaTeXStrings, Colors
import PyPlot
const plt = PyPlot
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

function plt_rectangle(ax, center, widths;
                           color="b", alpha=0.1, noFaceColor=false,
                           label="None")
    """
    Plots a rectangle with a given center and total widths
    arguments:  - center    - (2,) center
                - widths    - (2,) total widths from and to edges of rectangle
    """
    deltas = [widths[1], widths[2]]
    bottom_left = (center[1] - deltas[1]/2.0, center[2] - deltas[2]/2.0)
    rect = plt.matplotlib.patches.Rectangle((bottom_left[1],bottom_left[2]),
                                             deltas[1],deltas[2],
                                             linewidth=1,
                                             edgecolor=color,
                                             facecolor=color,
                                             alpha=alpha)
    ax.add_patch(rect)
    return ax
end




# ---------------------------------------------------
# ---------------------------------------------------
# ---------------------------------------------------
# Conversion utils for plotting
function rad2deg_arr!(v::Array{Float64,2})
    for ii = 1:size(v,1)
        for jj = 1:size(v,2)
            v[ii,jj] = rad2deg(v[ii,jj])
        end
    end
end
function quat_2_dcm(q::Vector{Float64})
    if length(q)!=4
        error("Quaternion input must be of length 4")
    end
    q1 = q[1]
    q2 = q[2]
    q3 = q[3]
    q4 = q[4]

    return [ 2*q1^2+2*q4^2-1.0      2*(q1*q2+q3*q4)     2*(q1*q3-q2*q4);
             2*(q1*q2-q3*q4)        2*q2^2+2*q4^2-1.0   2*(q2*q3+q1*q4);
             2*(q1*q3+q2*q4)        2*(q2*q3-q1*q4)     2*q3^2+2*q4^2-1.0 ]
end
function quat_2_rpy(q::Array{Float64,2})
    if size(q,1)!=4
        q = transpose(q)
    end
    if size(q,1)!=4
        error("Quaternion input must have four rows")
    end
    N = size(q,2)
    rpy = zeros(3,N)
    for k = 1:size(q,2)
        Rk = quat_2_dcm(q[:,k])

        rpy[1,k] =  atan( Rk[2,3], Rk[3,3] ) # roll
        rpy[2,k] = -asin( Rk[1,3] )          # pitch
        rpy[3,k] =  atan( Rk[1,2], Rk[1,1] ) # yaw
    end
    return rpy
end
# Conversion utils for plotting
# ---------------------------------------------------
# ---------------------------------------------------
# ---------------------------------------------------



function csm_plots_freeflyer(X_all, U_all, N, m::Astrobee)
    # integrate nonlinear dynamics
    t_grid = LinRange(0,m.tf,N)
    T_grid = LinRange(0,m.tf,250)
    # u      = prob.new_sol.control
    # f(t,y) = dynamics(t,y,u,t_grid,prob.pars.mdl_pars)
    # X      = rk4(f,T_grid,x0)

    fmt = CSMPlotFmt()

    # create first figure
    fig = plt.figure(num=10, figsize=fmt.dblwide)
    grd = gsp.GridSpec(2, 4, wspace=0.5, hspace=0.5)
    ax1 = plt.subplot(get(grd,(slice(0,2),slice(0,2))))
    ax2 = plt.subplot(get(grd,(0,2)))
    ax3 = plt.subplot(get(grd,(0,3)))
    ax4 = plt.subplot(get(grd,(1,2)))
    ax5 = plt.subplot(get(grd,(1,3)))

    csm_plot_freeflyer_trj(ax1,m,N,X_all,U_all,fmt)
    csm_plot_freeflyer_thrust(ax2,m,N,X_all,U_all,fmt)
    csm_plot_freeflyer_torque(ax3,m,N,X_all,U_all,fmt)
    csm_plot_freeflyer_attitude(ax4,m,N,X_all,U_all,fmt)
    csm_plot_freeflyer_attituderate(ax5,m,N,X_all,U_all,fmt)

    plt.show()
    plt.savefig("figs/astrobee/freeflyer_final_trj.png",bbox_inches="tight",dpi=300)

    # create the second figure
    # fig = plt.figure(figsize=fmt.figsize)
    # ax = plt.gca()
    # # csm_plot_freeflyer_alltrjs(ax,prob,fmt)

    # plt.show()
    # plt.savefig("figs/freeflyer_all_trjs.png",bbox_inches="tight",dpi=300)

    # csm_freeflyer_sd(prob,fmt,X)
    return nothing
end


function csm_plot_freeflyer_trj(ax,m,N,X_all,U_all,fmt::CSMPlotFmt)
    x   = X_all[end]
    u   = U_all[end]
    circle = fmt.circle
    scl = 0.3

    idx = [1,2]


    # plot discrete solution
    ax.plot(x[1,:],x[2,:],
            label="Solution",
            marker="o",color=fmt.col.blue,linestyle="-",
            markersize=fmt.markersize)

    # add thrust vectors
    udir = u[:,1]/norm(u[:,1])
    xs = [ x[1,1], x[1,1]+scl*udir[1] ]
    ys = [ x[2,1], x[2,1]+scl*udir[2] ]
    lines = Any[collect(zip(xs,ys))]
    for k = 1:(N-1)
        udir = u[:,k]/norm(u[:,k])
        xs = [ x[1,k], x[1,k]+scl*udir[1] ]
        ys = [ x[2,k], x[2,k]+scl*udir[2] ]
        push!(lines,collect(zip(xs,ys)))
    end
    horz_thrust_vecs = plt.matplotlib.collections.LineCollection(lines,
                                                    color=fmt.col.green,
                                                    label="Thrust Direction",
                                                    linewidth=fmt.lw)

    ## ----- obstacles ------ ##
    for obs_i = 1:length(m.obstacles)
        p_obs, obs_radius = m.obstacles[obs_i][1], m.obstacles[obs_i][2]
        plt_circle(ax, p_obs[idx], obs_radius; color=fmt.col.red, alpha=0.4)
    end
    ax.plot(Inf*[1,1],Inf*[1,1], "r-", label="Obstacle") # for legend
    ## -------- ISS ---------- ##
    lims_btm, lims_up = Vector([5.5,-2.5, 3.5]), Vector([12.2, 8.0, 6.5])
    obstacles, poly_obs = [], []
    keepin_zones, keepout_zones = get_ISS_zones()
    ax = plt_obstacles(ax, obstacles, keepin_zones, poly_obs, idx)
    # *********************************************


    ax.add_collection(horz_thrust_vecs)
    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(5,13)
    ax.set_ylim(-3,8)
    ax.set_xlabel(L"x_{\mathcal{I}}\ [m]",fontsize=fmt.fontsize)
    ax.set_ylabel(L"y_{\mathcal{I}}\ [m]",fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    ax.legend(fontsize=12,loc="upper left")
    ax.set_title("Final FreeFlyer Trajectory",fontsize=fmt.titlesize)

    return nothing
end

function csm_plot_freeflyer_thrust(ax,m,N,X_all,U_all,fmt::CSMPlotFmt)
    id_F = [1,2,3]
    scl = 1000; # scale factor [ N - > mN ]
    F = scl * U_all[end][id_F,:]
    tf  = model.tf
    T  = LinRange(0,tf,N-1)

    F_nrm_max = scl * 72e-3

    # compute control input norms for plotting
    F_nrm = zeros(N-1)
    for k = 1:(N-1)
        F_nrm[k] = norm(F[:,k])
    end

    ax.plot(T,F[1,:],label=L"x",
            color=fmt.col.red,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,F[2,:],label=L"y",
            color=fmt.col.blue,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,F[3,:],label=L"z",
            color=fmt.col.green,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,F_nrm,label=L"||F||",
            color=[0;0;0],
            linestyle="--",
            linewidth=fmt.lw)
    ax.plot([0;tf],[F_nrm_max;F_nrm_max],label=L"\|F\|_{2,max}",
            color=[1;0;0],
            linestyle="--",
            linewidth=fmt.lw)
    ylmin, ylmax = ax.get_ylim()
    ax.fill_between([0;tf],[F_nrm_max;F_nrm_max],[ylmax;ylmax],
                    facecolor=fmt.col.red,alpha=0.1)

    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(0,tf)
    ax.set_ylim(ylmin,ylmax)
    ax.set_xlabel("Time [s]",fontsize=fmt.fontsize)
    ax.set_ylabel("Thrust [mN]",fontsize=fmt.fontsize)
    # ax.legend(fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    # ax.legend(fontsize=fmt.fontsize)
    # ax.set_title(L"Final\ FreeFlyer\ Thrust",fontsize=fmt.titlesize)

    return nothing
end

function csm_plot_freeflyer_torque(ax,m,N,X_all,U_all,fmt::CSMPlotFmt)
    id_M = [4,5,6]
    scl = 1000 # scale factor [Nm -> mNm]
    M = scl * U_all[end][id_M,:]
    tf  = m.tf
    T  = LinRange(0,tf,N-1)

    M_nrm_max = scl * 2e-3

    # compute control input norms for plotting
    M_nrm = zeros(N-1)
    for k = 1:(N-1)
        M_nrm[k] = norm(M[:,k])
    end

    # Plot moment/torque trajectory
    ax.plot(T,M[1,:],label=L"x",
            color=fmt.col.red,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,M[2,:],label=L"y",
            color=fmt.col.blue,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,M[3,:],label=L"z",
            color=fmt.col.green,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,M_nrm,label=L"||M||",
            color=[0;0;0],
            linestyle="--",
            linewidth=fmt.lw)
    ax.plot([0;tf],[M_nrm_max;M_nrm_max],label=L"\|M\|_{2,max}",
            color=[1;0;0],
            linestyle="--",
            linewidth=fmt.lw)
    ylmin, ylmax = ax.get_ylim()
    ax.fill_between([0;tf],[M_nrm_max;M_nrm_max],[ylmax;ylmax],
                    facecolor=fmt.col.red,alpha=0.1)

    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(0,tf)
    ax.set_ylim(ylmin,ylmax)
    ax.set_xlabel("Time [s]",fontsize=fmt.fontsize)
    ax.set_ylabel("Torque [mNm]",fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    ax.legend(fontsize=10)
    # ax.set_title(L"Final\ FreeFlyer\ Controls",fontsize=fmt.titlesize)

    return nothing
end

function csm_plot_freeflyer_attitude(ax,m,N,X_all,U_all,fmt::CSMPlotFmt)
    id_q = [7,8,9,10]
    quat_d = X_all[end][id_q,:]
    quat_c = X_all[end][id_q,:]
    # convert quaternion to 3-2-1 Euler angles
    rpy_d  = quat_2_rpy(quat_d)
    rpy_c  = quat_2_rpy(quat_c)

    # convert to degrees
    rad2deg_arr!(rpy_d)
    rad2deg_arr!(rpy_c)

    tf  = m.tf
    T   = LinRange(0,tf,N)
    Ti  = LinRange(0,tf,size(X,2))

    ax.plot(T,rpy_d[1,:],
            color=fmt.col.red,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize)
    ax.plot(Ti,rpy_c[1,:],
            label="roll",
            color=fmt.col.red,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,rpy_d[2,:],
            color=fmt.col.green,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize)
    ax.plot(Ti,rpy_c[2,:],
            label="pitch",
            color=fmt.col.green,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,rpy_d[3,:],
            color=fmt.col.blue,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize)
    ax.plot(Ti,rpy_c[3,:],
            label="yaw",
            color=fmt.col.blue,
            linestyle="-",
            linewidth=fmt.lw)
    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(0,tf)
    ax.set_ylim(-180,180)
    ax.set_yticks([-180,-90,0,90,180])
    ax.set_xlabel("Time [s]",fontsize=fmt.fontsize)
    ax.set_ylabel("Attitude [deg]",fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    ax.legend(fontsize=8,loc="upper right")
    # ax.set_title(L"Attitude",fontsize=fmt.titlesize)

    return nothing
end

function csm_plot_freeflyer_attituderate(ax,m,N,X_all,U_all,fmt::CSMPlotFmt)
    id_w = [11,12,13]
    wB_d   = X_all[end][id_w,:]
    wB_c   = X_all[end][id_w,:]

    # convert to degrees
    rad2deg_arr!(wB_d)
    rad2deg_arr!(wB_c)

    tf  = m.tf
    T   = LinRange(0,tf,N)
    Ti  = LinRange(0,tf,size(X,2))

    w_nrm_max = rad2deg(1.)

    ax.plot(T,wB_d[1,:],
            color=fmt.col.red,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize,
            linewidth=fmt.lw)
    ax.plot(Ti,wB_c[1,:],
            label=L"\omega_x",
            color=fmt.col.red,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,wB_d[2,:],
            color=fmt.col.green,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize,
            linewidth=fmt.lw)
    ax.plot(Ti,wB_c[2,:],
            label=L"\omega_y",
            color=fmt.col.green,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,wB_d[3,:],
            color=fmt.col.blue,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize,
            linewidth=fmt.lw)
    ax.plot(Ti,wB_c[3,:],
            label=L"\omega_z",
            color=fmt.col.blue,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot([0;tf],[w_nrm_max;w_nrm_max],label=L"\|\omega\|_{\infty,max}",
            color=[1;0;0],linestyle="--",
            linewidth=fmt.lw)
    ylmin, ylmax = ax.get_ylim()
    ax.fill_between([0;tf],[w_nrm_max;w_nrm_max],[ylmax;ylmax],
                    facecolor=fmt.col.red,alpha=0.1)

    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(0,tf)
    ax.set_ylim(ylmin,ylmax)
    ax.set_xlabel("Time [s]",fontsize=fmt.fontsize)
    ax.set_ylabel("Angular Rate [deg/s]",fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    ax.legend(fontsize=10)
    # ax.set_title(L"Attitude Rate",fontsize=fmt.titlesize)

    return nothing
end

function csm_plot_freeflyer_alltrjs(ax,m,N,X_all,U_all,fmt::CSMPlotFmt)
    x_all = prob.all_trj.state
    circle = fmt.circle
    col1 = fmt.col.cyan
    col2 = fmt.col.magenta
    cols = zeros(3,prob.solved)
    for i = 1:3
        cols[i,:] = LinRange(col1[i],col2[i],prob.solved)
    end

    # plot obstacles
    for i = 1:pars.obsN
        H = I(2)/pars.obsiH[1:2,1:2,i]
        c = pars.obsC[1:2,i]

        obs = H * circle .+ c
        ax.plot(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.8,
                linewidth=1,linestyle="-")
        if i==1
            ax.fill_between(obs[1,:],obs[2,:],
                    color=fmt.col.red,alpha=0.1,
                    linewidth=1,linestyle="-")
            else
            ax.fill_between(obs[1,:],obs[2,:],
                    color=fmt.col.red,alpha=0.1,
                    linewidth=1,linestyle="-")
        end
    end
    # plot keepout zones
    for i = 1:pars.kozN
        obs = pars.koz[i]
        center = obs.center[1:2]
        widths = 2.0*[obs.dx;obs.dy]
        if i < 13
            # ax = plt_rectangle(ax, center, widths, color=fmt.col.red, alpha=0.1, label="None")
        else
            ax = plt_rectangle(ax, center, widths, color=fmt.col.green, alpha=0.1, label="None")
        end
    end

    # plot discrete solutions
    for iter = 1:prob.solved
        ax.plot(x_all[1,:,iter],x_all[2,:,iter],
                label="Iteration $(iter)",
                marker="o",color=cols[:,iter],linestyle="-",
                markersize=fmt.markersize)
    end

    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(5,13)
    ax.set_ylim(-3,8)
    ax.set_xlabel(L"x_{\mathcal{I}}\ [m]",fontsize=fmt.fontsize)
    ax.set_ylabel(L"y_{\mathcal{I}}\ [m]",fontsize=fmt.fontsize)
    ax.legend(fontsize=12,loc="upper left")
    ax.grid(alpha=fmt.gridalpha)
    ax.set_title(L"All\ FreeFlyer\ Trajectories",fontsize=fmt.titlesize)

    return nothing
end

function csm_freeflyer_sd(m,N,X_all,U_all,fmt::CSMPlotFmt,X)
    id_r = prob.pars.mdl_pars.id_r
    x  = prob.new_sol.state
    tf = prob.new_sol.tf
    N  = prob.pars.N
    Ni = size(X,2)
    T  = LinRange(0,tf,N)
    Ti = LinRange(0,tf,Ni)

    kozN = prob.pars.mdl_pars.kozN
    koz  = prob.pars.mdl_pars.koz

    # discrete
    sd = zeros(N)
    for k = 1:N
        temp = zeros(kozN)
        rk = x[id_r,k]
        for i = 1:kozN
            temp[i], = signed_distance(rk,koz[i])
        end
        sd[k] = minimum(temp)
    end

    # continuous
    SD = zeros(Ni)
    for k = 1:Ni
        temp = zeros(kozN)
        rk = X[id_r,k]
        for i = 1:kozN
            temp[i], = signed_distance(rk,koz[i])
        end
        SD[k] = minimum(temp)
    end

    fig = plt.figure(figsize=fmt.figsize)
    ax  = plt.gca()

    ax.plot([0;tf],[0;0],
        color=fmt.col.red,
        linestyle="--",
        linewidth=fmt.lw)
    ax.plot(Ti,SD,
        color=fmt.col.blue,
        linestyle="-",
        linewidth=fmt.lw)
    ax.plot(T,sd,
        color=fmt.col.blue,
        marker="o",
        markersize=fmt.markersize,
        linestyle="none")
    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    plt.xlim(0,tf)
    # plt.ylim(0,70)
    plt.grid(alpha=fmt.gridalpha)
    plt.title("Signed Distance",fontsize=fmt.titlesize)
    plt.xlabel("Time [s]",fontsize=fmt.fontsize)
    plt.ylabel(L"\min_i\;d_i(r)",fontsize=fmt.fontsize)
    plt.tight_layout()
    plt.show()
    return nothing
end
