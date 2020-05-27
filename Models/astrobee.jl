# ---------------------------------------------------------------------------------
# -   Model for the freeflyer example - Thomas Lew and Riccardo Bonalli 12/2019   -
# ---------------------------------------------------------------------------------

include("./polygonal_obstacles.jl")
include("./ISS.jl")

# Model freeflyer as a Julia class

export Astrobee

mutable struct Astrobee
    # State (r, v, q, ω) and control (F, T) dimensions
    x_dim
    u_dim

    # Dynamics and linearized dynamics
    f
    A
    B

    # Model constants
    model_radius # Radius of the sphere that contains the robot
    mass
    J
    Jinv

    # Problem settings
    dimLinearConstraintsU
    dimSecondOrderConeConstraintsU
    x_init
    x_final
    tf
    xMax
    xMin
    uMax
    uMin

    # Spherical obstacles (modeled by a center (x,y,z) and a radius r) and polygonal obstacles (ISS)
    obstacles
    poly_obstacles

    # GuSTO parameters
    Delta0
    omega0
    omegamax
    # threshold for constraints satisfaction : constraints <= epsilon
    epsilon
    epsilon_xf_constraint
    rho0
    rho1
    beta_succ
    beta_fail
    gamma_fail
    convergence_threshold
end

function Astrobee()
    x_dim = 13
    u_dim = 6

    model_radius = sqrt.(3)*0.5*0.05 # Each side of cubic robot is 5.0cm & inflate to sphere
    mass = 7.2
    J_norm = 0.1083
    J = J_norm*Matrix(1.0I,3,3)
    Jinv = inv(J)

    # problem 
    dimLinearConstraintsU = 0
    dimSecondOrderConeConstraintsU = 2
    s13     = sqrt(1. / 3.)
    x_init  = [ 7.2,-0.4,5.0,  3.5e-2,3.5e-2,1e-4,  s13,0.,s13,s13,     0,0,0]
    x_final = [11.3,6.0,4.5,  1e-4,1e-4,1e-4,  -0.5,0.5,-0.5,0.5,  0,0,0]
    tf = 100.
    myInf = 1.0e6 # Adopted to detect initial and final condition-free state variables
    xMax = [100.;100.;100. ; 0.4;0.4;0.4 ; 100.;100.;100.;100. ; 1.;1.;1]
    xMin = -xMax

    # constraints / limits
    hard_limit_accel = 0.01          # m/s^2
    hard_limit_alpha = 1. *3.14/180. # rad/s^2
    uMax = Array([mass*hard_limit_accel,mass*hard_limit_accel,mass*hard_limit_accel,  J_norm*hard_limit_alpha,J_norm*hard_limit_alpha,J_norm*hard_limit_alpha])
    uMin = -uMax


    # GuSTO Parameters
    Delta0 = 100.
    omega0 = 100.
    omegamax = 1.0e9
    epsilon = 1.0e-3
    epsilon_xf_constraint = 1e-4
    rho0 = 10.
    rho1 = 20.
    beta_succ = 2.
    beta_fail = 0.5
    gamma_fail = 5.
    convergence_threshold = 2.5


    # Spherical obstacles in the form [(x,y),r]
    obstacles = []
    obs = [[11.3,3.8,4.8], 0.3]
    push!(obstacles, obs)
    obs = [[8.5,-0.04,5.01], 0.3]
    push!(obstacles, obs)
    obs = [[11.2,1.84,5.01], 0.3]
    push!(obstacles, obs)

    # Polygonal obstacles
    print("Initializing the ISS.")
    keepin_zones, keepout_zones = get_ISS_zones()
    poly_obstacles = keepout_zones

    Astrobee(x_dim, u_dim,
               [], [], [],
               model_radius, mass, J, Jinv,
               dimLinearConstraintsU, dimSecondOrderConeConstraintsU, x_init, x_final, tf, 
               xMax, xMin, uMax, uMin,
               obstacles, poly_obstacles,
               Delta0, omega0, omegamax,
               epsilon, epsilon_xf_constraint,
               rho0, rho1,
               beta_succ, beta_fail,
               gamma_fail,
               convergence_threshold)
end



# Method that returns the GuSTO parameters (used for set up)

function get_initial_gusto_parameters(m::Astrobee)
    return m.Delta0, m.omega0, m.omegamax, m.epsilon, m.rho0, m.rho1, m.beta_succ, m.beta_fail, m.gamma_fail, m.convergence_threshold
end



# GuSTO is intialized by zero controls, and a straight-line in the state space

function initialize_trajectory(model::Astrobee, N::Int)
  x_dim,  u_dim   = model.x_dim, model.u_dim
  x_init, x_final = model.x_init, model.x_final
  
  X = hcat(range(x_init, stop=x_final, length=N)...)
  U = zeros(u_dim, N-1)

  return X, U
end



# Method that returns the convergence ratio between iterations (in percentage)
# The quantities X, U denote the actual solution over time, whereas Xp, Up denote the solution at the previous step over time

function convergence_metric(model::Astrobee, X, U, Xp, Up)
	x_dim = 3
    N = length(X[1,:])

    # Normalized maximum relative error between iterations
    max_num, max_den = -Inf, -Inf
    for k in 1:N
        val = norm(X[1:x_dim,k] - Xp[1:x_dim,k])
        max_num = val > max_num ? val : max_num

        val = norm(X[1:x_dim,k])
        max_den = val > max_den ? val : max_den
    end

    # Returning percentage error
    return max_num*100.0/max_den
end



# Method that returns the original cost

function true_cost(model::Astrobee, X, U, Xp, Up)
	x_dim, u_dim = model.x_dim, model.u_dim
    cost = 0.

    for k = 1:length(U[1,:])
        cost += sum(U[i,k]^2 for i = 1:u_dim) # This corresponds to ∫ ( || F(t) ||^2 + || T(t) ||^2 ) dt
    end

    return cost
end



# The following methods return the i-th coordinate at the k-th iteration of the various constraints and their linearized versions (when needed)
# These are returned in the form " g(t,x(t),u(t)) <= 0 "

# Method that gathers all the linear control constraints

function control_linear_constraints(model::Astrobee, X, U, Xp, Up, k, i)
    x_dim, u_dim = model.x_dim, model.u_dim
    uMin, uMax = model.uMin, model.uMax

    # Control bounds on Γ
    if i==1
      norm(U[1:3,k]) - uMax[1]
    elseif i==2
      norm(U[4:6,k]) - uMax[4]
    # if i >= 1 && i<=6
    elseif i >= 1 && i<=6
      # uMin <= u
      return uMin[i] - U[i,k]
    elseif i >= 7 && i<=12
      # u <= uMax
      return U[i-6,k] - uMax[i-6]
    else
      println("[astrobee.jl::control_linear_constraints] ERROR - too many constraints.")
    end
end

function control_second_order_cone_constraints(model::Astrobee, X, U, Xp, Up, k, i)
    u_dim = model.u_dim
    uMin, uMax = model.uMin, model.uMax

    if i == 1
      return (uMax[1], U[1:3,k])
    elseif i==2
      return (uMax[4], U[4:6,k])
    else
      println("[astrobee.jl::control_second_order_cone_constraints] ERROR - too many constraints.")
    end
end

# State bounds and trust-region constraints (these are all convex constraints)

function state_max_convex_constraints(model::Astrobee, X, U, Xp, Up, k, i)
    return ( X[i, k] - model.xMax[i] )
end



function state_min_convex_constraints(model::Astrobee, X, U, Xp, Up, k, i)
    return ( model.xMin[i] - X[i, k] )
end



function trust_region_max_constraints(model::Astrobee, X, U, Xp, Up, k, i, Delta)
    return ( (X[i, k] - Xp[i, k]) - Delta )
end



function trust_region_min_constraints(model::Astrobee, X, U, Xp, Up, k, i, Delta)
    return ( -Delta - (X[i, k] - Xp[i, k]) )
end



# Method that checks whether trust-region constraints are satisifed or not (recall that trust-region constraints are penalized)

function is_in_trust_region(model::Astrobee, X, U, Xp, Up, Delta)
    B_is_inside = true

    for k = 1:length(X[1,:])
        for i = 1:model.x_dim
            if trust_region_max_constraints(model, X, U, Xp, Up, k, i, Delta) > 0.
                B_is_inside = false
            end
            if trust_region_min_constraints(model, X, U, Xp, Up, k, i, Delta) > 0.
                B_is_inside = false
            end
        end
    end

    return B_is_inside
end



# Initial and final conditions on state variables

function state_initial_constraints(model::Astrobee, X, U, Xp, Up)
    return ( X[:,1] - model.x_init )
end



function state_final_constraints(model::Astrobee, X, U, Xp, Up)
    return ( X[:,end] - model.x_final )
end



# Methods that return the cylindrical obstacle-avoidance constraint and its lienarized version
# Here, a merely classical distance function is considered

function obstacle_constraint(model::Astrobee, X, U, Xp, Up, k, obs_i,
                                                 obs_type::String="sphere")
    #  obs_type    : Type of obstacles, can be 'sphere' or 'poly'
    if obs_type=="sphere"
      p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
      bot_radius        = model.model_radius
      total_radius      = obs_radius + bot_radius
      p_k  = X[1:3, k]
      
      dist = norm(p_k - p_obs, 2)
      constraint = -( dist - total_radius )
    elseif obs_type=="poly"
      obs            = model.poly_obstacles[obs_i]
      p_k            = X[1:3, k]
      dist_prev, pos = signed_distance_with_closest_point_on_surface(p_k, obs)
      constraint = -( dist_prev - model.model_radius )
    else
      print("[astrobee.jl::obstacle_constraint_convexified] Unknown obstacle type.")
    end

    return constraint
end



function obstacle_constraint_convexified(model::Astrobee, X, U, Xp, Up, k, obs_i,
                                                 obs_type::String="sphere")
    #  obs_type    : Type of obstacles, can be 'sphere' or 'poly'
    if obs_type=="sphere"
      p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
      bot_radius        = model.model_radius
      total_radius      = obs_radius + bot_radius
      p_k, p_kp         = X[1:3, k], Xp[1:3, k]
      
      dist_prev = norm(p_kp - p_obs, 2)
      n_prev    = (p_kp-p_obs) / dist_prev
      constraint = -( dist_prev - total_radius + sum(n_prev[i] * (p_k[i]-p_kp[i]) for i=1:3) )
    elseif obs_type=="poly"
      obs            = model.poly_obstacles[obs_i]
      p_k, p_kp      = X[1:3, k], Xp[1:3, k]
      dist_prev, pos = signed_distance_with_closest_point_on_surface(p_kp, obs)

      n_prev = (p_kp-obs.c[1:3]) / norm((p_kp-obs.c[1:3]),2)
      constraint = -( dist_prev - model.model_radius + sum(n_prev[i] * (p_k[i]-p_kp[i]) for i=1:3) )
    else
      print("[astrobee.jl::obstacle_constraint_convexified] Unknown obstacle type.")
    end
    
    return constraint
end



# The following methods return the dynamical constraints and their linearized versions
# These are returned as time-discretized versions of the constraints " x' - f(x,u) = 0 " or " x' - A(t)*(x - xp) - B(t)*(u - up) = 0 ", respectively



# Method that returns the dynamical constraints and their linearized versions all at once

function compute_dynamics(model, Xp, Up)
    N = length(Xp[1,:])

    f_all, A_all, B_all = [], [], []

    for k in 1:N-1
        x_k = Xp[:,k]
        u_k = Up[:,k]

        f_dyn_k, A_dyn_k, B_dyn_k = f_dyn(x_k, u_k, model), A_dyn(x_k, u_k, model), B_dyn(x_k, u_k, model)

        push!(f_all, f_dyn_k)
        push!(A_all, A_dyn_k)
        push!(B_all, B_dyn_k)
    end

    return f_all, A_all, B_all
end



# These methods return the dynamics and its linearizations with respect to the state (matrix A(t)) and the control (matrix B(t)), respectively

function f_dyn(x::Vector, u::Vector, model::Astrobee)
  x_dim = model.x_dim
  f = zeros(x_dim)

  r, v, ω = x[1:3], x[4:6], x[11:13]
  qw, qx, qy, qz = x[7:10]
  ωx, ωy, ωz = x[11:13]
  F, M = u[1:3], u[4:6]

  f[1:3] = v
  f[4:6] = F/model.mass

  f[7]  = 1/2*(-ωx*qx - ωy*qy - ωz*qz)
  f[8]  = 1/2*( ωx*qw - ωz*qy + ωy*qz)
  f[9]  = 1/2*( ωy*qw + ωz*qx - ωx*qz)
  f[10] = 1/2*( ωz*qw - ωy*qx + ωx*qy)
  f[11:13] = model.Jinv*(M - cross(ω,model.J*ω))

  return f
end



function A_dyn(x::Vector, u::Vector, model::Astrobee)
  x_dim = model.x_dim
  A = zeros(x_dim, x_dim)

  A[1:6,1:6] = kron([0 1; 0 0], Matrix(1.0I,3,3))

  Jxx, Jyy, Jzz = diag(model.J)
  qw, qx, qy, qz = x[7:10] 
  ωx, ωy, ωz = x[11:13]

  A[7,8] = -ωx/2
  A[7,9] = -ωy/2
  A[7,10] = -ωz/2
  A[7,11] = -qx/2
  A[7,12] = -qy/2
  A[7,13] = -qz/2

  A[8,7]  = ωx/2
  A[8,9]  = -ωz/2
  A[8,10] = ωy/2
  A[8,11] = qw/2
  A[8,12] = qz/2
  A[8,13] = -qy/2

  A[9,7]  = ωy/2
  A[9,8]  = ωz/2
  A[9,10] = -ωx/2
  A[9,11] = -qz/2
  A[9,12] = qw/2
  A[9,13] = qx/2

  A[10,7]  = ωz/2
  A[10,8]  = -ωy/2
  A[10,9]  = ωx/2
  A[10,11] = qy/2
  A[10,12] = -qx/2
  A[10,13] = qw/2

  A[11,12] =  (Jyy-Jzz)*ωz/Jxx
  A[11,13] =  (Jyy-Jzz)*ωy/Jxx
  A[12,11] = -(Jxx-Jzz)*ωz/Jyy
  A[12,13] = -(Jxx-Jzz)*ωx/Jyy
  A[13,11] =  (Jxx-Jyy)*ωy/Jzz
  A[13,12] =  (Jxx-Jyy)*ωx/Jzz

  return A
end



function B_dyn(x::Vector, u::Vector, model::Astrobee)
  x_dim, u_dim = model.x_dim, model.u_dim

  B = zeros(x_dim, u_dim)

  B[4:6,1:3] = Matrix(1.0I,3,3)/model.mass
  B[11:13,4:6] = model.Jinv

  return B
end