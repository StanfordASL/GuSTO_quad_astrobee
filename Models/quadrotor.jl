# ------------------------------------------------------------------
# -   Model for the quadrotor example - Riccardo Bonalli 12/2019   -
# ------------------------------------------------------------------



export Quadrotor



# Model quadrotor as a Julia class

mutable struct Quadrotor
    # State (r, v) and control (u, Γ) dimensions
    x_dim
    u_dim

    # Dynamics and linearized dynamics
    f
    A
    B

    # Model constants
    gravity
    θMax

    # Problem settings
    dimLinearConstraintsU
    dimSecondOrderConeConstraintsU
    x_init
    x_final
    tf
    xMin
    xMax
    uMin
    uMax

    # Cylindrical obstacles (modeled by a center (x,y) and a radius r) and polygon obstacles (not used in this example)
    obstacles
    poly_obstacles

    # GuSTO parameters
    Delta0
    omega0
    omegamax
    epsilon
    rho0
    rho1
    beta_succ
    beta_fail
    gamma_fail
    convergence_threshold
end



# The problem is set in the class constructor

function Quadrotor()
    x_dim = 6
    u_dim = 4

    gravity = [0;0;9.81]
    θMax = pi/3.0

    dimLinearConstraintsU = 9
    dimSecondOrderConeConstraintsU = 0
    x_init  = [0;0;0 ; 0;0;0]
    x_final = [2.5;6.;0 ; 0;0;0]
    tf = 2.
    myInf = 1.0e6 # Adopted to detect initial and final condition-free state variables
    xMin = [-0.1;-0.1;-myInf ; -myInf;-myInf;-myInf]
    xMax = [4;7;myInf ; myInf;myInf;myInf]
    uMin = 0.6
    uMax = 23.2

    Delta0 = 100.
    omega0 = 1.
    omegamax = 1.0e9
    epsilon = 1.0e-3
    rho0 = 10.0
    rho1 = 20.0
    beta_succ = 2.
    beta_fail = 0.5
    gamma_fail = 5.
    convergence_threshold = 5.0


    # Cylindrical obstacles in the form [(x,y),r]
    obstacles = []
    obs = [[1.0,2.0],1.0/2.5]
    push!(obstacles, obs)
    obs = [[2.0,5.0],1.0/2.5]
    push!(obstacles, obs)

    # Polygonal obstacles are not used in this example
    poly_obstacles = []

    Quadrotor(x_dim, u_dim,
             [], [], [],
             gravity, θMax,
             dimLinearConstraintsU, dimSecondOrderConeConstraintsU, x_init, x_final, tf, xMin, xMax, uMin, uMax,
             obstacles,
             poly_obstacles,
             Delta0,
             omega0,
             omegamax,
             epsilon,
             rho0,
             rho1,
             beta_succ,
             beta_fail,
             gamma_fail,
             convergence_threshold)
end



# Method that returns the GuSTO parameters (used for set up)

function get_initial_gusto_parameters(m::Quadrotor)
    return m.Delta0, m.omega0, m.omegamax, m.epsilon, m.rho0, m.rho1, m.beta_succ, m.beta_fail, m.gamma_fail, m.convergence_threshold
end



# GuSTO is intialized by zero controls and velocities, and a straight-line in position

function initialize_trajectory(model::Quadrotor, N::Int)
  x_dim,  u_dim   = model.x_dim, model.u_dim
  x_init, x_final = model.x_init, model.x_final
  dimR  = Int(model.x_dim/2)
  
  X = zeros(x_dim, N)
  x0 = x_init[1:dimR]
  x1 = x_final[1:dimR]
  X[1:dimR,:] = hcat(range(x0, stop=x1, length=N)...)
  for i = 1:N
    X[dimR+1:2*dimR,i] = zeros(dimR)
  end
  U = zeros(u_dim, N-1)

  return X, U
end



# Method that returns the convergence ratio between iterations (in percentage)
# The quantities X, U denote the actual solution over time, whereas Xp, Up denote the solution at the previous step over time

function convergence_metric(model::Quadrotor, X, U, Xp, Up)
    x_dim = model.x_dim
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

function true_cost(model::Quadrotor, X, U, Xp, Up)
    x_dim, u_dim = model.x_dim, model.u_dim
    cost = 0.

    for k = 1:length(U[1,:])
        cost += U[u_dim,k]^2 # This corresponds to ∫ Γ(t)^2 dt
    end

    return cost
end



# The following methods return the i-th coordinate at the k-th iteration of the various constraints and their linearized versions (when needed)
# These are returned in the form " g(t,x(t),u(t)) <= 0 "



# Method that gathers all the lienar control constraints

function control_linear_constraints(model::Quadrotor, X, U, Xp, Up, k, i)
    x_dim, u_dim = model.x_dim, model.u_dim
    uMin, uMax, θMax = model.uMin, model.uMax, model.θMax

    # Control bounds on Γ
    if i == 1
      return uMin - U[u_dim,k]
    elseif i == 2
      return U[u_dim,k] - uMax

    # Directional control constraint
    elseif i == 3
      return U[u_dim,k]*cos(θMax) - U[u_dim-1,k]

    # Control bounds on u via Γ
    elseif i == 4
      return U[1,k] - U[u_dim,k]
    elseif i == 5
      return  -U[u_dim,k] - U[1,k]
    elseif i == 6
      return U[2,k] - U[u_dim,k]
    elseif i == 7
      return  -U[u_dim,k] - U[2,k]
    elseif i == 8
      return U[3,k] - U[u_dim,k]
    elseif i == 9
      return  -U[u_dim,k] - U[3,k]
    else
      println("ERROR - TOO MANY LINEAR CONTROL CONSTRAINTS")
    end
end



# State bounds and trust-region constraints (these are all convex constraints)

function state_max_convex_constraints(model::Quadrotor, X, U, Xp, Up, k, i)
    return ( X[i, k] - model.xMax[i] )
end



function state_min_convex_constraints(model::Quadrotor, X, U, Xp, Up, k, i)
    return ( model.xMin[i] - X[i, k] )
end



function trust_region_max_constraints(model::Quadrotor, X, U, Xp, Up, k, i, Delta)
    return ( (X[i, k] - Xp[i, k]) - Delta )
end



function trust_region_min_constraints(model::Quadrotor, X, U, Xp, Up, k, i, Delta)
    return ( -Delta - (X[i, k] - Xp[i, k]) )
end



# Method that checks whether trust-region constraints are satisifed or not (recall that trust-region constraints are penalized)

function is_in_trust_region(model::Quadrotor, X, U, Xp, Up, Delta)
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

function state_initial_constraints(model::Quadrotor, X, U, Xp, Up)
    return (X[:,1] - model.x_init)
end



function state_final_constraints(model::Quadrotor, X, U, Xp, Up)
    return (X[:,end] - model.x_final)
end



# Methods that return the cylindrical obstacle-avoidance constraint and its lienarized version
# Here, a merely classical distance function is considered

function obstacle_constraint(model::Quadrotor, X, U, Xp, Up, k, obs_i)
    dimR  = Int(model.x_dim/2)
    p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
    p_k  = X[1:dimR-1,k]
    
    dist = norm(p_k - p_obs, 2)
    constraint = -( dist - obs_radius )

    return constraint
end



function obstacle_constraint_convexified(model::Quadrotor, X, U, Xp, Up, k, obs_i)
    dimR  = Int(model.x_dim/2)
    p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
    p_k  = X[1:dimR-1,k]
    p_kp = Xp[1:dimR-1,k]

    dist_prev = norm(p_kp - p_obs, 2)
    dir_prev = (p_kp - p_obs)/dist_prev
    constraint = -( dist_prev - obs_radius + sum(dir_prev[i] * (p_k[i] - p_kp[i]) for i=1:dimR-1) )

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

function f_dyn(x::Vector, u::Vector, model::Quadrotor)
  x_dim, u_dim = model.x_dim, model.u_dim
  dimR = Int(model.x_dim/2)

  g = model.gravity
  f = zeros(x_dim)

  r, v = x[1:dimR], x[dimR+1:2*dimR]
  w, Γ = u[1:u_dim-1], u[u_dim]

  f[1:dimR] = v
  f[dimR+1:2*dimR] = w - g

  return f
end



function A_dyn(x::Vector, u::Vector, model::Quadrotor)
  x_dim, u_dim = model.x_dim, model.u_dim
  dimR = Int(model.x_dim/2)
  A = zeros(x_dim,x_dim)

  for i = 1:dimR
    A[i,dimR+i] = 1.0
  end

  return A
end



function B_dyn(x::Vector, u::Vector, model::Quadrotor)
  x_dim, u_dim = model.x_dim, model.u_dim
  dimR = Int(model.x_dim/2)
  B = zeros(x_dim,u_dim)

  for i = 1:dimR
    B[dimR+i,i] = 1.0
  end

  return B
end