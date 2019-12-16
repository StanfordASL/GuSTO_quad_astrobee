# --------------------------------------------------------------------------
# -   GuSTO solver environment - Thomas Lew and Riccardo Bonalli 12/2019   -
# --------------------------------------------------------------------------



export GuSTOProblem



# GuSTO solver as a Julia class

mutable struct GuSTOProblem
    # Number of time-discretization steps and step-size, respectively
    N
    dt

    # Penalization weight ω and trsut-region constraints radius Δ, respectively
    omega
    Delta

    # Number of linear control constraints and second-order cone control constraints, respectively
    dimLinearConstraintsU
    dimSecondOrderConeConstraintsU

    # Model class
    solver_model

    # Current trajectory X and control U
    X
    U

    # Trajectory X and control U at the previous iteration
    Xp
    Up

    # The intial constraints are defined in the model class (used for set up)
    initial_constraint
end



# Standard constructor

function GuSTOProblem(model, N, dimLinearConstraintsU, dimSecondOrderConeConstraintsU, Xp, Up, solver=Ipopt.Optimizer)
    N     = N
    dt    = model.tf / (N-1)
    omega = model.omega0
    Delta = model.Delta0

    solver_model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    X = @variable(solver_model, X[1:model.x_dim,1:N  ])
    U = @variable(solver_model, U[1:model.u_dim,1:N-1])

    GuSTOProblem(N, dt,
                 omega, Delta, dimLinearConstraintsU, dimSecondOrderConeConstraintsU,
                 solver_model,
                 X, U, Xp, Up,
                 [])
end



# Methods that define the convex subproblem at each new SCP iteration

function reset_problem(scp_problem::GuSTOProblem, model, solver=Ipopt.Optimizer)
    scp_problem.solver_model = Model(with_optimizer(solver, print_level=0))
    N = scp_problem.N
    X = @variable(scp_problem.solver_model, X[1:model.x_dim,1:N  ])
    U = @variable(scp_problem.solver_model, U[1:model.u_dim,1:N-1])
    scp_problem.X       = X
    scp_problem.U       = U
end



function set_parameters(scp_problem::GuSTOProblem, model,
                        Xp, Up, omega, Delta)
    scp_problem.Xp = Xp
    scp_problem.Up = Up
    scp_problem.omega = omega
    scp_problem.Delta = Delta
end



# Method that compute the accuracy ration ρ (in percentage)

function accuracy_ratio(problem::GuSTOProblem, model, X, U, Xp, Up)
    N = length(X[1,:])

    num = 0.0
    den = 0.0

    # Contribution of the dynamics
    for k in 1:N-1
        x_k  = X[:, k]
        u_k  = U[:, k]
        x_kp = Xp[:,k]
        u_kp = Up[:,k]

        f_k  = f_dyn(x_k,  u_k,  model)
        A_k  = A_dyn(x_k,  u_k,  model)
        B_k  = B_dyn(x_k,  u_k,  model)
        f_kp = f_dyn(x_kp, u_kp, model)
        A_kp = A_dyn(x_kp, u_kp, model)
        B_kp = B_dyn(x_kp, u_kp, model)

        linearized = f_kp + A_kp*(x_k-x_kp) + B_kp*(u_k-u_kp)
        num += norm(f_k - linearized, 2)
        den += norm(linearized,       2)
    end

    # Contribution of the obstacles
    for k = 1:N
        for i = 1:length(model.obstacles)
            constraint            = obstacle_constraint(            model, X, U, Xp, Up, k, i)
            constraint_linearized = obstacle_constraint_convexified(model, X, U, Xp, Up, k, i)
            
            num += norm(constraint - constraint_linearized, 2)
            den += norm(constraint_linearized,              2)
        end
    end

    # Returning percentage ratio
    return num*100.0/den
end



# The following methods define the convex subproblem at each iteration in the "JuMP" framework



# Methods that define the cost

function define_cost(scp_problem::GuSTOProblem, model)
    total_cost =  add_true_cost(scp_problem, model)
    total_cost += add_penalties(scp_problem, model)
    @objective(scp_problem.solver_model, Min, total_cost)
end

function add_true_cost(scp_problem::GuSTOProblem, model)
    solver_model = scp_problem.solver_model
    X, U, Xp, Up = scp_problem.X, scp_problem.U, scp_problem.Xp, scp_problem.Up
    return true_cost(model, X, U, Xp, Up)
end



# This method adds to the cost the linearized state constraints as penalizations
# To simplify the formalism, penalizations are derived by introducing slack control variables λ, so that every scalar linear constraint g(t,x) <= 0 is rather penalized as
#
# ... + ∫ ω*( g(t,x) - λ )^2 dt
#
# with the additional linear control constraint " λ <= 0 "

function add_penalties(scp_problem::GuSTOProblem, model)
    solver_model = scp_problem.solver_model
    x_dim, u_dim = model.x_dim, model.u_dim
    omega, Delta = scp_problem.omega, scp_problem.Delta

    X, U, Xp, Up = scp_problem.X, scp_problem.U, scp_problem.Xp, scp_problem.Up
    N, dt        = length(X[1,:]), scp_problem.dt

    penalization = 0.
    
    # Contribution of trust-region constraints

    @variable(solver_model, lambdas_trust_max_convex_constraints[i=1:x_dim, k=1:N])
    @variable(solver_model, lambdas_trust_min_convex_constraints[i=1:x_dim, k=1:N])
    for k = 1:N
        for i = 1:x_dim
            lambda_max     = lambdas_trust_max_convex_constraints[i,k]
            constraint_max = trust_region_max_constraints(model, X, U, Xp, Up, k, i, Delta)
            lambda_min     = lambdas_trust_min_convex_constraints[i,k]
            constraint_min = trust_region_min_constraints(model, X, U, Xp, Up, k, i, Delta)

            @constraint(solver_model, lambda_max <= 0.)
            penalization += omega*(constraint_max-lambda_max)^2
            @constraint(solver_model, lambda_min <= 0.)
            penalization += omega*(constraint_min-lambda_min)^2
        end
    end
    
    # Contribution of obstacle-avoidance constraints

    # Non-polygonal obstacles

    Nb_obstacles = length(model.obstacles)
    if Nb_obstacles > 0
        @variable(solver_model, lambdas_obstacles[i=1:Nb_obstacles, k=1:N])
        for k = 1:N
            for i = 1:Nb_obstacles
                lambda     = lambdas_obstacles[i,k]
                constraint = obstacle_constraint_convexified(model, X, U, Xp, Up, k, i)

                @constraint(solver_model, lambda <= 0.)
                penalization += omega*(constraint-lambda)^2
            end
        end
    end

    # Polygonal obstacles

    Nb_poly_obstacles = length(model.poly_obstacles)
    if Nb_poly_obstacles > 0
        @variable(solver_model, lambdas_poly_obstacles[i=1:Nb_poly_obstacles, k=1:N])
        for k = 1:N
            for i = 1:Nb_poly_obstacles
                lambda     = lambdas_poly_obstacles[i,k]
                constraint = poly_obstacle_constraint_convexified(model, X, U, Xp, Up, k, i)

                @constraint(solver_model, lambda <= 0.)
                penalization += omega*(constraint-lambda)^2
            end
        end
    end

    # Contribution of state constraints that are different from trust-region constraints and obstacle-avoidance constraints

    @variable(solver_model, lambdas_state_max_convex_constraints[i=1:x_dim, k=1:N])
    @variable(solver_model, lambdas_state_min_convex_constraints[i=1:x_dim, k=1:N])
    for k = 1:N
        for i = 1:x_dim
            lambda_max     = lambdas_state_max_convex_constraints[i,k]
            constraint_max = state_max_convex_constraints(model, X, U, Xp, Up, k, i)
            lambda_min     = lambdas_state_min_convex_constraints[i,k]
            constraint_min = state_min_convex_constraints(model, X, U, Xp, Up, k, i)

            @constraint(solver_model, lambda_max <= 0.)
            penalization += omega*(constraint_max-lambda_max)^2
            @constraint(solver_model, lambda_min <= 0.)
            penalization += omega*(constraint_min-lambda_min)^2
        end
    end

    return penalization
end



# Method that checks whether penalized state constraints are hardly satisfied (up to the threshold ε)

function satisfies_state_inequality_constraints(scp_problem::GuSTOProblem, model, X, U, Xp, Up, Delta)
    B_satisfies_constraints = true
    x_dim = model.x_dim
    N = scp_problem.N
    epsilon = model.epsilon

    # Contribution of trust-region constraints

    for k = 1:N
        for i = 1:x_dim
            constraint_max = trust_region_max_constraints(model, X, U, Xp, Up, k, i, Delta)
            constraint_min = trust_region_min_constraints(model, X, U, Xp, Up, k, i, Delta)
            if constraint_max > epsilon || constraint_min > epsilon
                B_satisfies_constraints = false
            end
        end
    end

    # Contribution of obstacle-avoidance constraints

    for k = 1:N

        # Non-polygonal obstacles

        Nb_obstacles = length(model.obstacles)
        if Nb_obstacles > 0
            for i = 1:Nb_obstacles
                constraint = obstacle_constraint(model, X, U, [], [], k, i)
                if constraint > epsilon
                    B_satisfies_constraints = false
                end
            end
        end

        # Polygonal obstacles

        Nb_poly_obstacles = length(model.poly_obstacles)
        if Nb_poly_obstacles > 0
            for i = 1:Nb_poly_obstacles
                constraint = poly_obstacle_constraint(model, X, U, [], [], k, i)
                if constraint > epsilon
                    B_satisfies_constraints = false
                end
            end
        end
    end

    # Contribution of state constraints that are different from trust-region constraints and obstacle-avoidance constraints

    for k = 1:N
        for i = 1:x_dim
            constraint_max = state_max_convex_constraints(model, X, U, [], [], k, i)
            constraint_min = state_min_convex_constraints(model, X, U, [], [], k, i)
            if constraint_max > epsilon || constraint_min > epsilon
                B_satisfies_constraints = false
            end
        end
    end

    return B_satisfies_constraints
end



# These methods add the remaining constraints, such as linearized dyamical, intial/final conditions constraints and control constraints



function define_constraints(scp_problem::GuSTOProblem, model)
    add_initial_constraints(scp_problem, model)
    add_final_constraints(scp_problem, model)
    add_dynamics_constraints(scp_problem, model)
    add_control_constraints(scp_problem, model)
end



function add_initial_constraints(scp_problem::GuSTOProblem, model)
    solver_model = scp_problem.solver_model
    x_dim, u_dim = model.x_dim, model.u_dim
    omega, Delta = scp_problem.omega, scp_problem.Delta

    X, U, Xp, Up = scp_problem.X, scp_problem.U, scp_problem.Xp, scp_problem.Up

    constraint = state_initial_constraints(model, X, U, Xp, Up)
    scp_problem.initial_constraint = @constraint(solver_model, constraint .== 0.)
end



# To improve robustness, final conditions are imposed up to the threshold ε > 0 (originally defined to check satisfaction of state constraints)

function add_final_constraints(scp_problem::GuSTOProblem, model)
    solver_model = scp_problem.solver_model
    x_dim, u_dim = model.x_dim, model.u_dim
    omega, Delta = scp_problem.omega, scp_problem.Delta
    epsilon = model.epsilon

    X, U, Xp, Up = scp_problem.X, scp_problem.U, scp_problem.Xp, scp_problem.Up

    constraint = state_final_constraints(model, X, U, Xp, Up)
    @constraint(solver_model,  constraint - epsilon .<= 0.)
    @constraint(solver_model, -constraint - epsilon .<= 0.)
end



function add_dynamics_constraints(scp_problem::GuSTOProblem, model)
    solver_model = scp_problem.solver_model
    x_dim, u_dim = model.x_dim, model.u_dim
    omega, Delta = scp_problem.omega, scp_problem.Delta

    X, U, Xp, Up = scp_problem.X, scp_problem.U, scp_problem.Xp, scp_problem.Up
    N, dt        = length(X[1,:]), scp_problem.dt

    for k = 1:N-1
        X_knext  = X[:, k+1]
        X_k      = X[:, k]
        U_k      = U[:, k]
        X_kp     = Xp[:, k]
        U_kp     = Up[:, k]
        f_dyn_kp = model.f[k]
        A_kp     = model.A[k]
        B_kp     = model.B[k]

        # Standard forward Euler integration method
        constraint =  X_knext - ( X_k + dt * (  f_dyn_kp + 
                                                A_kp * (X_k-X_kp) + 
                                                B_kp * (U_k-U_kp)
                                             )
                                )
        @constraint(solver_model, constraint .== 0.)
    end
end



function add_control_constraints(scp_problem::GuSTOProblem, model)
    solver_model = scp_problem.solver_model
    x_dim, u_dim = model.x_dim, model.u_dim
    dimLinearConstraintsU, dimSecondOrderConeConstraintsU = scp_problem.dimLinearConstraintsU, scp_problem.dimSecondOrderConeConstraintsU

    X, U, Xp, Up = scp_problem.X, scp_problem.U, scp_problem.Xp, scp_problem.Up
    N            = length(X[1,:])

    for k = 1:N-1

        # Linear constraints
        if dimLinearConstraintsU > 0
            for i = 1:dimLinearConstraintsU
                constraint = control_linear_constraints(model, X, U, Xp, Up, k, i)
                @constraint(solver_model, constraint <= 0.)
            end
        end

        # Second order cone constraints
        if dimSecondOrderConeConstraintsU > 0
            for i = 1:dimSecondOrderConeConstraintsU
                t, x = control_second_order_cone_constraints(model, X, U, Xp, Up, k, i)
                @constraint(solver_model, [t; x] in SecondOrderCone())
            end
        end
    end
end
