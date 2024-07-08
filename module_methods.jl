using LatinHypercubeSampling
using Plots
using StatsPlots
using JuMP
using Ipopt
using PowerModels
using PowerModelsAnnex
using Suppressor
using RCall
using Polyhedra
using LinearAlgebra
using CSV
using DataFrames
using PowerModelsSecurityConstrained
using Statistics 
using Distributions
using DelimitedFiles
R"library(Rcpp)"
solver_args = Dict{Symbol,Any}()
solver_args[:tol] = 1e-6
solver_args[:print_level] = 0
solver = Ipopt.Optimizer
PowerModels.silence()

const TOLERANCE = 1e-4;

function update_gen_network(network, new_gen)
    ind = 1
    for (i, gen) in network["gen"]
        if network["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && network["gen"][i]["pmax"] > 0.0
            gen["pg"] += new_gen[ind]
            if gen["pg"] < gen["pmin"]
                gen["pg"] = gen["pmin"]
            elseif gen["pg"] > gen["pmax"]
                gen["pg"] = gen["pmax"]
            end

            ind += 1
        end
    end 
end

function update_bus_network(network, active_gens_buses ,new_bus)
    ind = 1
    for i in active_gens_buses
        network["bus"]["$i"]["vm"] += new_bus[ind]
        if network["bus"]["$i"]["vm"] < network["bus"]["$i"]["vmin"]
            network["bus"]["$i"]["vm"] = network["bus"]["$i"]["vmin"]
        elseif network["bus"]["$i"]["vm"] > network["bus"]["$i"]["vmax"]
            network["bus"]["$i"]["vm"] = network["bus"]["$i"]["vmax"]
        end

        ind += 1
    
    end 
end

function get_SetPoint(network_data,active_gens_buses)
    Pg_SetPoint = Dict{String, Float64}()
    for (i, gen) in network_data["gen"]
        if network_data["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && network_data["gen"][i]["pmax"] > 0.0
            Pg_SetPoint[i] = gen["pg"]
            #push!(Pg_SetPoint, (i, gen["pg"] ))
        end
    end 
    return Pg_SetPoint
end

function get_SetPoint_V(network_data, active_gens_buses)
    Pg_SetPoint = Dict{Int, Float64}()
    for (bus_key, bus_data) in network_data["bus"]
        bus_idx = parse(Int, bus_key)
        if bus_idx in active_gens_buses
            Pg_SetPoint[bus_idx] = bus_data["vm"]
        end
    end
    return Pg_SetPoint
end

function update_SetPoint(SetPoint, new_value, nb)
    newSetPoint = deepcopy(SetPoint)
    newSetPoint[nb] = new_value
    return newSetPoint
end


function comment_line_in_file(file_path::String, target_line::String)
    # Read the file
    lines = readlines(file_path)
    
    # Process the lines
    modified_lines = String[]
    for line in lines
        # Check if the line starts with the target string
        if startswith(line, target_line)
            # Comment out the line by adding '%'
            push!(modified_lines, "%" * line)
        else
            # Keep the line as is
            push!(modified_lines, line)
        end
    end
    
    # Write the modified lines back to the file
    open(file_path, "w") do file
        for line in modified_lines
            println(file, line)
        end
    end
end


function instantiate_system_QCRM(data_model)
    pm = instantiate_model(data_model, QCRMPowerModel, PowerModels.build_opf)
    vars = []
    slack_gen_idx = get_slack_idx(pm)
    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
        if i != slack_gen_idx && data_model["gen"]["$i"]["pmax"] > 0.0
            push!(vars, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
        end
    end
    gen_indexes = unique(map(x -> x["gen_bus"], values(pm.data["gen"])))
    for g in gen_indexes
        push!(vars, JuMP.variable_by_name(pm.model, string("0_vm[",g,"]")))
    end
    N = length(vars)
    header = get_header_names(vars)
    return pm, N, vars, header
end


function create_pol(N::Integer)
    """Create N-dim unit cube."""
    rSTDO = R"""
    library(volesti)
    P <- gen_cube($N,'H')
    values_b <- as.numeric(cbind(matrix(0,1,$N), matrix(1,1,$N))) 
    slot(P,'b') <- values_b
    """ # mondify cube from [-1 1] to [0 1]
    return rcopy(R"slot(P, 'A')"),rcopy(R"slot(P, 'b')")
end

function comp_vol()
    R"""
    v <- volume(P, settings = list("error" = 10^-12))
    """
    return rcopy(R"v")
end

function return_P()
    return rcopy(R"P")
end

#R"""
#library(volesti)
#P = gen_rand_hpoly(20, 60, generator = list('constants' = 'sphere'))"""

function create_scaled_pol(N::Integer,L::AbstractArray,U::AbstractArray)
    """Create polytope with increased upper bounds."""
        rSTDO = R"""
        library(volesti)
        P <- gen_cube($N,'H')
        values_b <- as.numeric(cbind(matrix(-c(unlist($L)),1,$N), matrix(c(unlist($U)),1,$N)))
        slot(P,'b') <- values_b
        #Hpolytope$new(P$A, P$b)
        """ # mondify cube from [-1 1] to [0 1]
        return rcopy(R"slot(P, 'A')"),rcopy(R"slot(P, 'b')"),rcopy(R"P")
end


function sample_pol(number_of_samples::Integer=1)
    samples = rcopy(R"sample_points(P, $number_of_samples)")
    return samples::Array{Float64,2}
end

function sample_pol_Walk(number_of_samples::Integer=1; random_walk::String="BRDHR")
    rSTDO = R"""
    library(volesti)
    # Create a list for the random walk parameter
    random_walk_list <- list("walk" = $random_walk)
    # Sample points
    samples <- sample_points(P, n = $number_of_samples, list("walk" = $random_walk))
    """
    samples = rcopy(R"samples")
    return samples
end

function get_pol()
    A = rcopy(R"slot(P, 'A')")
    b = rcopy(R"slot(P, 'b')")
    return A::Array{Float64,2}, b::Union{AbstractArray{<:Number},Number}
end


function add_to_pol(Al::AbstractArray,bl::Union{AbstractArray{<:Number},Number})
    R"""
    slot(P, 'A') <- rbind(slot(P, 'A'), $Al)
    slot(P, 'b') <- c(slot(P, 'b'), $bl)"""
    #nothing
    return rcopy(R"P")
end

function get_slack_idx(power_model::AbstractPowerModel)
    # Slack type == 3
    bus_idx =  [k for (k,v) in power_model.data["bus"] if v["bus_type"] == 3]
    bus_idx = parse(Int64,bus_idx[1])
    gen_idx = [k for (k,v) in power_model.data["gen"] if v["gen_bus"] == bus_idx]
    return parse(Int64, gen_idx[1])
end

function get_slack_idx(network)
    # Slack type == 3
    bus_idx =  [k for (k,v) in network["bus"] if v["bus_type"] == 3]
    bus_idx = parse(Int64,bus_idx[1])
    gen_idx = [k for (k,v) in network["gen"] if v["gen_bus"] == bus_idx]
    return  gen_idx[1]
end


function get_header_names(variables::AbstractArray)
    header = [] # used as header in x_opt results csv
    push!(header, JuMP.name.(variables))
    header = map((x) -> uppercase(replace(x, r"0_|\[|\]" => "")),header[1]) # prettify
    return [header]::AbstractArray
end

function gen_samples(n_samples,n_dimensions,level_min, level_max)
    scaling_list = [(level_min, level_max) for _ in 1:n_dimensions]
    plan, _ = LHCoptim(n_samples,n_dimensions,1000)
    scaled_plan = scaleLHC(plan,scaling_list)
    return scaled_plan
end

function gen_samples_vectors(n_samples, n_dimensions, level_min, level_max)
    scaling_list = [(level_min[i], level_max[i]) for i in eachindex(level_min)]
    plan, _ = LHCoptim(n_samples, n_dimensions, 1000)
    scaled_plan = scaleLHC(plan, scaling_list)
    return scaled_plan
end

function update_all_demand(network_data, sample_array)
    for i in 1:length(sample_array)
        network_data["load"]["$i"]["pd"] *= sample_array[i]
        network_data["load"]["$i"]["qd"] *= sample_array[i]
    end
    return network_data
end

function update_all_demand_reactive(network_data, sample_array)
    for i in 1:length(sample_array)
        network_data["load"]["$i"]["qd"] *= sample_array[i]
    end
    return network_data
end


function update_all_limits_reactive(network_data, sample_array)
    for i in 1:length(sample_array)
        network_data["gen"]["$i"]["qmax"] *= sample_array[i]
        network_data["gen"]["$i"]["qmin"] *= sample_array[i]
    end
    return network_data
end

function OPF_feasible_samples(sample_array, network_data)
    global feasible = 0
    global infeasible = 0 
    global non_feasible_index = []
    global feasible_index = []
    for i in 1:length(sample_array[:,1])
        network_data_tmp = deepcopy(network_data)
        demand_profile = sample_array[i,:]
        update_all_demand(network_data_tmp, demand_profile)
        results_QCOPF = solve_opf(network_data_tmp,QCLSPowerModel,optimizer_with_attributes(solver, "print_level" => 0))
        # Check feasibility
        feasible_iteration = 0
        if results_QCOPF["termination_status"] == LOCALLY_SOLVED
            feasible_iteration = 1 
            push!(feasible_index,i)
        else
            println(results_QCOPF["termination_status"])
            push!(non_feasible_index,i)
        end

        global feasible += feasible_iteration
        global infeasible += (1 - feasible_iteration)
    end
    return feasible, infeasible, non_feasible_index, feasible_index
end

function OPF_feasible_samples_reactive(sample_array, network_data)
    global feasible = 0
    global infeasible = 0 
    global non_feasible_index = []
    global feasible_index = []
    for i in 1:length(sample_array[:,1])
        network_data_tmp = deepcopy(network_data)
        demand_profile = sample_array[i,:]
        update_all_demand_reactive(network_data_tmp, demand_profile)
        results_QCOPF = solve_opf(network_data_tmp,QCLSPowerModel,optimizer_with_attributes(solver, "print_level" => 0))
        # Check feasibility
        feasible_iteration = 0
        if results_QCOPF["termination_status"] == LOCALLY_SOLVED
            feasible_iteration = 1 
            push!(feasible_index,i)
        else
            println(results_QCOPF["termination_status"])
            push!(non_feasible_index,i)
        end

        global feasible += feasible_iteration
        global infeasible += (1 - feasible_iteration)
    end
    return feasible, infeasible, non_feasible_index, feasible_index
end


function dim_and_limits_variables(network_data)
    ndim = 0
    bus_gen = []
    gen_index = []
    for j in eachindex(1:length(network_data["gen"]))
        push!(bus_gen, network_data["gen"]["$j"]["gen_bus"])
        if "$j" != get_slack_idx(network_data) && network_data["gen"]["$j"]["pmax"] > 0.0 && network_data["gen"]["$j"]["gen_status"] == 1
            push!(gen_index, network_data["gen"]["$j"]["index"])
            ndim += 1
        end
    end
    bus_gen = unique(bus_gen)
    ndim += length(bus_gen)

    min_lim = []
    max_lim = []
    for i in gen_index
        push!(min_lim, network_data["gen"]["$i"]["pmin"])
        push!(max_lim, network_data["gen"]["$i"]["pmax"])
    end
    for i in bus_gen
        push!(min_lim, network_data["bus"]["$i"]["vmin"])
        push!(max_lim, network_data["bus"]["$i"]["vmax"])
    end
    
    return ndim, bus_gen, gen_index, min_lim, max_lim
end 


#system_data = "/Users/lolacharles/Downloads/Thesis/Code/pglib_opf_case39_epri.m"
#network_data = PowerModels.parse_file(system_data)
#ndim, bus_gen, gen_index, min_lim, max_lim = dim_and_limits_variables(network_data)

function extract_number_and_type(variable_names::Vector{String})
    pg_numbers = Int[]
    vm_numbers = Int[]
    for variable_name in variable_names
        match_result = match(r"(PG|VM)(\d+)", variable_name)
        if match_result !== nothing
            variable_type = match_result.captures[1]
            variable_number = parse(Int, match_result.captures[2])
            if variable_type == "PG"
                push!(pg_numbers, variable_number)
            elseif variable_type == "VM"
                push!(vm_numbers, variable_number)
            end
        end
    end
    return pg_numbers, vm_numbers
end

function extract_number_and_type(variable_name::String)
    match_result = match(r"(PG|VM)(\d+)", variable_name)
    if match_result !== nothing
        variable_type = match_result.captures[1]
        variable_number = parse(Int, match_result.captures[2])
        return variable_type, variable_number
    else
        return nothing, nothing
    end
end

#extract_number_and_type(vcat(header[1]))

function run_qc_relax(pm::AbstractPowerModel, number_of_iterations::Integer, vol::Integer=0)
    """ Iteratively solve modified QC-AC-OPF
    inputs: power model and number of iterations"""
    # Initialize variables
    start = time();
    #info(logger,"QC relaxation tolerance set to $TOLERANCE .")
    optimal_setpoints = []
    vars = []
    pm_xhat = deepcopy(pm)
    resu = optimize_model!(pm_xhat, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    x_hat = []
    # build optimization problem
    # get variables: PG\{slack} and Vm (all generators)
    # slack: bus_type == 3
    #
    # to get all: vars = [pm.var[:nw][0][:pg].data; pm.var[:nw][0][:vm].data];
    #
    slack_gen_idx = get_slack_idx(pm)
    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
        if i != slack_gen_idx
            push!(vars, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
            push!(x_hat, value(JuMP.variable_by_name(pm_xhat.model, string("0_pg[",i,"]"))))
        end
    end
    gen_indexes = unique(map(x -> x["gen_bus"], values(pm.data["gen"])))
    for g in gen_indexes
        push!(vars, JuMP.variable_by_name(pm.model, string("0_vm[",g,"]")))
        push!(x_hat, value(JuMP.variable_by_name(pm_xhat.model, string("0_vm[",g,"]"))))
    end
    # extract header names
    header = get_header_names(vars) # call before normalization
    println(header)
    N = length(vars);
#TODO: Currently the variables of the optimization problem are scaled between
# [0;1] so we sample from a unit cube. if create_scaled_pol() is used the polytope
# polytope gets scaled and we can remove the scaling of the vars. This later
# might prove to be the proper solution as I could not make sure that scaling
# only vars is sufficient.
    nFactor = JuMP.upper_bound.(vars)
    # auxiliary variables to bridge AffExpr <-> NLconstraint
    #   See "Syntax notes" http://www.juliaopt.org/JuMP.jl/dev/nlp/#Nonlinear-Modeling-1)
    @variable(pm.model,aux[1:N])
    @constraint(pm.model,aux .== vars./nFactor)
    vars = vars./nFactor; # normalize after aux is defined

    # Create x_hat random sample
    create_pol(N) # Create initial unit cube
    #create_scaled_pol(N_dim_polytope, lower_bounds_vars , upper_bounds_vars) # Create initial unit cube

    #x_hat = sample_pol() # Initial sample from N unit cube
    # [["PG1", "PG2", "PG3", "PG5", "VM4", "VM1", "VM5", "VM3"]]
    #network_data_check["gen"]["1"]["pg"] = x_hat[1]*nFactor[1] 
    #network_data_check["gen"]["2"]["pg"] = x_hat[2]*nFactor[2]  
    #etwork_data_check["gen"]["3"]["pg"] = x_hat[3]*nFactor[3] 
    #network_data_check["gen"]["5"]["pg"] = x_hat[1]*nFactor[4] 
    
    #network_data_check["bus"]["4"]["vm"] = x_hat[5]*nFactor[5] 
    #network_data_check["bus"]["1"]["vm"] = x_hat[6]*nFactor[6] 
    #network_data_check["bus"]["5"]["vm"] = x_hat[7]*nFactor[7] 
    #network_data_check["bus"]["3"]["vm"] = x_hat[8]*nFactor[8] 
    #result_pf = solve_ac_pf(network_data_check, solver)
    #println(result_pf["termination_status"])
    # Build optimization problem
    # Alter the QCCR problem
    @variable(pm.model, r);
    @objective(pm.model, Min, r);
    @NLparameter(pm.model, x_hat_p[i = 1:N] == x_hat[i]) # updating the parameter has performance gain vs rebuilding the model
    @NLconstraint(pm.model, con_sphere, sqrt(sum((aux[i]-x_hat_p[i])^2 for i in 1:N))<= r);
    index = 0
    #debug(logger,"number of it: $number_of_iterations")
    println(JuMP.objective_function(pm))
    for k = 1:number_of_iterations # Julia for loop has it's own scope, x_hat not accessible
        # Execute optimization
        result = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
        # Calculate outputs
        println(vars)
        x_opt = JuMP.value.(vars);
        println(x_opt)
        println(JuMP.value(r))
        n_normT = transpose(x_opt[:] - x_hat); # This can be simply -> objective_value(pm.model) NO, in the article that is not abs() but l^1 norm
        #debug(logger,string("Objective_value: " ,JuMP.value(r)))
        push!(optimal_setpoints, (x_opt[:] .* nFactor)' )
        if !(isapprox(JuMP.value(r), 0; atol=TOLERANCE))
            println(index)
            index +=1
            # Update results
            #debug(logger,string("new A row: ", n_normT))
            #debug(logger,string("new b val: ", n_normT*x_opt))
            add_to_pol(n_normT, n_normT*x_opt)
            #debug(logger,string("Polytope values |get_pol() \n", get_pol() ))
        end
        try
            x_hat = sample_pol()
            println(x_hat)
            println(nFactor)
            #x_hat = x_hat .* nFactor
            # network_data_check["gen"]["1"]["pg"] = x_hat[1]*nFactor[1] 
            # network_data_check["gen"]["2"]["pg"] = x_hat[2]*nFactor[2]  
            # network_data_check["gen"]["3"]["pg"] = x_hat[3]*nFactor[3] 
            # network_data_check["gen"]["5"]["pg"] = x_hat[1]*nFactor[4] 
            
            # network_data_check["bus"]["4"]["vm"] = x_hat[5]*nFactor[5] 
            # network_data_check["bus"]["1"]["vm"] = x_hat[6]*nFactor[6] 
            # network_data_check["bus"]["5"]["vm"] = x_hat[7]*nFactor[7] 
            # network_data_check["bus"]["3"]["vm"] = x_hat[8]*nFactor[8] 
            # result_pf = solve_ac_pf(network_data_check, solver)
            # println(result_pf["termination_status"])
            #debug(logger,"New points sampled.")
        catch e
            #debug(logger,typeof(e))
            #debug(logger,e)
            return A,b
        end
        if vol > 0
            if k % vol == 0 || k == 1
                # save file
                A,b = get_pol()
                if k == 1
                    writedlm(string("polytope_A_",0,".csv"), A, ',')
                    writedlm(string("polytope_b_",0,".csv"), b, ',')
                else
                    writedlm(string("polytope_A_",k,".csv"), A, ',')
                    writedlm(string("polytope_b_",k,".csv"), b, ',')
                end
            end
        end
        # print("Volume: ", get_volume(N,A,b))
        JuMP.set_value.(x_hat_p, x_hat)
        #debug(logger,"New x_hat JuMP parameter values set.")
    end
    A,b = get_pol()
    #debug(logger,"Return polytope.")
    #debug(logger, string("Variable Names: ", header))
    return A,b,optimal_setpoints,header,nFactor
end


function run_qc_relax2(pm::AbstractPowerModel, number_of_iterations::Integer, vol::Integer=0)
    """ Iteratively solve modified QC-AC-OPF
    inputs: power model and number of iterations"""
    # Initialize variables
    start = time();
    #info(logger,"QC relaxation tolerance set to $TOLERANCE .")
    optimal_setpoints = []
    vars = []
    x_hat = []
    upper_bounds_vars = Vector{Float64}()
    lower_bounds_vars = Vector{Float64}()
    zero_upper_bounds_indices = Int[]

    pm_hat = deepcopy(pm)
    optimize_model!(pm_hat, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    # build optimization problem
    # get variables: PG\{slack} and Vm (all generators)
    # slack: bus_type == 3
    #
    # to get all: vars = [pm.var[:nw][0][:pg].data; pm.var[:nw][0][:vm].data];
    #
    slack_gen_idx = get_slack_idx(pm)
    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
        if i != slack_gen_idx
            push!(vars, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
            push!(upper_bounds_vars, network_data["gen"]["$i"]["pmax"])
            push!(lower_bounds_vars, network_data["gen"]["$i"]["pmin"])
            #push!(x_hat, value(JuMP.variable_by_name(pm_hat.model, string("0_pg[",i,"]"))))
            if upper_bound == 0
                push!(zero_upper_bounds_indices, length(upper_bounds_vars))  # Store index for deletion
            end
        end
    end
    gen_indexes = map(x -> x["gen_bus"], values(pm.data["gen"]))
    for g in gen_indexes
        push!(vars, JuMP.variable_by_name(pm.model, string("0_vm[",g,"]")))
        push!(upper_bounds_vars, network_data["bus"]["$g"]["vmax"])
        push!(lower_bounds_vars, network_data["bus"]["$g"]["vmin"])
        #push!(x_hat, value(JuMP.variable_by_name(pm_hat.model, string("0_vm[",g,"]"))))
    end
    # extract header names
    header = get_header_names(vars) # call before normalization
    deleteat!(upper_bounds_vars, zero_upper_bounds_indices)
    deleteat!(lower_bounds_vars, zero_upper_bounds_indices)
    N = length(vars);
    create_scaled_pol(N, lower_bounds_vars , upper_bounds_vars) # Create initial unit cube
    near_xopt = gen_samples_vectors(50, N, lower_bounds_vars , upper_bounds_vars)
    x_hat = near_xopt[1,:]
#TODO: Currently the variables of the optimization problem are scaled between
# [0;1] so we sample from a unit cube. if create_scaled_pol() is used the polytope
# polytope gets scaled and we can remove the scaling of the vars. This later
# might prove to be the proper solution as I could not make sure that scaling
# only vars is sufficient.
    nFactor = JuMP.upper_bound.(vars)
    # auxiliary variables to bridge AffExpr <-> NLconstraint
    #   See "Syntax notes" http://www.juliaopt.org/JuMP.jl/dev/nlp/#Nonlinear-Modeling-1)
    @variable(pm.model,aux[1:N])
    @constraint(pm.model,aux .== vars)
    #vars = vars./nFactor; # normalize after aux is defined

    # Create x_hat random sample
    #create_pol(N) # Create initial unit cube
    #x_hat = sample_pol() # Initial sample from N unit cube
    #x_hat .= nFactor
    volu = comp_vol()
    println("Volume : ",volu)
    
    # Build optimization problem
    # Alter the QCCR problem
    @variable(pm.model, r);
    @objective(pm.model, Min, r);
    @NLparameter(pm.model, x_hat_p[i = 1:N] == x_hat[i])                         # updating the parameter has performance gain vs rebuilding the model
    @NLconstraint(pm.model, con_sphere, sqrt(sum((aux[i]-x_hat_p[i])^2 for i in 1:N))<= r);

    #debug(logger,"number of it: $number_of_iterations")
    for k = 2:number_of_iterations 
        # Execute optimization
        result = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
        # Calculate outputs
        x_opt = JuMP.value.(vars)
        r_opt = JuMP.value(r)

        println("Xopt : ",x_opt)
        println("Xhat : ",x_hat)
        println("Optimal R :", r_opt)
        #n_normT = transpose(x_opt[:] - x_hat); # This can be simply -> objective_value(pm.model) NO, in the article that is not abs() but l^1 norm
        
        #debug(logger,string("Objective_value: " ,JuMP.value(r)))
        push!(optimal_setpoints, (x_opt[:] .* nFactor)' )
        if !(isapprox(JuMP.value(r), 0; atol=TOLERANCE))
            n_normT = transpose(x_opt[:] - x_hat)
            add_to_pol(n_normT, n_normT*x_opt)

            # Update results
            #debug(logger,string("new A row: ", n_normT))
            #debug(logger,string("new b val: ", n_normT*x_opt))
            println("Add a plane ---------------")
            n_normT = transpose(x_opt[:] - x_hat); # This can be simply -> objective_value(pm.model) NO, in the article that is not abs() but l^1 norm
            println("Norm :",n_normT)
            b = n_normT*x_opt 
            println("b : ",b)
            add_to_pol(n_normT, n_normT*x_opt)
            volu = comp_vol()
            println("Volume :  ", volu)
            #println("Volume after plane : ",R"volume(P)")
            #debug(logger,string("Polytope values |get_pol() \n", get_pol() ))
        end
        try
            #x_opt_minus_10_percent = 0.9 * x_opt
            #x_opt_plus_10_percent = 1.1 * x_opt
            #near_xopt = gen_samples_vectors(50, N, x_opt_minus_10_percent, x_opt_plus_10_percent)
            x_hat = near_xopt[k,:]
            #x_hat = sample_pol()
            #debug(logger,"New points sampled.")
        catch e
            #debug(logger,typeof(e))
            #debug(logger,e)
            return A,b
        end
        if vol > 0
            if k % vol == 0 || k == 1
                # save file
                A,b = get_pol()
                if k == 1
                    writedlm(string("polytope_A_",0,".csv"), A, ',')
                    writedlm(string("polytope_b_",0,".csv"), b, ',')
                else
                    writedlm(string("polytope_A_",k,".csv"), A, ',')
                    writedlm(string("polytope_b_",k,".csv"), b, ',')
                end
            end
        end
        # print("Volume: ", get_volume(N,A,b))
        JuMP.set_value.(x_hat_p, x_hat)
        #debug(logger,"New x_hat JuMP parameter values set.")
    end
    A,b = get_pol()
    #debug(logger,"Return polytope.")
    #debug(logger, string("Variable Names: ", header))
    return A,b
end

function randomize_values(n::Int, values::Vector{Float64})
    if isempty(values)
        throw(ArgumentError("Input vector cannot be empty"))
    end
    
    min_val = minimum(values)
    max_val = maximum(values)
    
    random_values = [rand(length(values)) .* (max_val - min_val) .+ min_val for _ in 1:n]
    
    return random_values
end

function infeasibility_certif_Constraints_check(network_data , Nb_HP, Nb_insamples)

    iter = 3
    volumes = []

    global optimal_setpoints = []
    global original_optimal = []

    data_tight_tmp, stats_tmp = solve_obbt_opf!(network_data, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0), max_iter=iter, model_constructor=QCLSPowerModel)

    ndim, bus_gen, gen_index, min_lim, max_lim = dim_and_limits_variables(data_tight_tmp)
    sample_ops = gen_samples_vectors(Nb_HP,ndim,min_lim,max_lim)
    pm = instantiate_model(data_tight_tmp, QCRMPowerModel, PowerModels.build_opf)
    create_scaled_pol(ndim, min_lim , max_lim)
    v = comp_vol()

    println("Initial volume : ", v)
    push!(volumes, v)
    vars = []
    global x_hat = []

    slack_gen_idx = get_slack_idx(pm)
    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
        if i != slack_gen_idx && data_tight_tmp["gen"]["$i"]["pmax"] > 0.0
            push!(vars, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
        end
    end
    gen_indexes = unique(map(x -> x["gen_bus"], values(pm.data["gen"])))
    for g in gen_indexes
        push!(vars, JuMP.variable_by_name(pm.model, string("0_vm[",g,"]")))
    end
    N = length(vars)
    # extract header names
    header = get_header_names(vars)
    x_hat = sample_ops[1,:]
    array_res = []
    @variable(pm.model, r);
    @variable(pm.model, aux[1:N])
    @constraint(pm.model, aux .== vars)
    @objective(pm.model, Min, r);
    @NLparameter(pm.model, x_hat_p[i = 1:N] == x_hat[i]) # updating the parameter has performance gain vs rebuilding the model
    @NLconstraint(pm.model, con_sphere, sqrt(sum((aux[i]-x_hat_p[i])^2 for i in 1:N))<= r)
    for j in 2:length(sample_ops[:,1])
        println("Iteration $(j-1)-----------------------------------")
        result = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
      
        x_opt = JuMP.value.(vars)
        r_opt = JuMP.value(r)
        push!(optimal_setpoints, x_opt)
        println("x_hat",x_hat)
        println("x_opt",x_opt)
        println(r_opt)
        n_normT = transpose(x_hat - x_opt[:] )
        if !(isapprox(JuMP.value(r), 0; atol=TOLERANCE)) && (r_opt > 0)
            # Update results
            add_to_pol(n_normT, n_normT*x_opt)
            v = comp_vol()
            
            println("Volume : ", v)
            push!(volumes, v)
        else
            push!(original_optimal, x_hat)

        end

        global x_hat = sample_ops[j,:]
        JuMP.set_value.(x_hat_p, x_hat)

    end
    nb_samples = Nb_insamples
    global in_samples = sample_pol(nb_samples)
    global nb_feasible = 0

    pm, N, vars, header = instantiate_system_QCRM(data_tight_tmp)
    pg_numbers, vm_numbers = extract_number_and_type(vcat(header[1]))
    pf_results = []
    for i in 1:nb_samples          #length(optimal_setpoints[:,1]) 
        #println("$i______________________")
        data_opf_verif = deepcopy(data_tight_tmp)
        
        for g in eachindex(pg_numbers)
            data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = in_samples[g,i] 
        end
        for v in eachindex(vm_numbers)
            data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = in_samples[length(pg_numbers)+v,i] 
        end
        PF_res1 = adjust_PVPQ(data_opf_verif, 8)
        
        if PF_res1["termination_status"] == LOCALLY_SOLVED
            global nb_feasible += 1
            push!(pf_results, PF_res1)
        else 
            println(PF_res1["termination_status"] )
        end
   
    end

    global over_array = []
    global under_array = []
    global pg_array = []
    global qg_array = []
    global sm_array = []
    global res_arr = []
    global total_over_array = []
    global total_under_array = []
    global total_pg_array = []
    global total_qg_array = []
    global total_sm_array = []

    for i in 1:length(pf_results)
        vm_vio_over, vm_vio_under = check_vm_violations(data_tight_tmp, pf_results[i]["solution"])
        push!(over_array, vm_vio_over)
        push!(under_array, vm_vio_under)
        pg_vio, qg_vio = check_pg_pq_violations(data_tight_tmp, pf_results[i]["solution"])
        sm_vio = check_flow_violations(data_tight_tmp, pf_results[i]["solution"])
        push!(pg_array, pg_vio)
        push!(qg_array, qg_vio)
        push!(sm_array, sm_vio)

    end

    threshold_qg = 10^-06
    threshold_sm = 10^-06
    
    nb_pg, index_pg = find_zeros_and_indexes(pg_array,threshold_qg)
    nb_qg, index_qg = find_zeros_and_indexes(qg_array,threshold_qg)
    nb_sm, index_sm = find_zeros_and_indexes(sm_array,10^-06)
    nb_ovi, index_ovi = find_zeros_and_indexes(over_array,10^-06)
    nb_uvi, index_uvi = find_zeros_and_indexes(under_array,10^-06)
    common_elements = [x for x in index_qg if (x in index_sm  && x in index_ovi && x in index_uvi && x in index_pg)]
    rest_of_the_indices = [x for x in eachindex(pf_results) if !(x in common_elements)]

    OPs = []
    for i in common_elements
        op = get_Full_OP(data_tight_tmp, pf_results[i]["solution"])
        push!(OPs, op)
    end

    OPs_notFeas = []
    for i in rest_of_the_indices
        op = get_Full_OP(data_tight_tmp, pf_results[i]["solution"])
        push!(OPs_notFeas, op)
    end

    return OPs, pf_results, data_tight_tmp, OPs_notFeas, common_elements, rest_of_the_indices
end





function comment_line_in_file(file_path::String, target_line::String)
    # Read the file
    lines = readlines(file_path)
    
    # Process the lines
    modified_lines = String[]
    for line in lines
        # Check if the line starts with the target string
        if startswith(line, target_line)
            # Comment out the line by adding '%'
            push!(modified_lines, "%" * line)
        else
            # Keep the line as is
            push!(modified_lines, line)
        end
    end
    
    # Write the modified lines back to the file
    open(file_path, "w") do file
        for line in modified_lines
            println(file, line)
        end
    end
end

#function update_gen_network(network, new_gen)
#    ind = 1
#    for (i, gen) in network["gen"]
#        if network["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && network["gen"][i]["pmax"] > 0.0
#            gen["pg"] += new_gen[ind]
#            ind += 1
#        end
#    end 
#end

function get_SetPoint(network_data)
    Pg_SetPoint = Dict{String, Float64}()
    for (i, gen) in network_data["gen"]
        if network_data["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && network_data["gen"][i]["pmax"] > 0.0
            Pg_SetPoint[i] = gen["pg"]
            #push!(Pg_SetPoint, (i, gen["pg"] ))
        end
    end 
    return Pg_SetPoint
end

function update_SetPoint(SetPoint, new_value, nb)
    newSetPoint = deepcopy(SetPoint)
    newSetPoint[nb] = new_value
    return newSetPoint
end

machine_oneDoneQ() = OneDOneQMachine(
           R = 0.0,
           Xd = 1.3125,
           Xq = 1.2578,
           Xd_p = 0.1813,
           Xq_p = 0.25,
           Td0_p = 5.89,
           Tq0_p = 0.6,
       )

machine1() = RoundRotorQuadratic(
R = 0.0043 ,
Td0_p = 0.5871,
Td0_pp = 0.0248,
Tq0_p = 0.1351,
Tq0_pp = 0.0267,
Xd = 1.670,
Xq = 1.600,
Xd_p = 0.265,
Xq_p = 0.460,
Xd_pp = 0.205,
Xl = 0.150,
Se = (0.091, 0.400 )
)

avr_type1() = AVRTypeI(
           Ka = 20.0,
           Ke = 0.01,
           Kf = 0.063,
           Ta = 0.2,
           Te = 0.314,
           Tf = 0.35,
           Tr = 0.001,
           Va_lim = (min = -5.0, max = 5.0),
           Ae = 0.0039, #1st ceiling coefficient
           Be = 1.555, #2nd ceiling coefficient
       )

tg_none() = TGFixed(efficiency = 1.0)       
pss_none() = PSSFixed(V_pss = 0.0)
shaft_no_damping() = SingleMass(
           H = 3.01,
           D = 0.0,
       )

TGOV1() = SteamTurbineGov1(
    R = 0.05,
    T1 = 0.49,
    valve_position_limits = (min = 33 , max = 0.4),
    T2 = 2.1,
    T3 = 7,
    D_T = 0,
    DB_h = 0,
    DB_l = 0,
    T_rate = 0,
)       


tg_type1() = TGTypeI(
0.02, #R
0.1, #Ts
0.45, #Tc
0.0, #T3
12.0, #T4
50.0, #T5
(min = 0.0, max = 1.2), #P_lims
)


function create_system(network)
    export_matpower("file.m", network)
    file_path = "file.m"
    # The target line to detect and comment out
    target_line = "mpc.multinetwork"
    # Call the function to process the file
    comment_line_in_file(file_path, target_line)
    target_line = "mpc.multiinfrastructure"
    # Call the function to process the file
    comment_line_in_file(file_path, target_line)
    sys = System("file.m")

    return sys
end

function add_dyn_system_basic(sys)
    for g in get_components(Generator, sys)

        #Create the dynamic generator
        case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine_oneDoneQ(),
            shaft = shaft_no_damping(),
            avr = avr_type1(),
            prime_mover = tg_type1(), #tg_none(),
            pss = pss_none(),
            )
        #Attach the dynamic generator to the system
        add_component!(sys, case_gen, g)
    end
    return sys
end


function add_dyn_system_basic_withoutTG(sys)
    for g in get_components(Generator, sys)

        #Create the dynamic generator
        case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine1(),
            shaft = shaft_no_damping(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
            )
        #Attach the dynamic generator to the system
        add_component!(sys, case_gen, g)
    end
    return sys
end

machine1() = RoundRotorQuadratic(
    R = 0.0043 ,
    Td0_p = 0.5871,
    Td0_pp = 0.0248,
    Tq0_p = 0.1351,
    Tq0_pp = 0.0267,
    Xd = 1.670,
    Xq = 1.600,
    Xd_p = 0.265,
    Xq_p = 0.460,
    Xd_pp = 0.205,
    Xl = 0.150,
    Se = (0.091, 0.400 )
)

machine2() = RoundRotorQuadratic(
    R = 0.0035 ,
    Td0_p = 1.100,
    Td0_pp = 0.0277,
    Tq0_p = 0.1086 ,
    Tq0_pp = 0.0351,
    Xd = 1.180 ,
    Xq = 1.050 ,
    Xd_p = 0.220,
    Xq_p = 0.380 ,
    Xd_pp = 0.145 ,
    Xl = 0.075 ,
    Se = (0.0933, 0.4044  )
)

machine3() = RoundRotorQuadratic(
    R = 0.0 ,
    Td0_p = 11.600 ,
    Td0_pp = 0.058,
    Tq0_p = 0.159 ,
    Tq0_pp = 0.201,
    Xd = 2.373,
    Xq = 1.172 ,
    Xd_p = 0.343,
    Xq_p = 1.172,
    Xd_pp = 0.231 ,
    Xl = 0.150,
    Se = (0.091, 0.400 )
)

machine4() = RoundRotorQuadratic(
    R = 0.0025 ,
    Td0_p = 8.000,
    Td0_pp = 0.0525,
    Tq0_p = 0.008 ,
    Tq0_pp = 0.0151,
    Xd = 1.769 ,
    Xq = 0.855 ,
    Xd_p = 0.304 ,
    Xq_p = 0.5795,
    Xd_pp = 0.2035 ,
    Xl = 0.1045,
    Se = (0.304, 0.666)
)

avr_none() = AVRFixed(0.0)

machine_genrou() = RoundRotorExponential(;
    R = 0.0,
    Td0_p = 8.0,
    Td0_pp = 0.03,
    Tq0_p = 0.4,
    Tq0_pp = 0.05,
    Xd = 1.8,
    Xq = 1.7,
    Xd_p = 0.3,
    Xq_p = 0.55,
    Xd_pp = 0.25,
    Xl = 0.2,
    Se = (0.0, 0.0),
)

shaft_no_damping_genrou() = SingleMass(
           H = 6.5,
           D = 0.0,
       )

avr_type1() = AVRTypeI(
           Ka = 20.0,
           Ke = 0.01,
           Kf = 0.063,
           Ta = 0.2,
           Te = 0.314,
           Tf = 0.35,
           Tr = 0.001,
           Va_lim = (min = -5.0, max = 5.0),
           Ae = 0.0039, #1st ceiling coefficient
           Be = 1.555, #2nd ceiling coefficient
       )

sexs() = SEXS(
    Ta_Tb = 0.1,
    Tb = 10,
    K = 100,
    Te = 0.1,
    V_lim = (min = 0.0, max = 5.0),
)

tg_none() = TGFixed(efficiency = 1.0)       
pss_none() = PSSFixed(V_pss = 0.0)

shaft1() = SingleMass(H = 2.656 , D = 2.0)
shaft2() = SingleMass(H = 4.985 , D = 2.0)
shaft3() = SingleMass(H = 1.520 , D = 0.0)
shaft4() = SingleMass(H = 1.20 , D = 0.0)

tg_type1() = TGTypeI(
    0.02, #R
    0.1, #Ts
    0.45, #Tc
    0.0, #T3
    12.0, #T4
    50.0, #T5
    (min = 0.0, max = 1.2), #P_lims
)

function construct_dyn_14bus(sys_14)
    for g in get_components(Generator, sys_14)
        if get_number(get_bus(g)) == 1
            case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine1(),
            shaft = shaft1(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
        )
            #Attach the dynamic generator to the system
            add_component!(sys_14, case_gen, g)

        end 

        if get_number(get_bus(g)) == 2
            case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine2(),
            shaft = shaft2(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
        )
            #Attach the dynamic generator to the system
            add_component!(sys_14, case_gen, g)

        end 

        if get_number(get_bus(g)) == 3
            case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine3(),
            shaft = shaft3(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
        )
            #Attach the dynamic generator to the system
            add_component!(sys_14, case_gen, g)

        end 

        if get_number(get_bus(g)) == 6  || get_number(get_bus(g)) == 8
            case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine4(),
            shaft = shaft4(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
        )
            #Attach the dynamic generator to the system
            add_component!(sys_14, case_gen, g)
        end 
    end
end


function get_OP_from_System(network)
    slack_gen_idx = get_slack_idx(network)
    OP_x = []
    for (g, gen) in network["gen"]
        if g != slack_gen_idx && gen["pmax"] > 0.0
            push!(OP_x, gen["pg"] )
        end
    end
    gen_indexes = unique(map(x -> x["gen_bus"], values(network["gen"])))
    for b in gen_indexes
        push!(OP_x, network["bus"]["$b"]["vm"])
    end
    return OP_x
end


function N_1_step(network_data, OPs_Feas, cont_out)

    global over_array = []
    global under_array = []
    global pg_array = []
    global qg_array = []
    global sm_array = []
    global res_arr = []
    global total_over_array = []
    global total_under_array = []
    global total_pg_array = []
    global total_qg_array = []
    global total_sm_array = []
    global flag = false 
    global N_1_feasible = []
    threshold = 10-6



    for i in eachindex(OPs_Feas)#common_elements
        OP = OPs_Feas[i] 
        data_opf_verif = deepcopy(network_data)
        
        for g in eachindex(pg_numbers)
            data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = OP[g]
        end
        for v in eachindex(vm_numbers)
            data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = OP[length(pg_numbers)+v]
        end
        global network_data = deepcopy(data_opf_verif)

        network_data["area_gens"] = Dict()
        network_data["gen_contingencies"] = []
        network_data["branch_contingencies"] = []
        contingencies = []
        contingencies_gen = []

        for (i, branch) in network_data["branch"]
            if i in cont_out
            continue
            else
                push!(contingencies, (idx=parse(Int,i), label="LINE-$(i)", type="branch"))
            end
        end

        network_data["branch_contingencies"] = contingencies
        network_data["gen_contingencies"] = contingencies_gen
        global multinetwork = build_c1_scopf_multinetwork_modif(network_data)
        if multinetwork["per_unit"] == true
            for (n, network) in multinetwork["nw"]
                multinetwork["nw"]["$n"]["per_unit"] = true
            end
        end
        
        result = run_c1_scopf_modif(multinetwork, ACPPowerModel, Ipopt.Optimizer)
        for l in 0:(length(multinetwork["nw"])-1)
            update_data!(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
        end
        result = run_c1_scopf_modif(multinetwork, ACPPowerModel, Ipopt.Optimizer)

        for l in 0:(length(multinetwork["nw"])-1)
            update_data!(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
        end
        global result = run_c1_scopf_modif(multinetwork, ACPPowerModel, Ipopt.Optimizer)
        
        

        for l in 0:(length(multinetwork["nw"])-1)
            pg_vio, qg_vio = check_pg_pq_violations(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
            sm_vio = check_flow_violations(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
            vm_over_tmp, vm_under_tmp = check_vm_violations(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
            push!(over_array, vm_over_tmp)
            push!(under_array, vm_under_tmp)
            push!(pg_array, pg_vio)
            push!(qg_array, qg_vio)
            push!(sm_array, sm_vio)
            
        end

        push!(total_over_array,over_array)
        push!(total_under_array,under_array)
        push!(total_pg_array,pg_array)
        push!(total_qg_array, qg_array)
        push!(total_sm_array, sm_array)

        if (all_below_threshold(over_array, threshold) &&
            all_below_threshold(under_array, threshold) &&
            all_below_threshold(pg_array, threshold) &&
            all_below_threshold(qg_array, threshold) &&
            all_below_threshold(sm_array, threshold))
            push!(N_1_feasible, i)
        end
        global over_array = []
        global under_array = []
        global pg_array = []
        global qg_array = []
        global sm_array = []

    end
    return total_sm_array, total_qg_array, total_over_array, total_under_array , total_pg_array, N_1_feasible, result["solution"]["nw"], multinetwork["nw"]
end


function SSS_eval(data_tight_HP, global_OPs, pg_numbers, vm_numbers, dir_repo)

    global total_damp = []
    global total_dist = []
    global total_eigen = []

    for i in eachindex(global_OPs)   
        data_build = deepcopy(data_tight_HP)
        for g in eachindex(pg_numbers)
            data_build["gen"]["$(pg_numbers[g])"]["pg"] = global_OPs[i][g] 
        end
        for v in eachindex(vm_numbers)
            data_build["bus"]["$(vm_numbers[v])"]["vm"] = global_OPs[i][length(pg_numbers)+v] 
        end

        syst = create_system(data_build)
        add_dyn_system_basic(syst)
        damping_tmp, dist_tmp, eigen = small_signal_module_sys(dir_repo,syst, i, in_file= false)
        push!(total_damp, damping_tmp)
        push!(total_dist, dist_tmp)
        push!(total_eigen, eigen)
    
    end
    return total_damp, total_dist, total_eigen
end 


function closest_to_zero_indices(arr, N::Int)
    # Create an array of tuples with absolute value and original index
    abs_with_index = [(abs(val), idx) for (idx, val) in enumerate(arr)]
    
    # Sort the array by absolute value
    sorted_abs_with_index = sort(abs_with_index, by = x -> x[1])
    
    # Get the first N indices from the sorted array
    indices = [sorted_abs_with_index[i][2] for i in 1:N]
    
    return indices
end


function DW_step(data_tight, index_res, cls_OP, pf_results_prev, distance, alpha)

    file = open("/Users/lolacharles/Downloads/Thesis/Code/Final/DW_file.txt", "w")

    global file_nb = 0
    global damping_array = []
    global dist_arr = []
    global plot_damping = []
    global plot_dist = []
    global DW_OPS = []

    gen_keys = collect(keys(data_tight["gen"]))
    gen_slack = get_slack_idx(data_tight)
    filtered_gen_keys = filter(key -> key != gen_slack && data_tight["gen"][key]["pmax"] > 0.0, gen_keys)
    # Get the list of active generators
    active_gens = [data_tight["gen"][key] for key in filtered_gen_keys]

    # Compute the length of active generators
    lgt_active_gen = length(active_gens)

    for i in cls_OP 
        data_build = deepcopy(data_tight)
        update_data!(data_build, pf_results_prev[index_res[i]]["solution"])
        sys_studied = create_system(data_build)
        add_dyn_system_basic(sys_studied)
        damping_tmp, dist_tmp, eigenval = small_signal_module_sys(dir_repo,sys_studied, 100)
        println(file, "Initial damping :", damping_tmp)
        println(file, "Initial dist :", dist_tmp)   
        println(file, "Initial eigenv :",eigenval)   
        global initial_DR = damping_tmp
        global damping_array_global = []
        global dist_array_global = []
        for DW in 1:15
            for (g, gen) in data_build["gen"]
                if data_build["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_build["gen"][g]["pmax"] > 0.0
                    network_data_copy = deepcopy(data_build)
                    global updated_value = data_build["gen"]["$g"]["pg"] + 0.2 #* data_build["gen"]["$g"]["pmax"]
                    if network_data_copy["gen"]["$g"]["pmax"] < updated_value 
                        network_data_copy["gen"]["$g"]["pg"] = network_data_copy["gen"]["$g"]["pmax"]
                        println(file, "Out of the bounds ")
                    else
                        network_data_copy["gen"]["$g"]["pg"] = updated_value
                    end
                    
                    sys_studied = create_system(network_data_copy)
                    add_dyn_system_basic(sys_studied)
        
                    global file_nb += 1 
                    damping_tmp, dist_tmp = small_signal_module_sys(dir_repo,sys_studied, file_nb)
                    push!(damping_array, damping_tmp)
                    push!(dist_arr, dist_tmp)
        
                end
            end 
        
        
            for (g, gen) in data_build["gen"]
                if data_build["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_build["gen"][g]["pmax"] > 0.0
                    network_data_copy = deepcopy(data_build)
                    global updated_value = data_build["gen"]["$g"]["pg"] - 0.2 #* data_build["gen"]["$g"]["pmax"]
                    if network_data_copy["gen"]["$g"]["pmin"] > updated_value 
                        network_data_copy["gen"]["$g"]["pg"] = network_data_copy["gen"]["$g"]["pmin"]
                        println(file, "Out of the bounds ")
                    else
                        network_data_copy["gen"]["$g"]["pg"] = updated_value
                    end
                    
                    sys_studied = create_system(network_data_copy)
                    add_dyn_system_basic(sys_studied)
        
                    global file_nb += 1 
                    damping_tmp, dist_tmp = small_signal_module_sys(dir_repo,sys_studied, file_nb)
                    push!(damping_array, damping_tmp)
                    push!(dist_arr, dist_tmp)
        
                end
            end 
            index_gen_disturb = argmin(abs.(damping_array))
            sorted_index = sortperm(abs.(damping_array))
            min_damp = minimum(damping_array)
            
            if abs(min_damp) > distance[1]
                step = alpha[1]
            elseif distance[2] < abs(min_damp)  < distance[1]
                step = alpha[2]
            elseif distance[3] < abs(min_damp) < distance[2]
                step = alpha[3] #0.15
            elseif abs(min_damp) < distance[3]
                step = alpha[4] #0.15    
            end

            println(file, "$DW __________________--")
            sign_DW = 1
            if index_gen_disturb > (lgt_active_gen)
                sign_DW = -1
                index_gen_disturb -= lgt_active_gen
            end    
            

            if data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pg"] <= data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pmin"] || data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pg"] >= data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pmax"]
                index_gen_disturb = sorted_index[2]
                if index_gen_disturb > (lgt_active_gen)
                    sign_DW = -1
                    index_gen_disturb -= lgt_active_gen
                end    

                if data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pg"] <= data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pmin"] || data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pg"] >= data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pmax"]
                    index_gen_disturb = sorted_index[3]
                    if index_gen_disturb > (lgt_active_gen)
                        sign_DW = -1
                        index_gen_disturb -= lgt_active_gen
                    end    
                    
                end

            end

            sP = get_SetPoint(data_build)
            println(file, collect(values(sP)))
            index_gen_disturb = filtered_gen_keys[index_gen_disturb]
            newOP_PG = data_build["gen"][index_gen_disturb]["pg"] + sign_DW * 0.2 #* data_build["gen"][index_gen_disturb]["pmax"]
            new_SP = update_SetPoint(sP, newOP_PG, index_gen_disturb)
            println(file,collect(values(new_SP)))
            println(collect(values(new_SP)))
            println(collect(values(sP)))

            grad = (collect(values(new_SP)) - collect(values(sP)))/norm(collect(values(new_SP)) - collect(values(sP)))
            println(file, grad)
            new =  step*data_build["gen"][index_gen_disturb]["pg"].*(grad)#.*collect(values(new_SP)))
            println(file, new)
            println(file,data_build["gen"] )
            update_gen_network(data_build, new)
            println(file,data_build["gen"] )
            sys_studied = create_system(data_build)
            add_dyn_system_basic(sys_studied)
            damping_tmp, dist_tmp, eigenval = small_signal_module_sys(dir_repo,sys_studied, 100)
            println(file,data_build["gen"] )
            println(file, damping_tmp)
            println(file, dist_tmp)   
            println(file, eigenval)        
            push!(plot_damping, damping_tmp)
            push!(plot_dist, dist_tmp)
            push!(damping_array_global , damping_array)
            global damping_array = []
            push!(dist_array_global , dist_arr)
            global dist_arr = []
            OP_tmp = get_Full_OP(data_build,data_build)
            push!(DW_OPS, OP_tmp)
            push!(DW_OPS, (damping_tmp, dist_tmp))


            if initial_DR < damping_tmp
                break
            end
        end

        println(file, dist_array_global)
        println(file, damping_array_global)
        dist_array_global = []
        damping_array_global = []
        
    end
    close(file)

    return DW_OPS 

end


function sample_MVND(OP_HP, network_basic, data_tight, Nb_ops)
    mean_vec = mean(OP_HP, dims=1)

    cov_matrix = cov(OP_HP)
    eps = 1e-6  # Small value to add to the diagonal
    cov_matrix += eps * I
    # Bias the sampling by scaling the covariance matrix
    biased_cov_matrix = 0.25 * cov_matrix

    # Ensure the biased covariance matrix is positive definite
    while !isposdef(biased_cov_matrix)
        eps *= 10
        biased_cov_matrix = bias_factor * (cov_matrix + eps * I)
    end

    # Create a Multivariate Normal Distribution with the biased covariance matrix


    # Create a Multivariate Normal Distribution with the biased covariance matrix
    mvn = MvNormal(mean_vec[1], cov_matrix)


    # Draw N3 samples from the biased distribution
    biased_samples = rand(mvn, Nb_ops)
    pm, N, vars, header = instantiate_system_QCRM(data_tight)
    pg_numbers, vm_numbers = extract_number_and_type(vcat(header[1]))
    pf_results = []
    for i in eachindex(biased_samples[1,:])
        data_opf_verif = deepcopy(data_tight)
        
        for g in eachindex(pg_numbers)
            data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = biased_samples[g,i] #optimal_setpoints[i][g] #in_samples[g,i]
        end
        for v in eachindex(vm_numbers)
            data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = biased_samples[length(pg_numbers)+v,i] #optimal_setpoints[i][length(pg_numbers)+v] #in_samples[length(pg_numbers)+v,i]
        end

        PF_res1 = adjust_PVPQ(data_opf_verif, 8)

        if PF_res1["termination_status"] == LOCALLY_SOLVED
            global nb_feasible += 1
            push!(pf_results, PF_res1)
        else 
            println(PF_res1["termination_status"] )
        end


    end

    global over_array = []
    global under_array = []
    global pg_array = []
    global qg_array = []
    global sm_array = []
    global res_arr = []
    global total_over_array = []
    global total_under_array = []
    global total_pg_array = []
    global total_qg_array = []
    global total_sm_array = []

    for i in 1:length(pf_results)
        vm_vio_over, vm_vio_under = check_vm_violations(network_basic, pf_results[i]["solution"])
        push!(over_array, vm_vio_over)
        push!(under_array, vm_vio_under)
        pg_vio, qg_vio = check_pg_pq_violations(network_basic, pf_results[i]["solution"])
        sm_vio = check_flow_violations(network_basic, pf_results[i]["solution"])
        push!(pg_array, pg_vio)
        push!(qg_array, qg_vio)
        push!(sm_array, sm_vio)

    end

    threshold_qg = 10^-06
    threshold_sm = 10^-06

    nb_pg, index_pg = find_zeros_and_indexes(pg_array,threshold_qg)
    nb_qg, index_qg = find_zeros_and_indexes(qg_array,threshold_qg)
    nb_sm, index_sm = find_zeros_and_indexes(sm_array,10^-06)
    nb_ovi, index_ovi = find_zeros_and_indexes(over_array,10^-06)
    nb_uvi, index_uvi = find_zeros_and_indexes(under_array,10^-06)
    common_elements = [x for x in index_qg if (x in index_sm  && x in index_ovi && x in index_uvi && x in index_pg)]
    rest_of_the_indices = [x for x in eachindex(pf_results) if !(x in common_elements)]

    OPs = []
    for i in common_elements
        op = get_Full_OP(data_tight, pf_results[i]["solution"])
        push!(OPs, op)
    end

    OPs_notFeas = []
    for i in rest_of_the_indices
        op = get_Full_OP(data_tight, pf_results[i]["solution"])
        push!(OPs_notFeas, op)
    end

    return OPs, OPs_notFeas

end

Demand_3 = 
[[0.8666666666666666, 0.7, 0.9777777777777777],
[0.8111111111111111, 1.2, 0.9222222222222222],
[1.0333333333333332, 0.8666666666666666, 1.2],
[0.9777777777777777, 1.1444444444444444, 1.1444444444444444],
[0.9222222222222222, 0.8111111111111111, 0.7],
[1.1444444444444444, 0.7555555555555555, 0.8666666666666666],
[0.7, 0.9777777777777777, 1.0888888888888888],
[0.7555555555555555, 0.9222222222222222, 0.8111111111111111],
[1.2, 1.0333333333333332, 1.0333333333333332],
[1.0888888888888888, 1.0888888888888888, 0.7555555555555555]]

Demand_5 = 
[[1.0333333333333332, 0.7, 0.9777777777777777],
 [1.1444444444444444, 0.9222222222222222, 1.2],
 [0.9222222222222222, 1.2, 0.8111111111111111],
 [0.9777777777777777, 1.1444444444444444, 1.0888888888888888],
 [0.7, 1.0888888888888888, 1.0333333333333332],
 [0.7555555555555555, 0.7555555555555555, 0.8666666666666666],
 [0.8111111111111111, 0.9777777777777777, 0.7],
 [1.0888888888888888, 0.8111111111111111, 0.7555555555555555],
 [0.8666666666666666, 0.8666666666666666, 1.1444444444444444],
 [1.2, 1.0333333333333332, 0.9222222222222222]]

 Demand_14 = 
 [[0.7555555555555555, 0.7, 0.8111111111111111, 0.8111111111111111, 0.7555555555555555, 0.8111111111111111, 0.9777777777777777, 1.0888888888888888, 1.0333333333333332, 1.0888888888888888, 0.7555555555555555],
 [1.0333333333333332, 0.8666666666666666, 0.7, 0.8666666666666666, 0.9222222222222222, 1.2, 0.8666666666666666, 0.7555555555555555, 1.1444444444444444, 0.9777777777777777, 1.2],
 [0.9222222222222222, 0.7555555555555555, 1.2, 1.0333333333333332, 0.8666666666666666, 1.0333333333333332, 0.7, 0.9777777777777777, 0.7555555555555555, 1.1444444444444444, 0.9222222222222222],
 [1.1444444444444444, 1.1444444444444444, 0.9222222222222222, 0.7555555555555555, 0.9777777777777777, 0.7555555555555555, 0.7555555555555555, 1.2, 0.8666666666666666, 0.8666666666666666, 1.0888888888888888],
 [0.9777777777777777, 0.8111111111111111, 0.8666666666666666, 1.2, 1.2, 0.7, 0.9222222222222222, 0.8111111111111111, 0.9777777777777777, 0.7555555555555555, 1.0333333333333332],
 [0.7, 1.2, 0.9777777777777777, 1.0888888888888888, 0.8111111111111111, 0.9777777777777777, 0.8111111111111111, 0.8666666666666666, 0.9222222222222222, 0.7, 0.8111111111111111],
 [1.2, 0.9222222222222222, 1.1444444444444444, 0.9222222222222222, 0.7, 0.9222222222222222, 1.2, 0.9222222222222222, 1.0888888888888888, 0.8111111111111111, 0.8666666666666666],
 [0.8111111111111111, 1.0333333333333332, 1.0888888888888888, 1.1444444444444444, 1.1444444444444444, 1.0888888888888888, 1.0333333333333332, 1.1444444444444444, 1.2, 1.0333333333333332, 1.1444444444444444],
 [1.0888888888888888, 0.9777777777777777, 0.7555555555555555, 0.9777777777777777, 1.0333333333333332, 1.1444444444444444, 1.1444444444444444, 1.0333333333333332, 0.7, 0.9222222222222222, 0.7],
 [0.8666666666666666, 1.0888888888888888, 1.0333333333333332, 0.7, 1.0888888888888888, 0.8666666666666666, 1.0888888888888888, 0.7, 0.8111111111111111, 1.2, 0.9777777777777777]]


 function read_file_to_vectors(filename::String)
    vectors = []  # Array to hold the vectors
    open(filename, "r") do file
        for line in eachline(file)
            # Split the line by whitespace and convert to numbers
            vector = parse.(Float64, split(line))
            push!(vectors, vector)
        end
    end
    return vectors
end

file_txt = "/Users/lolacharles/Downloads/Thesis/Code/Final/Load_demand_39.txt"
Demand_39 = read_file_to_vectors(file_txt)