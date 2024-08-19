function variable_gen_power_modif(pm::AbstractPowerModel; kwargs...)
    variable_gen_power_real_modif(pm; kwargs...)
    variable_gen_power_imaginary_modif(pm; kwargs...)
end


"variable: `pg[j]` for `j` in `gen`"
function variable_gen_power_real_modif(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pg = var(pm, nw)[:pg] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pg",
        start = comp_start_value(ref(pm, nw, :gen, i), "pg_start")
    )
    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_lower_bound(pg[i], gen["pmin"])
            JuMP.set_upper_bound(pg[i], gen["pmax"])
        end
    else
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_lower_bound(pg[i], 0.0)
            if gen["pmax"] == 0.0
                JuMP.set_upper_bound(pg[i], gen["pmax"])
            end
        end    
    end

    report && sol_component_value(pm, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end

"variable: `qq[j]` for `j` in `gen`"
function variable_gen_power_imaginary_modif(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qg = var(pm, nw)[:qg] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_qg",
        start = comp_start_value(ref(pm, nw, :gen, i), "qg_start")
    )

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_lower_bound(qg[i], gen["qmin"])
            JuMP.set_upper_bound(qg[i], gen["qmax"])
        end
    end

    report && sol_component_value(pm, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end


function solve_ac_opf_modif(file, optimizer; kwargs...)
    return solve_opf_modif(file, ACPPowerModel, optimizer; kwargs...)
end

function solve_opf_modif(file, model_type::Type, optimizer; kwargs...)
    return solve_model(file, model_type, optimizer, build_opf_modif; kwargs...)
end

function build_opf_modif(pm::AbstractPowerModel)
    variable_bus_voltage(pm, bounded = false)
    variable_gen_power_modif(pm, bounded = false)
    variable_branch_power(pm , bounded = false)
    variable_dcline_power(pm , bounded = false)
   
    objective_min_fuel_and_flow_cost(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end


function check_flow_violations(network, solution)
    
    update_data!(network, solution)
    branch_res = calc_branch_flow_ac(network)
    nb_branch = length(network["branch"])
    sm_vio = NaN
    
    if haskey(branch_res, "branch")
        sm_vio = 0.0
        for (i,branch) in branch_res["branch"]
            if isnan(branch["qf"]) &&  isnan(branch["qf"]) && isnan(branch["qt"]) && isnan(branch["qt"])
                continue 
            else
                s_fr = abs(branch["pf"])
                s_to = abs(branch["pt"])
                
                if !isnan(branch["qf"]) && !isnan(branch["qt"])
                    s_fr = sqrt(branch["pf"]^2 + branch["qf"]^2)
                    s_to = sqrt(branch["pt"]^2 + branch["qt"]^2)
                end

                # note true model is rate_c
                #vio_flag = false
                rating = network["branch"]["$i"]["rate_c"]

                if s_fr > rating
                    sm_vio += ((s_fr - rating)/rating)*100
                    #vio_flag = true
                end
                if s_to > rating
                    sm_vio += ((s_to - rating)/rating)*100
                    #vio_flag = true
                end
                #if vio_flag
                #    info(_LOGGER, "$(i), $(branch["f_bus"]), $(branch["t_bus"]): $(s_fr) / $(s_to) <= $(branch["rate_c"])")
                #end
            end
        end
    end
    sm_vio = sm_vio/nb_branch
    return sm_vio
end

function check_flow_violations_bis(network, solution)
    #if !haskey(network, "branch")
    #    update_data!(network, solution)
    #    branch_res = calc_branch_flow_ac(network)
    #    nb_branch = length(network["branch"])
    #    sm_vio = NaN
    #end
    if haskey(network, "branch")
        sm_vio = []
        s_fr_arr = []
        s_to_arr = []
        for (i,branch) in network["branch"]
            if isnan(branch["qf"]) &&  isnan(branch["qf"]) && isnan(branch["qt"]) && isnan(branch["qt"])
                continue 
            else
                s_fr = abs(branch["pf"])
                s_to = abs(branch["pt"])
                
                if !isnan(branch["qf"]) && !isnan(branch["qt"])
                    s_fr = sqrt(branch["pf"]^2 + branch["qf"]^2)
                    s_to = sqrt(branch["pt"]^2 + branch["qt"]^2)
                end
                push!(s_fr_arr, s_fr)
                push!(s_to_arr, s_to)
                # note true model is rate_c
                #vio_flag = false
                rating = network["branch"]["$i"]["rate_c"]

                if s_fr > rating
                    #sm_vio += ((s_fr - rating)/rating)*100
                    push!(sm_vio,i)
                    push!(sm_vio, ((s_fr - rating)/rating)*100)
                    #vio_flag = true
                end
                if s_to > rating
                    push!(sm_vio,i)
                    #sm_vio += ((s_to - rating)/rating)*100
                    push!(sm_vio, ((s_to - rating)/rating)*100)
                    #vio_flag = true
                end
                #if vio_flag
                #    info(_LOGGER, "$(i), $(branch["f_bus"]), $(branch["t_bus"]): $(s_fr) / $(s_to) <= $(branch["rate_c"])")
                #end
            end
        end
    end
    #sm_vio = sm_vio/nb_branch
    return sm_vio, s_fr_arr, s_to_arr
end

function check_pg_pq_violations(network, solution)
    pg_vio = 0.0
    qg_vio = 0.0
    nb_gen = length(network["gen"])
    gen_qg_lim = []
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0  && network["bus"]["$(gen["gen_bus"])"]["bus_type"] != 4
            gen_sol = solution["gen"][i]

            if gen_sol["pg"] < gen["pmin"]
                if gen["pmax"] > 0
                    pg_vio += ((gen["pmin"] - gen_sol["pg"])/gen["pmax"])*100
                else
                    pg_vio += (gen["pmin"] - gen_sol["pg"])
                end
            end
            if gen_sol["pg"] > gen["pmax"]
                if gen["pmax"] > 0
                    pg_vio += ((gen_sol["pg"] - gen["pmax"])/gen["pmax"])*100
                else
                    pg_vio += (gen["pmin"] - gen_sol["pg"])
                end
            end

            if gen_sol["qg"] < gen["qmin"]
                push!(gen_qg_lim, i)
                #println(abs(gen["qmin"] - gen_sol["qg"]))
                qg_vio += (abs(gen["qmin"] - gen_sol["qg"])/gen["qmax"])*100
            end
            if gen_sol["qg"] > gen["qmax"]
                push!(gen_qg_lim, i)
                #println(gen_sol["qg"] - gen["qmax"])
                qg_vio += ((gen_sol["qg"] - gen["qmax"])/gen["qmax"])*100
            end
        end
    end
    pg_vio = pg_vio/nb_gen
    qg_vio = qg_vio/nb_gen
    return pg_vio, qg_vio, gen_qg_lim

end


function check_vm_violations(network, solution; vm_digits = 6)
    vm_vio_over = 0.0
    vm_vio_under = 0.0
    nb_vm = length(network["bus"])

    for (i,bus) in network["bus"]
        if bus["bus_type"] != 4
            bus_sol = solution["bus"][i]

            # helps to account for minor errors in equality constraints
            sol_val = round(bus_sol["vm"], digits=vm_digits)

            if sol_val < bus["vmin"]
                vm_vio_under += bus["vmin"] - sol_val
    
            end
            if sol_val > bus["vmax"]
                vm_vio_over += sol_val - bus["vmax"]
                
            end
            
        end
    end
    vm_vio_over = vm_vio_over/nb_vm
    vm_vio_under = vm_vio_under/nb_vm
    return vm_vio_over, vm_vio_under 
end

function check_vm_violations_bis(network, solution; vm_digits = 6)
    vm_vio_over = 0.0
    vm_vio_under = 0.0
    nb_vm = length(network["bus"])

    for (i,bus) in network["bus"]
        if bus["bus_type"] != 4
            bus_sol = solution["bus"][i]

            # helps to account for minor errors in equality constraints
            sol_val = round(bus_sol["vm"], digits=vm_digits)

            #vio_flag = false
            if sol_val < bus["vmin"]
                vm_vio_under += bus["vmin"] - sol_val
                #vio_flag = true
            end
            if sol_val > bus["vmax"]
                vm_vio_over += sol_val - bus["vmax"]
                #vio_flag = true
            end
            #if vio_flag
            #    info(_LOGGER, "$(i): $(bus["vmin"]) - $(sol_val) - $(bus["vmax"])")
            #end
        end
    end
    vm_vio_over = vm_vio_over
    vm_vio_under = vm_vio_under
    return vm_vio_over, vm_vio_under 
end

function plot_filtered_eigenvalues(total_eigen, min_index, max_index, real_limit::Float64)
    # Filter eigenvalues based on real part limit
    filtered_eigenvalues_min = filter(x -> real(x) >= real_limit, total_eigen[min_index])
    filtered_eigenvalues_max = filter(x -> real(x) >= real_limit, total_eigen[max_index])
    # Extract real and imaginary parts
    scatter(filtered_eigenvalues_min)
    scatter!(filtered_eigenvalues_max)
    
    # Function to calculate damping ratio lines
    function damping_ratio_lines(damping_ratio)
        real_parts = real_limit:0.1:0
        imag_parts = real_parts .* sqrt(1 - damping_ratio^2) / damping_ratio
        return real_parts, imag_parts
    end
    
    # Calculate damping ratio lines for 3% and 5%
    real_3, imag_3 = damping_ratio_lines(0.03)
    real_5, imag_5 = damping_ratio_lines(0.05)
    
    # Plot damping ratio lines
    plot!(real_3, imag_3, linestyle=:dash, label="3% Damping Ratio", color=:purple)
    plot!(real_3, -imag_3, linestyle=:dash, color=:purple)
    plot!(real_5, imag_5, linestyle=:dash, label="5% Damping Ratio", color=:orange)
    plot!(real_5, -imag_5, linestyle=:dash, color=:orange)
    
    # Set plot properties
    #xlabel!(L"\mathrm{Re}(\lambda)")  # LaTeX-formatted label for the real part of lambda
    #ylabel!(L"\mathrm{Im}(\lambda)")  # LaTeX-formatted label for the imaginary part of lambda
    title!("Filtered Eigenvalues and Damping Ratio Lines")
    
    # Save plot as PNG
    savefig("/Users/lolacharles/Downloads/Thesis/Code/filtered_eigenvalues_plot.png")
end

function build_c1_scopf_multinetwork_modif(network::Dict{String,<:Any})

contingencies = length(network["gen_contingencies"]) + length(network["branch_contingencies"])

#info(_LOGGER, "building scopf multi-network with $(contingencies+1) networks")

if contingencies > 0
    mn_data = replicate(network, contingencies)
    base_network = mn_data["nw"]["0"] = deepcopy(mn_data["nw"]["1"])

    for (n, network) in mn_data["nw"]
        if n == "0"
            continue
        end

        for (i,bus) in network["bus"]
            if haskey(bus, "evhi")
                bus["vmax"] = bus["evhi"]
            end
            if haskey(bus, "evlo")
                bus["vmin"] = bus["evlo"]
            end
        end

        for (i,branch) in network["branch"]
            if haskey(branch, "rate_c")
                branch["rate_a"] = branch["rate_c"]
            end
        end
    end

    network_id = 1
    for cont in base_network["gen_contingencies"]
        cont_nw = mn_data["nw"]["$(network_id)"]
        cont_nw["name"] = cont.label
        cont_gen = cont_nw["gen"]["$(cont.idx)"]
        cont_gen["gen_status"] = 0

        gen_buses = Set{Int}()
        for (i,gen) in cont_nw["gen"]
            if gen["gen_status"] != 0
                push!(gen_buses, gen["gen_bus"])
            end
        end
        cont_nw["gen_buses"] = gen_buses

        network["response_gens"] = Set()
        gen_bus = cont_nw["bus"]["$(cont_gen["gen_bus"])"]
        cont_nw["response_gens"] = []
        #cont_nw["response_gens"] = cont_nw["area_gens"][gen_bus["area"]]

        network_id += 1
    end
    for cont in base_network["branch_contingencies"]
        cont_nw = mn_data["nw"]["$(network_id)"]
        cont_nw["name"] = cont.label
        cont_branch = cont_nw["branch"]["$(cont.idx)"]
        cont_branch["br_status"] = 0

        gen_buses = Set{Int}()
        for (i,gen) in cont_nw["gen"]
            if gen["gen_status"] != 0
                push!(gen_buses, gen["gen_bus"])
            end
        end
        cont_nw["gen_buses"] = gen_buses

        fr_bus = cont_nw["bus"]["$(cont_branch["f_bus"])"]
        to_bus = cont_nw["bus"]["$(cont_branch["t_bus"])"]

        cont_nw["response_gens"] = Set()
        if haskey(cont_nw["area_gens"], fr_bus["area"])
            cont_nw["response_gens"] = cont_nw["area_gens"][fr_bus["area"]]
        end
        if haskey(network["area_gens"], to_bus["area"])
            cont_nw["response_gens"] = union(cont_nw["response_gens"], cont_nw["area_gens"][to_bus["area"]])
        end

        network_id += 1
    end

else
    mn_data = replicate(network, 1)
    mn_data["nw"]["0"] = mn_data["nw"]["1"]
    delete!(mn_data["nw"], "1")
end

return mn_data
end


function run_c1_scopf_modif(file, model_constructor, solver; kwargs...)
    return solve_model(file, model_constructor, solver, build_c1_scopf_modif; multinetwork=true, kwargs...)
end

# enables support for v[1], required for objective_variable_pg_cost when pg is an expression
Base.getindex(v::JuMP.GenericAffExpr, i::Int64) = v


""
function build_c1_scopf_modif(pm::AbstractPowerModel)
    # base-case network id is 0

    variable_bus_voltage(pm,bounded = false, nw=0)
    variable_gen_power(pm, bounded = false,nw=0)
    
    variable_branch_power(pm, bounded = false,nw=0)

    constraint_model_voltage(pm, nw=0)

    for i in ids(pm, :ref_buses, nw=0)
        constraint_theta_ref(pm, i, nw=0)
        constraint_voltage_magnitude_setpoint(pm, i, nw=0)

        # if multiple generators, fix power generation degeneracies
        if length(ref(pm,nw=0, :bus_gens, i)) > 1
            for j in collect(ref(pm, nw=0, :bus_gens, i))[2:end]
                constraint_gen_setpoint_active(pm, j, nw=0)
                constraint_gen_setpoint_reactive(pm, j, nw=0)
            end
        end
    end
     
    tolerance = 1e-6
    for i in ids(pm, :bus, nw=0)
        constraint_power_balance(pm, i, nw=0)
         # PV Bus Constraints
         if length(ref(pm,nw=0, :bus_gens, i)) > 0 && !(i in ids(pm, nw=0,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            for j in ref(pm, :bus_gens, i)
                constraint_gen_setpoint_active(pm, j, nw=0)
                
                gen = ref(pm, nw = 0, :gen, j)
                if gen["qg"] > gen["qmax"] || isapprox(gen["qg"], gen["qmax"], atol=tolerance)
                    constraint_gen_setpoint_reactive(pm, 0, gen["index"], gen["qmax"])
                elseif gen["qmin"] > gen["qg"] || isapprox(gen["qg"], gen["qmin"], atol=tolerance)
                    constraint_gen_setpoint_reactive(pm, 0, gen["index"], gen["qmin"])
                else
                    constraint_voltage_magnitude_setpoint(pm, i, nw=0)
                end

            end
        end
    end

    for i in ids(pm, :branch, nw=0)
        constraint_ohms_yt_from(pm, i, nw=0)
        constraint_ohms_yt_to(pm, i, nw=0)

    end


    contigency_ids = [id for id in nw_ids(pm) if id != 0]
    for nw in contigency_ids
        variable_bus_voltage(pm,bounded = false, nw=nw)
        variable_gen_power(pm, bounded = false,nw=nw)
        variable_branch_power(pm, bounded = false,nw=nw)

       


        constraint_model_voltage(pm, nw=nw)

        for i in ids(pm, :ref_buses, nw=nw)
            constraint_theta_ref(pm, i, nw=nw)
            constraint_voltage_magnitude_setpoint(pm, i, nw=nw)

            #if multiple generators, fix power generation degeneracies
            if length(ref(pm,nw=nw, :bus_gens, i)) > 1
               for j in collect(ref(pm,nw=nw, :bus_gens, i))[2:end]
                   constraint_gen_setpoint_active(pm, j, nw=nw)
                   constraint_gen_setpoint_reactive(pm, j, nw=nw)
               end
            end

        end

        gen_buses = ref(pm, :gen_buses, nw=nw)
        for i in ids(pm, :bus, nw=nw)
            constraint_power_balance(pm, i, nw=nw)

            if length(ref(pm,nw=nw, :bus_gens, i)) > 0 && !(i in ids(pm, nw=nw,:ref_buses))
                # this assumes inactive generators are filtered out of bus_gens
    
                #constraint_voltage_magnitude_setpoint(pm, i, nw=nw)
                for j in ref(pm, :bus_gens, i)
                
                    constraint_gen_setpoint_active(pm, j, nw=nw)
                    
                    gen = ref(pm, nw = nw, :gen, j)
                    if gen["qg"] > gen["qmax"] || isapprox(gen["qg"], gen["qmax"], atol=tolerance)
                        constraint_gen_setpoint_reactive(pm, nw, gen["index"], gen["qmax"])
                    elseif gen["qmin"] > gen["qg"] || isapprox(gen["qg"], gen["qmin"], atol=tolerance)
                        constraint_gen_setpoint_reactive(pm, nw, gen["index"], gen["qmin"])
                    else
                        constraint_voltage_magnitude_setpoint(pm, i, nw=nw)
                    end

    
                    
                end
            end
        end

        for i in ids(pm, :branch, nw=nw)
            constraint_ohms_yt_from(pm, i, nw=nw)
            constraint_ohms_yt_to(pm, i, nw=nw)

    
        end
    end


    ##### Setup Objective #####
    
end


function update_active_reactive_power_data!(network::Dict{String,<:Any}, data::Dict{String,<:Any}; branch_flow=false)
    for (i,bus) in data["bus"]
        nw_bus = network["bus"][i]
        nw_bus["va"] = bus["va"]
    end

    for (i,gen) in data["gen"]
        nw_gen = network["gen"][i]
        nw_gen["pg"] = gen["pg"]
        nw_gen["qg"] = gen["qg"]
    end

    if branch_flow
        for (i,branch) in data["branch"]
            nw_branch = network["branch"][i]
            nw_branch["pf"] = branch["pf"]
            nw_branch["pt"] = branch["pt"]
        end
    end
end
