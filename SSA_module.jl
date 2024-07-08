using Pkg
using SIIPExamples
using PowerSystems
const PSY = PowerSystems
using TimeSeries
using PowerSimulationsDynamics
using Sundials
using Plots
using DataFrames
PSD = PowerSimulationsDynamics


#Compute the distance from the imaginary axis 
function dist_to_imaxis(eigenvalues)
    #Input = eigenvalues of the system -> Vector{Float64}
    real_eigen = real.(eigenvalues)
    imag_eigen = imag.(eigenvalues)
    if maximum(real_eigen) > 10^-8
        println("Not small-signal stable")
    else -10^-8 < maximum(real_eigen) < 10^-8
        zero_real_indices = findall(real_eigen .== 0)
        zero_real_and_imag_indices = filter(i -> imag_eigen[i] == 0, zero_real_indices)
        if isempty(zero_real_and_imag_indices)
            println("Warning : 0 but stable")
        else
            println("Warning : there is an eigenvalue at the origin")
        end
        real_eigen = [eig_val for eig_val in real_eigen if eig_val <= -10^-8]
    end
    dist = maximum(real_eigen)
    return dist
end

#Compute the minimum damping of the system
function min_damping(eigenvalues; damping_limit = 2)
    filtered_eigenval = [eig_val for eig_val in eigenvalues if !(-10^-8 <= real.(eig_val) <= 10^-8)]
    real_part_eigen = real.(filtered_eigenval)
    #real_part_eigen = [eig_val for eig_val in real_part_eigen if !(-10^-8 <= eig_val <= 10^-8)]
    im_part_eigen = imag.(filtered_eigenval)
    damping = []
    for i in eachindex(real_part_eigen)
        temp_damp = -real_part_eigen[i]/sqrt(real_part_eigen[i].^2+im_part_eigen[i].^2)
        if abs(temp_damp) <= damping_limit
            push!(damping,temp_damp)
        end
    end
    return minimum(unique(damping))
end



#Compute the two small-signal indices
function small_signal_module(dir_repo,file_stat,file_dyn)
    #Input : dir_repo -> folder path (String)
    #        file_stat -> .raw file (String)
    #        file_dyn -> .dyn file (String)
    output_file_path = joinpath(dir_repo, "output_SSanalysis_results.txt")
    output_file = open(output_file_path, "w")
    sys = System(file_stat, file_dyn)
    time_span = (0.0, 30.0)
    sim_studied = Simulation(ResidualModel, sys , pwd(), time_span)
    res_studied = small_signal_analysis(sim_studied)
    dist = dist_to_imaxis(res_studied.eigenvalues)
    damping = min_damping(res_studied.eigenvalues)
    println(output_file, "Distance: $dist, Damping: $damping")
    close(output_file)

end

function small_signal_module_m(dir_repo,sys,file_dyn)
    #Input : dir_repo -> folder path (String)
    #        file_stat -> .raw file (String)
    #        file_dyn -> .dyn file (String)
    output_file_path = joinpath(dir_repo, "output_SSanalysis_results.txt")
    output_file = open(output_file_path, "w")
    add_dyn_injectors!(sys, file_dyn)
    time_span = (0.0, 30.0)
    sim_studied = Simulation(ResidualModel, sys , pwd(), time_span)
    res_studied = small_signal_analysis(sim_studied)
    dist = dist_to_imaxis(res_studied.eigenvalues)
    damping = min_damping(res_studied.eigenvalues)
    println(output_file, "Distance: $dist, Damping: $damping")
    close(output_file)

end

function small_signal_module_sys(dir_repo,sys, nb; in_file=false)
    #Input : dir_repo -> folder path (String)
    #        file_stat -> .raw file (String)
    #        file_dyn -> .dyn file (String)
    
    time_span = (0.0, 30.0)
    sim_studied = Simulation(ResidualModel, sys , dir_repo, time_span)
    res_studied = small_signal_analysis(sim_studied)
    dist = dist_to_imaxis(res_studied.eigenvalues)
    damping = min_damping(res_studied.eigenvalues)
    eigenval = res_studied.eigenvalues
    if in_file == true
        output_file_path = joinpath(dir_repo, "output_SSanalysis_results_$nb.txt")
        output_file = open(output_file_path, "w")
        println(output_file, "Distance: $dist, Damping: $damping")
        close(output_file)
    end 
    return damping, dist, eigenval

end


#Path of the .raw file 
#stat = "/Users/lolacharles/Downloads/Dynamics Networks 4 PSSe/IEEE14_jconto/IEEE14_v33.raw"
#Path of the .dyn file
#dyn = "/Users/lolacharles/Downloads/Dynamics Networks 4 PSSe/IEEE14_jconto/ieee14.dyr"
#sys = System(stat)
#println(sys)
#bus_dict_gen = parse_dyr_components(dyn)
#ndyn = "/Users/lolacharles/Downloads/Special course/Small-signal_stability/SSA_module/Code_SSAmodule/11BUS_KUNDUR_TGOV.dyr"
#add_dyn_injectors!(sys, dyn)
#small_signal_module("/Users/lolacharles/Downloads/Special course/Small-signal_stability/SSA_module/Code_SSAmodule", stat, dyn )


function inertia_change_sys(dir_repo,file_stat,pourcentage,figure_name; construct_method =nothing, file_dyn=nothing,scatter_ind = false)
    #Input : dir_repo -> folder path (String)
    #        file_stat -> .raw file (String)
    #        file_dyn -> .dyn file (String), optional 
    #        pourcentage -> range of pourcentages for the values modification (Array{Float64})
    #        figure_name -> name of the plot constructed (String)
    #        construct_method -> method to construct the dynamic part of a system without a .dyr file, optional
    #        scatter_ind -> True = scatter plot, False = non-scatter plot
    output_file_path = joinpath(dir_repo, "output_Inertiaanalysis_results.txt")
    # Open the output file in write mode
    output_file = open(output_file_path, "w")
    x_axis_values = pourcentage*100
    damping_values = Float64[]
    distance_values = Float64[]
    for i in eachindex(pourcentage)
        if !isnothing(file_dyn)
            sys = System(file_stat, file_dyn)
        elseif !isnothing(construct_method)
            sys = construct_method(dir_repo,file_stat)
        else
            println("System couldn't be constructed")
        end

        for g in get_components(DynamicGenerator, sys)
            shaft = get_shaft(g)
            value = get_H(shaft)
            set_H!(shaft, value*pourcentage[i])
        end
        res = solve_powerflow(sys)
        print(res)
        time_span = (0.0, 30.0)
        sim_studied = PSD.Simulation(ResidualModel, sys , pwd(), time_span)
        res_studied = small_signal_analysis(sim_studied)
        dist = dist_to_imaxis(res_studied.eigenvalues)
        damping_arr = min_damping(res_studied.eigenvalues)
        damping_min = minimum(damping_arr)
        push!(damping_values, damping_min)
        push!(distance_values, dist)
        value =pourcentage[i]
        println(output_file, "Loading: $value, Dist: $dist, Damping: $damping_min") 
    end
    close(output_file)
    if scatter_ind
        scatter(x_axis_values, damping_values, marker = :circle, label = "Damping", xlabel = "Change of the system inertia [%]", ylabel = "Damping",xticks=(x_axis_values),xtickfontsize=6)
        # Create a twin y-axis for distance
        final_plot = scatter!(twinx(), x_axis_values, distance_values, marker = :square, markercolor = :green, label = "Distance", ylabel = "Distance",legend=:top,grid=true)
    else
        # Create a plot for damping
        plot(x_axis_values, damping_values, marker = :circle, label = "Damping", xlabel = "Change of the system inertia [%]", ylabel = "Damping",xticks=(x_axis_values),xtickfontsize=6)
        # Create a twin y-axis for distance
        final_plot = plot!(twinx(), x_axis_values, distance_values, marker = :square, linecolor = :green,markercolor = :green, label = "Distance", ylabel = "Distance",legend=:top,grid=true)
    end
    # Show the plot
    display(final_plot)
    savefig(final_plot, figure_name)
end


function loading_change_sys(dir_repo,file_stat,pourcentage,figure_name; construct_method =nothing, file_dyn=nothing,scatter_ind = false)
    #Input : dir_repo -> folder path (String)
    #        file_stat -> .raw file (String)
    #        file_dyn -> .dyn file (String), optional 
    #        pourcentage -> range of pourcentages for the values modification (Array{Float64})
    #        figure_name -> name of the plot constructed (String)
    #        construct_method -> method to construct the dynamic part of a system without a .dyr file, optional
    #        scatter_ind -> True = scatter plot, False = non-scatter plot
    output_file_path = joinpath(dir_repo, "output_analysis_results.txt")
    # Open the output file in write mode
    output_file = open(output_file_path, "w")
    x_axis_values = pourcentage*100
    damping_values = Float64[]
    distance_values = Float64[]
    for i in eachindex(pourcentage)
        if !isnothing(file_dyn)
            sys = System(file_stat, file_dyn)
        elseif !isnothing(construct_method)
            sys = construct_method(dir_repo,file_stat)
        else
            println("System couldn't be constructed")
        end
        #sys = System(file_stat, file_dyn)
        for c in get_components(PowerLoad, sys)
            value = get_active_power(c)
            set_active_power!(c,value*pourcentage[i])
        end
        for gen in get_components(ThermalStandard, sys)
            value = get_active_power(gen)
            set_active_power!(gen,value*pourcentage[i])
        end
        res = solve_powerflow(sys)
        print(res)
        time_span = (0.0, 30.0)
        sim_studied = PSD.Simulation(ResidualModel, sys , pwd(), time_span)
        res_studied = small_signal_analysis(sim_studied)
        dist = dist_to_imaxis(res_studied.eigenvalues)
        damping = min_damping(res_studied.eigenvalues)
        push!(damping_values, damping)
        push!(distance_values, dist)
        value =pourcentage[i]
        println(output_file, "Loading: $value, Dist: $dist, Damping: $damping") 
    end
    close(output_file)
    if scatter_ind
        scatter(x_axis_values, damping_values, marker = :circle, label = "Damping",legend=:topleft, xlabel = "Loading of the system [%]", ylabel = "Damping",xticks=(x_axis_values),xtickfontsize=6)
        # Create a twin y-axis for distance
        final_plot = scatter!(twinx(), x_axis_values, distance_values, marker = :square, markercolor = :green, label = "Distance", ylabel = "Distance",legend=:top,grid=true)
    else
        # Create a plot for damping
        plot(x_axis_values, damping_values, marker = :circle, label = "Damping",legend=:left, xlabel = "Loading of the system [%]", ylabel = "Damping",xticks=(x_axis_values),xtickfontsize=6)
        # Create a twin y-axis for distance
        final_plot = plot!(twinx(), x_axis_values, distance_values, marker = :square, linecolor = :green,markercolor = :green, label = "Distance", ylabel = "Distance",legend=:top,grid=true)

    end
    # Show the plot
    display(final_plot)
    savefig(final_plot, figure_name)
end

function study_variations_dynamic(dir_repo, file_stat, data_dict, number_variations, x_axis_nb, figure_name;file_path_dyn = nothing,construct_method = nothing , scatter_ind = false)
    #Input : dir_repo -> folder path (String)
    #        file_stat -> .raw file (String)
    #        file_dyn -> .dyn file (String), optional 
    #        data_dict -> dictionary where the modifications are stored(Dict{})
    #        number_variations -> number of modifications
    #        x_axis_nb -> set the modification whose values are going to be the values on the x-axis(Int)
    #        figure_name -> name of the plot constructed (String)
    #        construct_method -> method to construct the dynamic part of a system without a .dyr file, optional
    #        scatter_ind -> True = scatter plot, False = non-scatter plot
    output_file_path = joinpath(dir_repo, "output_analysisdynamic_results.txt")
    output_file = open(output_file_path, "w")
    x_axis_values = Float64[]
    damping_values = Float64[]
    distance_values = Float64[]
    if !isnothing(file_path_dyn)
        sys = System(file_stat, file_path_dyn)
    elseif !isnothing(construct_method)
        sys = construct_method(dir_repo,file_stat)
    else
        println("System couldn't be constructed")
    end
    for i in 1:number_variations
        for (index, key) in enumerate(keys(data_dict))
            new_value_array = data_dict[key]
            if index == x_axis_nb && i == 1
                push!(x_axis_values, new_value_array...)
            end
            generator, parameter = key
            dyn_gen = get_component(DynamicGenerator, sys, generator)
            shaft = get_shaft(dyn_gen)
            if parameter == "Inertia"
                set_H!(shaft, new_value_array[i])
            elseif parameter == "Damping"
                set_D!(shaft, new_value_array[i])
            end
            # You probably only need to reinitialize the simulation if changing
            # the parameter updates the operating condition. If not you can skip it
            if index == length(keys(data_dict))
                time_span = (0.0, 30.0)
                sim_studied = Simulation(ResidualModel, sys, mktempdir(), time_span)
                res_studied = small_signal_analysis(sim_studied)
                dist = dist_to_imaxis(res_studied.eigenvalues)
                damping = min_damping(res_studied.eigenvalues)
                push!(damping_values, damping)
                push!(distance_values, dist)
                value = new_value_array[i]
                println(output_file, "$generator $parameter: $value, Dist: $dist, Damping: $damping") 
            end
        end
    end
    close(output_file)
    if scatter_ind
        scatter(x_axis_values, damping_values, marker = :circle, label = "Damping", xlabel = "Load 1 Magnitude [MW]", ylabel = "Damping")
        # Create a twin y-axis for distance
        final_plot = scatter!(twinx(), x_axis_values, distance_values, marker = :square, markercolor = :green, label = "Distance", ylabel = "Distance",legend=:top,grid=true)
    else
        # Create a plot for damping
        plot(x_axis_values, damping_values, marker = :circle, label = "Damping", xlabel = "Inertia Generators 3 & 4", ylabel = "Damping",xticks=(x_axis_values),xtickfontsize=5, legend=:topright)
        # Create a twin y-axis for distance
        final_plot = plot!(twinx(), x_axis_values, distance_values, marker = :square, linecolor = :green,markercolor = :green, label = "Distance", ylabel = "Distance",legend=:top,grid=true)
    end
    # Show the plot
    display(final_plot)
    savefig(final_plot, figure_name)
end   

#Method to construct the WSCC9 model with inverter 
function construct_WSCC_9_inverter(dir_name,file_name)
    #Input : dir_repo -> folder path (String)
    #        file_name -> .raw file (String)
    WSCC_9_sys = System(joinpath(dir_name, file_name))

    slack_bus = [b for b in get_components(Bus, WSCC_9_sys) if b.bustype == BusTypes.REF][1]

    #Dynamic components creation
    machine_oneDoneQ() = OneDOneQMachine(
        0.0, #R
        0.146, #Xd
        0.0969, #Xq
        0.0608, #Xd_p
        0.0969, #Xq_p
        8.96, #Td0_p
        0.31, #Tq0_p
    )

    machine_twoDoneQ() = OneDOneQMachine(
        0.0, #R
        0.8958, #Xd
        0.8645, #Xq
        0.1198, #Xd_p
        0.1969, #Xq_p
        6.0, #Td0_p
        0.535, #Tq0_p
    )

    machine_threeDoneQ() = OneDOneQMachine(
        0.0, #R
        1.3125, #Xd
        1.2578, #Xq
        0.1813, #Xd_p
        0.25, #Xq_p
        5.89, #Td0_p
        0.6, #Tq0_p
    )
    #shaft_special() = FiveMassShaft(23.64,0,0,0,0,0.125,0,0,0,0,0,0,0,0,0,0,0,0)
    shaft_one() = SingleMass(
        23.64, #H 
        0.125, #D
    )

    shaft_two() = SingleMass(
        6.4, #H 
        0.033, #D
    )

    shaft_three() = SingleMass(
        3.01, #H 
        0.016, #D
    )

    avr_type1() = AVRTypeI(
        20.0, #Ka - Gain
        1.0, #Ke
        0.063, #Kf
        0.2, #Ta
        0.314, #Te
        0.35, #Tf
        0.001, #Tr
        (min = -5.0, max = 5.0),
        0.0039, #Ae - 1st ceiling coefficient
        1.555, #Be - 2nd ceiling coefficient
    )

    #No TG
    tg_none() = TGFixed(1.0) #efficiency

    #No PSS
    pss_none() = PSSFixed(0.0) #Vs

    converter_high_power = AverageConverter(rated_voltage = 13.80, rated_current = 100.0)
    outer_cont = OuterControl(
            VirtualInertia(Ta = 2.0, kd = 400.0, kω = 20.0),
            ReactivePowerDroop(kq = 0.2, ωf = 1000.0),
        )
    inner_cont = VoltageModeControl(
            kpv = 0.59,     #Voltage controller proportional gain
            kiv = 736.0,    #Voltage controller integral gain
            kffv = 0.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
            rv = 0.0,       #Virtual resistance in pu
            lv = 0.2,       #Virtual inductance in pu
            kpc = 1.27,     #Current controller proportional gain
            kic = 14.3,     #Current controller integral gain
            kffi = 0.0,     #Binary variable enabling the current feed-forward in output of current controllers
            ωad = 50.0,     #Active damping low pass filter cut-off frequency
            kad = 0.2,      #Active damping gain
        )   
    dc_source_lv = FixedDCSource(voltage = 600.0)   
    pll = KauraPLL(
            ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
            kp_pll = 0.084,  #PLL proportional gain
            ki_pll = 4.69,   #PLL integral gain
        )
    filt = LCLFilter(lf = 0.08, rf = 0.003, cf = 0.074, lg = 0.2, rg = 0.01)


    for g in get_components(Generator, WSCC_9_sys)
        #Find the generator at bus 102
        if get_number(get_bus(g)) == 1
            #Create the dynamic generator
            case_gen = DynamicGenerator(
                get_name(g),
                1.0, # ω_ref,
                machine_oneDoneQ(), #machine
                shaft_one(), #shaft
                avr_type1(), #avr
                tg_none(), #tg
                pss_none(), #pss
            )
            #Attach the dynamic generator to the system by specifying the dynamic and static part
            add_component!(WSCC_9_sys, case_gen, g)

        elseif get_number(get_bus(g)) == 2
            #Create the dynamic generator
            case_gen = DynamicGenerator(
                get_name(g),
                1.0, # ω_ref,
                machine_twoDoneQ(), #machine
                shaft_two(), #shaft
                avr_type1(), #avr
                tg_none(), #tg
                pss_none(), #pss
            )
            #Attach the dynamic generator to the system by specifying the dynamic and static part
            add_component!(WSCC_9_sys, case_gen, g)
        
        elseif get_number(get_bus(g)) == 3
            #Create the dynamic inverter
            case_inv = DynamicInverter(
                name = get_name(g),
                ω_ref = 1.0,
                converter = converter_high_power,
                outer_control = outer_cont,
                inner_control = inner_cont,
                dc_source = dc_source_lv,
                freq_estimator = pll,
                filter = filt,
            )
            #Attach the dynamic inverter to the system
            add_component!(WSCC_9_sys, case_inv, g)
        end
    end
    return WSCC_9_sys
end

#Method to construct the WSCC9 model 
function construct_WSCC_9(dir_name,file_name)
    #Input : dir_repo -> folder path (String)
    #        file_name -> .raw file (String)
    WSCC_9_sys = System(joinpath(dir_name, file_name))


    slack_bus = [b for b in get_components(Bus, WSCC_9_sys) if b.bustype == BusTypes.REF][1]

    #Dynamic components creation
    machine_oneDoneQ() = OneDOneQMachine(
        0.0, #R
        0.146, #Xd
        0.0969, #Xq
        0.0608, #Xd_p
        0.0969, #Xq_p
        8.96, #Td0_p
        0.31, #Tq0_p
    )

    machine_twoDoneQ() = OneDOneQMachine(
        0.0, #R
        0.8958, #Xd
        0.8645, #Xq
        0.1198, #Xd_p
        0.1969, #Xq_p
        6.0, #Td0_p
        0.535, #Tq0_p
    )

    machine_threeDoneQ() = OneDOneQMachine(
        0.0, #R
        1.3125, #Xd
        1.2578, #Xq
        0.1813, #Xd_p
        0.25, #Xq_p
        5.89, #Td0_p
        0.6, #Tq0_p
    )

    shaft_one() = SingleMass(
        23.64, #H 
        0.125, #D
    )

    shaft_two() = SingleMass(
        6.4, #H 
        0.033, #D
    )

    shaft_three() = SingleMass(
        3.01, #H 
        0.016, #D
    )


    avr_type1() = AVRTypeI(
        20.0, #Ka - Gain
        1.0, #Ke
        0.063, #Kf
        0.2, #Ta
        0.314, #Te
        0.35, #Tf
        0.001, #Tr
        (min = -5.0, max = 5.0),
        0.0039, #Ae - 1st ceiling coefficient
        1.555, #Be - 2nd ceiling coefficient
    )

    #No TG
    tg_none() = TGFixed(1.0) #efficiency

    #No PSS
    pss_none() = PSSFixed(0.0) #Vs

    for g in get_components(Generator, WSCC_9_sys)
        print(get_number(get_bus(g)))
        #Find the generator at bus 102
        if get_number(get_bus(g)) == 1
            #Create the dynamic generator
            case_gen = DynamicGenerator(
                get_name(g),
                1.0, # ω_ref,
                machine_oneDoneQ(), #machine
                shaft_one(), #shaft
                avr_type1(), #avr
                tg_none(), #tg
                pss_none(), #pss
            )
            #Attach the dynamic generator to the system by specifying the dynamic and static part
            add_component!(WSCC_9_sys, case_gen, g)

        elseif get_number(get_bus(g)) == 2
            #Create the dynamic generator
            case_gen = DynamicGenerator(
                get_name(g),
                1.0, # ω_ref,
                machine_twoDoneQ(), #machine
                shaft_two(), #shaft
                avr_type1(), #avr
                tg_none(), #tg
                pss_none(), #pss
            )
            #Attach the dynamic generator to the system by specifying the dynamic and static part
            add_component!(WSCC_9_sys, case_gen, g)
        
        elseif get_number(get_bus(g)) == 3
            #Create the dynamic generator
            case_gen = DynamicGenerator(
                get_name(g),
                1.0, # ω_ref,
                machine_threeDoneQ(), #machine
                shaft_three(), #shaft
                avr_type1(), #avr
                tg_none(), #tg
                pss_none(), #pss
            )
            #Attach the dynamic generator to the system by specifying the dynamic and static part
            add_component!(WSCC_9_sys, case_gen, g)
        end
    end
    return WSCC_9_sys
end



#Method to construct the WSCC9 model 
function construct_WSCC_9_sys(system_data)
    #Input : dir_repo -> folder path (String)
    #        file_name -> .raw file (String)
    WSCC_9_sys = System(system_data)


    slack_bus = [b for b in get_components(Bus, WSCC_9_sys) if b.bustype == BusTypes.REF][1]

    #Dynamic components creation
    machine_oneDoneQ() = OneDOneQMachine(
        0.0, #R
        0.146, #Xd
        0.0969, #Xq
        0.0608, #Xd_p
        0.0969, #Xq_p
        8.96, #Td0_p
        0.31, #Tq0_p
    )

    machine_twoDoneQ() = OneDOneQMachine(
        0.0, #R
        0.8958, #Xd
        0.8645, #Xq
        0.1198, #Xd_p
        0.1969, #Xq_p
        6.0, #Td0_p
        0.535, #Tq0_p
    )

    machine_threeDoneQ() = OneDOneQMachine(
        0.0, #R
        1.3125, #Xd
        1.2578, #Xq
        0.1813, #Xd_p
        0.25, #Xq_p
        5.89, #Td0_p
        0.6, #Tq0_p
    )

    shaft_one() = SingleMass(
        23.64, #H 
        0.125, #D
    )

    shaft_two() = SingleMass(
        6.4, #H 
        0.033, #D
    )

    shaft_three() = SingleMass(
        3.01, #H 
        0.016, #D
    )


    avr_type1() = AVRTypeI(
        20.0, #Ka - Gain
        1.0, #Ke
        0.063, #Kf
        0.2, #Ta
        0.314, #Te
        0.35, #Tf
        0.001, #Tr
        (min = -5.0, max = 5.0),
        0.0039, #Ae - 1st ceiling coefficient
        1.555, #Be - 2nd ceiling coefficient
    )
    avr_none() = AVRFixed(0.0)
    #No TG
    tg_none() = TGFixed(1.0) #efficiency

    #No PSS
    pss_none() = PSSFixed(0.0) #Vs

    for g in get_components(Generator, WSCC_9_sys)
        print(get_number(get_bus(g)))
        #Find the generator at bus 102
        if get_number(get_bus(g)) == 1
            #Create the dynamic generator
            case_gen = DynamicGenerator(
                get_name(g),
                1.0, # ω_ref,
                machine_oneDoneQ(), #machine
                shaft_one(), #shaft
                avr_none(), #avr
                tg_none(), #tg
                pss_none(), #pss
            )
            #Attach the dynamic generator to the system by specifying the dynamic and static part
            add_component!(WSCC_9_sys, case_gen, g)

        elseif get_number(get_bus(g)) == 2
            #Create the dynamic generator
            case_gen = DynamicGenerator(
                get_name(g),
                1.0, # ω_ref,
                machine_twoDoneQ(), #machine
                shaft_two(), #shaft
                avr_none(), #avr
                tg_none(), #tg
                pss_none(), #pss
            )
            #Attach the dynamic generator to the system by specifying the dynamic and static part
            add_component!(WSCC_9_sys, case_gen, g)
        
        elseif get_number(get_bus(g)) == 3
            #Create the dynamic generator
            case_gen = DynamicGenerator(
                get_name(g),
                1.0, # ω_ref,
                machine_threeDoneQ(), #machine
                shaft_three(), #shaft
                avr_none(),#avr_type1(), #avr
                tg_none(), #tg
                pss_none(), #pss
            )
            #Attach the dynamic generator to the system by specifying the dynamic and static part
            add_component!(WSCC_9_sys, case_gen, g)
        end
    end
    return WSCC_9_sys
end