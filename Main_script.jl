include("/Users/lolacharles/Downloads/Thesis/Code/Exported/module_methods.jl")
include("/Users/lolacharles/Downloads/Thesis/Code/Dummy/New_module.jl")
include("/Users/lolacharles/Downloads/Thesis/Code/N-1 tests/N_1_methods.jl")
include("/Users/lolacharles/Downloads/Special course/Small-signal_stability/SSA_module/Code_SSAmodule/SSA_module.jl")


#system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case3_lmbd.m"
#system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case5_pjm.m"
#system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case14_ieee.m"
system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case39_epri.m"
#system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case57_ieee.m"

network_data = PowerModels.parse_file(system_data)
pm, N, vars, header = instantiate_system_QCRM(network_data)
pg_numbers, vm_numbers = extract_number_and_type(vcat(header[1]))

#nb_load = length(network_data["load"])
#sample_array = gen_samples(100, nb_load, 0.7, 1.2)
#feasible, infeasible, non_feasible_index, feasible_index = OPF_feasible_samples(sample_array, network_data)
#sample_array_filtered = [sample_array[i,:] for i in eachindex(sample_array[:,1]) if !(i in non_feasible_index)]
#demand_profile = sample_array_filtered[2] 
#update_all_demand_reactive(network_data, demand_profile)
#writedlm("/Users/lolacharles/Downloads/Thesis/Code/Final/Load_demand_39.txt", sample_array_filtered )
#-----------------Start of the methodology----------------------------------------
global total_OPs = []
global total_Infeas_OPs = []


for demand in Demand_5
    update_all_demand_reactive(network_data, demand)
    #Parameters 
    Nb_HP = 40
    Nb_insamples = 10000

    #OBBT + Hyperplanes 
    OPs_HP, pf_results_HP, data_tight_HP, OPs_notFeas_HP, common_elements_HP, rest_of_the_indices_HP = infeasibility_certif_Constraints_check(network_data, Nb_HP, Nb_insamples)
    push!(total_OPs, OPs_HP)
end

Nb_HP = 40
Nb_insamples = 1000
update_all_demand(network_data, Demand_39[1])
@time OPs_HP, pf_results_HP, data_tight_HP, OPs_notFeas_HP, common_elements_HP, rest_of_the_indices_HP = infeasibility_certif_Constraints_check(network_data, Nb_HP, Nb_insamples)


df_DW_f = DataFrame(hcat(OPs_HP...)', Symbol.(header[1]))
df_DW_f.Feas = ones(nrow(df_DW_f))

df_DW_i = DataFrame(hcat(OPs_notFeas_HP...)', Symbol.(header[1]))
df_DW_i.Feas = zeros(nrow(df_DW_i))

# Check if column names match
if names(df_DW_f) != names(df_DW_i)
    throw(ArgumentError("The DataFrames have different column names"))
end

# Combine the DataFrames vertically
df_DW = vcat(df_DW_f, df_DW_i)

CSV.write("/Users/lolacharles/Downloads/Thesis/Code/Final/Datasets/DS_14bus_ACOPF_Demand1_10k.csv", df_DW)




Nb_MVD = 10000
MND_ops, MND_ops_Infeas = sample_MVND(OPs_HP, network_data, data_tight_HP, Nb_MVD)

df_DW_f = DataFrame(hcat(OPs_HP...)', Symbol.(header[1]))
df_DW_f.Feas = ones(nrow(df_DW_f))

df_DW_f_MVND = DataFrame(hcat(MND_ops...)', Symbol.(header[1]))
df_DW_f_MVND.Feas = ones(nrow(df_DW_f_MVND))

df_DW_i = DataFrame(hcat(OPs_notFeas_HP...)', Symbol.(header[1]))
df_DW_i.Feas = zeros(nrow(df_DW_i))

df_DW_i_MVND = DataFrame(hcat(MND_ops_Infeas...)', Symbol.(header[1]))
df_DW_i_MVND.Feas = zeros(nrow(df_DW_i_MVND))

# Check if column names match
if names(df_DW_f) != names(df_DW_i)
    throw(ArgumentError("The DataFrames have different column names"))
end

# Combine the DataFrames vertically
df_DW = vcat(df_DW_f, df_DW_i, df_DW_f_MVND, df_DW_i_MVND)

CSV.write("/Users/lolacharles/Downloads/Thesis/Code/Final/Datasets/DS_14bus_ACOPF_MVND_10k_Demand1.csv", df_DW)



# N-1 step
cont_out = ["1","3"]
cont_out = ["1", "13", "14"]
cont_out = ["5","14", "20","27", "32", "33", "34", "37", "39", "41","46" ]
df = CSV.read("/Users/lolacharles/Downloads/Thesis/Code/Final/DS_39bus_ACOPF.csv", DataFrame)
first_eight_columns = df[:, 1:19]
vector_of_arrays = [collect(row) for row in eachrow(first_eight_columns)]
#update_all_demand_reactive(network_data, Demand_39[1])
@time total_sm_array, total_qg_array, total_over_array, total_under_array , total_pg_array , Nb_N_1, resul, net = N_1_step(network_data, OPs_HP, cont_out)

threshold = 1e-6
indices_below_threshold = []
for (i, sub_vector) in enumerate(total_sm_array)
    #if mean(sub_vector) < 0.5

    if (all_below_threshold(total_sm_array[i], threshold) &&
        all_below_threshold(total_under_array[i], threshold) &&
        all_below_threshold(total_over_array[i], threshold) #&&
        #all_below_threshold(total_qg_array[i], threshold) &&
        #all_below_threshold(total_pg_array[i], threshold)
        )
        push!(indices_below_threshold, i)
    end
end

fli = filter(arr -> all_below_threshold(arr, 0.001), total_over_array)


sm_index = [mean(arr) for arr in total_sm_array]
over_index = [mean(arr) for arr in total_over_array]
under_index = [mean(arr) for arr in total_under_array]

# Create DataFrames
df_DW_f = DataFrame(hcat(OPs_HP...)', Symbol.(header[1]))
df_DW_f.Feas = ones(nrow(df_DW_f))

df_DW_i = DataFrame(hcat(OPs_notFeas_HP...)', Symbol.(header[1]))
df_DW_i.Feas = zeros(nrow(df_DW_i))

# Add N_1 column based on feasibility
df_DW_f.SM_index = sm_index[1:length(df_DW_f.Feas)]
df_DW_i.SM_index = -1 * ones(nrow(df_DW_i))

df_DW_f.Over_index = over_index[1:length(df_DW_f.Feas)]
df_DW_i.Over_index = -1 * ones(nrow(df_DW_i))

df_DW_f.Under_index = under_index[1:length(df_DW_f.Feas)]
df_DW_i.Under_index= -1 * ones(nrow(df_DW_i))


df_DW_f.N_1 = under_index[1:length(df_DW_f.Feas)]
df_DW_i.N_1 = -1 * ones(nrow(df_DW_i))

# Combine the DataFrames
df_DW = vcat(df_DW_f, df_DW_i)

df_DW.N_1 = zeros(nrow(df_DW))

# Set N_1 to 1 for feasible rows
for idx in indices_below_threshold
    df_DW.N_1[idx] = 1
end

# Check if column names match
if names(df_DW_f) != names(df_DW_i)
    throw(ArgumentError("The DataFrames have different column names"))
end

CSV.write("/Users/lolacharles/Downloads/Thesis/Code/Final/Datasets/DS_14bus_ACOPF_N_1_Demand1_10k.csv", df_DW)




# SSS step : 

global_OPs = vcat(OPs_HP, OPs_notFeas_HP)
index_res = vcat(common_elements_HP, rest_of_the_indices_HP)
dir_repo = "/Users/lolacharles/Downloads/Thesis/Code/HP_AND_SSS"
@time total_damp, total_dist, total_eigen = SSS_eval(data_tight_HP, global_OPs, pg_numbers, vm_numbers, dir_repo)
writedlm("/Users/lolacharles/Downloads/Thesis/Code/Final/Datasets/Total_damp_14_1000bus.txt", total_damp)

# Create DataFrames
df_DW_f = DataFrame(hcat(OPs_HP...)', Symbol.(header[1]))
df_DW_f.Feas = ones(nrow(df_DW_f))

df_DW_i = DataFrame(hcat(OPs_notFeas_HP...)', Symbol.(header[1]))
df_DW_i.Feas = zeros(nrow(df_DW_i))

# Add N_1 column based on feasibility
df_DW_f.SM_index = sm_index[1:length(df_DW_f.Feas)]
df_DW_i.SM_index = -1 * ones(nrow(df_DW_i))

df_DW_f.Over_index = over_index[1:length(df_DW_f.Feas)]
df_DW_i.Over_index = -1 * ones(nrow(df_DW_i))

df_DW_f.Under_index = under_index[1:length(df_DW_f.Feas)]
df_DW_i.Under_index= -1 * ones(nrow(df_DW_i))


df_DW_f.N_1 = under_index[1:length(df_DW_f.Feas)]
df_DW_i.N_1 = -1 * ones(nrow(df_DW_i))

# Combine the DataFrames
df_DW = vcat(df_DW_f, df_DW_i)

df_DW.N_1 = zeros(nrow(df_DW))

# Set N_1 to 1 for feasible rows
for idx in indices_below_threshold
    df_DW.N_1[idx] = 1
end

df_DW.SSS = total_damp[1:length(df_DW.Feas)]

# Check if column names match
if names(df_DW_f) != names(df_DW_i)
    throw(ArgumentError("The DataFrames have different column names"))
end

CSV.write("/Users/lolacharles/Downloads/Thesis/Code/Final/DS_14bus_ACOPF_SSS.csv", df_DW)


closet_OPS = closest_to_zero_indices(total_damp,10)

function closest_to_zero_indices_no_abs(arr, N::Int)
    # Create an array of tuples with absolute value and original index
    abs_with_index = [(val, idx) for (idx, val) in enumerate(arr)]
    
    # Sort the array by absolute value
    sorted_abs_with_index = sort(abs_with_index, by = x -> x[1])
    
    # Get the first N indices from the sorted array
    indices = [sorted_abs_with_index[i][2] for i in 1:N]
    
    return indices
end

closet_OPS = closest_to_zero_indices_no_abs(total_damp,10)
cls_OP = global_OPs[closet_OPS]

distance = [0.045, 0.02, 0.01]
#alpha = [0.4, 0.4, 0.05, 0.01]
#distance = [0.032, 0.02, 0.01]
alpha = [0.1, 0.01, 0.1, 0.01]
dw_ops = DW_step(data_tight_HP, index_res, closet_OPS,  pf_results_HP, distance, alpha)




closet_OPS_bis = [689, 100]
closest_to_zero_indices(total_damp,1)
cls_OP_bis = global_OPs[ 689]

closet_OPS_bis = closet_OPS[1:20]

distance = [0.045, 0.02, 0.01]
alpha = [0.4, 0.4, 0.05, 0.01]
dw_ops_bis = DW_step(data_tight_HP, index_res, closet_OPS_bis,  pf_results_HP, distance, alpha)

