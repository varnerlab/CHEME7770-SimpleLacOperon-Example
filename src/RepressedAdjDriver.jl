include("Include.jl")

# Script to solve the balance equations -
time_start = 0.0
time_stop = 240.0
time_step_size = 0.1
number_of_timesteps = length(time_start:time_step_size:time_stop)
number_of_states = 15

# Load the data dictionary (default parameter values) -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# Which parameter do we want to change?
# parameter_name_mapping_array = [
#   "n_gene_lacA_gene_lacI"	;		# 1
#   "K_gene_lacA_gene_lacI"	;		# 2
#   "n_gene_lacI_gene_system"	;	# 3
#   "K_gene_lacI_gene_system"	;	# 4
#   "n_gene_lacY_gene_lacI"	;		# 5
#   "K_gene_lacY_gene_lacI"	;		# 6
#   "n_gene_lacZ_gene_lacI"	;		# 7
#   "K_gene_lacZ_gene_lacI"	;		# 8
#
#   "W_gene_lacA_RNAP"	;					# 1 9
#   "W_gene_lacA_gene_lacI"	;			# 2	10
#   "W_gene_lacI_RNAP"	;					# 3 11
#   "W_gene_lacI_gene_system"	;		# 4	12
#   "W_gene_lacY_RNAP"	;					# 5	13
#   "W_gene_lacY_gene_lacI"	;			# 6	14
#   "W_gene_lacZ_RNAP"	;					# 7	15
#   "W_gene_lacZ_gene_lacI"	;			# 8	16
#   "W_gene_system_RNAP"	;				# 9	17
#
#   "rnapII_concentration"	; 		# 1 18
#   "ribosome_concentration"	;		# 2 19
#   "degradation_constant_mRNA"	;	# 3 20
#   "degradation_constant_protein"	;	# 4 21
#   "kcat_transcription"	;				# 5 22
#   "kcat_translation"	;					# 6 23
#   "maximum_specific_growth_rate"	;	# 7 24
#   "saturation_constant_transcription"	;	# 8 25
#   "saturation_constant_translation"	;	# 9 26
# ]


parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
average_scaled_sensitivity_array = zeros(number_of_states,1)
for (parameter_index,parameter_value) in enumerate(parameter_name_mapping_array)

  local_data_dictionary = deepcopy(data_dictionary)

  @show parameter_index

  # Solve the adj balances -
  (T,X) = repressed_adj_simulation(time_start,time_stop,time_step_size,parameter_index,local_data_dictionary)

  # split -
  state_array = X[:,1:15]
  sensitivity_array = X[:,16:end]
  scaled_sensitivity_array = scale_sensitivity_array(T,state_array,sensitivity_array,parameter_index,local_data_dictionary)

  # time average -
  average_sensitivity_col = time_average_array(T,scaled_sensitivity_array)

  # grab -
  average_scaled_sensitivity_array = [average_scaled_sensitivity_array average_sensitivity_col]
end

average_scaled_sensitivity_array = average_scaled_sensitivity_array[:,2:end]
