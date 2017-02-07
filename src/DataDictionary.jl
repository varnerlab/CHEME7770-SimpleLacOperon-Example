# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2017-02-02T08:53:57.086
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start::Float64,time_stop::Float64,time_step_size::Float64)

	# stoichiometric_matrix and dilution_matrix -
	stoichiometric_matrix = readdlm("./Network.dat")
	dilution_matrix = readdlm("./Dilution.dat")
	degradation_matrix = readdlm("./Degradation.dat")

	# array of gene lengths -
	gene_coding_length_array = [
		609		;	# 1	gene_lacA
		1080	;	# 2	gene_lacI
		1251	;	# 3	gene_lacY
		3075	;	# 4	gene_lacZ
		1200	;	# 5	gene_system
	]

	# array of mRNA coding lengths -
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	# 6	1	mRNA_gene_lacA
		gene_coding_length_array[2]	;	# 7	2	mRNA_gene_lacI
		gene_coding_length_array[3]	;	# 8	3	mRNA_gene_lacY
		gene_coding_length_array[4]	;	# 9	4	mRNA_gene_lacZ
		gene_coding_length_array[5]	;	# 10	5	mRNA_gene_system
	]

	# array of mRNA coding lengths -
	protein_coding_length_array = [
		203																				;	# 11	1	protein_gene_lacA
		360																				;	# 12	2	protein_gene_lacI
		417																				;	# 13	3	protein_gene_lacY
		1025																			;	# 14	4	protein_gene_lacZ
		round((0.33)*mRNA_coding_length_array[5])	;	# 15	5	protein_gene_system
	]

	# Default: HL-60, we've updated these numbers from bionumbers -
	# ------------------------------------------------------------------------------------------#
	# constants (from bionumbers)       units
	# ------------------------------------------------------------------------------------------#
	cell_diameter = 1.1                 # mum
	number_of_rnapII = 4600            	# copies/cells
	number_of_ribosome = 50000         	# copies/cells
	mRNA_half_life_TF = 0.083           # hrs
	protein_half_life = 70              # hrs
	doubling_time_cell = 0.33           # hrs
	max_translation_rate = 16.5         # aa/sec
	max_transcription_rate = 60.0       # nt/sec
	average_transcript_length = 1200   	# nt
	average_protein_length = 300       	# aa
	fraction_nucleus = 0.0             	# dimensionless
	av_number = 6.02e23                 # number/mol
	avg_gene_number = 2                 # number of copies of a gene
	polysome_number = 10								# number of ribsomoses per transcript
	# ------------------------------------------------------------------------------------------#
	#
	# ------------------------------------------------------------------------------------------#
	# Calculate constants using bionumber values
	# ------------------------------------------------------------------------------------------#
	# Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/V)*1e9                   # nM
	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/V)*1e9               # nM

	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                       # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(0.5)                    # hr^-1

	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)      # hr^-1
	kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)             # hr^-1

	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                      # hr^-1

	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/av_number)*(1/V)*1e9                  # nM

	# How fast do my cells die?
	death_rate_constant = 0.2*maximum_specific_growth_rate                            # hr^-1

	# Saturation constants for translation and trascription -
	saturation_transcription = 4600*(1/av_number)*(1/V)*1e9                           # nM
	saturation_translation = 100000*(1/av_number)*(1/V)*1e9                           # nM
	# -------------------------------------------------------------------------------------------#

	# initial condition array -
	initial_condition_array = [
		avg_gene_concentration	;	# 1	gene_lacA
		avg_gene_concentration	;	# 2	gene_lacI
		avg_gene_concentration	;	# 3	gene_lacY
		avg_gene_concentration	;	# 4	gene_lacZ
		avg_gene_concentration	;	# 5	gene_system

		0.0	;	# 6	mRNA_gene_lacA
		0.0	;	# 7	mRNA_gene_lacI
		0.0	;	# 8	mRNA_gene_lacY
		0.0	;	# 9	mRNA_gene_lacZ
		0.0	;	# 10	mRNA_gene_system

		0.0	;	# 11	protein_gene_lacA
		0.0	;	# 12	protein_gene_lacI
		0.0	;	# 13	protein_gene_lacY
		0.0	;	# 14	protein_gene_lacZ
		0.0	;	# 15	protein_gene_system
	]

	binding_parameter_dictionary = Dict{AbstractString,Float64}()
	binding_parameter_dictionary["n_gene_lacA_gene_lacI"] = 1.0
	binding_parameter_dictionary["K_gene_lacA_gene_lacI"] = 120.0
	binding_parameter_dictionary["n_gene_lacI_gene_system"] = 1.0
	binding_parameter_dictionary["K_gene_lacI_gene_system"] = 120.0
	binding_parameter_dictionary["n_gene_lacY_gene_lacI"] = 1.0
	binding_parameter_dictionary["K_gene_lacY_gene_lacI"] = 120.0
	binding_parameter_dictionary["n_gene_lacZ_gene_lacI"] = 1.0
	binding_parameter_dictionary["K_gene_lacZ_gene_lacI"] = 120.0

	# Alias the control function parameters -
	control_parameter_dictionary = Dict{AbstractString,Float64}()
	control_parameter_dictionary["W_gene_lacA_RNAP"] = 0.01
	control_parameter_dictionary["W_gene_lacA_gene_lacI"] = 10.0
	control_parameter_dictionary["W_gene_lacI_RNAP"] = 0.0
	control_parameter_dictionary["W_gene_lacI_gene_system"] = 0.01
	control_parameter_dictionary["W_gene_lacY_RNAP"] = 0.01
	control_parameter_dictionary["W_gene_lacY_gene_lacI"] = 10.0
	control_parameter_dictionary["W_gene_lacZ_RNAP"] = 0.01
	control_parameter_dictionary["W_gene_lacZ_gene_lacI"] = 10.0
	control_parameter_dictionary["W_gene_system_RNAP"] = 0.0

	# Parameter name index array -
	parameter_name_mapping_array = [
		"n_gene_lacA_gene_lacI"	;	# 1
		"K_gene_lacA_gene_lacI"	;	# 2
		"n_gene_lacI_gene_system"	;	# 3
		"K_gene_lacI_gene_system"	;	# 4
		"n_gene_lacY_gene_lacI"	;	# 5
		"K_gene_lacY_gene_lacI"	;	# 6
		"n_gene_lacZ_gene_lacI"	;	# 7
		"K_gene_lacZ_gene_lacI"	;	# 8
		"W_gene_lacA_RNAP"	;	# 9
		"W_gene_lacA_gene_lacI"	;	# 10
		"W_gene_lacI_RNAP"	;	# 11
		"W_gene_lacI_gene_system"	;	# 12
		"W_gene_lacY_RNAP"	;	# 13
		"W_gene_lacY_gene_lacI"	;	# 14
		"W_gene_lacZ_RNAP"	;	# 15
		"W_gene_lacZ_gene_lacI"	;	# 16
		"W_gene_system_RNAP"	;	# 17
	]


	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["average_transcript_length"] = average_transcript_length
	data_dictionary["average_protein_length"] = average_protein_length
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	data_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
	data_dictionary["ribosome_concentration"] = ribosome_concentration # muM
	data_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
	data_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
	data_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
	data_dictionary["kcat_translation"] = kcat_translation  # hr^-1
	data_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
	data_dictionary["death_rate_constant"] = death_rate_constant
	data_dictionary["avg_gene_concentration"] = avg_gene_concentration
	data_dictionary["saturation_constant_transcription"] = saturation_transcription
	data_dictionary["saturation_constant_translation"] = saturation_translation

	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["dilution_matrix"] = dilution_matrix
	data_dictionary["degradation_matrix"] = degradation_matrix

	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	data_dictionary["parameter_name_mapping_array"] = parameter_name_mapping_array

	# Signaling data added by JV -
	data_dictionary["lactate_abundance"] = 0.0
	data_dictionary["lactate_order"] = 2.0
  data_dictionary["lactate_affinity"] = 2.0
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
