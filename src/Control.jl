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
# Function: Control
# Description: Calculate the transcriptional control array at time t
# Generated on: 2017-02-02T08:53:57.656
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function Control(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# initialize the control -
	control_array = zeros(5)

	# Alias the species -
	gene_lacA = x[1]
	gene_lacI = x[2]
	gene_lacY = x[3]
	gene_lacZ = x[4]
	gene_system = x[5]
	mRNA_gene_lacA = x[6]
	mRNA_gene_lacI = x[7]
	mRNA_gene_lacY = x[8]
	mRNA_gene_lacZ = x[9]
	mRNA_gene_system = x[10]
	protein_gene_lacA = x[11]
	protein_gene_lacI = x[12]
	protein_gene_lacY = x[13]
	protein_gene_lacZ = x[14]
	protein_gene_system = x[15]

	# Alias the binding parameters -
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_gene_lacA_gene_lacI = binding_parameter_dictionary["n_gene_lacA_gene_lacI"]
	K_gene_lacA_gene_lacI = binding_parameter_dictionary["K_gene_lacA_gene_lacI"]
	n_gene_lacI_gene_system = binding_parameter_dictionary["n_gene_lacI_gene_system"]
	K_gene_lacI_gene_system = binding_parameter_dictionary["K_gene_lacI_gene_system"]
	n_gene_lacY_gene_lacI = binding_parameter_dictionary["n_gene_lacY_gene_lacI"]
	K_gene_lacY_gene_lacI = binding_parameter_dictionary["K_gene_lacY_gene_lacI"]
	n_gene_lacZ_gene_lacI = binding_parameter_dictionary["n_gene_lacZ_gene_lacI"]
	K_gene_lacZ_gene_lacI = binding_parameter_dictionary["K_gene_lacZ_gene_lacI"]

	# Alias the control function parameters -
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_gene_lacA_RNAP = control_parameter_dictionary["W_gene_lacA_RNAP"]
	W_gene_lacA_gene_lacI = control_parameter_dictionary["W_gene_lacA_gene_lacI"]
	W_gene_lacI_RNAP = control_parameter_dictionary["W_gene_lacI_RNAP"]
	W_gene_lacI_gene_system = control_parameter_dictionary["W_gene_lacI_gene_system"]
	W_gene_lacY_RNAP = control_parameter_dictionary["W_gene_lacY_RNAP"]
	W_gene_lacY_gene_lacI = control_parameter_dictionary["W_gene_lacY_gene_lacI"]
	W_gene_lacZ_RNAP = control_parameter_dictionary["W_gene_lacZ_RNAP"]
	W_gene_lacZ_gene_lacI = control_parameter_dictionary["W_gene_lacZ_gene_lacI"]
	W_gene_system_RNAP = control_parameter_dictionary["W_gene_system_RNAP"]

	# We need to call to the signaling function to calculate any "active" species -
	protein_gene_lacI = calculate_active_species(t,x,data_dictionary)
	
	# Transfer function target:gene_lacA actor:gene_lacI
	actor_set_gene_lacA_gene_lacI = [
		protein_gene_lacI
	]
	actor = prod(actor_set_gene_lacA_gene_lacI)
	b_gene_lacA_gene_lacI = (actor^(n_gene_lacA_gene_lacI))/(K_gene_lacA_gene_lacI^(n_gene_lacA_gene_lacI)+actor^(n_gene_lacA_gene_lacI))

	# Control function for gene_lacA -
	control_array[1] = (W_gene_lacA_RNAP)/(1+W_gene_lacA_RNAP+W_gene_lacA_gene_lacI*b_gene_lacA_gene_lacI)

	# Transfer function target:gene_lacI actor:gene_system
	actor_set_gene_lacI_gene_system = [
		protein_gene_system
	]
	actor = prod(actor_set_gene_lacI_gene_system)
	b_gene_lacI_gene_system = (actor^(n_gene_lacI_gene_system))/(K_gene_lacI_gene_system^(n_gene_lacI_gene_system)+actor^(n_gene_lacI_gene_system))

	# Control function for gene_lacI -
	control_array[2] = (W_gene_lacI_RNAP+W_gene_lacI_gene_system*b_gene_lacI_gene_system)/(1+W_gene_lacI_RNAP+W_gene_lacI_gene_system*b_gene_lacI_gene_system)

	# Transfer function target:gene_lacY actor:gene_lacI
	actor_set_gene_lacY_gene_lacI = [
		protein_gene_lacI
	]
	actor = prod(actor_set_gene_lacY_gene_lacI)
	b_gene_lacY_gene_lacI = (actor^(n_gene_lacY_gene_lacI))/(K_gene_lacY_gene_lacI^(n_gene_lacY_gene_lacI)+actor^(n_gene_lacY_gene_lacI))

	# Control function for gene_lacY -
	control_array[3] = (W_gene_lacY_RNAP)/(1+W_gene_lacY_RNAP+W_gene_lacY_gene_lacI*b_gene_lacY_gene_lacI)

	# Transfer function target:gene_lacZ actor:gene_lacI
	actor_set_gene_lacZ_gene_lacI = [
		protein_gene_lacI
	]
	actor = prod(actor_set_gene_lacZ_gene_lacI)
	b_gene_lacZ_gene_lacI = (actor^(n_gene_lacZ_gene_lacI))/(K_gene_lacZ_gene_lacI^(n_gene_lacZ_gene_lacI)+actor^(n_gene_lacZ_gene_lacI))

	# Control function for gene_lacZ -
	control_array[4] = (W_gene_lacZ_RNAP)/(1+W_gene_lacZ_RNAP+W_gene_lacZ_gene_lacI*b_gene_lacZ_gene_lacI)

	# Control function for gene_system -
	control_array[5] = (W_gene_system_RNAP)/(1+W_gene_system_RNAP)

	# return -
	return control_array
end
