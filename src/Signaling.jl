function calculate_active_species(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

  # get data from the data_dictionary -
  lactate_abundance = data_dictionary["lactate_abundance"]
  lactate_order = data_dictionary["lactate_order"]
  Kd = data_dictionary["lactate_affinity"]

  # Calculate the correction -
  phi = (lactate_abundance^(lactate_order))/(lactate_abundance^(lactate_order)+Kd^(lactate_order))

  # Apply the correction -
  protein_gene_lacI = x[12]
  protein_gene_lacI = protein_gene_lacI*(1-phi)

  # return the "active" species -
  return protein_gene_lacI
end
