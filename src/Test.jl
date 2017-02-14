include("Include.jl")

data_dictionary = DataDictionary(0.0,0.0,0.0)
(T,X) = repressed_simulation(0,240.0,0.1,data_dictionary)

#a_row = calculate_jacobian_row(T[20],vec(X[20,:]),6,data_dictionary)
#b_row = calculate_bmatrix_row(T[600],transpose(X[600,:]),7,data_dictionary)
BM = calculate_bmatrix(T[600],transpose(X[600,:]),data_dictionary)
#JM = calculate_jacobian(T[600],transpose(X[600,:]),data_dictionary)
#JMFD = finite_diff_jacobian(T[600],transpose(X[600,:]),data_dictionary)
