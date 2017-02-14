function time_average_array(time_array,data_array)

  # what is the delta T?
  delta_time = (time_array[end] - time_array[1])

  # initialize -
  average_array = Float64[]

  # what is the size of the array?
  (number_of_timesteps,number_of_states) = size(data_array)
  for state_index = 1:number_of_states

    # grab the data column -
    data_col = data_array[:,state_index]

    # average -
    average_value = (1/delta_time)*trapz(time_array,data_col)

    # push -
    push!(average_array,average_value)
  end

  return average_array

end

function scale_sensitivity_array(time_array,state_array,sensitivity_array,parameter_index,data_dictionary)

  # what is small?
  epsilon = 1e-6

  # initialize -
  (number_of_timesteps,number_of_states) = size(state_array)
  scaled_sensitivity_array = zeros(number_of_timesteps,number_of_states)

  # What is the nominal parameter value?
  parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
  parameter_name = parameter_name_mapping_array[parameter_index]

  # build the total parameter dictionary -
  binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  misc_parameter_dictionary = data_dictionary["misc_parameter_dictionary"]
  total_parameter_dictionary = merge(binding_parameter_dictionary,control_parameter_dictionary,misc_parameter_dictionary)

  # Grab the default value -
  default_parameter_value = total_parameter_dictionary[parameter_name]

  # main loop -
  for (time_index,time_value) in enumerate(time_array)

    for state_index = 1:number_of_states

      state_value = state_array[time_index,state_index]
      if (state_value<epsilon)
        state_value = epsilon
      end

      old_value = sensitivity_array[time_index,state_index]
      new_value = old_value*(default_parameter_value/state_value)
      scaled_sensitivity_array[time_index,state_index] = new_value
    end
  end

  return scaled_sensitivity_array
end

function calculate_jacobian(time,state_array,data_dictionary)

  # what is the size of the system?
  number_of_states = length(state_array)

  # calculate each row of the jacobian -
  jacobian_array = zeros(1,number_of_states)
  for (state_index,state_value) in enumerate(state_array)

    jacobian_row = calculate_jacobian_row(time,state_array,state_index,data_dictionary)
    jacobian_array = [jacobian_array  ; transpose(jacobian_row)]
  end

  jacobian_array = jacobian_array[2:end,:]
  return jacobian_array
end

function calculate_bmatrix(time,state_array,data_dictionary)

  # what is the size of the system?
  parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
  number_of_parameters = length(parameter_name_mapping_array)
  number_of_states = length(state_array)

  # calculate each row of the jacobian -
  b_array = zeros(number_of_parameters,1)
  for (state_index,state_value) in enumerate(state_array)

    bmatrtix_row = calculate_bmatrix_row(time,state_array,state_index,data_dictionary)
    b_array = [b_array bmatrtix_row]
  end

  b_array = b_array[:,2:end]
  return transpose(b_array)
end

function calculate_bmatrix_row(time,state_array,balance_index,data_dictionary)

  # define some constants -
  const epsilon = 1e-6
  const delta = 0.01

  # parameter name dictionary -
  parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
  binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  misc_parameter_dictionary = data_dictionary["misc_parameter_dictionary"]
  number_of_parameters = length(parameter_name_mapping_array)

  # create a mega dictionary -
  total_parameter_dictionary = merge(binding_parameter_dictionary,control_parameter_dictionary,misc_parameter_dictionary)

  # create delta parameter array -
  parameter_delta_array = Float64[]
  for (parameter_index,parameter_name) in enumerate(parameter_name_mapping_array)

    # Grab the default value -
    default_parameter_value = total_parameter_dictionary[parameter_name]

    # state perturbation -
    peturbed_parameter_value = default_parameter_value*(delta);

    #@show (peturbed_parameter_value,default_parameter_value)

    # check -
    if (peturbed_parameter_value<epsilon)
      peturbed_parameter_value = epsilon
    end

    # capture -
    push!(parameter_delta_array,(peturbed_parameter_value))
  end

  # Create the diag array -
  diag_delta_array = diagm(vec(parameter_delta_array))

  # Create bVec -
  f_nominal = Balances(time,vec(state_array),data_dictionary)

  # estimate the perturbed balances -
  rhs_delta_array = Float64[]
  for (parameter_index,parameter_name) in enumerate(parameter_name_mapping_array)

    # copy -
    local_data_dictionary = deepcopy(data_dictionary)

    # Grab the default value -
    default_parameter_value = total_parameter_dictionary[parameter_name]

    # update the state -
    perturbed_parameter_array = zeros(number_of_parameters)
    for local_index = 1:number_of_parameters
      if (parameter_index == local_index)
        perturbed_parameter_array[parameter_index] = default_parameter_value*(1+delta);
      else
        local_parameter_name = parameter_name_mapping_array[local_index]
        perturbed_parameter_array[local_index] = total_parameter_dictionary[local_parameter_name]
      end
    end

    if (parameter_index<=8)

      # we are in the binding section -
      local_data_dictionary["binding_parameter_dictionary"][parameter_name] = perturbed_parameter_array[parameter_index]
    elseif (parameter_index>8 && parameter_index<=17)

      # we are in the control section -
      local_data_dictionary["control_parameter_dictionary"][parameter_name] = perturbed_parameter_array[parameter_index]
    else
      # we are in the bar parameter section -
      local_data_dictionary[parameter_name] = perturbed_parameter_array[parameter_index]
    end


    # calculate the perturbed balances -
    f_perturbed = Balances(time,vec(state_array),local_data_dictionary)
    f_perturbed = f_perturbed[balance_index] - f_nominal[balance_index]

    # capture -
    push!(rhs_delta_array,f_perturbed)
  end

  # calculate the bmatrix row -
  bmatrix_row = diag_delta_array\rhs_delta_array

  # return -
  return bmatrix_row
end

function finite_diff_jacobian(time,state_array,data_dictionary)

  # define some constants -
  const epsilon = 1e-6
  const delta = 0.05
  number_of_states = length(state_array)

  # initialize -
  jacobian_array = zeros(number_of_states,number_of_states)

  # nominal -
  f_nominal = Balances(time,vec(state_array),data_dictionary)

  #@show state_array

  for row_index = 1:number_of_states
    for col_index = 1:number_of_states

      perturbed_state_array = zeros(number_of_states)
      for perturbation_index = 1:number_of_states

        if (col_index == perturbation_index)
          perturbed_state_array[col_index] = state_array[col_index]*(1+delta)
        else
          perturbed_state_array[perturbation_index] = state_array[perturbation_index]
        end
      end

      #@show perturbed_state_array

      # calculate the balances -
      f_perturbed = Balances(time,vec(perturbed_state_array),data_dictionary)

      # calculate the entry -
      perturbation_size = state_array[col_index]*delta
      if (perturbation_size<epsilon)
        perturbation_size = epsilon
      end
      jacobian_array[row_index,col_index] = (f_perturbed[row_index] - f_nominal[row_index])/(perturbation_size)
    end
  end

  return jacobian_array
end

function calculate_jacobian_row(time,state_array,balance_index,data_dictionary)

  # define some constants -
  const epsilon = 1e-6
  const delta = 0.001
  number_of_states = length(state_array)

  # Create the delta state array -
  state_delta_array = Float64[]
  for (state_index,state_value) in enumerate(state_array)

    # state perturbation -
    peturbed_state = state_value*(delta);

    # check -
    if (peturbed_state<epsilon)
      peturbed_state = epsilon
    end

    # capture -
    push!(state_delta_array,peturbed_state)
  end

  # Create the diag array -
  diag_delta_array = diagm(vec(state_delta_array))

  # Create bVec -
  f_nominal = Balances(time,vec(state_array),data_dictionary)

  # estimate the perturbed balances -
  rhs_delta_array = Float64[]
  for (state_index,state_value) in enumerate(state_array)

    # update the state -
    perturbed_state_array = zeros(number_of_states)

    for local_index = 1:number_of_states
      if (state_index == local_index)
        perturbed_state_array[state_index] = state_value*(1+delta);
      else
        perturbed_state_array[local_index] = state_array[local_index]
      end
    end


    # calculate the perturbed balances -
    f_perturbed = Balances(time,vec(perturbed_state_array),data_dictionary)
    f_perturbed_delta = f_perturbed[balance_index] - f_nominal[balance_index]

    # capture -
    push!(rhs_delta_array,f_perturbed_delta)
  end

  # calculate the jacobian row -
  jacobian_row = diag_delta_array\rhs_delta_array

  # return -
  return jacobian_row
end

function trapz{Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end
    r = zero(zero(Tx) + zero(Ty))
    if n == 1; return r; end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    #= correction -h^2/12 * (f'(b) - f'(a))
    ha = x[2] - x[1]
    he = x[end] - x[end-1]
    ra = (y[2] - y[1]) / ha
    re = (y[end] - y[end-1]) / he
    r/2 - ha*he/12 * (re - ra)
    =#
    return r/2
end
