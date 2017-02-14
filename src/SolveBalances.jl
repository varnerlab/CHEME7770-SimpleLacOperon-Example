# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# School of Chemical Engineering Purdue University
# W. Lafayette IN 46907 USA

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
# SolveBalances: Solves model equations from TSTART to TSTOP given parameters in data_dictionary.
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# TSTART  - Time start
# TSTOP  - Time stop
# Ts - Time step
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# TSIM - Simulation time vector
# X - Simulation state array (NTIME x NSPECIES)
# ----------------------------------------------------------------------------------- #
function SolveAdjBalances(TSTART,TSTOP,Ts,parameter_index,data_dictionary)

  # Get required stuff from DataFile struct -
  TSIM = collect(TSTART:Ts:TSTOP);
  initial_condition_vector = data_dictionary["initial_condition_array"];

  # Call the ODE solver -
  fbalances(t,y) = AdjBalances(t,y,parameter_index,data_dictionary);
  (t,y) = ode45(fbalances,initial_condition_vector,TSIM;points=:specified,minstep=1e-1000);

  # Map -
  number_of_timesteps = length(t)
  number_of_states = length(initial_condition_vector)
  X = zeros(number_of_timesteps,number_of_states)
  for state_index = 1:number_of_states
    tmp = map(y->y[state_index],y)
    for time_index = 1:number_of_timesteps
      X[time_index,state_index] = tmp[time_index]
    end
  end

  # # Check for smalls -
  epsilon = 1e-4
  idx_n = find(abs(X).<epsilon)
  X[idx_n] = 0.0

  # return time and state -
  return (t,X);
end

function SolveBalances(TSTART,TSTOP,Ts,data_dictionary)

  # Get required stuff from DataFile struct -
  TSIM = collect(TSTART:Ts:TSTOP);
  initial_condition_vector = data_dictionary["initial_condition_array"];

  # Call the ODE solver -
  fbalances(t,y) = Balances(t,y,data_dictionary);
  (t,y) = ode23s(fbalances,initial_condition_vector,TSIM;points=:specified);

  # Map -
  number_of_timesteps = length(t)
  number_of_states = length(initial_condition_vector)
  X = zeros(number_of_timesteps,number_of_states)
  for state_index = 1:number_of_states
    tmp = map(y->y[state_index],y)
    for time_index = 1:number_of_timesteps
      X[time_index,state_index] = tmp[time_index]
    end
  end

  # Check for negatives -
  idx_n = find(X.<0)
  X[idx_n] = 0.0

  # return time and state -
  return (t,X);
end
