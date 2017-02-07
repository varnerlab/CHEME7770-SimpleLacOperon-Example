## CHEME 7770 Example: Effective model of the lac-operon logic
This repository contains the model code for an effective model of the lac operon implemented in the [Julia](http://julialang.org) programming language.

### Installation and Requirements
You can download this repository as a zip file, or `clone`/`pull` it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/CHEME7770-SimpleLacOperon-Example.git

or

	$ git clone https://github.com/varnerlab/CHEME7770-SimpleLacOperon-Example.git

The model code was machine generated using the [Gene Regulatory Network in Julia (JuGRN)](https://github.com/varnerlab/JuGRN-Generator) code generation system. The model code uses several [Julia](http://julialang.org) packages:

Package | Description | Command
--- | --- | ---
ODE | Contains the ``ode23s`` subroutine to solve the model equations | Pkg.add("ODE")
PyPlot | Used to make figures (assume you have Python installed) | Pkg.add("PyPlot")

### Solve the model equations?
The model equations can be solved by executing one of the predefined driver routines.
``Induced.jl`` solves the model equations where we add lactose to a nominal system,
while ``Repressed.jl`` solves the case where we express a lactose insensitive lacI.
The model equation, and simulation code is are contained in the [src](https://github.com/varnerlab/CHEME7770-SimpleLacOperon-Example/tree/master/src) directory.
