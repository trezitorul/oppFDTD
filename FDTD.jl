#FDTD.jl
#This is the file that will contain the main FDTD loop and initialization.
#Additionally this file will contain the Field Data Structure (FDS) array.
#FDS is structured as follows. (X Coord, Y Coord, Z Coord, t step*, FieldQuantityIndex**, FieldIndex)
#*tstep, the tstep is a cycling coordinate, given the current time increment, everything above the current index is in the past
#This way only part of FDS is rewritten everytime the code runs. Instead of rewriting the entire matrix.
#**Field Quantity Index: E field component would be represented by 0 and H by 1 and Dipole orientation by 2...etc
#Field is the dimension of the FieldQuantity being referenced. IE 0 might be Ex
#FDS will also contain a time component which allows it to store historical data as 
#Required by the algorithms implemented by the operators system.
#include("2DOperatorsTM.jl")


module OppFDTD 

export simulate


include("operators.jl")#Contains the definition and some associated functions of operators.
require("operators.jl")
include("plottingUtilities.jl")#Contains plotting utilities
#include("fileIO.jl")#Not Implemented Yet
#Contains the functions for easily creating and performing 2D Yee algorithm simulations
include("2DYeeUtils.jl")
#Simulate is the core of the FDTD engine, this is where the FDS structure (which holds all of the E and H fields) is updated for the next timestep.
#steps is the number of timesteps to execute before terminating operation.
#sampleRate is how often FDS should be written to disk for later animation
#FDS is the Field Data Structure, it holds all of the information about the field
#simulator, is the operator that updates the FDS for each timestep
#OutputOperation, this is a function that can sample FDS during a run, it can be a fileIO function, plotting
function simulate(steps::Int64, sampleRate::Int64, FDS::Array{Float64}, simulator, outputOperation)
	for t=1:steps#The integer time step. The way the algorithm increments including the time step size is set in the simulator operator
		if mod(t,sampleRate)==0#Checks to see if it is time to save an FDS snapshot to disk
			#plotEField(FDS)#Saves FDS to disk
			outputOperation(FDS)
		end
#		println(FDS)
		for i=1:size(FDS,1)#The x spacial coordinate
			for j=1:size(FDS,2)#The y spacial coordinate
				for k=1:size(FDS,3)#The z spacial coordinate
					#println("simulate")
 					#println((i,j,k,t))
					applyOperator!(i,j,k,t,simulator,FDS)#Applys the simulator operator to FDS, the operator type definition is found in operators.jl
					#The simulator operator is found in a separate file that is included in a particular simulation run.

				end
			end
		end
	end		
end

end