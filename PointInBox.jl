#PointInBox.jl
#This is a simulation of a hard source radiating in a box that is surrounded by a PEC
#The hard source is a point. It simply varies the electric field according to a function in time.
#The dimensions of the box are 10cm^3 

dx=.01
dy=.01
dz=.01
dt=.001
include(TEMOperators)#Correct later!
using OppFDTD

function OppBC(i,j,k,t)
	#Implement operator BC 

function MatBC(i,j,k,t)
	return (epp, mu, magLoss, sigma)

function initNullFDS()
	return "Not done yet"#Totally empty FDS

simulator = initTEOperator(OppBC, MatBC)

FDS = initNullFDS()

simTimeSteps = 100000

simulate(simTimeSteps, 1000, FDS, simulator)

#simulator.Transforms appendSource#Pseudocode


