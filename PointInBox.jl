#3DPointInBox.jl
#This is a simulation of a hard source radiating in a box that is surrounded by a PEC
#The hard source is a point. It simply varies the electric field according to a function in time.
#The dimensions of the box are 10cm^3

dx=.01
dy=.01
dz=.01
dt=.001

using OppFDTD#This allows for access to the simulator
using 2DTE#This brings in the 2DTE operators

function createYee2DOperator(OppBC, MatBC)
	OH=operator(2DTE.HStep, MatBC)
	OE=operator(2DTE.EStep, MatBC)
	I=operator(2DTE.Identity)
	Wall=operator(2DTE, PECWallBC)
	return operator((I,OE, OH, Wall))
#OppBC for Yee's algorithm is as follows. Even steps are integer values, Odd Steps are half steps. 
#This applies for both space and time domains.
#By looking at
function OppBC(i,j,k,t)
	ret=0
	ret=generateYeeGrid(i,j,k,t,2,3,1)
	if isWall(i,j,k,t)==True
		ret=4
	end
end

function isSource(i,j,k,t)
	

function generateYeeGrid(i,j,k,t,OEVal, OHVal, Identity)
	istep=mod(i,2)&&size(FDS,1)>1
	jstep=mod(j,2)&&size(FDS,2)>1
	kstep=mod(j,2)&&size(FDS,3)>1
	ret = Identity
	if istep==1&&jstep==1&&kstep==1&&mod(t+1,2)#If i,j,k are at half steps and t even
		ret=OEVal#Compute the HStep on these points
	elseif istep==0&&jstep==0&&kstep==0&&mod(t,2)#If i,j,k are at integer steps and t is odd 
		ret=OHVal#Compute the EStep on these points (i,j,k are )
	end
 	return ret
end

function isWall(i,j,k,t)
	iwall=(i==1||i==size(FDS,1))&&size(FDS,1)!=1#These check to see if we are on a wall and in a dimension that is being used
	jwall=(j==1||j==size(FDS,2))&&size(FDS,2)!=1#If we did not check dimension size then we would put wall over entire domain.
	kwall=(k==1||k==size(FDS,3))&&size(FDS,3)!=1#Since over the unused dimension every grid point is a wall point
	if (iwall||jwall||kwall)
		return 1
	else
		return 0
	end 
end

function setSources(i,j,k,t)


function PECWallBC(i,j,k,t)
	if i==1||i==size(FDS,1)
		return 1
	elseif j==1||j==size(FDS,2)
		return 2
	elseif k==1||j==size(FDS,3)
		return 3
	end	

function MatBC(i,j,k,t)
	return (epp, mu, magLoss, sigma, dx, dy, dz, dt)

function initNullFDS()
	return "Not done yet"#Totally empty FDS

simulator = initTEOperator(OppBC, MatBC)

FDS = initNullFDS()

simTimeSteps = 100000

simulate(simTimeSteps, 1000, FDS, simulator)

#simulator.Transforms appendSource#Pseudocode


