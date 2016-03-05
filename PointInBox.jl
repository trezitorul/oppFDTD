#3DPointInBox.jl
#This is a simulation of a hard source radiating in a box that is surrounded by a PEC
#The hard source is a point. It simply varies the electric field according to a function in time.
#The dimensions of the box are 10cm^3
srcFreq=10000
wavelength=2*pi/srcFreq
dx=wavelength/10#m
dy=wavelength/10
dz=wavelength/10
dt=1/((3*10^8)*sqrt(1/(dx^2)+1/(dy^2))) #magic timestep
xSteps=500
ySteps=500 
zSteps=1
tSteps=1000
tDim = 2
sampleRate=10
sourceBC=(100,srcFreq,(1,2))

include("FDTD.jl")
using OppFDTD#This allows for access to the simulator

function sineSource(i,j,k,t,FDS,BC)#BC = (I,w, normal)S#Replace with imported hardSource file
	modT=mod(t,2)+1
	for i=1:length(BC[3])
		FDS[i,j,k,modT,BC[3][i]]=BC[1]*sin(BC[2]*t)
	end
end

function sourceLocation(i,j,k,t,srcNum)
	if i==1+floor((size(FDS,1)+1)/2)&&j==floor((size(FDS,2)+1)/2)&&mod(t,2)==0
	 	return srcNum
	end
end

source=createNewHardPointSource(sineSource,sourceBC)

#Sets the Material Boundary Conditions
function MatBC(i,j,k,t)
	epp=8.854187817*10.0^(-12)
	mu=4*pi*10.0^(-7)
#	epp=1
#	mu=1
	magLoss=0
	sigma=0
	return (epp, mu, magLoss, sigma, dx, dy, dz, dt)
end




function OppBC(i,j,k,t)
	return compositeBC(i,j,k,t,[EmptyYeeGrid(i,j,k,t,2,3,1), PECWallBC(i,j,k,t),sourceLocation(i,j,k,t,5)])
end

masterOperator=createYee2DTEOperator(OppBC, MatBC)
addOperator(masterOperator, source)

#visualizeOppBC()
FDS=zeros(Float64,xSteps, ySteps, zSteps, tDim,3)

simulate(tSteps, sampleRate, FDS, simulator, plotEField)



