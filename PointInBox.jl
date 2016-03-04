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


include("2DOperatorsTE.jl")
using FlatTE#This allows for access to the simulator

function sineSource(i,j,k,t,FDS,BC)#BC = (I,w, normal)S
	modT=2
	FDS[i,j,k,modT,2]=BC[1]*sin(BC[2]*t)
	FDS[i,j,k,modT,1]=BC[1]*sin(BC[2]*t)
end

function createYee2DOperator(OppBC, MatBC)
	OH=FlatTE.operator((FlatTE.HStep), MatBC)
	OE=FlatTE.operator((FlatTE.EStep), MatBC)
	I=FlatTE.operator((FlatTE.Identity),0)
	Wall=FlatTE.operator((FlatTE.EMWall), PECWallBC)
	sineSrc=FlatTE.operator((sineSource), (1,srcFreq,3))#BC for this source 
	return FlatTE.operator((I,OE, OH, Wall, sineSrc),OppBC)
end
#OppBC for Yee's algorithm is as follows. Even steps are integer values, Odd Steps are half steps. 
#This applies for both space and time domains.
#By looking at
function OppBC(i,j,k,t)
	ret=1
	ret=generateYeeGrid(i,j,k,t,2,3,1)
	if isWall(i,j,k,t)==true
		ret=4
	end
	if i==1+floor((size(FDS,1)+1)/2)&&j==floor((size(FDS,2)+1)/2)&&mod(t,2)==0
	 	ret = 5
	 end
	 return ret 
end
	
function generateYeeGrid(i,j,k,t,OEVal, OHVal, Identity)
	istep=mod(i,2)==1
	jstep=mod(j,2)==1
	ret = Identity
	if (istep!=jstep)&&mod(t,2)==0#If i,j,k are at half steps and t even
#	if (istep||jstep)&&mod(t+1,2)==1#If i,j,k are at half steps and t even
		ret=OEVal#Compute the HStep on these points
	elseif (istep==jstep)&&mod(t,2)==1#If i,j,k are at integer steps and t is odd 
#	elseif (istep==false||jstep==false)&&mod(t,2)==1#If i,j,k are at integer steps and t is odd 
		ret=OHVal#Compute the EStep on these points (i,j,k are )
	end
 	return ret
end

function visualizeOppBC()
	println([OppBC(i,j,1,1) for i=1:xSteps, j=1:ySteps, k=1:zSteps])
	println([OppBC(i,j,1,2) for i=1:xSteps, j=1:ySteps, k=1:zSteps])
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


function PECWallBC(i,j,k,t)
	if i==1||i==size(FDS,1)
		return 1
	elseif j==1||j==size(FDS,2)
		return 2
	elseif k==1||j==size(FDS,3)
		return 3
	end	
end

function MatBC(i,j,k,t)
	epp=8.854187817*10.0^(-12)
	mu=4*pi*10.0^(-7)
#	epp=1
#	mu=1
	magLoss=0
	sigma=0
	return (epp, mu, magLoss, sigma, dx, dy, dz, dt)
end

function initNullFDS()
	return "Not done yet"#Totally empty FDS
end 

simulator = createYee2DOperator(OppBC, MatBC)

#visualizeOppBC()

FDS=zeros(Float64,xSteps, ySteps, zSteps, tDim,3)
#FDS[500,501,1,2,1]=10000000
#FDS[25,26,1,2,2]=1
#FDS[10,10,1,2,2]=10000000
simTimeSteps = 100000
#println(FDS)
FlatTE.simulate(tSteps, sampleRate, FDS, simulator)



