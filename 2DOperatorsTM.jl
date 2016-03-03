#2DOperators.jl

#TEMOperators.jl
#This file contains the operator definitions for a 2D TEM simulation using YEE's algorithm
#For YEE'S algorithm the structure of the FDS is as follows (x,y,z=1,t,xyzdim)
#BC is structured as follows for OH BC(i,j,k=1,t)=(epsilon, mu, magfieldLoss, sigma)
#This is built around Yee's FDTD algorithm and time stepping.
#HSTEP increments the current i,j,k,t location by

#Decided to add a comment
function HStep(i,j,k,t,FDS,BC)
	sizeFDS=size(FDS,4)
	modT=mod(t,sizeFDS)+1#Finds the current index we should be on.
	modTn=mod(t-1,sizeFDS)+1#Finds the 1/2 index in the past

	Hxn=FDS(i,j,k,modT,1)
	Hyn=FDS(i,j,k,modT,2)
	Hzn=FDS(i,j,k,modT,3)

	Exy=FDS(i,j+1,k,modTn,1)
	Exny=FDS(i,j-1,k,modTn,1)
	Exz=FDS(i,j,k,modTn,1)
	Exnz=FDS(i,j,k,modTn),1)
	Eyx=FDS(i+1,j,k,modTn,2)
	Eynx=FDS(i-1,j,k,modTn,2)
	Eyz=FDS(i,j,k,modTn,2)
	Eynz=FDS(i,j,k,modTn,2)
	Ezx=FDS(i+1,j,k,modTn,3)
	Eznx=FDS(i-1,j,k,modTn,3)
	Ezy=FDS(i,j+1,k,modTn,3)
	Ezny=FDS(i,j-1,k,modTn,3)

	mloss=BC(i,j,k,t)[3]
	mu=BC(i,j,k,t)[2]
	A=((1-(mloss*dt/(2*mu)))/((1+mloss*dt)/(2*mu)))
	B=((dt/mu)/(1+(mloss*dt)/(2*mu)))

# Changes here are due to 2D . No variation WRT z direction.
	Hx=A*Hxn+B*(((Eyz-Eynz)/dz)*0-((Ezy-Ezny)/dy))
	Hy=A*Hyn+B*(((Ezx-Eznx)/dx)-0*((Exz-Exnz)/dz))
	Hz=A*Hzn+B*(((Exy-Exny)/dy)-((Eyx-Eynx)/dx))

# Changes here are due to TM mode 
	FDS(i,j,k,modT,1)=Hx
	FDS(i,j,k,modT,2)=Hy
	FDS(i,j,k,modT,3)=Hz*0
end

function EStep(i,j,k,t,FDS,BC)
	sizeFDS=size(FDS,4)
	modT=mod(t,sizeFDS)+1#Finds the current index we should be on.
	moTn=mod(t-1,sizeFDS)+1#Finds the opposite time step which is 1/2 timestep in the past.

	Exn=FDS(i,j,k,modT,1)#The current time
	Eyn=FDS(i,j,k,modT,2)
	Ezn=FDS(i,j,k,modT,3)

	Hxy=FDS(i,j+1,k,modTn,1)
	Hxny=FDS(i,j-1,k,modTn,1)
	Hxz=FDS(i,j,k,modTn,1)
	Hxnz=FDS(i,j,k,modTn,1)
	Hyx=FDS(i+1,j,k,modTn,2)
	Hynx=FDS(i-1,j,k,modTn,2)
	Hyz=FDS(i,j,k,modTn,2)
	Hynz=FDS(i,j,k,modTn,2)
	Hzx=FDS(i+1,j,k,modTn,3)
	Hznx=FDS(i-1,j,k,modTn,3)
	Hzy=FDS(i,j+1,k,modTn,3)
	Hzny=FDS(i,j-1,k,modTn,3)

	sigma=BC(i,j,k,t)[4]
	ep=BC(i,j,k,t)[1]

	A=((1-(sigma*dt/(2*ep)))/((1+sigma*dt)/(2*ep)))
	B=((dt/ep)/(1+(sigma*dt)/(2*ep)))

# Changes here are due to 2D . No variation WRT z direction.
	Ex=A*Exn+B*(((Hzy-Ezny)/dy)-0*((Eyz-Eynz)/dz))
	Ey=A*Eyn+B*(((Hxz-Hxnz)/dz)*0-((Hzx-Hznx)/dx))
	Ez=A*Ezn+B*(((Hyx-Hynx)/dx)-((Hxy-Hxny)/dy))


# Changes here are due to TM mode
	FDS(i,j,k,modT,1)=Ex*0
	FDS(i,j,k,modT,2)=Ey*0
	FDS(i,j,k,modT,3)=Ez
end

#EM Wall can either be a PMC or PEC depending on the boundary condition.
#If the wall is executed on odd time steps (1/2 time steps) and on even spacial steps then it will act as PMC
#If the wall is executed on even time steps (integer time steps) and on odd spacial steps then it will be a PEC
function EMWall(i,j,k,t,FDS,BC)
	modT=mod(t,size(FDS,4))+1
		if i==1 || i==size(FDS,1)
			FDS(i,j,k,modT,2)=0
			FDS(i,j,k,modT,3)=0
		elseif j==1 || j==size(FDS,2)
			FDS(i,j,k,modT,1)=0
			FDS(i,j,k,modT,3)=0
		elseif k==1 || k==size(FDS,3)
			FDS(i,j,k,modT,1)=0
			FDS(i,j,k,modT,2)=0
end

#This operation does nothing. It is used if you want to avoid making any changes on the FDS
function Identity(i,j,k,t, FDS, BC)
	return 0
end

function initTEOperator(OppBC,MatBC)
	OH=operator(HStep,MatBC)
	OE=operator(EStep,MatBC)
	PEC=operator(PECBoundary,0)
	return operator((OE, OH, EMWall, Identity), OppBC)
end
