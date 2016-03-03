#2DOperators.jl

#TEMOperators.jl
#This file contains the operator definitions for a 2D TEM simulation using YEE's algorithm
#For YEE'S algorithm the structure of the FDS is as follows (x,y,z=1,t,xyzdim)
#BC is structured as follows for OH BC(i,j,k=1,t)=(epsilon, mu, magfieldLoss, sigma)
#This is built around Yee's FDTD algorithm and time stepping.
#HSTEP increments the current i,j,k,t location by

#This function computes the HField at the next time step and updates the FDS to reflect this change.
#i,j,k,t are all integer arguments, these are the space-time coordinates of the location in FDS being acted on.
#FDS is the field data structure. It has the following structure FDS[i,j,k,t,parameterxyzdim] 
#BC is the BC set in the operator calling HSTEP. HSTEP expects BC_ijk=(epsilon,mu,magloss,sigma)
function HStep(i,j,k,t,FDS,BC)
	sizeFDS=size(FDS,4)#This is the number of time steps that FDS contains (Yee's algorithm only requires 2 previous steps t-1/2 and t-1)
	modT=mod(t,sizeFDS)+1#The current cyclical timestep to access in FDS (for yee will be either 1 or 2)
	modTn=mod(t-1,sizeFDS)+1#This gives the index of the t-1/2 data set (For Yee this will be the electric field values)

	Hzn=FDS(i,j,k,modT,3)#The n-1/2 Hz value (is currently stored in FDS in the modT location)

	#The following are the E values at the n timestep (IE the current timestep -1/2)
	#The naming convention works as follows: Exy is the Ex value at one spacial step forward in the y direction. 
	#(recall indices are 1/2 steps)
	#Exny is the same as Exy it is just one index down in y
	Exy=FDS(i,j+2,k,modTn,1)
	Exny=FDS(i,j-2,k,modTn,1)
	Eyx=FDS(i+2,j,k,modTn,2)
	Eynx=FDS(i-2,j,k,modTn,2)

	#BC=(epp, mu, magloss, sigma)
	mloss=BC(i,j,k,t)[3]
	mu=BC(i,j,k,t)[2]

	#These are the coefficients that are found in Yee's standard algorithm
	#See Taflove pg. 72 for the algorithm and coefficients.
	A=((1-(mloss*dt/(2*mu)))/((1+mloss*dt)/(2*mu)))
	B=((dt/mu)/(1+(mloss*dt)/(2*mu)))

	#The current n+1/2 Hz value. In TE only Hz changes.
	Hz=A*Hzn+B*(((Exy-Exny)/dy)-((Eyx-Eynx)/dx))

	FDS(i,j,k,modT,3)=Hz#inplace overwrite of the old Hz value stored in FDS
end
#This function computes the n+1 Efield values and updates them inplace in the FDS.
#i,j,k,t are the spacetime coordinates for the FDS index being updated. 
#FDS is the field data structure. It has the following structure FDS[i,j,k,t,parameterxyzdim] 
#BC is the BC set in the operator calling HSTEP. HSTEP expects BC_ijk=(epsilon,mu,magloss,sigma)
function EStep(i,j,k,t,FDS,BC)
	sizeFDS=size(FDS,4)#This is the number of time steps that FDS contains (Yee's algorithm only requires 2 previous steps t-1/2 and t-1)
	modT=mod(t,sizeFDS)+1#The current cyclical timestep to access in FDS (for yee will be either 1 or 2)
	modTn=mod(t-1,sizeFDS)+1#This gives the index of the t-1/2 data set (For Yee this will be the electric field values)

	Exn=FDS(i,j,k,modT,1)#The n-1 Ex and Ey values currently stored in FDS
	Eyn=FDS(i,j,k,modT,2)
	
	#The following are the H values at the n timestep (IE the current timestep -1/2)
	#The naming convention works as follows: Hxy is the Hx value at one spacial step forward in the y direction. 
	#(recall indices are 1/2 steps)
	#Hxny is the same as Hxy it is just one index down in y instead of one index up.
	Hzx=FDS(i+2,j,k,modTn,3)
	Hznx=FDS(i-2,j,k,modTn,3)
	Hzy=FDS(i,j+2,k,modTn,3)
	Hzny=FDS(i,j-2,k,modTn,3)

	#BC=(epp, mu, magloss, sigma)
	sigma=BC(i,j,k,t)[4]#Just creating some named variables to avoid errors
	ep=BC(i,j,k,t)[1]

	#These are the coefficients that are found in Yee's standard algorithm
	#See Taflove pg. 72 for the algorithm and coefficients.
	A=((1-(sigma*dt/(2*ep)))/((1+sigma*dt)/(2*ep)))
	B=((dt/ep)/(1+(sigma*dt)/(2*ep)))

 	#The updated n+1 timestep Ex and Ey values at i,j,k,t  
	Ex=A*Exn+B*((Hzy-Hzny)/dy)
	Ey=A*Eyn+B*(-(Hzx-Hznx)/dx)

	#This updates the i,j,k,t E field values in FDS.
	FDS(i,j,k,modT,1)=Ex
	FDS(i,j,k,modT,2)=Ey
end

#This function creates and returns a composite operator that is used to compute TE wave propagation
#It requires OppBC an indexed function or object which represents the functional Boundary conditions. 
#E.G. When is a particular operator applied as compared to another.
#While MatBC is also an indexed function which takes in i,j,k,t and returns the material properties at that index
#It returns a tuple of the form (epp, mu, magloss, sigma) 
#NOTICE, this initializer just creates a master operator, but it does not contain the EMWall 
function initTEOperator(OppBC,MatBC)
	OH=operator(HStep,MatBC)#This is the H Field calculating operator
	OE=operator(EStep,MatBC)#This is the E Field calculating operator, they apply the Steps given above.
	I=operator(I, 0)
	return operator([I, OE, OH], OppBC)
end
