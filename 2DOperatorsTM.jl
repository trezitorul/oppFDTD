#2DOperatorsTE.jl

#TEOperators.jl
#This file contains the operator definitions for a 2D TE simulation using YEE's algorithm
#For YEE'S algorithm the structure of the FDS is as follows (x,y,z=1,t,xyzdim)
#BC is structured as follows for OH BC(i,j,k=1,t)=(epsilon, mu, magfieldLoss, sigma)
#This is built around Yee's FDTD algorithm and time stepping.
#HSTEP increments the current i,j,k,t location by

#This function computes the HField at the next time step and updates the FDS to reflect this change.
#i,j,k,t are all integer arguments, these are the space-time coordinates of the location in FDS being acted on.
#FDS is the field data structure. It has the following structure FDS[i,j,k,t,parameterxyzdim] 
#BC is the BC set in the operator calling HSTEP. HSTEP expects BC_ijk=(epsilon,mu,magloss,sigma, dx, dy, dz, dt)

module FlatTM
export HStep, EStep
function HStep(i,j,k,t,FDS,BC)
	sizeFDS=size(FDS,4)#This is the number of time steps that FDS contains (Yee's algorithm only requires 2 previous steps t-1/2 and t-1)
	#modT=mod(t-1,sizeFDS)+1#The current cyclical timestep to access in FDS (for yee will be either 1 or 2)
	#modTn=mod(t,sizeFDS)+1#This gives the index of the t-1/2 data set (For Yee this will be the electric field values)
	modT=1	
	modTn=2
#	println("HSTEP")
#	println(modT)
#	println(modTn)
#	println((i,j,k,t))
	Hxn=FDS[i,j,k,modT,1]#The n-1/2 Hz value (is currently stored in FDS in the modT location)
	Hyn=FDS[i,j,k,modT,2]
	#The following are the E values at the n timestep (IE the current timestep -1/2)
	#The naming convention works as follows: Exy is the Ex value at 1/2 spacial step forward in the y direction. 
	#(recall indices are 1/2 steps)
	#Exny is the same as Exy it is just one index down in y
	Ezy=FDS[i,j+1,k,modTn,3]
	Ezny=FDS[i,j-1,k,modTn,3]
	Ezx=FDS[i+1,j,k,modTn,3]
	Eznx=FDS[i-1,j,k,modTn,3]

	#BC=(epp, mu, magloss, sigma,dx, dy, dz, dt)
	mloss=BC(i,j,k,t)[3]
	mu=BC(i,j,k,t)[2]
	dx=BC(i,j,k,t)[5]
	dy=BC(i,j,k,t)[6]
	dz=BC(i,j,k,t)[7]
	dt=BC(i,j,k,t)[8]
	#These are the coefficients that are found in Yee's standard algorithm
	#See Taflove pg. 72 for the algorithm and coefficients.
	A=((1-((mloss*dt)/(2*mu)))/((1+((mloss*dt)/(2*mu)))))
	B=(dt/mu)/(1+((mloss*dt)/(2*mu)))
#	println("A,B")
#	println(A)
#	println(B)
#	println("---------------")
#	println(Exy)
#	println(Exny)
#	println(Eyx)
#	println(Eynx)
	#The current n+1/2 Hz value. In TE only Hz changes.
	Hx=A*Hxn+B*((Ezny-Ezy)/dy)
	Hy=A*Hyn+B*((Ezx-Eznx)/dx)
#	println(A)
#	println(B)
	FDS[i,j,k,modT,1]=Hx#inplace overwrite of the old Hz value stored in FDS
	FDS[i,j,k,modT,2]=Hy
end

#This function computes the n+1 Efield values and updates them inplace in the FDS.
#i,j,k,t are the spacetime coordinates for the FDS index being updated. 
#FDS is the field data structure. It has the following structure FDS[i,j,k,t,parameterxyzdim] 
#BC is the BC set in the operator calling HSTEP. HSTEP expects BC_ijk=(epsilon,mu,magloss,sigma)
function EStep(i,j,k,t,FDS,BC)
 
	sizeFDS=size(FDS,4)#This is the number of time steps that FDS contains (Yee's algorithm only requires 2 previous steps t-1/2 and t-1)
#	modT=mod(t-1,sizeFDS)+1#The current cyclical timestep to access in FDS (for yee will be either 1 or 2)
#	modTn=mod(t,sizeFDS)+1#This gives the index of the t-1/2 data set (For Yee this will be the electric field values)
	modT=2	
	modTn=1
#	println("ESTEP")
#	println(modT)
#	println(modTn)
#	println((i,j,k,t))
	Ezn=FDS[i,j,k,modT,3]#The n-1 Ex and Ey values currently stored in FDS
	
	#The following are the H values at the n timestep (IE the current timestep -1/2)
	#The naming convention works as follows: Hxy is the Hx value at one spacial step forward in the y direction. 
	#(recall indices are 1/2 steps)
	#Hxny is the same as Hxy it is just one index down in y instead of one index up.
	Hyx=FDS[i+1,j,k,modTn,2]
	Hynx=FDS[i-1,j,k,modTn,2]
	Hxy=FDS[i,j+1,k,modTn,1]
	Hxny=FDS[i,j-1,k,modTn,1]

	#BC=(epp, mu, magloss, sigma, dx, dy, dz, dt)
	sigma=BC(i,j,k,t)[4]#Just creating some named variables to avoid errors
	ep=BC(i,j,k,t)[1]
	dx=BC(i,j,k,t)[5]
	dy=BC(i,j,k,t)[6]
	dz=BC(i,j,k,t)[7]
	dt=BC(i,j,k,t)[8]
	#These are the coefficients that are found in Yee's standard algorithm
	#See Taflove pg. 72 for the algorithm and coefficients.
	A=((1-((sigma*dt)/(2*ep)))/(1+((sigma*dt)/(2*ep))))
	B=((dt/ep)/(1+(sigma*dt)/(2*ep)))

 	#The updated n+1 timestep Ex and Ey values at i,j,k,t  
	Ez=A*Ezn+B*(((Hyx-Hynx)/dx)+((Hxny-Hxy)/dy))
	#This updates the i,j,k,t E field values in FDS.
	FDS[i,j,k,modT,3]=Ez
 
end
end