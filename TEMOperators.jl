#TEMOperators.jl
#This file contains the operator definitions for a 2D TEM simulation using YEE's algorithm
#For YEE'S algorithm the structure of the FDS is as follows (x,y,z,t,{E=1 H=2}, xyzdim)
#BC is structured as follows for OH BC(i,j,k,t)=(epsilon, mu, magfieldLoss, sigma)
dx=.01
dy=.01
dz=.01
dt=.001

#Gets the parameter
function getParam(i,j,k,t,)

function OH(i,j,k,t,FDS,BC)=
	modT=mod(t,size(FDS,4))#Finds the current index we should be on.
	sizeFDS=size(FDS,4)

	Hxn=FDS(i,j,k,mod(t-2,sizeFDS),2,1)
	Hyn=FDS(i,j,k,mod(t-2,sizeFDS),2,2)
	Hzn=FDS(i,j,k,mod(t-2,sizeFDS),2,3)

	Eyz=FDS(i,j,k+1,mod(t-1,sizeFDS),1,2)
	Eynz=FDS(i,j,k-1, mod(t-1,sizeFDS),1,2)
	Ezy=FDS(i,j+1,k, mod(t-1,sizeFDS),1,3)
	Ezny=FDS(i,j-1,k, mod(t-1,sizeFDS),1,3)
	Ezx=FDS(i+1,j,k, mod(t-1,sizeFDS),1,3)
	Eznx=FDS(i-1,j-1,k, mod(t-1,sizeFDS),1,3)
	Eyz=FDS(i,j,k+1,mod(t-1,sizeFDS),1,1)
	Eynz=FDS(i,j,k-1, mod(t-1,sizeFDS),1,1)
	Exz=FDS(i,j,k+1,mod(t-1,sizeFDS),1,1)
	Exnz=FDS(i,j,k-1, mod(t-1,sizeFDS),1,1)
	Exy=FDS(i,j+1,k ,mod(t-1,sizeFDS),1,1)
	Exny=FDS(i,j-1,k, mod(t-1,sizeFDS),1,1)
	Eyx=FDS(i+1,j,k, mod(t-1,sizeFDS),1,2)
	Eynx=FDS(i-1,j-1,k, mod(t-1,sizeFDS),1,2)

	mloss=BC(i,j,k,t)[3]
	mu=BC(i,j,k,t)[2]
	A=((1-(mloss*dt/(2*mu)))/((1+mloss*dt)/(2*mu)))
	B=((dt/mu)/(1+(mloss*dt)/(2*mu)))

	Hx=A*Hxn+B*(((Eyz-Eynz)/dz)-((Ezy-Ezny)/dy))
	Hy=A*Hyn+B*(((Ezx-Eznx)/dx)-((Exz-Exnz)/dz))
	Hz=A*Hzn + B*(((Exy-Ezny)/dy)-((Eyx-Eynx)/dz))
 

	FDS(i,j,k,modT,2,1)=Hx
	FDS(i,j,k,modT,2,2)=Hy
	FDS(i,j,k,modT,2,3)=Hz

function OE(i,j,k,t,FDS,BC)=
	modT=mod(t,size(FDS,4))#Finds the current index we should be on.
	sizeFDS=size(FDS,4)

	Exn=FDS(i,j,k,mod(t-2,sizeFDS),1,1)
	Eyn=FDS(i,j,k,mod(t-2,sizeFDS),1,2)
	Ezn=FDS(i,j,k,mod(t-2,sizeFDS),1,3)

	Hyz=FDS(i,j,k+1,mod(t-1,sizeFDS),2,2)
	Hynz=FDS(i,j,k-1, mod(t-1,sizeFDS),2,2)
	Hzy=FDS(i,j+1,k, mod(t-1,sizeFDS),2,3)
	Hzny=FDS(i,j-1,k, mod(t-1,sizeFDS),2,3)
	Hzx=FDS(i+1,j,k, mod(t-1,sizeFDS),2,3)
	Hznx=FDS(i-1,j-1,k, mod(t-1,sizeFDS),2,3)
	Hyz=FDS(i,j,k+1,mod(t-1,sizeFDS),2,1)
	Hynz=FDS(i,j,k-1, mod(t-1,sizeFDS),2,1)
	Hxz=FDS(i,j,k+1,mod(t-1,sizeFDS),2,1)
	Hxnz=FDS(i,j,k-1, mod(t-1,sizeFDS),2,1)
	Hxy=FDS(i,j+1,k ,mod(t-1,sizeFDS),2,1)
	Hxny=FDS(i,j-1,k, mod(t-1,sizeFDS),2,1)
	Hyx=FDS(i+1,j,k, mod(t-1,sizeFDS),2,2)
	Hynx=FDS(i-1,j-1,k, mod(t-1,sizeFDS),2,2)

	sigma=BC(i,j,k,t)[4]
	ep=BC(i,j,k,t)[1]
 
	A=((1-(sigma*dt/(2*ep)))/((1+sigma*dt)/(2*ep)))
	B=((dt/ep)/(1+(sigma*dt)/(2*ep)))

	Ex=A*Exn+B*(((Hzy-Ezny)/dy)-((Eyz-Eynz)/dz))
	Ey=A*Eyn+B*((-(Hzx-Hznx)/dx)+((Hxz-Hxnz)/dz))
	Ez=A*Ezn + B*((-(Hxy-Hzny)/dy)+((Hyx-Hynx)/dz))

	FDS(i,j,k,modT,1,1)=Ex
	FDS(i,j,k,modT,1,2)=Ey
	FDS(i,j,k,modT,1,3)=Ez

function masterOperator(i,j,k,FDS,BC)
	