#operator.jl
#This file contains the operator type definitions and base functions
export applyOperator!, EMWall, Identity, hardSource, operator

type operator
#each concrete subtype of operator should implement the following methods
#BC::Function, returns value of transform type to be used. 
# These are BC for which Transform to be used in space and time 
#Transforms::function applied using elements of inputted FDS
#transform(i,j,k,t, FDS): performs in place update of FDS according to BC 
#and transform 
	Transforms#This contains the operations to be done on FDS
	BC #A boundary condition, can be either a matrix or function anything indexable. It can be used in different ways
	#If Transforms contains operators then BC will contain indices of which operators should be applied at which timestep index 
	#In addition to which spacial location
	#If Transforms is simply a function that operatates directly on FDS then BC will hold data relevant to the transformation required such 
	#as epsilon and mu.

	operator(Transforms,BC)=new(Transforms,BC)#Constructor for the operator. 
end

function applyOperator!(i,j,k,t,opp::operator,FDS)#Applys the operators transforms, if Transforms contains functions then they are applied
#However, operators will frequently contain other operators with their own BC IE refractive index...etc. 
	if typeof(opp.Transforms)==Function
		applyOperator!(i,j,k,t,opp.Transforms, FDS,opp.BC)#Applies the operator in Transforms based on what the BC conditions require.
	else
		applyOperator!(i,j,k,t,opp.Transforms[opp.BC(i,j,k,t)],FDS)
	end
end

function applyOperator!(i,j,k,t,opp::Function, FDS,BC)#If we get to a base class of operator then the base class will 
#	println("----------")
#	println((i,j,k,t))
#	println(opp)

	opp(i,j,k,t,FDS,BC)#Apply the operator Transform to the FDS at index i,j,k,t
	#Once we reach a base operator (IE an operator that has Transforms = a function that updates FDS) then it simply operates on the data
	#This terminates the hierarchical recursive operator structure. 
end

#EM Wall can either be a PMC or PEC depending on the boundary condition.
#If the wall is executed on odd time steps (1/2 time steps) and on even spacial steps then it will act as PMC
#If the wall is executed on even time steps (integer time steps) and on odd spacial steps then it will be a PEC
#Both the PEC and PML essentially set the parallel components of the field equal to 0. 
#BC gives the index of the normal direction. This function then sets all of the parallel components equal to zero. 
#The location and operation of this function is dictated by the master operator BC.
function EMWall(i,j,k,t,FDS,BC)
	modT=mod(t,size(FDS,4))+1
 	for a=1:size(FDS,5)
 		if a!=BC(i,j,k,t)
			FDS[i,j,k,modT,a]=0
		end
	end
end

#This is the identity operation. It does not modify the FDS in anyway.
function Identity(i,j,k,t, FDS, BC)
	return 0
end
 