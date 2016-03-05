#2DYeeUtils.jl
#This file contains many functions that are common to 2DYee algorthim runs.
#Additionally, this file contains the TE and TM modules that can be imported for execution.
include("2DOperatorsTE.jl")
require("2DOperatorsTE.jl")
import FlatTE
include("2DOperatorsTM.jl")
require("2DOperatorsTM.jl")
import FlatTM

export createNewHardPointSource, EmptyYeeGridBC, YeeGrid, EMWallBC, PECOppBC, PMCOppBC,createYee2DTEOperator, createYee2DTMOperator
export updateWallBC, addOperator, compositeBC
function YeeGrid(i,j,k,t,OEVal, OHVal, Identity)
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

function createNewHardPointSource(srcFxn,BC)
	return operator([srcFxn],BC)
end
#This function generates an empty box, the Transforms it uses are stored in an array. This base object,
#can then be easily extended with more operators such as walls or sources in the simulation file
function EmptyYeeGridBC(i,j,k,t)
	ret=1
	ret=YeeGrid(i,j,k,t,2,3,1)
	return ret 
end

function EMWallBC(i,j,k,t)
	if (i==1||i==size(FDS,1))&&size(FDS,1)!=1
		return 1
	elseif (j==1||j==size(FDS,2))&&size(FDS,2)!=1
		return 2
	elseif k==1||k==size(FDS,3)&&size(FDS,3)!=1
		return 3
	end	
end

#Creates the OperatorBC to build a PEC wall around the entire domain.
function PECOppBC(i,j,k,t)
	if mod(t,2)==0&&isWall(i,j,k,t)
		return 4
	end
end

#Creates the OperatorBC to build a PMC wall around the entire domain
function PMCOppBC(i,j,k,t)
	if mod(t,2)==1&&isWall(i,j,k,t)
		return 4
	end
end
#Creates the bare minimum YEE 2D TEOperator which is in an empty box with walls the kill both E and H fields
#The emptyBox Opp is simply the yeeGridon all points. This opp, nees to be updated with 
function createYee2DTEOperator(OppBC,MatBC)
	OH=operator((FlatTE.HStep), MatBC)
	OE=operator((FlatTE.EStep), MatBC)
	I=operator((Identity),0)
	Wall=operator((EMWall), EMWallBC)
	return operator((I,OE,OH,Wall),OppBC)
end 

#This returns the operator with prebuilt BC for an empty box Surrounded by a 
function createYee2DTMOperator(OppBC, MatBC)
	OH=operator([FlatTM.HStep], MatBC)
	OE=operator([FlatTM.EStep], MatBC)
	I=operator([Identity],0)
	Wall=operator([EMWall], NormWallBC)
	return operator([I,OE,OH,Wall],OppBC)
end

#Upgrades the masterOperators wall operator with the given BC, useful if you want to change wall types, on a whim.
function updateWallBC(masterOperator::operator, BC)
	masterOperator.Transforms[4].BC=BC
	return masterOperator
end

#This adds an operator to the master operator transform list. In addition to this new operator addition requires the BC 
#To be upgraded.
function addOperator(masterOperator::operator, operator,updatedBC)
	push!(masterOperator.Transforms,operator)
	master.Operator.BC=updatedBC
end

function addOperator(masterOperator::operator, operator)
	push!(masterOperator.Transforms,operator)
end
#This function is used to create Composite Operators. For instance, one can start with an empty box and then
#by executing this composite add in source terms along with other terms.
function compositeBC(i,j,k,t,listOfBC)
	ret = 1
	for i=1:length(listOfBC)
		ret = listOfBC(i,j,k,t)
	end
	return ret
end

function visualizeBC()
	println("The HField Operators BC")
	println([OppBC(i,j,1,1) for i=1:xSteps, j=1:ySteps, k=1:zSteps])
	println("The EField Operators BC")
	println([OppBC(i,j,1,2) for i=1:xSteps, j=1:ySteps, k=1:zSteps])
end