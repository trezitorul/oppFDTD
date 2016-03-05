#plottingUtilities.jl
using PyPlot
PyPlot.ion()
clf()
#plt=plot(rand(10,10))
function getZMat(FDS,Param,dim)
 	Z=[sqrt(FDS[i,j,1,2,1]^2+FDS[i,j,1,2,2]^2)  for i=1:2:size(FDS,1),j=2:2:size(FDS,2)]
 	println(size(Z))
end
function plotEField(FDS)
	X=[x for x=1:2:size(FDS,1), y=1:2:size(FDS,2)]
	Y=[y for x=1:2:size(FDS,1), y=1:2:size(FDS,2)]
	getZMat(FDS,1,1)

	Z=[FDS[i,j,1,2,1]^2+FDS[i,j,1,2,2]^2  for i=1:2:size(FDS,1),j=2:2:size(FDS,2)]
	Z2=[FDS[i,j,1,2,1]^2+FDS[i,j,1,2,2]^2  for i=2:2:size(FDS,1),j=1:2:size(FDS,2)]
	#Z=[Z1[i,j+1]+Z2[i+1,j] for i=1:size(Z1,1)-1, j=1:size(Z1,2)-1]
#	println(Z)
#	println(size(Z))
#	print(size(Z))S
#	print(size(Y))
#	clf()
#	println(FDS)
#	size(Z)
		
#	pcolor(Z)
	imshow(Z, cmap="bwr",interpolation="gaussian", vmin=-.1, vmax=.1)
#	pcolor(Z,cmap=ColorMap("RdPu"))
	PyPlot.draw()
end
function plotHField(FDS)
	X=[x for x=1:2:size(FDS,1), y=1:2:size(FDS,2)]
	Y=[y for x=1:2:size(FDS,1), y=1:2:size(FDS,2)]
	getZMat(FDS,1,1)
	Z=[FDS[i,j,1,1,1]^2+FDS[i,j,1,1,2]^2 for i=1:size(FDS,1),j=1:size(FDS,2)]
#	Z=[FDS[i,j,1,1,1]^2+FDS[i,j,1,1,2]^2  for i=1:size(FDS,1),j=1:size(FDS,2)]
#	Z2=[FDS[i,j,1,2,1]^2+FDS[i,j,1,2,2]^2  for i=2:2:size(FDS,1),j=1:2:size(FDS,2)]
	#Z=[Z1[i,j+1]+Z2[i+1,j] for i=1:size(Z1,1)-1, j=1:size(Z1,2)-1]
#	println(Z)
#	println(size(Z))
#	print(size(Z))S
#	print(size(Y))
#	clf()
#	println(FDS)
#	size(Z)
		
#	pcolor(Z)
	imshow(Z, cmap="bwr",interpolation="gaussian", vmin=-.000001, vmax=.000001)
#	pcolor(Z,cmap=ColorMap("RdPu"))
	PyPlot.draw()
end