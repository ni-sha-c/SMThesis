module PlotStuff

export plotvel
using LsqFit
using PyPlot
using toolbox
function plotvel(f)


	f1 = open(f)
	v = Float64[]
	z = Float64[]
	dens = Float64[]
	count = 0

	for l in enumerate(eachline(f1))

		data = readdlm(IOBuffer(l[2]))
		push!(v, convert(Float64, 100*data[3]))
		push!(z, convert(Float64, data[1]))
		push!(dens, convert(Float64, data[2]))
		if(data[2] > 10.)
				count = count + 1
		end

	end
	figure()
	vplot = Float64[]
	zplot = Float64[]
	for i=1:length(dens)
		if(dens[i]>10.0)
				push!(vplot, v[i])
				push!(zplot, z[i])
		end
	end

	plot(zplot-zplot[1], vplot, "o-")
	hold(true)
#Curve fitting.

	model(x,p) = p[3] + p[1]*(x-p[2]).^2.0

	fit = curve_fit(model, zplot[2:end-1], vplot[2:end-1], [1.0, 1.0, 1.0,1.0])

	dudy1 = 2.0*fit.param[1]*(zplot[1]-fit.param[2])
	dudy2 = 2.0*fit.param[1]*(zplot[end]-fit.param[2])
	vfit = model(zplot,fit.param)
	

	Vs = mean([vfit[end], vfit[1]])
	Ls1 = vfit[1]/dudy1 
	Ls2 = vfit[end]/dudy2
	plot(zplot-zplot[1],vfit, "o-")
	
	
	
	
	
	return zplot, vplot, vfit, fit.param, Ls1, Ls2





end
function plottemp(f)
	f1 = open(f)
	T = Float64[]
	for l in enumerate(eachline(f1))
		if(l[1]>70)
			data = readdlm(IOBuffer(l[2]))
			if(typeof(data[1])!=SubString{ASCIIString})
				push!(T, convert(Float64, data[1]))
			end
		end
	end
	#figure()
	Tmax = maximum(T)
	Tmin = minimum(T)
	Tavg = mean(T)

	plot(smooth(T[1:50:end]))

	return Tmax,Tmin,Tavg
end
end
