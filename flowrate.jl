module Flowrate

export flowrate
using LsqFit
using PyPlot
function flowrate(z,v, h, g, n)
#Same convention and units as the function, velocity.
	z = (z-z[1])
	g    = g*0.06649 #units: angstrom/ps^2
	rho  = n/(99.9363*71.38*h)*14.512e0*1660.468 #units: kg/m^3

	uavg = mean(v)
	q    = uavg*h #units : angstrom^2/ps

	chn  = q/g/rho*1.e5 #units: m (nm)^3 s/kg 


	return chn
end
function flowrate_bulk(z,v, h, g, n)
#Same convention and units as the function, velocity.
	z = (z-z[1])
	g    = g*0.06649 #units: angstrom/ps^2
	rho  = n/(99.9363*71.38*h)*14.512e0*1660.468 #units: kg/m^3
	v = v[5:end-4]
	uavg = mean(v)
	h = length(v)
	q    = uavg*h #units : angstrom^2/ps

	chn  = q/g/rho*1.e5 #units: m (nm)^3 s/kg 


	return chn
end
function velocity(z, v, h, g, n)
# g = g'/1.e-4 eV/Angstrom.
# v = v' Angstrom/ps
# z = z' Angstrom (actual coordinate)
# h = h' Angstrom (measured film thickness)
# n = number of atoms

#Nondimensionalize
#	h0   = 50.e0 #units: angstrom
#	U0   = 50.e0 #units: m/s
#	t0   = 100.e0  #units: ps
#	rho0 = 1000.e0 #units: kg/m^3

	
	z = (z-z[1])
	g    = g*0.06649 #units: angstrom/ps^2
	rho  = n/(99.5e0*72.3e0*h)*14.512e0*1660.468 #units: kg/m^3

	v = v/rho/g*1.e8 # Angstrom^2 m s/kg

	model(x,p) = p[1]*x.^2.0 + p[2]*x + p[3]

	fit = curve_fit(model, z[2:end-1], v[2:end-1], [1.0, 1.0, 1.0])
	vfit = model(z,fit.param)*g*rho/1.e6
	v = v*g*rho/1.e6
	res1 = fit.param[1]
	res2 = fit.param[2]
	res3 = fit.param[3]

	mu = -0.5/res1 #in Pa-s
	Ls = res3/res2 #in Angstrom
	#plot(z-z[1],vfit/maximum(vfit),"*-")
	plot(z[5:end-5]-z[5],v[5:end-5]/maximum(v),"o-")
	return mu,Ls

end
function velocity_scaled(z, v, h, g, n)
# g = g'/1.e-4 eV/Angstrom.
# v = v' Angstrom/ps
# z = z' Angstrom (actual coordinate)
# h = h' Angstrom (measured film thickness)
# n = number of atoms

#Nondimensionalize
#	h0   = 50.e0 #units: angstrom
#	U0   = 50.e0 #units: m/s
#	t0   = 100.e0  #units: ps
#	rho0 = 1000.e0 #units: kg/m^3

	
	z = (z-z[1])
	g    = g*0.06649 #units: angstrom/ps^2
	rho  = n/(99.5e0*72.3e0*h)*14.512e0*1660.468 #units: kg/m^3

	v = v/rho/g*1.e8 # Angstrom^2 m s/kg

	model(x,p) = p[1]*x.^2.0 + p[2]*x + p[3]

	fit = curve_fit(model, z[2:end-1], v[2:end-1], [1.0, 1.0, 1.0])
	vfit = model(z,fit.param)

	res1 = fit.param[1]
	res2 = fit.param[2]
	res3 = fit.param[3]

	plot(z-z[1],vfit,"--",linewidth=2.5)
	plot(z-z[1],v,"r",linewidth=2.5)

end
function velocity_plot(z, v, h, g, n, Ls, mu)
# g = g'/1.e-4 eV/Angstrom.
# v = v' Angstrom/ps
# z = z' Angstrom (actual coordinate)
# h = h' Angstrom (measured film thickness)
# n = number of atoms

#Nondimensionalize
#	h0   = 50.e0 #units: angstrom
#	U0   = 50.e0 #units: m/s
#	t0   = 100.e0  #units: ps
#	rho0 = 1000.e0 #units: kg/m^3

	
	z = z - z[1]
	h = z[end]
	g    = g*0.06649 #units: angstrom/ps^2
	rho  = n/(99.5e0*72.3e0*h)*14.512e0*1660.468 #units: kg/m^3

	v = v/rho/g*1.e8 # Angstrom^2 m s/kg

	model(x,p) = -0.5/p[1]*((x).^2.0 - h.*(x) - p[2]*h)
	z1 = linspace(0.e0, z[end],100) 
	vfit = model(z1,[mu,Ls])

	plot(z1,vfit*1.e-10,"b--",linewidth=2.5)
	plot(z,v*1.e-10,"r-",linewidth=2.5)

end

function velocity_couette(z, v, h, u)

	model(x,p) = p[1]*x + p[2]

	fit = curve_fit(model, z[2:end-1]-z[1], v[2:end-1], [1.0, 1.0])
	vfit = model(z,fit.param)

	res1 = fit.param[1]
	res2 = fit.param[2]
	
	Ls1 = 0.5*(2.0*u/res1 - h)
	Ls2 = 0.5*(-u*h/res2 - h)
	return Ls1

end


function mu_real(mustar)
	h0   = 50.e0 #units: angstrom
	U0   = 50.e0 #units: m/s
	t0   = 100.e0  #units: ps
	rho0 = 1000.e0 #units: kg/m^3

	mureal = mustar*rho0*h0^2.e0/t0*1.e-5 #units: mPa s
	return mureal
end
end
