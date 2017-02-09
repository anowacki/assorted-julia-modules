module PREM
# Module PREM provides a number of routines for evaluating PREM
# as given by Dziewonski & Anderson, PEPI, 1981.

__precompile__()

export
    eta,
    pressure,
    Qkappa,
    Qmu,
    rho,
    vc,
    vp,
    vph,
    vpv,
    vs,
    vsh,
    vsv

const NewtonG = 6.67428e-11

# Radius of the Earth in the model
const a = 6371.0

const rad  = [1221.5, 3480.0, 3630.0, 5600.0, 5701.0, 5771.0, 5971.0, 6151.0, 
               6291.0, 6346.0, 6356.0, 6368.0, 6371.0]
const rho0 = [13.0885, 12.5815,  7.9565,  7.9565,  7.9565,  5.3197, 11.2494, 
                7.1089,  2.6910,  2.6910,  2.9000,  2.6000, 1.0200]
const rho1 = [ 0.0000, -1.2638, -6.4761, -6.4761, -6.4761, -1.4836, -8.0298, 
               -3.8045,  0.6924,  0.6924,  0.0000,  0.0000,  0.0000]
const rho2 = [-8.8381, -3.6426,  5.5283,  5.5283,  5.5283,  0.0000,  0.0000, 
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000]
const rho3 = [ 0.0000, -5.5281, -3.0807, -3.0807, -3.0807,  0.0000,  0.0000, 
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000]
const vp0  = [11.2622, 11.0487, 15.3891, 24.9520, 29.2766, 19.0957, 39.7027, 
               20.3926,  4.1875,  4.1875,  6.8000,  5.8000,  1.4500]
const vp1 =  [ 0.0000, -4.0362, -5.3181,-40.4673,-23.6027, -9.8672,-32.6166, 
              -12.2569,  3.9382,  3.9382,  0.0000,  0.0000,  0.0000]
const vp2 =  [-6.3640,  4.8023,  5.5242, 51.4832,  5.5242,  0.0000,  0.0000, 
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000]
const vp3 =  [ 0.0000,-13.5732, -2.5514,-26.6419, -2.5514,  0.0000,  0.0000, 
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000]
const vs0 =  [ 3.6678,  0.0000,  6.9254, 11.1671, 22.3459,  9.9839, 22.3512, 
                8.9496,  2.1518,  2.1519,  3.9000,  3.2000,  0.0000]
const vs1 =  [ 0.0000,  0.0000,  1.4672,-13.7818,-17.2473, -4.9324,-18.5956, 
               -4.4597,  2.3481,  2.3481,  0.0000,  0.0000,  0.0000]
const vs2 =  [-4.4475,  0.0000, -2.0834, 17.4575, -2.0834,  0.0000,  0.0000, 
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000]
const vs3 =  [ 0.0000,  0.0000,  0.9783, -9.2777,  0.9783,  0.0000,  0.0000, 
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000]
const vpv0 = [11.2622, 11.0487, 15.3891, 24.9520, 29.2766, 19.0957, 39.7027, 
               20.3926,  0.8317,  0.8317,  6.8000,  5.8000,  1.4500]
const vpv1 = [ 0.0000, -4.0362, -5.3181,-40.4673,-23.6027, -9.8672,-32.6166, 
              -12.2569,  7.2180,  7.2180,  0.0000,  0.0000,  0.0000]
const vph0 = [11.2622, 11.0487, 15.3891, 24.9520, 29.2766, 19.0957, 39.7027, 
               20.3926,  3.5908,  3.5908,  6.8000,  5.8000,  1.4500]
const vph1 = [ 0.0000, -4.0362, -5.3181,-40.4673,-23.6027, -9.8672,-32.6166, 
              -12.2569,  4.6172,  4.6172,  0.0000,  0.0000,  0.0000]
const vsv0 = [ 3.6678,  0.0000,  6.9254, 11.1671, 22.3459,  9.9839, 22.3512, 
                8.9496,  5.8582,  5.8582,  3.9000,  3.2000,  0.0000]
const vsv1 = [ 0.0000,  0.0000,  1.4672,-13.7818,-17.2473, -4.9324,-18.5956, 
               -4.4597, -1.4678, -1.4678,  0.0000,  0.0000,  0.0000]
const vsh0 = [ 3.6678,  0.0000,  6.9254, 11.1671, 22.3459,  9.9839, 22.3512, 
                8.9496, -1.0839, -1.0839,  3.9000,  3.2000,  0.0000]
const vsh1 = [ 0.0000,  0.0000,  1.4672,-13.7818,-17.2473, -4.9324,-18.5956, 
               -4.4597,  5.7176,  5.7176,  0.0000,  0.0000,  0.0000]
const Qmu_model =  [84.6, -1., 312., 312., 312., 143., 143., 143., 80., 600., 
                     600., 600., -1.]
const Qkappa_model = [1327.7, 57823., 57823., 57823., 57823., 57823., 57823., 
                       57823., 57823., 57823., 57823., 57823., 57823.]
const eta0 = [1., 1., 1., 1., 1., 1., 1., 1., 3.3687, 3.3687, 1., 1., 1.]
const eta1 = [1., 1., 1., 1., 1., 1., 1., 1.,-2.4778,-2.4778, 1., 1., 1.]

# Dimensionalised coefficients in SI units (rad in m, rho in kg/m^-3, V in m/s)
const Rad = rad.*1.e3
const Rho0 = rho0.*1.e3
const Rho1 = rho1./a
const Rho2 = rho2.*1.e-3/a^2
const Rho3 = rho3.*1.e-6/a^3

# Return values in g/cm^3 or km/s
"`rho(r)`:  Return the density at radius `r` km in g/cm^3"
rho(r) = poly(r, rho0, rho1, rho2, rho3)
"`vp(r)`:  Return the isotropic P-wave velocity at radius `r` km in km/s"
vp(r)  = poly(r, vp0, vp1, vp2, vp3)
"`vs(r)`:  Return the isotropic S-wave velocity at radius `r` km in km/s"
vs(r)  = poly(r, vs0, vs1, vs2, vs3)
"`vc(r)`:  Return the isotropic bulk sound velocity at radius `r` km in km/s"
vc(r) = sqrt(vp(r).^2 - 4/3*vs(r).^2)
# Note anisotropic velocities are only different in 0 and 1 terms
"`vpv(r)`:  Return Vpv at radius `r` km in km/s"
vpv(r) = poly(r, vpv0, vpv1, vp2, vp3)
"`vph(r)`:  Return Vph at radius `r` km in km/s"
vph(r) = poly(r, vph0, vph1, vp2, vp3)
"`vsv(r)`:  Return Vsv at radius `r` km in km/s"
vsv(r) = poly(r, vsv0, vsv1, vs2, vs3)
"`vsh(r)`:  Return Vsh at radius `r` km in km/s"
vsh(r) = poly(r, vsh0, vsh1, vs2, vs3)
"`Qmu(r)`:  Return Qmu at radius `r` km"
Qmu(r) = Qmu_model[findlayer(r)]
"`Qkappa(r)`:  Return Qkappa at radius `r` km."
Qkappa(r) = Qkappa_model[findlayer(r)]
"`eta(r)`:  Return anisotropic parameter eta at radius `r` km."
eta(r) = poly(r, eta0, eta1, zeros(eta0), zeros(eta0))

function findlayer(r)
	r < 0. || r > a && error("PREM: radius must be in range 0 to $a km")
	for i = 1:length(rad)
		r < rad[i] && return i
	end
	return length(rad)
end

function poly(r, a0, a1, a2, a3)
	# Evaluate a set of arrays of polynomial coeffients at a certain radius
	i = findlayer(r)
	a0[i] + a1[i]*(r/a) + a2[i]*(r/a)^2 + a3[i]*(r/a)^3
end

function mass(r)
	# Compute the mass within radius r km
	# Mass of Earth ≈ 5.97e24 kg
	l = findlayer(r)
	r *= 1.e3
	M = 0.
	for i = 1:l
		if i == 1
			R0 = 0.
		else
			R0 = rad[i-1]*1.e3
		end
		if i < l
			R = rad[i]*1.e3
		else
			R = r
		end
		dM = Rho0[i]*(R^3 - R0^3)/3 + Rho1[i]*(R^4 - R0^4)/4 + 
			 Rho2[i]*(R^5 - R0^5)/5 + Rho3[i]*(R^6 - R0^6)/6
		M += dM
	end
	return 4*π*M
end

# Mass of the Earth
const Me = mass(a)

# Mass between the surface and radius r
surface_mass(r) = Me - mass(r)

# Acceleration due to gravity at radius r km
g(r) = (r == 0) ? 0. : NewtonG*mass(r)/(r*1.e3)^2

integration_func(r) = 1.e3*rho(r)*g(r)
"`pressure(r)`:  Return the pressure in Pa at radius `r` km"
function pressure(r)
	# Compute the lithostatic pressure at radius r km by integrating downward
	return 1.e3*quadgk(integration_func, r, a)[1]
end

end # Module
