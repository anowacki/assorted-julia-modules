"""
AK135 evaluates the seismic velocities, density and attenuation for the
1D Earth model AK135.
"""
module AK135

export
	rho,
	vp,
	vs

"Number of radial points in AK135 model (number of layers + 1)"
const N = 136
"Radius of Earth in AK135 model"
const R = 6371.

const rad = [
       0.,   50.,  101.,  152.,  202.,  253.,  304.,  354.,
     405.,  456.,  507.,  557.,  659.,  710.,  760.,  811.,
     862.,  912.,  963., 1014., 1064., 1115., 1166., 1217.,
    1217., 1267., 1317., 1368., 1418., 1468., 1519., 1569.,
    1670., 1720., 1770., 1821., 1871., 1921., 1972., 2022.,
    2072., 2123., 2173., 2223., 2274., 2324., 2374., 2425.,
    2475., 2525., 2576., 2626., 2676., 2727., 2777., 2827.,
    2878., 2928., 2978., 3029., 3079., 3129., 3180., 3230.,
    3280., 3331., 3381., 3431., 3479., 3479., 3531., 3581.,
    3631., 3631., 3681., 3731., 3779., 3829., 3878., 3928.,
    3977., 4027., 4076., 4126., 4175., 4225., 4274., 4324.,
    4373., 4423., 4472., 4522., 4571., 4621., 4670., 4720.,
    4769., 4819., 4868., 4918., 4967., 5017., 5066., 5116.,
    5165., 5215., 5264., 5314., 5363., 5413., 5462., 5512.,
    5561., 5611., 5661., 5711., 5711., 5761., 5811., 5861.,
    5911., 5961., 5961., 6011., 6061., 6111., 6161., 6161.,
    6206., 6251., 6293., 6336., 6336., 6351., 6351., 6371.] # / km

# Indices of nodes in each major layer
const irad_inner_core = 1:24
const irad_outer_core = 25:69
const irad_mantle = 70:136

const rho_model = [
    13.0122, 13.0117, 13.0100, 13.0074, 13.0036, 12.9988, 12.9929, 12.9859,
    12.9779, 12.9688, 12.9586, 12.9474, 12.9217, 12.9070, 12.8917, 12.8751,
    12.8574, 12.8387, 12.8188, 12.7980, 12.7760, 12.7530, 12.7289, 12.7037,
    12.1391, 12.1133, 12.0867, 12.0593, 12.0311, 12.0001, 11.9722, 11.9414,
    11.8772, 11.8437, 11.8092, 11.7737, 11.7373, 11.6998, 11.6612, 11.6216,
    11.5809, 11.5391, 11.4962, 11.4521, 11.4069, 11.3604, 11.3127, 11.2639,
    11.2137, 11.1623, 11.1095, 11.0555, 11.0001, 10.9434, 10.8852, 10.8257,
    10.7647, 10.7023, 10.6385, 10.5731, 10.5062, 10.4378, 10.3679, 10.2964,
    10.2233, 10.1485, 10.0722,  9.9942,  9.9145,  5.7721,  5.7458,  5.7196,
     5.6934,  5.4387,  5.4176,  5.3962,  5.3748,  5.3531,  5.3313,  5.3092,
     5.2870,  5.2646,  5.2420,  5.2192,  5.1963,  5.1732,  5.1499,  5.1264,
     5.1027,  5.0789,  5.0548,  5.0306,  5.0062,  4.9817,  4.9570,  4.9321,
     4.9069,  4.8817,  4.8562,  4.8307,  4.8050,  4.7790,  4.7528,  4.7266,
     4.7001,  4.6735,  4.6467,  4.6198,  4.5926,  4.5654,  4.5162,  4.4650,
     4.4118,  4.3565,  4.2986,  4.2387,  3.9201,  3.9206,  3.9218,  3.9233,
     3.9273,  3.9317,  3.5068,  3.4577,  3.4110,  3.3663,  3.3243,  3.3243,
     3.3711,  3.4268,  3.3450,  3.3200,  2.9200,  2.9200,  2.7200,  2.7200] # g/cm^3

const vp_model = [
    11.2622, 11.2618, 11.2606, 11.2586, 11.2557, 11.2521, 11.2477, 11.2424,
    11.2364, 11.2295, 11.2219, 11.2134, 11.1941, 11.1830, 11.1715, 11.1590,
    11.1457, 11.1316, 11.1166, 11.0983, 11.0850, 11.0718, 11.0585, 11.0427,
    10.2890, 10.2854, 10.2745, 10.2565, 10.2329, 10.2049, 10.1739, 10.1415,
    10.0768, 10.0439, 10.0103,  9.9761,  9.9410,  9.9051,  9.8682,  9.8304,
     9.7914,  9.7513,  9.7100,  9.6673,  9.6232,  9.5777,  9.5306,  9.4814,
     9.4297,  9.3760,  9.3205,  9.2634,  9.2042,  9.1426,  9.0792,  9.0138,
     8.9461,  8.8761,  8.8036,  8.7283,  8.6496,  8.5692,  8.4861,  8.4001,
     8.3122,  8.2213,  8.1283,  8.0382,  8.0000, 13.6601, 13.6570, 13.6533,
    13.6498, 13.6498, 13.5899, 13.5311, 13.4741, 13.4156, 13.3584, 13.3017,
    13.2465, 13.1895, 13.1337, 13.0786, 13.0226, 12.9663, 12.9093, 12.8524,
    12.7956, 12.7384, 12.6807, 12.6226, 12.5638, 12.5030, 12.4427, 12.3813,
    12.3181, 12.2558, 12.1912, 12.1247, 12.0571, 11.9891, 11.9208, 11.8491,
    11.7768, 11.7020, 11.6265, 11.5493, 11.4704, 11.3897, 11.3068, 11.2228,
    11.1355, 11.0553, 10.9222, 10.7909, 10.2000, 10.0320,  9.8640,  9.6962,
     9.5280,  9.3601,  9.0302,  8.8476,  8.6650,  8.4822,  8.3007,  8.3007,
     8.1750,  8.0505,  8.0450,  8.0400,  6.5000,  6.5000,  5.8000,  5.8000] # km/s

const vs_model = [
    3.6678, 3.6675, 3.6667, 3.6653, 3.6633, 3.6608, 3.6577, 3.6540,
    3.6498, 3.6450, 3.6396, 3.6337, 3.6202, 3.6130, 3.6044, 3.5957,
    3.5864, 3.5765, 3.5661, 3.5551, 3.5435, 3.5314, 3.5187, 3.5043,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 7.2817, 7.2700, 7.2593,
    7.2485, 7.2485, 7.2253, 7.2031, 7.1804, 7.1584, 7.1368, 7.1144,
    7.0932, 7.0722, 7.0504, 7.0286, 7.0069, 6.9852, 6.9625, 6.9416,
    6.9194, 6.8972, 6.8743, 6.8517, 6.8289, 6.8056, 6.7820, 6.7579,
    6.7323, 6.7070, 6.6813, 6.6554, 6.6285, 6.6009, 6.5728, 6.5431,
    6.5131, 6.4822, 6.4514, 6.4182, 6.3860, 6.3519, 6.3164, 6.2799,
    6.2424, 6.2100, 6.0898, 5.9607, 5.6104, 5.5047, 5.3989, 5.2922,
    5.1864, 5.0806, 4.8702, 4.7832, 4.6964, 4.6094, 4.5184, 4.5184,
    4.5090, 4.5000, 4.4900, 4.4800, 3.8500, 3.8500, 3.4600, 3.4600] # km/s

const Qkappa_model = [
      601.27,   601.32,   601.46,   601.70,   602.05,   602.49,   603.04,   603.69,
      604.44,   605.28,   606.26,   607.31,   609.74,   611.18,   612.62,   614.21,
      615.93,   617.78,   619.71,   621.50,   624.08,   626.87,   629.89,   633.26,
    57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00,
    57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00,
    57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00,
    57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00,
    57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00, 57822.00,
    57822.00, 57822.00, 57822.00, 57822.00, 57822.00,   723.12,   725.11,   726.87,
      722.73,   933.21,   940.88,   952.00,   960.36,   968.46,   976.81,   985.63,
      990.77,   999.44,  1008.79,  1018.38,  1032.14,  1042.07,  1048.09,  1058.03,
     1064.23,  1070.38,  1085.97,  1097.16,  1108.58,  1120.09,  1127.02,  1134.01,
     1141.32,  1148.76,  1156.04,  1163.16,  1170.53,  1178.19,  1186.06,  1193.99,
     1202.04,  1210.02,  1217.91,  1226.52,  1234.54,  1243.02,  1251.69,  1260.68,
     1269.44,  1277.93,  1311.17,  1350.54,   428.69,   425.51,   422.55,   419.94,
      417.32,   413.66,   377.93,   366.34,   355.85,   346.37,   338.47,   200.97,
      188.72,   182.57,   182.03,   182.03,   972.77,   972.77,  1368.02,  1368.02]

const Qmu_model = [
     85.03,  85.03,  85.03,  85.03,  85.03,  85.03,  85.03,  85.03,
     85.03,  85.03,  85.03,  85.03,  85.03,  85.03,  85.03,  85.03,
     85.03,  85.03,  85.03,  85.03,  85.03,  85.03,  85.03,  85.03,
      0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
      0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
      0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
      0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
      0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,
      0.00,   0.00,   0.00,   0.00,   0.00, 273.97, 273.97, 273.97,
    271.74, 350.88, 354.61, 359.71, 363.64, 367.65, 371.75, 375.94,
    378.79, 383.14, 387.60, 392.16, 398.41, 403.23, 406.50, 411.52,
    414.94, 418.41, 425.53, 431.03, 436.68, 442.48, 446.43, 450.45,
    454.55, 458.72, 462.96, 467.29, 471.70, 476.19, 480.77, 485.44,
    490.20, 495.05, 500.00, 505.05, 510.20, 515.46, 520.83, 526.32,
    531.91, 537.63, 543.48, 549.45, 172.93, 170.82, 168.78, 166.80,
    164.87, 162.50, 146.57, 142.76, 139.38, 136.38, 133.72,  79.40,
     76.55,  76.06,  75.60,  75.60, 403.93, 403.93, 599.99, 599.99]

# Generate functions to evaluate properties at given radius
for (v, doc_name, doc_unit) in zip((:rho, :vp, :vs), ("rho", "Vp", "Vs"), ("g/cm^3", "km/s", "km/s"))
	func_name, var_name = symbol("$v"), symbol("$(v)_model")
	@eval begin
		@doc """
		`$($func_name)(r) -> $($doc_name)`
		
		Return the value of $($doc_name) in the AK135 model at radius `r` km ($($doc_unit)).
		""" ->
		($func_name)(r) = interp(r, $var_name)
	end
end
for (v, doc_name) in zip((:Qkappa, :Qmu), ("Qkappa", "Qmu"))
	func_name, var_name = symbol("$v"), symbol("$(v)_model")
	@eval begin
		@doc """
		`$($func_name)(r) -> $($doc_name)`
		
		Return the value of $($doc_name) in the AK135 model at radius `r` km.
		""" ->
		($func_name)(r) = $var_name[findlayer(r)]
	end
end

"""
`interp(r, v) -> V`

Return the interpolated value `V` from the array `v` at radius `r` km,
providing the array `v` is the same length as the number of layers in
AK135 (136, confusingly).
"""
function interp(r, v)
	length(v) == N || error("AK135.interp: Size of `v` is not $(N)")
	# i is the index of the node *below* this radius
	i = findlayer(r)
	# Exactly at a node or at the surface
	if rad[i] == r || i == N
		return v[i]
	# Discontinuity: choose the lower layer
	elseif rad[i] == rad[i+1]
		return v[i]
	else
		return (v[i]*(rad[i+1] - r) + v[i+1]*(r - rad[i]))/(rad[i+1] - rad[i])
	end
end

"""
`findlayer(r) -> i`

Return the node index `i` at or above which the radius `r` occurs.
"""
function findlayer(r)
	rad[1] <= r <= rad[end] ||
		error("AK135.findlayer: Radius must be in range $(rad[1]) to $(rad[end]) km")
	@inbounds for i = 1:N
		r < rad[i] && return i
	end
	N
end

end # module