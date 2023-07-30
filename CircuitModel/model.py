# contains ODE solvers
import scipy

def diff(t, x, params):
	"""
	@param[in]	t 	Time parameter (float)
	@param[in] 	x	Initial conditions (i0, i1, i2, Vatp)
	@param[in]	params	Global 
	"""
	## Constants
	# Potentials
	Vcho = 0.0
	Vfat = 0.0
	# Resistances
	Rbox = 0.0
	Roxphos = 0.0
	Retc = 0.0
	Rcvant = 0.0
	Rl = 0.0
	# System capacitance
	Cpcr = 0.0
	## Non-linear components
	Rg = 0.02065*x[3]**3 + 0.128*x[3] - 0.017
	Rcvant = ((8-x[3])**2) / ( ((8-x[3])**2) + 1)

	f_ = []
	# loop currents
	f.append( Rg + Rbox )
	# Rg
	f.append( )
	# tau
	f.append( R*Cpcr )

# Circuit model
def CircuitModel():
	# user-defined function for solving diffeq, timespan (in seconds), initial conditions
    sol = scipy.integrate.solve_ivp(diff, [0, 1], xo)

CircuitModel()
