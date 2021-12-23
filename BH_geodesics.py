"""
BH_geodesic.py
==============

Author:
	Stefano Schmidt 2021

Description:
	Script to compute the geodesic motion around a Black Hole (or a planet), both in Schwartschild solution and in the Newtonian approximation

Usage:
	python BH_geodesic.py example.ini

You need to have installed numpy, scipy and matplotlib

Here's an example of ini file that simulates an hyperbolic encounter:

[HYPERBOLIC]
L= 4.02
phi_0= 0.
r_0 = 1.5
r_dot_0 = -1.02
GR = 1
t_max= 2e3
max_step = 1e11
relative_units = 1

Here's an example to simulate an elliptic orbit with precession:

[ELLIPTIC]
L= 10.02
phi_0= 0.
r_0 = 1.5
r_dot_0 = -0.2
GR = 1
t_max= 3e5
max_step = 1e11
relative_units = 1

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import configparser, sys, os

############################################################################################################

def to_cartesian(r, phi):
	return r*np.cos(phi), r*np.sin(phi)

def V(r, L, GR = False):
		#Mass M of the body is set to 1
		#All the units are in units of M
	r = r+1e-5
	E = 0.5* np.square(L/r) - 1./r
	
	if GR:
		E = E - np.square(L)/r**3

	return E

def E(r, r_dot, L, GR = False):
	return 0.5*r_dot**2 + V(r, L, GR)

def dV_dr(r, L, GR = False):
		#Mass M of the body is set to 1
		#All the units are in units of M
	r = r+1e-5
	E_dot = - np.square(L)/r**3 + 1./r**2
	
	if GR:
		E_dot = E_dot + 3*np.square(L)/r**4

	return E_dot


def LHS(t, x, L, GR):
		#x is [r, phi, r_dot]
	if x.ndim ==1: x = x[None,:]
	assert x.shape[1] == 3
	return np.stack([x[:,2], L/x[:,0]**2, -dV_dr(x[:,0], L, GR)], axis = 1)

############################################################################################################

def integrate_geodesics(vals):
	"Perform the integration of with a given set of parameters"
	#########
	# Solving for the dynamics (plus some output to screen)
	#########

	print("##### Integrating {} #####".format(vals['name']))
	print("\tIntial values: ")
	for k_, v_ in vals.items():
		print('\t\t{}: {}'.format(k_, vals[k_]))
	print("\t\tInitial energy: {}".format(E(vals['r_0'], vals['r_dot_0'], vals['L'], vals['GR'])))

	sol = solve_ivp(LHS, [0,vals['t_max']], [vals['r_0'], vals['phi_0'], vals['r_dot_0']], args =  (vals['L'], vals['GR']),
			t_eval = np.linspace(0, vals['t_max'], 10000), max_step = vals['max_step'])
	print("\tSolver message: {}".format(sol['message']))

	sol_newton = None
	if vals['compare']:
		sol_newton = solve_ivp(LHS, [0,vals['t_max']], [vals['r_0'], vals['phi_0'], vals['r_dot_0']], args =  (vals['L'], False),
			t_eval = np.linspace(0, vals['t_max'], 10000), max_step = vals['max_step'])
	
	return sol, sol_newton


def plot(vals, sol, sol_newton):
	"Plot the result of the integration"
	#########
	# Plotting
	#########

	t, dyn = sol['t'], sol['y'].T
	x, y = to_cartesian(dyn[:,0], dyn[:,1])
	label = 'GR' if vals['GR'] else 'Newton'

	plt.figure()
	plt.title("Radius")
	plt.plot(t, dyn[:,0], label = label)
	plt.xlabel("t/M")
	plt.ylabel("r/M")
	plt.legend()

	plt.figure()
	plt.title("Angle")
	plt.plot(t, dyn[:,1], label = label)
	plt.xlabel("t/M")
	plt.ylabel(r"$\phi$")
	plt.legend()

	plt.figure()
	plt.title("Radial velocity")
	plt.plot(t, dyn[:,2], label = label)
	plt.xlabel("t/M")
	plt.ylabel(r"$\dot{r}$")
	plt.legend()

	plt.figure()
	plt.title("Trajectory")
	plt.plot(x, y, label = label)
	plt.xlabel("x/M")
	plt.ylabel("y/M")
	if vals['compare']:
		x_N, y_N = to_cartesian(sol_newton['y'][0,:], sol_newton['y'][1,:])
		plt.plot(x_N, y_N, label = 'Newton')
	plt.legend()
	plt.plot(0, 0, 'o', lw = .5, c = 'red')
	plt.plot(x[0], y[0], 'x', lw = 4, c = 'red')
	circle=plt.Circle((0,0),2, facecolor = 'k', edgecolor = 'k', lw = 2)
	plt.gca().add_patch(circle)
	unit_x = plt.xlim()[1]-plt.xlim()[0]
	unit_y = plt.ylim()[1]-plt.ylim()[0]
	if unit_x>unit_y:
		scale = unit_x/unit_y
		plt.ylim([plt.ylim()[0]*scale, plt.ylim()[1]*scale])
	else:
		scale = unit_y/unit_x
		plt.xlim([plt.xlim()[0]*scale, plt.xlim()[1]*scale])

	if vals['save_trajectory']: plt.savefig('{}/{}_trajectory.jpeg'.format(vals['folder'], vals['name']))
	if vals['show']: plt.show()

	return


############################################################################################################

#########
# Setting initial conditions: parsing ini file
#########
default_optional_vals = {
	'compare': 0,
	'show': 1,
	'save_trajectory': 0,
	'relative_units': 0,
	'max_step': 1e4
}

if len(sys.argv)> 1:
	if sys.argv[1] == '--help':
		print(__doc__)
		quit()
	if not os.path.isfile(sys.argv[1]):
		raise ValueError("The given input file '{}' doesn't exist".format(sys.argv[1]))
	else:
		print("Reading ini from file: {} ".format(sys.argv[1]))
	config = configparser.ConfigParser()
	config.optionxform = str
	config.read(sys.argv[1])
	folder = os.path.dirname(sys.argv[1])
	if folder == '': folder = '.'
		
		#loops on ini sections. Each is a different integration
	vals_list = []
	for sec in config.sections():
		vals = dict(config[sec])
			#converting to float
		for k_, v_ in vals.items():
			vals[k_] = float(v_)

		for k_, v_ in default_optional_vals.items():
			if not (k_ in vals): vals[k_] =  v_
		
		if vals['relative_units']:
			vals['r_dot_0']= vals['r_dot_0']*np.sqrt(2*np.abs(V(vals['r_0'], vals['L'], vals['GR'])))
		
		if vals['compare']:
			vals['GR'] = 1
		
		vals['folder'] = folder
		vals['name'] = sec
		vals_list.append(vals)
	
else:
		#setting some default values
	vals = {
		'name': 'test',
		'folder': '.',
		'L': 4.02,
		'phi_0': 0.,
		'r_0': 169.68,
		'GR': True,
		't_max': 1e4,
	}
	vals['r_dot_0']= -.2*np.sqrt(2*np.abs(V(vals['r_0'], vals['L'], vals['GR'])))
	for k_, v_ in default_optional_vals.items(): vals[k_] =  v_
	vals_list = [vals]

for vals in vals_list:
	sol, sol_newton = integrate_geodesics(vals)

	plot(vals, sol, sol_newton)



























