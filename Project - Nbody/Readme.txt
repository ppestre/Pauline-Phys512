Final Project

Commented code for each of the questions can be found in nbody_questions.py.

Part 1:
	A single particle at rest remains at rest. See 'single.gif'.
	If I zoom way in on the plots, I can see motion on ~10^-14, 
	likely from floating point error in calculating the acceleration.

Part 2:
	Two particles on the x-axis, start moving in +-y direction and continue in approximately circular orbit.
	See 'cicular_orbit.gif' for a snapshot of each time step, or 'cicular orbit.png' for a plot after several orbits.
	The orbits are a bit wobbly, so there might be some rounding error, or some finetuning of the initial velocity needed.
	(v_i approx from Virial theorem)

Part 3: 
	For periodic boundary conditions, see 'circular_BC_many'. Didn't get around to non-periodic :(
	Energy is not very well conserved. There are a lot of approximations in this method, 
	and poor accounting for particles that disappear off the edge. 

Part 4: 
	Start with a 1/k^3 grid in kspace. Multiply by gaussian noise, and inverse FT to get starting density. 
	Code still needs so tweaks to handle this well.
	
