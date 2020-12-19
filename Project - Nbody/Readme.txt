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
	See gifs by number of particles and steps. Tried a few, the big ones took a while.
	Potential energy is conserved to within floating point error. 
	Kinetic energy is conserved to within 0.3% (in one example). I expect some slight rounding error from approximating 
	a point to have the average acceleration value of its grid cell, even if it's near the edge.

Part 4: 
	Start with a 1/k^3 grid in kspace. Multiply by gaussian noise, and inverse FT to get starting density. 
	The universe looks a bit lumpy! Could use some tweaks to make it a bit more efficient, but it's somethin.
	Seeing some weird edge behaviour?
	
