Sarah Dennis 
(sarahdennis@brandeis.edu, sarahdennis98@gmail.com)

PhD research at Brandeis University, Dept of Mathematics.
Collaboration with Thomas G. Fai (tfai@brandeis.edu)

Comparison of Stokes and Reynolds models of low Reynolds number fluids in cornered domains

WORK IN PROGRESS (Jan 2025) 

#---------------------------------------------------------------------------------------------------------------------

Part 1:
Methods for bi-harmonic Navier-Stokes equation

Appropriate for low Reynolds number flows in 2D piecewise linear domains

-- stokes_run.py: Choose domain geometry from stokes_examples. Choose run method from stokes_control.

-- stokes_examples.py: wrapper examples of stokes_heights, variations on 2D sliders bearings. 

-- stokes_control: helpers linking to stokes_solver: iteration, file saving, plotting

-- stokes_heights: boundary condition handling, wrapper for domain

-- stokes_solver: iterative finite difference method for stream-velocity biharmonic equation

-- stokes_pressure: 2D finite difference and contour integral for velocity-pressure equation

-- stokes_convergence: grid convergence tests

-- ./examples: saved stream function solutions to specific examples, load_run and load_plot methods in stokes_control 

#---------------------------------------------------------------------------------------------------------------------

Part 2:
Methods for Reynolds equation

Appropriate for lubrication flows in 2D piecewise linear domains

-- reyn_run.py: Choose domain geometry from reyn_examples. Choose run method from reyn_control.

-- reyn_examples.py: wrapper examples of reyn_heights, variations on 2D slider bearings.

-- reyn_control: solution helpers linking to reyn_pressures_finDiff and reyn_pressures_pwlinear

-- reyn_heights: piecewise linear, sinusoidal, and random height constructors, wrapper for domain

-- reyn_pressures_finDiff: second order finite difference method to Reynolds equation for pressure

-- reyn_pressures_pwlinear: analytic solution to Reynolds equation with piecewise linear height

-- reyn_velocity: double integration for Reynolds pressure-velocity equation

-- reyn_convergence: grid convergence tests

#---------------------------------------------------------------------------------------------------------------------

Part 3: 
Additional/Shared methods

-- domain.py: inhereted by both reyn_heights and stokes_heights

-- graphics.py: All basic graphing methods. Referenced by *_control and *_convergence

-- read_write.py: For stokes_control and stokes_solver iterative method, file names assigned in stokes_examples
