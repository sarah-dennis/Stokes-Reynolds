Sarah Dennis

PhD research at Brandeis University, Dept of Mathematics

Comparison of Stokes and Reynolds models of low reynolds number fluids in cornered domains

WORK IN PROGRESS (Nov 2024) 
email me: sarahdennis@brandeis.edu, sarahdennis98@gmail.com

Part 1:
Methods for Biharmonic Navier Stokes equations

Appropriate for low Reynolds number flows in 2D piecewise linear domains

-- stokes_run.py: Choose domain geometry from stokes_examples. Choose run method from stokes_control. 
-- stokes_examples.py: wrapper examples of stokes_heights, variations on 2D textured sliders. 
-- stokes_control: solution helpers linking to stokes_solver
-- stokes_heights: boundary condition handling, wrapper for domain
-- stokes_solver: iterative finite difference method for stream-velocity biharmonic equation
-- stokes_pressure: 2D finite difference and contor integral for velocity-pressure poisson equation
-- stokes_convergence: grid convergence tests
-- /examples:  

Part 2:
Methods for Reynolds equation

Appropriate for lubrication flows in 2D piecewise linear domains

-- reyn_run.py: Choose domain geomoetry from reyn_examples. Choose run method from reyn_control
-- reyn_examples.py: wrapper examples of reyn_heights, variations on 2D textured sliders.
-- reyn_control: solution helpers linking to reyn_pressures_finDiff and reyn_pressures_pwlinear
-- reyn_heights: piewise linear handling, wrapper for domain
-- reyn_pressures_finDiff: second order finite difference method to Reynolds equation for pressure
-- reyn_pressures_pwlinear: analytic solution to Reynolds equation with piecewise linear height
-- reyn_velocity: double integration for Reynolds pressure-velocity equation
-- reyn_convergence: grid convergence tests

Part 3: 
Additional/Shared:
-- domain.py: inhereted by both reyn_heights and stokes_heights
-- graphics.py: All basic graphing methods. Referenced by *_control and *_convergence
-- error.py: All basic error methods. Referenced by *_convergence
-- read_write.py: For stokes_control and stokes_solver iterative method, file names assigned in stokes_examples
