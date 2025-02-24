Sarah Dennis 
(sarahdennis@brandeis.edu, sarahdennis98@gmail.com)

PhD research at Brandeis University, Dept of Mathematics.
in c4ollaboration with Thomas G. Fai (tfai@brandeis.edu)

Comparison of Comparison of lubrication theory and Stokes flow models in step bearings with flow separation

WORK IN PROGRESS (updated Feb 2025) 
#---------------------------------------------------------------------------------------------------------------------

Part 1:
Methods for Reynolds equation in lubrication theory 
Solves for the pressure p(x) and velocity field (u,v) for incompressible fluid at zero Reynolds number
    in a 2D domain of length L:=xf-x0 and height H:=h(x)-y0 > 0.
    No-slip condition at fluid-surface interfaces. Relative surface velocity U:=u(x,y0)-u(x,h), u(x,h):=0, V:=0.
    Prescribed pressure drop dP:= p(xf,y)-p(x0,y), p(xf,0):=0

-- reyn_run.py: Run this file for the Reynolds solver.

                I. Set values for U, dP, and the grid scale N
                
                    N = 100
                    U = 1
                    dP = 0
                
                II. Choose domain geometry from reyn_examples. (Or uncomment an existing example)
            
                    # backward facing step with H/h=2 and L/l = 4
                    Example = examples.BFS
                    args = [H=2,l=4]     
                
                III. Initialize the solver
                   
                    solver = control.Reynolds_Solver(Example)
                   
                IV. Choose the solver method from reyn_control (Or uncomment an existing method)
                
                    solver.fd_solve(N, plot=plots_on)           # finite difference solution
                    solver.pwl_solve(N, plot=plots_on)          # analytic solution for piecewise-linear h(x)
                    solver.fd_adj_solve(N, plot=plots_on)       # finite difference solution to *adjusted Reynolds equation*
                    
-- reyn_examples.py: wrapper examples of reyn_heights, variations on 2D slider bearings.

-- reyn_control: solution helpers linking to reyn_pressures_finDiff and reyn_pressures_pwlinear

-- reyn_heights: piecewise linear, sinusoidal, and random height constructors, wrapper for domain

-- reyn_pressure: Pressure super class for various solution methods. Adjusted Reynolds pressure computed.

-- reyn_pressures_finDiff: second order finite difference method to Reynolds equation for pressure

-- reyn_pressures_pwlinear: analytic solution to Reynolds equation with piecewise linear height

-- reyn_velocity: double integration for Reynolds pressure-velocity equation

-- reyn_convergence: grid convergence tests

#---------------------------------------------------------------------------------------------------------------------

Part 2:
Methods for bi-harmonic Navier-Stokes equation.
Solves for the pressure p(x,y) and velocity field (u,v) for incompressible fluid at low Reynolsd number
    in a 2D domain of length L:=xf-x0 and height H:=h(x)-y0 > 0.
    No-slip condition at fluid-surface interfaces. Relative surface velocity U:=u(x,y0)-u(x,h), u(x,h):=0, V:=0.
    Prescribed flux Q:= sum_y u(x,y) dx
    Fully developed laminar flow profile (dP/dx ~ Q) at x0 and xf


-- stokes_run.py: Run this file for the Stokes solver. 

                I. Choose domain geometry from stokes_examples. (Or uncomment an existing example)
                
                    # backward facing step with H/h=2 and L/l = 4,  Re=0, Q=2, U=0
                    example = examples.BFS_H2L4_Re0_Q2_U0 
                    
                II. Initialize the solver from stokes_control
                
                    solver = control.Stokes_Solver(example)    
                    
                III. a) load or run an existing solution saved in ./examples    
                      
                    # ./examples/BFS_H2L4_Re0_Q2_U0/BFS_H2L4_Re0_Q2_U0_N160.csv
                    
                    N=160
                    solver.load_run(N) # runs iterative solution until convergence 
                    solver.load_plot(N)
                    
                III. b) scale up an existing solution to a new grid size
                
                    # ./examples/BFS_H2L4_Re0_Q2_U0/BFS_H2L4_Re0_Q2_U0_N80.csv
                   
                    N_old = 80
                    N_new = 100 #if this grid size already exists, load_scale will overwrite it
                    
                    solver.load_scale(N_old, N_new)
                    solver.load_run(N_new)

                    
                III. c) run a new solution (may take a while)
                
                    N = 100 #if this grid size already exists, new_run will overwrite it
                    solver.new_run(N)

                

-- stokes_examples.py: wrapper examples of stokes_heights, variations on 2D sliders bearings. 

-- stokes_control: helpers linking to stokes_solver: iteration, file saving, plotting

-- stokes_heights: boundary condition handling, wrapper for domain

-- stokes_solver: iterative finite difference method for stream-velocity biharmonic equation

-- stokes_pressure: 2D finite difference and contour integral for velocity-pressure equation

-- stokes_convergence: grid convergence tests

-- ./examples: saved stream function solutions to specific examples, load_run and load_plot methods in stokes_control 

#---------------------------------------------------------------------------------------------------------------------

Part 3: 
Additional/Shared methods

-- domain.py: inhereted by both reyn_heights and stokes_heights

-- graphics.py: All basic graphing methods. Referenced by *_control and *_convergence

-- read_write.py: For stokes_control and stokes_solver iterative method, file names assigned in stokes_examples
