Sarah Dennis 
(sarahdennis@brandeis.edu, sarahdennis98@gmail.com)

PhD research at Brandeis University, Dept of Mathematics.
In collaboration with Thomas G. Fai (tfai@brandeis.edu)

COMPARISON OF EXTENDED LUBRICATION THEORIES FOR STOKES FLOW

last updated Oct 2025

#---------------------------------------------------------------------------------------------------------------------

Part 1:
Methods for Reynolds equation in lubrication theory, and models in extended or perturbed lubrication theory.
Solves for the pressure p and velocity (u,v) for incompressible fluid at zero Reynolds number

-- reyn_run.py: Run this file for the Reynolds solver.

I. Choose domain geometry (set Example and args, see reyn_examples.py for details)

	Example = reyn_examples.BFS
	H = 2 #inlet height
	h = 1 #outlet height
	l = 2 #inlet length
	L = 4 #total length
	args = [h,H,l,L]

II. Set the grid scale N and the Fixed(U,dP) or Mixed(U,Q) boundary condition
                
	N = 100  # number of grid points in unit length

	U = 1    # u(x,0) = U, u(x,h)=v(x,0)=v(x,h)=0

	# mixed Neumann-Dirichlet pressure BC 
	Q = 1    # flux {dp/dx (x0,y) ~ Q, p(xL,y)=0}
	BC = bc.Mixed(U, Q)  

	# fixed Neumann pressure BC 
	dP = 10  # pressure drop {p(x0,y)=dP, p(xL,y)=0}
	BC = bc.Fixed(U,dP)
	
                
III. Initialize the solver
                   
	solver = control.Reynolds_Solver(Example, BC, args)
                   
IV. Choose the solver method(s) (uncomment an existing method, see reyn_control.py for details)
                
	solver.fd_solve(N, plot=plots_on)           		# finite difference solution
	solver.pwl_solve(N, plot=plots_on)          		# analytic solution for piecewise-linear h(x)
	solver.fd_adj_TG_solve(N, plot=plots_on)    		# finite difference solution to Takeuchi-Gu (2018) adjusted Reynolds equation
	solver.fd_adj_solve(N, plot=plots_on)       		# finite difference solution to velocity-adjusted Reynolds equation
	solver.fd_pert_solve(N, order=4, plot=plots_on)      # finite difference solution to perturbed Reynolds equation, (choose order = 0, 2, or 4)

OTHER IMPORTANT FILES FOR REYNOLDS SOLVER:
-- reyn_control.py : solver calls and plotting functions, main functionality for reyn_run.py
-- reyn_examples.py : all example geometries, see how args are defined for each example

#---------------------------------------------------------------------------------------------------------------------

Part 2:
	Methods for bi-harmonic Navier-Stokes equation.
	Solves for the pressure p(x,y) and velocity (u,v) for incompressible fluid at low Reynolds number


-- stokes_run.py: Run this file for the Stokes solver. 

I. Choose domain geometry (set Example and args, see stokes_examples.py for details)
                
	Example = stokes_examples.Logistic
	H = 2        # inlet height
	h = 1        # outlet height
	L = 16       # total length
	delta = 32   # slope: -delta*(H-h)/4
	args = [H, h, L, delta]

II. set boundary conditions & Reynolds number  (only Mixed(U,Q) boundary condition for biharmonic solver)

	U = 0      # u(x,0) = U, u(x,h)=v(x,0)=v(x,h)=0
	Q = 1      # flux {dp/dx (x0,y) ~ Q, p(xL,y)=0}
	Re = 0     # Re ~ U*Ly

II. Initialize the solver from stokes_control
                
	solver = control.Stokes_Solver(Example, args, U, Q, Re, max_iters=500000)    
                    
III. a) load or run an existing solution saved in ./examples    
                                       
	N=160
	solver.load_run(N)    # runs iterative solution until convergence or max_iters
	solver.load_plot(N)	  
                    
III. b) scale up an existing solution to a new grid size
                                   
	N_old = 80	  
	N_new = 100   # if this grid size already exists, load_scale will overwrite it
  
	solver.load_scale(N_old, N_new)
	solver.load_run(N_new)

                    
III. c) run a new solution
                
	N = 100       # if this grid size already exists, new_run will overwrite it
	solver.new_run(N)

OTHER IMPORTANT FILES FOR STOKES SOLVER:

OTHER IMPORTANT FILES FOR REYNOLDS SOLVER:
-- stokes_control.py : solver calls and plotting functions, main functionality for stokes_run.py
-- stokes_examples.py : all example geometries, see how args are defined for each example
                