# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:05 2024

@author: sarah
"""
import numpy as np
import graphics

import time
import domain as dm
import reyn_velocity as rv
import reyn_pressure as rp 
import reyn_perturbed as rpert
import reyn_boundary as bc

from reyn_heights import PWC_Height, PWL_Height, make_PWC, make_PWL

i_test_scale=2


lenx = 4
leny = 2
x_start = 6
y_start = 0
x_stop= x_start + lenx
y_stop = y_start + leny

# colorbar min max
vel_max = 5
p_min= 60
p_max = 110
        
log_linthresh=1e-8  
        
class Reynolds_Solver: 
    def __init__(self, Example, BC, args=None):
        self.Example = Example #initialize height = Example(args) in the solver
        self.args = args
        
        self.BC = BC        
        
        


    def fd_solve(self, N, plot=True, scaled=False, zoom=False,inc=False, uv=False):
        solver_title = "Reynolds Finite Difference"
        
        height = self.Example(self.args, N)

        t0 = time.time()
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)

        reyn_pressure = rp.FinDiff_ReynPressure(height, self.BC)
        tf = time.time()
    
       
        print('fd time: ', tf-t0)
        velocity = None
        if plot: 
            velocity = rv.ReynVelocity(height, self.BC, ps=reyn_pressure.ps_1D)
            self.p_plot(height, reyn_pressure, velocity.Q, solver_title, scaled, zoom)
            self.v_plot(self.BC, height, velocity, reyn_pressure, solver_title, scaled, zoom, inc, uv)

        return reyn_pressure, velocity, tf-t0

    def pwc_schur_solve(self, N, plot=True, scaled=False, zoom=False, inc=False, uv=False):
        solver_title = "Reynolds PWC Schur Complement"
        
        height = self.Example(self.args, N)

        if not isinstance(height, PWC_Height):
            height = make_PWC(height)
    
                
        
        
        t0 = time.time()
        pressure = rp.PwcSchur_ReynPressure(height, self.BC)
        tf = time.time()
        print('pwc schur time: ', tf-t0)
        
       
        velocity = None
        if plot: 
             height.hxs = np.zeros(height.Nx)
             velocity = rv.ReynVelocity(height, self.BC, ps=pressure.ps_1D)
     
             self.p_plot(height, pressure, velocity.Q, solver_title, scaled, zoom)
             self.v_plot(self.BC, height, velocity, pressure, solver_title, scaled, zoom, inc, uv)

        return pressure, velocity, tf-t0
    
    def pwl_schur_solve(self, N, plot=True, scaled=False, zoom=False,inc=False, uv=False):
        solver_title = "Reynolds PWL Schur Complement" #" Piecewise Linear"
        
        height = self.Example(self.args,N)
        
        if not (isinstance(self.BC, bc.Mixed)):
            raise TypeError('Only prescribed flux for pwc schur complement')
        
        if not (isinstance(height, PWL_Height) or isinstance(height, PWC_Height)) :
            height = make_PWL(height)

        t0 = time.time()
        pressure = rp.PwlSchur_ReynPressure(height, self.BC)
        tf = time.time()
        print('pwl schur time: ', tf-t0)
        
        
        
    
        velocity = None
        if plot:
            height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
            velocity = rv.ReynVelocity(height, self.BC, ps=pressure.ps_1D)
            self.p_plot(height, pressure, velocity.Q, solver_title, scaled, zoom)
            self.v_plot(self.BC, height, velocity, pressure, solver_title, scaled, zoom, inc, uv)

        
        return pressure, velocity, tf-t0
    
    def pwc_schur_parallel_solve(self, N, plot=True, scaled=False, zoom=False, inc=False, uv=False):
        solver_title = "Reynolds Schur Complement"
        
        height = self.Example(self.args, N)

        if not isinstance(height, PWC_Height):
            height = make_PWC(height)
        
        t0 = time.time()
                
        height.hxs = np.zeros(height.Nx)

        pressure = rp.PwcSchur_parallel_ReynPressure(height, self.BC)
        tf = time.time()
        print('schur parallel time: ', tf-t0)
        
     
        velocity = rv.ReynVelocity(height, self.BC, ps=pressure.ps_1D)
        if plot:        
             
             self.p_plot(height, pressure, velocity.Q, solver_title, scaled, zoom)
             self.v_plot(self.BC, height, velocity, pressure, solver_title, scaled, zoom, inc, uv)

        return pressure, velocity, tf-t0
       
    def pwl_gmres_solve(self, N, plot=True, scaled=False, zoom=False,inc=False, uv=False):
        solver_title = "Reynolds GMRes" #" Piecewise Linear"
        
        height = self.Example(self.args,N)
        
        t0 = time.time()
       

        if not (isinstance(height, PWL_Height) or isinstance(height, PWC_Height)) :
            height = make_PWL(height)
       
        pressure = rp.PwlGMRes_ReynPressure(height, self.BC)
        tf = time.time()
        
        print('pwl gmres time: ', tf-t0)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        velocity = rv.ReynVelocity(height, self.BC, ps=pressure.ps_1D)
    
        
        if plot:
            self.p_plot(height, pressure, velocity.Q, solver_title, scaled, zoom)
            self.v_plot(self.BC, height, velocity, pressure, solver_title, scaled, zoom, inc, uv)

        
        return pressure, velocity, tf-t0
     
    def fd_adj_solve(self, N, plot=True, scaled=False, zoom=False, inc=False, uv=False):
        solver_title = "VA-ELT"
        t0 = time.time()
        height = self.Example(self.args, N)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        height.h2xs = dm.center_second_diff(height.hs, height.Nx, height.dx)
        height.h3xs = dm.center_third_diff(height.hs, height.Nx, height.dx)
        
        adj_pressure = rp.VelAdj_ReynPressure(height, self.BC)
        tf = time.time()
        adj_velocity = rv.VelAdj_ReynVelocity(height,self.BC, adj_pressure)

        if plot:
            self.p_plot(height, adj_pressure , adj_velocity.Q, solver_title, scaled, zoom)
            self.v_plot(self.BC, height, adj_velocity, adj_pressure, solver_title, scaled, zoom, inc, uv)
       
        return adj_pressure, adj_velocity, tf-t0
    
    def fd_adj_TG_solve(self, N, plot=True, scaled=False, zoom=False, inc=False, uv=False):
        solver_title = "T.G.-ELT"
        t0 = time.time()
        height = self.Example(self.args, N)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        height.h2xs = dm.center_second_diff(height.hs, height.Nx, height.dx)
        height.h3xs = dm.center_third_diff(height.hs, height.Nx, height.dx)
        
        try:
            adj_pressure = rp.TGAdj_ReynPressure(height, self.BC)
            if isinstance(self.BC, bc.Mixed):
                raise Exception(f"adj.-TG solver prescribed Q={self.BC.Q:.1f} for P_Reyn; Q for P_adj + P_reyn will differ")
        except Exception as e:
            print(e)
        tf = time.time()        
        adj_velocity = rv.TGAdj_ReynVelocity(height, self.BC, adj_pressure)

        
        if plot:
            self.p_plot(height, adj_pressure , adj_velocity.Q, solver_title, scaled, zoom)
            self.v_plot(self.BC, height, adj_velocity, adj_pressure, solver_title, scaled, zoom, inc, uv)
       
        return  adj_pressure, adj_velocity, tf-t0


    def fd_pert_solve(self, N, order,  plot=True, scaled=False, zoom=False, inc=False, uv=False, get_all = False):
        height = self.Example(self.args, N)
        t0 = time.time()
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        height.h2xs = dm.center_second_diff(height.hs, height.Nx, height.dx)
        height.h3xs = dm.center_third_diff(height.hs, height.Nx, height.dx)
        
        try:
            reyn_pressure = rp.FinDiff_ReynPressure(height, self.BC)

            if isinstance(self.BC, bc.Fixed):
                raise Exception(f"pert. solver prescribed dP={self.BC.dP:.1f} for P_0; dP for P_k>0 will differ")
        except Exception as e:
            print(e)

        
    
        reyn_velocity = rv.ReynVelocity(height, self.BC, reyn_pressure.ps_1D)                   
        
        pert = rpert.PerturbedReynSol(height, self.BC, order, reyn_pressure, reyn_velocity)
        tf = time.time()
        solver_title = ""
        
        if plot:
    
            if order > 3:
                solver_title4 = solver_title + " $\epsilon^4$ PLT"
                self.p_plot(height, pert.pert4_pressure, pert.pert4_velocity.Q, solver_title4, scaled, zoom)
                self.v_plot(self.BC, height, pert.pert4_velocity, pert.pert4_pressure, solver_title4, scaled, zoom, inc, uv)
        
            # elif order > 1:
            solver_title2 = solver_title + " $\epsilon^2$ PLT"
            self.p_plot(height, pert.pert2_pressure, pert.pert2_velocity.Q, solver_title2, scaled, zoom)
            self.v_plot(self.BC, height, pert.pert2_velocity, pert.pert2_pressure, solver_title2, scaled, zoom, inc, uv)
           
            
        if get_all:
            return pert
        else:
            if order < 3:
                return pert.pert2_pressure, pert.pert2_velocity, tf-t0
            
            else:
                return pert.pert4_pressure, pert.pert4_velocity, tf-t0

    
    def p_plot(self, height, pressure, flux, solver_title, scaled=False, zoom=False):
        if pressure.ps_2D is None:
            pressure.make_2D_ps(height)
            
        # dps = [pressure.ps_1D[i] - pressure.ps_1D[i+1] for i in range(height.Nx-1)]
        paramstr = "$Q=%.2f$, $U=%.1f$, $\Delta P=%.2f$"%(flux, self.BC.U, -pressure.dP)
        p_title = solver_title +'\n' + paramstr
        p_labels = ["$p$", "$x$","$y$"]
           
        graphics.plot_2D(pressure.ps_1D, height.xs, p_title, p_labels)
        
        if scaled:
            x_scale = height.xs[-1]-height.xs[0]
            y_scale = min(height.hs)
            p_scale = flux*x_scale/y_scale #*visc
        
            graphics.plot_contour_mesh(pressure.ps_2D/p_scale, height.xs/x_scale, height.ys/y_scale, p_title, p_labels, vmin=p_min/p_scale, vmax=p_max/p_scale, log_cmap=False)
        
        else:
         
            graphics.plot_contour_mesh(pressure.ps_2D, height.xs, height.ys, p_title, p_labels, vmin=p_min, vmax=p_max, log_cmap=False)
        
        if zoom:
            if scaled:
                xs_zoom, ys_zoom = graphics.grid_zoom_1D(height.xs, height.ys, height, x_start, x_stop, y_start, y_stop)
                p_zoom = graphics.grid_zoom_2D(pressure.ps_2D, height, x_start, x_stop, y_start, y_stop)
                graphics.plot_contour_mesh(p_zoom/p_scale, xs_zoom/x_scale, ys_zoom/y_scale, p_title, p_labels, vmin=p_min, vmax=p_max, log_cmap=False,n_contours=100)
            
            else:
                xs_zoom, ys_zoom = graphics.grid_zoom_1D(height.xs, height.ys, height, x_start, x_stop, y_start, y_stop)
                p_zoom = graphics.grid_zoom_2D(pressure.ps_2D, height, x_start, x_stop, y_start, y_stop)
                graphics.plot_contour_mesh(p_zoom, xs_zoom, ys_zoom, p_title, p_labels, vmin=p_min, vmax=p_max, log_cmap=False,n_contours=50)
            
      
    
    def v_plot(self, BC, height, velocity, pressure, solver_title, scaled=False, zoom=False,  inc=False, uv=False):
        dP = pressure.dP 
        paramstr = "$Q=%.2f$, $U=%.2f$, $\Delta   P=%.2f$"%(velocity.Q, self.BC.U, -dP)
        v_title = solver_title + '\n' + paramstr
        v_ax_labels =  ['$|(  u,  v)|_2$','$  x$', '$  y$'] 
        if scaled:
            x_scale = height.xs[-1]-height.xs[0]
            y_scale = min(height.hs)
            u_scale = velocity.Q/y_scale
            v_scale = velocity.Q/x_scale
           
            uv_mag = np.sqrt((velocity.u/u_scale)**2 + (velocity.v/v_scale)**2)
            graphics.plot_stream_heat(velocity.u/u_scale, velocity.v/y_scale, height.xs/x_scale, height.ys/y_scale, uv_mag, v_title, v_ax_labels, vmin=0, vmax=vel_max/velocity.Q, log_cmap=False)

        else:
           
            uv_mag = np.sqrt((velocity.u)**2 + (velocity.v)**2)
            graphics.plot_stream_heat(velocity.u, velocity.v, height.xs, height.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=vel_max, log_cmap=False)

        if uv:
            
            i_test = int(i_test_scale*height.Nx)
            graphics.plot_2D(height.ys,velocity.u[:,i_test],  f'$u({height.xs[i_test]:.2f}),y)$',[f'$u({height.xs[i_test]:.2f},y)$','$y$'])
            graphics.plot_2D(height.ys,velocity.v[:,i_test],  f'$v({height.xs[i_test]:.2f}),y)$',[f'$v({height.xs[i_test]:.2f},y)$','$y$'])
           
            # graphics.plot_quiver(velocity.u, velocity.v, height.xs, height.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=vel_max)

            graphics.plot_contour_mesh(velocity.u, height.xs, height.ys, 'u', ['$u$', '$x$', '$y$'], -3, 3)
            graphics.plot_contour_mesh(velocity.v, height.xs, height.ys, 'v', ['$v$', '$x$', '$y$'], -3, 3)
           
        if zoom:
            xs_zoom, ys_zoom = graphics.grid_zoom_1D(height.xs, height.ys, height, x_start, x_stop, y_start, y_stop)
            u_2D_zoom = graphics.grid_zoom_2D(velocity.u, height, x_start, x_stop, y_start, y_stop)
            v_2D_zoom = graphics.grid_zoom_2D(velocity.v, height, x_start, x_stop, y_start, y_stop)
            uv_mag_zoom = graphics.grid_zoom_2D(uv_mag, height, x_start, x_stop, y_start, y_stop)
            graphics.plot_stream_heat(u_2D_zoom, v_2D_zoom, xs_zoom, ys_zoom, uv_mag_zoom, v_title, v_ax_labels, vmin=0, vmax=vel_max, log_cmap=False)

        if inc:
            # inc = velocity.make_inc(height)
            
            # graphics.plot_contour_mesh(inc, height.xs, height.ys, 'incompressibility', ['$u_x+v_y$', '$x$', '$y$'], -1, 1)
            
            # graphics.plot_contour_mesh(uv_mag, height.xs, height.ys, v_title, v_ax_labels, vmin=0, vmax=vel_max, log_cmap=False)
            
            qs = velocity.get_flux(height)
            # qs = velocity.get_adj_flux(BC,height, pressure) #adj
            # qs = velocity.get_reyn_flux(BC, height, pressure) #reyn
            graphics.plot_2D(qs, height.xs, 'flux $\mathcal{Q} = q(x) =\int_0^{h(x)} u(x,y) dy$', ['$x$', '$q(x)=\mathcal{Q}$'])

    def sinusoid_exact_sol(self, N, plot=True, scaled=False, zoom=False, inc=False, uv=False):
        height = self.Example(self.args, N)
       
        ps = np.zeros(height.Nx)
        
        h_recip_dx = [height.h_recip_deriv_fun(x) for x in height.xs]
        for i in range(height.Nx):
            h = height.hs[i] 
            ps[i] = -6*self.BC.U * (h + height.H)/((height.k * height.H)**2 * (2+height.delta**2)) * h_recip_dx[i]
        return ps