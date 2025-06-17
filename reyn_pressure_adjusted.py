# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:16:58 2025

@author: sarah
"""
import numpy as np
import domain as dm
import graphics

# def make_adj_ps(height, reyn_ps):
#     ps_adj = np.zeros((height.Ny, height.Nx))

#     hs = height.hs
#     hxs = height.hxs
#     dx = height.dx
    
#     #--------------------------------------------------------------------------
#     dP_Es = np.zeros(height.Nx)
#     dP_Ws = np.zeros(height.Nx)
#     dPs = np.zeros(height.Nx)
        
#     i=0
#     dP_Es[i] = dm.center_first(dx, reyn_ps[i : i+3])
#     dPs[i] =  dm.right_first(dx, reyn_ps[i : i+3])
#     dP_Ws[i] = dPs[i]
    
#     i=1
#     dP_Es[i] = dm.center_first(dx, reyn_ps[i : i+3])
#     dPs[i] = dm.center_first(dx, reyn_ps[i-1 : i+2])
#     dP_Ws[i] = dm.right_first(dx, reyn_ps[i-1 : i+2])
    
#     i=height.Nx-1
#     dP_Ws[i] = dm.center_first(dx, reyn_ps[i-2 : i+1])
#     dPs[i] = dm.left_first(dx, reyn_ps[i-2 : i+1])
#     dP_Es[i] = dPs[i]
    
#     i=height.Nx-2
#     dP_Ws[i] = dm.center_first(dx, reyn_ps[i-2 : i+1]) 
#     dPs[i] = dm.center_first(dx, reyn_ps[i-1 : i+2]) 
#     dP_Es[i] = dm.left_first(dx, reyn_ps[i-1 : i+2]) 
    
#     for i in range(2,height.Nx-2): 
#         dP_Es[i] = dm.center_first(dx, reyn_ps[i : i+3])
#         dPs[i] = dP_Es[i-1]
#         dP_Ws[i] = dP_Es[i-2] #dm.center_first(dx, reyn_ps[i-2 : i+1])
    
#     for i in height.i_peaks[1:-1]:
#         dP_Es[i-1] = dm.left_first(dx, reyn_ps[i-3 : i])
#         dP_Ws[i+1] = dm.right_first(dx, reyn_ps[i+1 : i+4])
#    #---------------------------------------------------------------------------
#     for j in range(height.Ny):

#         y = height.ys[j]
        
#         i = 0 #inlet bc
#         if y <= height.hs[i]:
#             adj_E = (hs[i+1]-y)*dP_Es[i]/2 + height.U*height.visc/hs[i+1]
#             adj_W = (hs[i]-y)*dP_Ws[i]/2 + height.U*height.visc/hs[i]
#             adj = y*(adj_E - adj_W)/(height.dx) #- ((hs[i]/2)*dPs[i] - height.U*height.visc/hs[i])*hxs[i]

#             ps_adj[j,i] = reyn_ps[i] + adj#*(hs[i]-y)/hs[i]
#         else:
#             ps_adj[j,i]=None
        
            
#         i = height.Nx-1 # outlet bc
#         if y <= height.hs[i]:
#             adj_E = (hs[i]-y)*dP_Es[i]/2 + height.U*height.visc/hs[i]
#             adj_W = (hs[i-1]-y)*dP_Ws[i]/2 + height.U*height.visc/hs[i-1]
#             adj = y* (adj_E - adj_W)/(height.dx) #- ((hs[i]/2)*dPs[i] - height.U*height.visc/hs[i])*hxs[i]
#             ps_adj[j,i] = reyn_ps[i] + adj#*(hs[i]-y)/hs[i]
#         else:
#             ps_adj[j,i]=None 
            
    
#         for i in range(1,height.Nx-1): 
                    
#             if y > hs[i]:
#                 ps_adj[j,i] = None
                
#             else:  
                
#                 #find dP_E = dp/px @ x_i+1 and dP_W = dp/dx # x_i-1
#                 h_E = height.hs[i+1]
#                 h=height.hs[i]
#                 h_W = height.hs[i-1]
                
#                 dP_E = dP_Es[i]
#                 dP=dPs[i]
#                 dP_W = dP_Ws[i]
                
#                 adj_E = (h_E-y) * dP_E/2 - height.U*height.visc/h_E
               
#                 adj_W = (h_W-y) * dP_W/2 - height.U*height.visc/h_W

#                 adj = y*(adj_E - adj_W)/(2*height.dx) 

#                 ps_adj[j,i] = reyn_ps[i] + adj
                
                
                
#                 if y + height.dy > hs[i]:
                    
#                     print(adj - ((h/2)*dP - height.U*height.visc/h)*hxs[i])


#     for i in height.i_peaks[1:-1]:
#         h=height.hs[i]

#         for j in range(height.Ny):
            
#             y=height.ys[j]

#             if y<= h and y <=height.hs[i-3] and y<= height.hs[i+3]:
                
#                 ps_adj[j,i-2:i+3] = dm.avg_3x(ps_adj[j,i-3:i+4])

#             elif y > h and y <= height.hs[i+3]:
#                 ps_adj[j,i+2] = ps_adj[j,i+3] 
#                 ps_adj[j,i+1] = ps_adj[j,i+2] 
#                 ps_adj[j,i] = ps_adj[j,i+1] 
                
#             elif y > h and y<= height.hs[i-3]:
#                 ps_adj[j,i-2] = ps_adj[j,i-3] 
#                 ps_adj[j,i-1] = ps_adj[j,i-2] 
#                 ps_adj[j,i] = ps_adj[j,i-1] 
               

#     # graphics.plot_2D(ps_adj[60,:], height.xs, 'p(x,y0)', ['x','p'])
    
#     return ps_adj 


def make_adj_ps(height, reyn_ps):
    ps_adj = np.zeros((height.Ny, height.Nx))

    errs = np.zeros(height.Nx)
    
    hs = height.hs
    hxs = height.hxs
    dx = height.dx
    
    #--------------------------------------------------------------------------
    pxs = np.zeros(height.Nx)
    pxxs = np.zeros(height.Nx)
        
    i=0
    
    pxs[i] =  dm.right_first(dx, reyn_ps[i : i+3])
    pxxs[i] = dm.right_second(dx, reyn_ps[i : i+4])
    
 
    i=height.Nx-1
    pxs[i] = dm.left_first(dx, reyn_ps[i-2 : i+1])
    pxxs[i] = dm.left_second(dx, reyn_ps[i-3 : i+1])
    
    
    for i in range(1,height.Nx-1): 
        pxs[i] = dm.center_first(dx, reyn_ps[i-1 : i+2])
        pxxs[i] = dm.center_second(dx, reyn_ps[i-1 : i+2])
       
    for i in height.i_peaks[1:-1]:
        pxs[i-1:i+2] = dm.avg_2x(pxs[i-2 : i+3])
        pxxs[i-1:i+2] = dm.avg_2x(pxxs[i-2 : i+3])
   #---------------------------------------------------------------------------
   
    for i in range(height.Nx): 
        h = hs[i]
        hx = hxs[i]
        px = pxs[i]
        pxx = pxxs[i]
        
        f1 = ((pxx*h + px*hx)/2 - height.U*height.visc/(h**2)*hx)
        f2 = ((h/2)*px - height.U*height.visc/h)*hx
        
        for j in range(height.Ny):
            y = height.ys[j]
            if y > h:
                ps_adj[j,i] = None
                continue
            else:  
                adj = -pxx*(y**2)/2 + f1*y  -f2

                ps_adj[j,i] = reyn_ps[i] + adj
                            
                if y + height.dy > h:
                    
                    errs[i] = abs(ps_adj[j,i]-reyn_ps[i])
    print('max err p(x,h)=p_re(x): ', np.max(errs))
    return ps_adj 