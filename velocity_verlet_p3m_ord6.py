import numpy as np

import yukawa_p3m_ord6 

def update(pos, vel, acc, dt, PBC, kappa, N, d, Lv, Lmax_v, Lmin_v, G, rc, acc_fact, U_fact,Z,G_k,kx_v,ky_v,kz_v,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz,p,mx_max,my_max,mz_max,acc_s_r,acc_fft,rho_r,E_x_p,E_y_p,E_z_p):
    
    vel = vel + 0.5*acc*dt
    
    pos = pos + vel*dt
    
    # periodic boundary condition
    if PBC == 1:
        for i in np.arange(N):
            for p in np.arange(d):
        
                if pos[i,p] > Lmax_v[p]:
                    pos[i,p] = pos[i,p] - Lv[p]
                if pos[i,p] < Lmin_v[p]:
                    pos[i,p] = pos[i,p] + Lv[p]
        
    
    U, acc = yukawa_p3m_ord6.force_pot(kappa, G, N, rc, Lv, acc_fact, U_fact, pos, acc, Z,G_k,kx_v,ky_v,kz_v,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz,p,mx_max,my_max,mz_max,acc_s_r,acc_fft,rho_r,E_x_p,E_y_p,E_z_p)
    
    vel = vel + 0.5*acc*dt
    
    return pos, vel, acc, U