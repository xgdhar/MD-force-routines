import numpy as np

import velocity_verlet_p3m_ord6 as velocity_verlet

def vscale(pos, vel, acc, T_desired, PBC, it, dt, N, d, Lv, Lmax_v, Lmin_v, kappa, G, rc, af, uf, kf,Z,G_k,kx_v,ky_v,kz_v,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz,p,mx_max,my_max,mz_max,acc_s_r,acc_fft,rho_r,E_x_p,E_y_p,E_z_p):

    pos, vel, acc, U = velocity_verlet.update(pos, vel, acc, dt, PBC, kappa, N, d, Lv, Lmax_v, Lmin_v, G, rc, af, uf, Z,G_k,kx_v,ky_v,kz_v,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz,p,mx_max,my_max,mz_max,acc_s_r,acc_fft,rho_r,E_x_p,E_y_p,E_z_p)

    K = kf*np.ndarray.sum(vel**2)
    
    T = K/kf/float(N)
    
    print dt*it, T
    
    if it <= 1999:
        
        fact = np.sqrt(T_desired/T)
        vel = vel*fact
        
    else:
        
        fact = np.sqrt((20.0*T_desired/T-1.0)/20.0)
        vel = vel*fact
        
    return pos, vel, acc, U