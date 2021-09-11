import numpy as np

import yukawa_pm_ord6 as y_pm_o6
import yukawa_pp as y_pp

def force_pot(kappa, G, N, rc, Lv, af, uf, pos, acc, Z,G_k,kx_v,ky_v,kz_v,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz,p,mx_max,my_max,mz_max,acc_s_r,acc_fft,rho_r,E_x_p,E_y_p,E_z_p):
    
    acc_s_r.fill(0.0)
    acc_fft.fill(0.0)
    
    rho_r.fill(0.0)
    E_x_p.fill(0.0)
    E_y_p.fill(0.0)
    E_z_p.fill(0.0)


    U_short, acc_s_r = y_pp.particle_particle(kappa,G,rc,Lv,pos,acc_s_r)
    
    U_fft, acc_fft = y_pm_o6.particle_mesh_fft_r(kappa,G,pos,Z,N,G_k,kx_v,ky_v,kz_v,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz,p,mx_max,my_max,mz_max,rho_r,acc_fft,E_x_p,E_y_p,E_z_p)

    acc = af*(acc_s_r + acc_fft)
    
    U_self = -N*G/np.sqrt(np.pi)
    U = uf*(U_short + U_fft + U_self)
    #U = uf*(U_short + U_long)
    
    print 'U_short, U_long, U_self =', [U_short, U_fft, U_self]
    
    return U, acc