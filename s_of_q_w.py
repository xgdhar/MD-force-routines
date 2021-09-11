import numpy as np

def n_of_q_w(n_q_t,N,dt,T):
    
    n_qx_w = (np.fft.fft(n_q_t[:,0]))*dt
    n_qy_w = (np.fft.fft(n_q_t[:,1]))*dt
    n_qz_w = (np.fft.fft(n_q_t[:,2]))*dt

    S_qx_w = (np.absolute(n_qx_w))**2/(N*T)
    S_qy_w = (np.absolute(n_qy_w))**2/(N*T)
    S_qz_w = (np.absolute(n_qz_w))**2/(N*T)
    
    n_w = np.size(S_qx_w)
    S_q_x_y_z_w = np.zeros((n_w,3))
    
    S_q_x_y_z_w[:,0] = S_qx_w
    S_q_x_y_z_w[:,1] = S_qy_w
    S_q_x_y_z_w[:,2] = S_qz_w
    
    return S_q_x_y_z_w