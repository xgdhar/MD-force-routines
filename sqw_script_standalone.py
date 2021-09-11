#This is a script to compute S(q,w) wih particle positions from MD runs of a Yukawa system
#Author: Gautham Dharuman, CMSE and ECE, Michigan State University

import numpy as np
import matplotlib.pyplot as plt

def n_of_q_w(n_q_t,N,dt,T):
    
    """function to compute S(q,w) from n(q,t) array of a single MD run"""
    
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

def s_of_q_w_avg_mat(n_q_t_array,Nq,sqw_mat,N_avg,N,dt,T):
    
    """
    function to compute S(q,w)s from multiple n(q,t) arrays of different MD runs 
    and returning the average of those S(q,w)s
    """
    
    for j in range(N_avg):
        
        n_q_t = n_q_t_array[j]
    
        for i in range(Nq):
        
            n_qi_t = n_q_t[:,i,:]
        
            s_q_w = n_of_q_w(n_qi_t,N,dt,T)
            s_q_w_avg = np.mean(s_q_w,axis=(1))
        
            sqw_mat[j,i,:] = s_q_w_avg
        
    return sqw_mat
    
if __name__ == '__main__':

    N = 10000 #number of particles
    dt = 0.05 #time step in yukawa units (inverse plasma frequency)
    Npd = 1 #no. of steps
    q_max = 2.5 #max. value in k-space (in units of inverse ion-sphere radius)
    N_avg = 1 #no. of MD runs over which the averaging is done
    T = Npd * dt
    print "Total time (T) = ", T
    
    L = (4*np.pi*N/3)**(1./3) #edge length of the cubic box in units of ion-sphere radius
    dq = 2*np.pi/L #smallest k for the box size
    Nq = int(q_max/dq) #no. of k modes
    
    dw = 2*np.pi/T #smallest temporal frequency
    w = dw*np.arange(Npd) #temporal modes
    
    #Array for loading the n_q_t data saved in each MD run
    n_q_t_array = np.zeros((N_avg,Npd,Nq,3),dtype='complex128')
    
    #Computing n(q,t) and S(q,w) for a single MD run
    n_q_t = np.zeros((Npd,Nq,3),dtype='complex128') #3 denotes 3D system
    
    #Random initial positions
    pos = L*np.random.random((N,3))
    
    #Loop over time steps
    for it in range(Npd):
        
        """
        update particle positions: pos
        """
        #Spatial Fourier transform
        for iq in range(Nq):
            q_p = (iq+1)*dq
            n_q_t[it,iq,0] = np.sum(np.exp(-1j*q_p*pos[:,0]))
            n_q_t[it,iq,1] = np.sum(np.exp(-1j*q_p*pos[:,1]))
            n_q_t[it,iq,2] = np.sum(np.exp(-1j*q_p*pos[:,2]))
            
    #np.save('n_q_t.npy',n_q_t)
    #n_q_t_array[1] = np.load('n_q_t_id.npy')
    
    n_q_t_array[0,:] = n_q_t.copy() #here averaging over only one MD run; need to loop over if averaging over many MD runs
    
    s_q_w_array = np.zeros((N_avg,Nq,Npd))
    s_q_w_array = s_of_q_w_avg_mat(n_q_t_array,Nq,s_q_w_array,N_avg,N,dt,T)
    s_q_w_avg_mat = np.mean(s_q_w_array,axis=(0)) # Averaging over S(qx,w), S(qy,w), S(qz,w)
    
    
    #Plotting the first 6 spatial modes of S(q,w)
    l0, = plt.plot(w,s_q_w_avg_mat[0,:],'g',alpha=1.0)
    l1, = plt.plot(w,s_q_w_avg_mat[1,:],'k',alpha=1.0)
    l2, = plt.plot(w,s_q_w_avg_mat[2,:],'k',alpha=0.8)
    l3, = plt.plot(w,s_q_w_avg_mat[3,:],'k',alpha=0.6)
    l4, = plt.plot(w,s_q_w_avg_mat[4,:],'k',alpha=0.4)
    l5, = plt.plot(w,s_q_w_avg_mat[5,:],'k',alpha=0.2)
    
    plt.ylim(0,3.)
    plt.xlim(0,2.5)
    plt.legend([l0,l1,l2,l3,l4,l5],['$q = \\Delta q$','$q = 2\\Delta q$','$q = 3\\Delta q$','$q = 4\\Delta q$','$q = 5\\Delta q$','$q = 6\\Delta q$',],loc='best',fontsize=15)
    plt.xlabel('$\\omega/\\omega_p$',fontsize=20)
    plt.ylabel('$S(q=ka,\\omega/\\omega_p)$',fontsize=20)
    plt.title('$\\Gamma = 2$, $\\kappa = 0.1$     ($\\Delta q = 2\\pi/L = 0.18$,   $N = 10^{4}$)',fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    plt.show()