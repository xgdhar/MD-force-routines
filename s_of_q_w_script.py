# Script for computing S(q,w) averaged over data from 'N_avg' number of MD runs with different
# initial conditions
import numpy as np
import matplotlib.pyplot as plt

N = 10000
dt = 0.05
Neq = 2000
Npd = 16000
q_max = 2.5 #max. value in q-space of interest

N_avg = 20 #no. of MD runs over which the averaging is done

T = Npd * dt
print "Total time (T) = ", T

L = (4*np.pi*N/3)**(1./3)
dq = 2*np.pi/L
Nq = int(q_max/dq)

dw = 2*np.pi/T
w = dw*np.arange(Npd)

# Array for loading the n_q_t data saved in each MD run
n_q_t_array = np.zeros((N_avg,Npd,3*Nq,3),dtype='complex128')

# Load the n_q_t data from each MD run using the following example
# n_q_t_array[1] = np.load('n_q_t_id.npy')

# importing the library that computes S(q,w) for one run
import s_of_q_w as sqw
# importing the library that computes S(q,w) averaged over N_avg runs
import s_of_q_w_avg as sqw_avg

s_q_w_array = np.zeros((N_avg,Nq,Npd))
s_q_w_array = sqw_avg.s_of_q_w_avg_mat(n_q_t_array,Nq,s_q_w_array,N_avg,N,dt,T)
s_q_w_avg_mat = np.mean(s_q_w_array,axis=(0)) # Averaging over S(qx,w), S(qy,w), S(qz,w)

# Plotting the first 6 spatial modes of S(q,w)
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
