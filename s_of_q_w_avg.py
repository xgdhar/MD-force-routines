
import numpy as np
import s_of_q_w as sqw

def s_of_q_w_avg_mat(n_q_t_array,Nq,sqw_mat,N_avg,N,dt,T):
    
    for j in range(N_avg):
        
        n_q_t = n_q_t_array[j]
    
        for i in range(Nq):
        
            n_qi_t = n_q_t[:,i*3,:]
        
            s_q_w = sqw.n_of_q_w(n_qi_t,N,dt,T)
            s_q_w_avg = np.mean(s_q_w,axis=(1))
        
            sqw_mat[j,i,:] = s_q_w_avg
        
    return sqw_mat