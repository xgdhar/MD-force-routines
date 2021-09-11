# Following is the main code that executes molecular dynamics simulation for a Yukawa plasma 
# using the efficient Particle-Particle-Particle-Mesh algorithm for force computation.
# The code constitutes a number of functions that are in separate files

# Code developer: Gautham Dharuman
# Dept. of Computational Mathematics, Science, and Engineering,
# Dept. of Electrical and Computer Engineering, 
# Michigan State University

# Importing Numpy module
import numpy as np
import time

t1 = time.time()

# Importing MD modules
import read
import initialize_pos_vel
import velocity_verlet_p3m_ord6 as velocity_verlet
import thermostat_p3m_ord6 as thermostat
import yukawa_gf_opt
import yukawa_p3m_ord6 as yukawa_p3m


# Reading MD conditions from input file
c = np.loadtxt('yukawa_MD_p3m.in', delimiter=' ', usecols=(0,), unpack=True)
Gamma = c[0]
kappa = c[1]
N = int(c[2])
dt = c[3]
Neq = int(c[4])
Nt = int(c[5])
PBC = int(c[6])
snap_int = int(c[7])
init = c[8]
write_output = c[9]
write_xyz = c[10]
seed_int = int(c[11])

# Other MD parameters
T_desired = 1.0/Gamma                 # desired temperature
L = (4.0*np.pi*N/3.0)**(1.0/3.0)      # box length
Lx = L
Ly = L
Lz = L
Lv = np.array([L, L, L])              # box length vector
d = np.count_nonzero(Lv)              # no. of dimensions
Lmax_v = np.array([L, L, L]) 
Lmin_v = np.array([0.0, 0.0, 0.0])

#Ewald parameters
G = 0.46
G_ew = G
rc = 7.8
#nmax = 7
#nmax_v = np.array([nmax, nmax, nmax])

#P3M parameters
Mx = 64
My = 64
Mz = 64
hx = Lx/Mx
hy = Ly/My
hz = Lz/Mz
p = 6
mx_max = 3
my_max = 3
mz_max = 3

t2 = time.time()

G_k, kx_v, ky_v, kz_v, A_pm = yukawa_gf_opt.gf_opt(kappa,G_ew,p,mx_max,my_max,mz_max,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz)

t3 = time.time()

# pre-factors as a result of using 'reduced' units
af = 1.0/3.0                          # acceleration factor
uf = 1.0                              # potential energy factor
kf = 1.5                              # kinetic energy factor

print '\n\n----------- Molecular Dynamics Simulation of Yukawa System ----------------------'
print 'Gamma = ', Gamma
print 'kappa = ', kappa
print 'Temperature = ', T_desired
print 'No. of particles = ', N
print 'Box length along x axis = ', Lv[0]
print 'Box length along y axis = ', Lv[1]
print 'Box length along z axis = ', Lv[2]
print 'No. of non-zero box dimensions = ', d
print 'time step = ',dt
print 'No. of equilibration steps = ', Neq
print 'No. of post-equilibration steps = ', Nt
print 'snapshot interval = ', snap_int
print 'Periodic boundary condition{1=yes, 0=no} =', PBC
print 'grid_size * Ewald_parameter (h * alpha) = ', hx*G_ew

# Diagnostic variables
#Tv_eq = np.zeros(Neq, dtype=float) # temperature during equilibration
#tv = dt*np.arange(Nt, dtype=float)  # time vector for production phase
#Ev = np.zeros(Nt, dtype=float)     # total energy during production
#Kv = np.zeros_like(Ev)              # total kinetic energy during production
#Uv = np.zeros_like(Ev)              # total potential energy during production
#Tv_pd = np.zeros_like(Ev)           # temperature during production
#rel_err_Ev = np.zeros_like(Ev)      # relative energy error during production

# Particle positions and velocities array
pos = np.zeros((N,d))
vel = np.zeros_like(pos)
acc = np.zeros_like(pos)
Z = np.ones(N)

acc_s_r = np.zeros_like(pos)
acc_fft = np.zeros_like(pos)

rho_r = np.zeros((Mz,My,Mx))
E_x_p = np.zeros(N)
E_y_p = np.zeros(N)
E_z_p = np.zeros(N)


# F(k,t): Spatial Fourier transform of density fluctutations
dq = 2.*np.pi/L
print 'smallest interval in Fourier space for S(q,w): dq = ', dq
q_max = 2.5
Nq = 3*int(q_max/dq)

n_q_t = np.zeros((Nt,Nq,3),dtype='complex128') #

# initializing the q vector
qv = np.zeros(Nq)

for iqv in range(0,Nq,3):

    iq = iqv/3.
    qv[iqv] = (iq+1.)*dq
    qv[iqv+1] = (iq+1.)*np.sqrt(2.)*dq
    qv[iqv+2] = (iq+1.)*np.sqrt(3.)*dq

#array for temperature, total energy, kinetic energy, potential energy
t_Tp_E_K_U2 = np.zeros((1,5))

# Initializing particle positions and velocities
if init == 1:
    
    print '\nReading initial particle positions and velocities from file...'
    
    f_input = 'init.out'           # name of input file
    pos, vel = read.initL(pos, vel, f_input)
    
else:
    
    print '\nAssigning random initial positions and velocities...'
    
    # initial particle positions uniformly distributed in the box
    # initial particle velocities with Maxwell-Boltzmann distribution
    pos, vel = initialize_pos_vel.initial(seed_int,N,Lv,pos,vel,T_desired) 
    
    #pos = initialize.positions(N,Lv,pos)          #initial particle positions uniformly distributed in the box
    #vel = initialize.velocities(N,vel,T_desired)  #initial particle velocities with Maxwell-Boltzmann distribution

    
# Visualizing initial particle positions and velocities
#visualize.initial_pos_vel(pos,vel,T_desired,Gamma,kappa,N)

# Calculating initial forces and potential energy
#U, acc = force.minimum_image(kappa, N, d, Lv, af, uf, pos, acc)
#U, acc = yukawa_ewald.force_pot(kappa, G, rc, nmax_v, Lv, af, uf, pos, acc)

t4 = time.time()

U, acc = yukawa_p3m.force_pot(kappa, G, N, rc, Lv, af, uf, pos, acc, Z,G_k,kx_v,ky_v,kz_v,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz,p,mx_max,my_max,mz_max,acc_s_r,acc_fft,rho_r,E_x_p,E_y_p,E_z_p)

#---------------

print '\n------------- Equilibration -------------'
print 'time - temperature'

for it in range(Neq):
    
    #pos, vel, acc, U = thermostat.vscale(pos, vel, acc, T_desired, PBC, it, dt, kappa, N, d, Lv, Lmax_v, Lmin_v, G, rc, nmax_v, af, uf, kf)
    pos, vel, acc, U = thermostat.vscale(pos, vel, acc, T_desired, PBC, it, dt, N, d, Lv, Lmax_v, Lmin_v, kappa, G, rc, af, uf, kf,Z,G_k,kx_v,ky_v,kz_v,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz,p,mx_max,my_max,mz_max,acc_s_r,acc_fft,rho_r,E_x_p,E_y_p,E_z_p)
#---------------

t5 = time.time()

print '\n------------- Production -------------'

# Opening files for writing particle positions, velcoities and forces
f_output = open('pos_vel_acc.out','w')
f_output_E = open('time_temp_totalEnergy_kinEnergy_potEnergy.out','w')
f_xyz = open('pos_vel_acc.xyz','w')

print 'time - total energy - kinetic energy - potential energy'

for it in range(Nt):
    
    pos, vel, acc, U = velocity_verlet.update(pos, vel, acc, dt, PBC, kappa, N, d, Lv, Lmax_v, Lmin_v, G, rc, af, uf, Z,G_k,kx_v,ky_v,kz_v,Mx,My,Mz,hx,hy,hz,Lx,Ly,Lz,p,mx_max,my_max,mz_max,acc_s_r,acc_fft,rho_r,E_x_p,E_y_p,E_z_p)

    K = kf*np.ndarray.sum(vel**2)
    E = K + U
    Tp = K/kf/float(N)

    print dt*it, E, K, U
    
    #Ev[it] = E
    #Kv[it] = K
    #Uv[it] = U
    #Tv_pd[it] = K/kf/float(N)
    #rel_err_Ev[it] = (E - Ev[0])/Ev[0]
    
    t_Tp_E_K_U = np.array([dt*it, Tp, E, K, U])
    t_Tp_E_K_U2[:] = t_Tp_E_K_U
    
    # Spatial Fourier transform
    for iqv in range(Nq):
        q_p = qv[iqv]
        n_q_t[it,iqv,0] = np.sum(np.exp(-1j*q_p*pos[:,0]))
        n_q_t[it,iqv,1] = np.sum(np.exp(-1j*q_p*pos[:,1]))
        n_q_t[it,iqv,2] = np.sum(np.exp(-1j*q_p*pos[:,2]))
    
    # writing particle positions and velocities to file
    if write_output == 1:
        if np.mod(it+1,snap_int) == 0:
            irp = np.hstack((pos,vel,acc))
            np.savetxt(f_output,irp)
            np.savetxt(f_output_E,t_Tp_E_K_U2)
            
            if write_xyz == 1:
                f_xyz.writelines('{0:d}\n'.format(N))
                f_xyz.writelines('x y z vx vy vz ax ay az\n')
                np.savetxt(f_xyz,irp)

np.save('n_q_t_1',n_q_t)

# closing output files        
f_output.close()
f_output_E.close()
f_xyz.close()

# saving last positions, velocities and accelerations
irp2 = np.hstack((pos,vel,acc))
np.savetxt('pos_vel_acc_final.out',irp2)

t6 = time.time()

print 'Time for importing required libraries = ', t2-t1
print 'Time for computing converged Greens function = ', t3-t2
print 'Time for initialization = ', t4-t3
print 'Time for equilibration = ', t5-t4
print 'Time for production = ', t6-t5
print 'Total elapsed time = ', t6-t1

#---------------

# Visualizing energy components from production phase
#visualize.energy_components(tv,Ev,Kv,Uv,rel_err_Ev,Gamma,kappa,N,dt)