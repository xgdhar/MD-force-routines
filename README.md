# Sarkas: pure-Python Yukawa-PPPM molecular dynamics code
# Hello!!

For the code to work you need:
1. The library 'pyfftw' for the code. It is a python wrapper to FFTW 3.3 or a higher version. So, FFTW 3.3 or higher has to be installed already. You can find details about the 'pyfftw' library here: https://github.com/hgomersall/pyFFTW

2. The library 'numba' which makes the code high performant (execution times are comparable to that of C or Fortran). 'numba' can be found here:
http://numba.pydata.org

I have set the input file  'yukawa_MD_p3m.in' to perform MD for 10000 particles with other parameter values for the case corresponding to Fig. 10 of the attached paper. Fig. 10a of the paper shows the kernel-smoothed dynamic structure factor S(q,w) which was obtained by averaging the data from 20 different MD runs (differing in their initial conditions for positions and velocities) followed by kernel smoothing of the noisy S(q,w) of MD.

Some info about the code files:

1) The input file 'yukawa_MD_p3m.in' has the following inputs:
a. 'Gamma' is the dimensionless Coulomb coupling parameter that determines the extent of coupling in the plasma
b. 'kappa' is the dimensionless inverse screening parameter that determines the extent of screening in the plasma (smaller the value, smaller is the screening)
c. 'Number of particles' 
d. 'dt' is the time step in units of inverse ion plasma frequency
e. 'Number of equilibration steps' is the number of steps for which the system is thermostatted
f. 'Number of post-equilibration steps' is the number of steps after equilibration. The data collected during this phase is used for any post-processing required to compute the observables of interest
g. 'Periodic boundary condition (1=yes, 0=no)': this is always 1. The code cannot perform non periodic boundary conditions.
h. 'snapshot interval' is the frequency with which quantities are written to output files
i. 'initial conditions (1=file, 0=random)' is used to specify whether the code has to read the initial conditions from a file or assign random initial positions and velocities
j. 'write to output (pos_vel_acc.out) file (1=yes, 0=no)' for writing to output file 'pos_vel_acc.out'
k. 'write to Ovito (pos_vel_acc.xyz) file (1=yes, 0=no)' for writing to output file 'pos_vel_acc.xyz'
l. 'seed for random number generator for initial positions and velocities' is the seed number for generating different initial conditions

2) 'Yukawa_MD_p3m_ord6.py' is the main code that imports the required functions to perform a MD run

3) 'read.py' is for reading in initial conditions from a file if it is required

4) 'initialize_pos_vel.py' is for assigning random initial particle positions and velocities using the seed number specified in the input file

4) 'velocity_verlet_p3m_ord6.py' is the velocity Verlet integrator

5)  'thermostat_p3m_ord6.py' is the thermostat

6)  'yukawa_gf_opt.py' computes the optimal Green function required by the PPPM algorithm

7)  'yukawa_p3m_ord6.py' is the code that computes the forces using the PPPM algorithm. This code imports two functions: one for the PP part (particle-particle) named 'yukawa_pp.py' and the other for the PM part (particle-mesh) named 'yukawa_pm_ord6.py'. Parameters required for PP and PM parts are set in the main file. They vary depending on the error tolerance for the PPPM algorithm.

The main MD code writes the following files as outputs:
- 'pos_vel_acc.out' has the particle positions, velocities and accelerations written to it at a frequency given by the value of 'snapshot interval' in the input file
- 'time_temp_totalEnergy_kinEnergy_potEnergy.out' has the temperature, total energy, kinetic energy and potential written to it at the same frequency.
- 'pos_vel_acc.xyz' is the xyz version of 'pos_vel_acc.out' for visualizing using Ovito software.
- 'pos_vel_acc_final.out' has the particle positions, velocities and accelerations corresponding to the last step of the post-equilibration phase
- 'n_q_t_a.npy' is the numpy array that has the spatial Fourier transform of the particle positions of every step (after equilibration) corresponding to MD run labeled as 'a' ('a' for the case of Fig. 10 runs from 1 to 20). The set of 'n_q_t_a.npy' for 'a' from 1 to 'N_avg' is processed by the the script 's_of_q_w_script.py' to generate the averaged S(q,w) (averaging over 20 runs still produces noisy curves; you would notice that the noise increases for larger 'q'). The script imports two functions 's_of_q_w.py' and 's_of_q_w_avg.py' that are required for computing S(q,w) from a single run and for averaging over many runs. Note that there is a manual element in the script: the set of 'n_q_t_a.npy' for 'a' from 1 to 'N_avg' must be manually loaded into the 'n_q_t_array'. 
