#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<cuda.h>
#include<time.h>

__global__
void updateKernel(float kappa, float G, float rc, float pi, int Ncell, int* bin_count, float* bin_atom_ln, int bin_atom_len, float* nl_list_ln, int cell_len)
{

  int ic = blockIdx.x + blockIdx.y * gridDim.x;
  int thId = threadIdx.x + blockDim.x * blockIdx.x + blockIdx.y * gridDim.x * blockDim.x;
  int ip = thId - ic*27;
  //num_threads - num_of_atoms_in_cutoff_sphere = 10;

  int jn;
  float Zc_ip, x_ip, y_ip, z_ip; 
  float Zc_jp, x_jp, y_jp, z_jp;
  int c_ipart, c_jpart; 
  float n_x_sh, n_y_sh, n_z_sh;
  float dx, dy, dz, r, f1, f2, f3, fr;
  float G_r_k_p, G_r_k_m;
 
  float sqt_pi = 1.7724538509055159;

  if(ic < Ncell){
    if(ip < (ic*37+bin_count[ic])){

  //for(int ic = 0; ic < Ncell; ic++){
    // for(int ip = 0; ip < bin_count[ic]; ip++){
       
        //c_ipart = bin_atom_ln[ic*bin_atom_len+ip*8];
        //Zc_ip = bin_atom_ln[ic*bin_atom_len+ip*8+1];
        //x_ip = bin_atom_ln[ic*bin_atom_len+ip*8+2];
        //y_ip = bin_atom_ln[ic*bin_atom_len+ip*8+3];
        //z_ip = bin_atom_ln[ic*bin_atom_len+ip*8+4];

        c_ipart = bin_atom_ln[ip*8];
        Zc_ip = bin_atom_ln[ip*8+1];
        x_ip = bin_atom_ln[ip*8+2];
        y_ip = bin_atom_ln[ip*8+3];
        z_ip = bin_atom_ln[ip*8+4];

        for(int jc = 0; jc < 27; jc++){
           jn = (int) nl_list_ln[ic*cell_len+jc*4];
           for(int jp = 0; jp < bin_count[jn]; jp++){

              c_jpart = bin_atom_ln[jn*bin_atom_len+jp*8];
              Zc_jp = bin_atom_ln[jn*bin_atom_len+jp*8+1];
              x_jp = bin_atom_ln[jn*bin_atom_len+jp*8+2];
              y_jp = bin_atom_ln[jn*bin_atom_len+jp*8+3];
              z_jp = bin_atom_ln[jn*bin_atom_len+jp*8+4];
              
              n_x_sh = nl_list_ln[ic*cell_len+jc*4+1];         
              n_y_sh = nl_list_ln[ic*cell_len+jc*4+2];         
              n_z_sh = nl_list_ln[ic*cell_len+jc*4+3];         
              //if(ic == 0) {printf("c = %d, c_i = %d, c_j = %d\n", ic, c_ipart, c_jpart);}
              if(c_ipart != c_jpart){
              
                //if((ic == 0) && (ip == 0)) {printf("c = %d, c_i = %d, c_j = %d\n", ic, c_ipart, c_jpart);}
                dx = x_ip - (x_jp + n_x_sh);
                dy = y_ip - (y_jp + n_y_sh);
                dz = z_ip - (z_jp + n_z_sh);
                r = sqrt(dx*dx + dy*dy + dz*dz);

                if(r < rc){

                  G_r_k_p = G*r + 0.5*kappa/G;
                  G_r_k_m = G*r - 0.5*kappa/G; 
                  f1 = (0.5/(r*r)) * exp(kappa*r) * erfc(G*r + 0.5*kappa/G) * (1.0 - kappa*r);
                  f2 = (0.5/(r*r)) * exp(-kappa*r) * erfc(G*r - 0.5*kappa/G) * (1.0 + kappa*r);
                  f3 = (G/(sqt_pi*r)) * (exp(-G_r_k_p*G_r_k_p) * exp(kappa*r) + exp(-G_r_k_m*G_r_k_m) * exp(-kappa*r) );
                  fr = Zc_ip*Zc_jp*(f1 + f2 + f3);           
               
                  //bin_atom_ln[ic*bin_atom_len+ip*8+5] = bin_atom_ln[ic*bin_atom_len+ip*8+5] + (fr*dx/r);   
                  //bin_atom_ln[ic*bin_atom_len+ip*8+6] = bin_atom_ln[ic*bin_atom_len+ip*8+6] + (fr*dy/r);   
                  //bin_atom_ln[ic*bin_atom_len+ip*8+7] = bin_atom_ln[ic*bin_atom_len+ip*8+7] + (fr*dz/r);   
                  
                  bin_atom_ln[ip*8+5] = bin_atom_ln[ip*8+5] + (fr*dx/r);   
                  bin_atom_ln[ip*8+6] = bin_atom_ln[ip*8+6] + (fr*dy/r);   
                  bin_atom_ln[ip*8+7] = bin_atom_ln[ip*8+7] + (fr*dz/r);   
                  
                  //bin_atom_ln[jn*bin_atom_len+jp*8+5] = bin_atom_ln[jn*bin_atom_len+jp*8+5] - (fr*dx/r);   
                  //bin_atom_ln[jn*bin_atom_len+jp*8+6] = bin_atom_ln[jn*bin_atom_len+jp*8+6] - (fr*dy/r);   
                  //bin_atom_ln[jn*bin_atom_len+jp*8+7] = bin_atom_ln[jn*bin_atom_len+jp*8+7] - (fr*dz/r);   

                }

              }    


           }


        }

    }

  }



}


void update(float kappa, float G, float rc, float pi, int Ncell, int* bin_count, float* bin_atom_ln, int bin_atom_len, float* nl_list_ln, int cell_len)
{

 // int cN, cp, cxN, cyN, czN, cxsh, cysh, czsh;
  //float rsh_x, rsh_y, rsh_z;

  int jn;
  float Zc_ip, x_ip, y_ip, z_ip; 
  float Zc_jp, x_jp, y_jp, z_jp;
  int c_ipart, c_jpart; 
  float n_x_sh, n_y_sh, n_z_sh;
  float dx, dy, dz, r, f1, f2, f3, fr;


  for(int ic = 0; ic < Ncell; ic++){
     for(int ip = 0; ip < bin_count[ic]; ip++){
       
        c_ipart = bin_atom_ln[ic*bin_atom_len+ip*8];
        Zc_ip = bin_atom_ln[ic*bin_atom_len+ip*8+1];
        x_ip = bin_atom_ln[ic*bin_atom_len+ip*8+2];
        y_ip = bin_atom_ln[ic*bin_atom_len+ip*8+3];
        z_ip = bin_atom_ln[ic*bin_atom_len+ip*8+4];

        for(int jc = 0; jc < 27; jc++){
           jn = (int) nl_list_ln[ic*cell_len+jc*4];
           for(int jp = 0; jp < bin_count[jn]; jp++){

              c_jpart = bin_atom_ln[jn*bin_atom_len+jp*8];
              Zc_jp = bin_atom_ln[jn*bin_atom_len+jp*8+1];
              x_jp = bin_atom_ln[jn*bin_atom_len+jp*8+2];
              y_jp = bin_atom_ln[jn*bin_atom_len+jp*8+3];
              z_jp = bin_atom_ln[jn*bin_atom_len+jp*8+4];
              
              n_x_sh = nl_list_ln[ic*cell_len+jc*4+1];         
              n_y_sh = nl_list_ln[ic*cell_len+jc*4+2];         
              n_z_sh = nl_list_ln[ic*cell_len+jc*4+3];         
              //if(ic == 0) {printf("c = %d, c_i = %d, c_j = %d\n", ic, c_ipart, c_jpart);}
              if(c_ipart < c_jpart){
              
                //if(ic == 0) {printf("c = %d, c_i = %d, c_j = %d\n", ic, c_ipart, c_jpart);}
                dx = x_ip - (x_jp + n_x_sh);
                dy = y_ip - (y_jp + n_y_sh);
                dz = z_ip - (z_jp + n_z_sh);
                r = sqrt(dx*dx + dy*dy + dz*dz);

                if(r < rc){
                  f1 = (0.5/(r*r)) * exp(kappa*r) * erfc(G*r + 0.5*kappa/G) * (1.0 - kappa*r);
                  f2 = (0.5/(r*r)) * exp(-kappa*r) * erfc(G*r - 0.5*kappa/G) * (1.0 + kappa*r);
                  f3 = (G/(sqrt(pi)*r)) * (exp(-pow((G*r + 0.5*kappa/G),2)) * exp(kappa*r) + exp(-pow((G*r - 0.5*kappa/G),2)) * exp(-kappa*r) );
                  fr = Zc_ip*Zc_jp*(f1 + f2 + f3);           
               
                  bin_atom_ln[ic*bin_atom_len+ip*8+5] = bin_atom_ln[ic*bin_atom_len+ip*8+5] + (fr*dx/r);   
                  bin_atom_ln[ic*bin_atom_len+ip*8+6] = bin_atom_ln[ic*bin_atom_len+ip*8+6] + (fr*dy/r);   
                  bin_atom_ln[ic*bin_atom_len+ip*8+7] = bin_atom_ln[ic*bin_atom_len+ip*8+7] + (fr*dz/r);   
                  
                  //bin_atom_ln[ip*8+5] = bin_atom_ln[ip*8+5] + (fr*dx/r);   
                  //bin_atom_ln[ip*8+6] = bin_atom_ln[ip*8+6] + (fr*dy/r);   
                  //bin_atom_ln[ip*8+7] = bin_atom_ln[ip*8+7] + (fr*dz/r);   
                  bin_atom_ln[jn*bin_atom_len+jp*8+5] = bin_atom_ln[jn*bin_atom_len+jp*8+5] - (fr*dx/r);   
                  bin_atom_ln[jn*bin_atom_len+jp*8+6] = bin_atom_ln[jn*bin_atom_len+jp*8+6] - (fr*dy/r);   
                  bin_atom_ln[jn*bin_atom_len+jp*8+7] = bin_atom_ln[jn*bin_atom_len+jp*8+7] - (fr*dz/r);   

                }

              }    


           }


        }

     }

  }

}

int main()
{

  float kappa, G;
  float L, rc;
  float rcx, rcy, rcz;
 
  float diff;
  struct timespec start, end; 
  float diff_c;
  struct timespec start_c, end_c; 

  int ipart;
  int i, j;
  int c, cx, cy, cz;


  int N = 10000000;

  float const pi = 3.141592653589793;
  float const emp = -50;

  float **pos = (float **)malloc(N*sizeof(float *));
  for(i = 0; i < N; i++){
     pos[i] = (float *)malloc(3*sizeof(float));
  }

  float *Z = (float *)malloc(N*sizeof(float));
  for(i = 0; i < N; i++){
     Z[i] = 1.0;
  }
  
  FILE *file;
  file = fopen("pos_1e7.txt", "r");

  for(i = 0; i < N; i++){
     for(j = 0; j < 3; j++){

        if(!fscanf(file, "%f", &pos[i][j]))
          break;
     }
  }

  fclose(file);

  kappa = 0.1;
  G = 0.5;
  rc = 3.0;

  L = pow(4.0 * pi * N/3.0, 1.0/3.0);
  printf("L = %f\n", L);

  float Lx = L;
  float Ly = L;
  float Lz = L;

  int Lxd = (int) floor(Lx/rc);
  int Lyd = (int) floor(Ly/rc);
  int Lzd = (int) floor(Lz/rc);

  int Ncell = Lxd*Lyd*Lzd;

  printf("%d %d %d %d\n", Lxd, Lyd, Lzd, Ncell);

  rcx = Lx/Lxd;
  rcy = Ly/Lyd;
  rcz = Lz/Lzd;

  printf("%f %f %f\n",rcx, rcy, rcz);

  int bin_atom_c = (int) (1* pow(rc,3) + 10.0);
  printf("%d\n",bin_atom_c);
  const int atm_len = 8;
  int bin_atom_len = atm_len * bin_atom_c;
  float bin_atom[Ncell][atm_len*bin_atom_c];
  //float bin_atom_ln[Ncell*bin_atom_len];
  float *bin_atom_ln = (float *)malloc(Ncell*bin_atom_len*sizeof(float));

  for(i = 0; i < Ncell; i++){
     for(j = 0; j < bin_atom_c; j++){
        bin_atom[i][j] = emp;
        bin_atom_ln[i*bin_atom_len + j] = emp;
     }
  }

  int *bin_count = (int *)malloc(Ncell*sizeof(int));
  int bcount, bin_idx;
  for(i = 0; i < Ncell; i++){
     bin_count[i] = 0;
  }

  printf("bin_count[1] = %d\n", bin_count[1]);

  for(ipart = 0; ipart < N; ipart++){
     
     cx = (int) floor(pos[ipart][0]/rcx);
     cy = (int) floor(pos[ipart][1]/rcy);
     cz = (int) floor(pos[ipart][2]/rcz);

     c = cx + cy*Lxd + cz*Lxd*Lyd;

     //if(c == 1){printf("c = %d, ipart = %d, bin_count[%d] = %d\n", c, ipart, c, bin_count[c]);}

     bcount = bin_count[c];
     bin_idx = atm_len*bcount;
     bin_atom[c][bin_idx] = ipart; 
     //if(c == 1){printf("c = %d, ipart = %d, bin_count[%d] = %d, bin_idx = %d, bin_atom[%d][bin_idx] = %f\n", c, ipart, c, bin_count[c], bin_idx, c, bin_atom[c][bin_idx]);}
     bin_atom[c][bin_idx+1] = Z[ipart]; 
     bin_atom[c][bin_idx+2] = pos[ipart][0]; 
     bin_atom[c][bin_idx+3] = pos[ipart][1]; 
     bin_atom[c][bin_idx+4] = pos[ipart][2]; 
     bin_atom[c][bin_idx+5] = 0.0; 
     bin_atom[c][bin_idx+6] = 0.0; 
     bin_atom[c][bin_idx+7] = 0.0;
    
     bin_atom_ln[c*bin_atom_len+bin_idx] = ipart; 
     //if(c == 1){printf("c = %d, ipart = %d, bin_count[%d] = %d, bin_idx = %d, bin_atom[%d][bin_idx] = %f\n", c, ipart, c, bin_count[c], bin_idx, c, bin_atom[c][bin_idx]);}
     bin_atom_ln[c*bin_atom_len+bin_idx+1] = Z[ipart]; 
     bin_atom_ln[c*bin_atom_len+bin_idx+2] = pos[ipart][0]; 
     bin_atom_ln[c*bin_atom_len+bin_idx+3] = pos[ipart][1]; 
     bin_atom_ln[c*bin_atom_len+bin_idx+4] = pos[ipart][2]; 
     bin_atom_ln[c*bin_atom_len+bin_idx+5] = 0.0; 
     bin_atom_ln[c*bin_atom_len+bin_idx+6] = 0.0; 
     bin_atom_ln[c*bin_atom_len+bin_idx+7] = 0.0;
     
     bin_count[c] += 1; 

  }

  int chk = 0;
  printf("bin_count[%d] = %d\n", chk, bin_count[chk]);
  
  for(i = 0; i < bin_count[chk]; i++){
     //printf("c = %d, ipart = %f, Z = %f, x = %f, y = %f, z = %f\n",chk, bin_atom[chk][i*8], bin_atom[chk][i*8+1], bin_atom[chk][i*8+2], bin_atom[chk][i*8+3], bin_atom[chk][i*8+4], bin_atom[chk][i*8+5], bin_atom[chk][i*8+6], bin_atom[chk][i*8+7]);
  }
  printf("-------\n");
  for(i = 0; i < bin_count[chk]; i++){
     //printf("c = %d, ipart = %f, Z = %f, x = %f, y = %f, z = %f\n",chk, bin_atom_ln[chk*bin_atom_len + i*8], bin_atom_ln[chk*bin_atom_len + i*8+1], bin_atom_ln[chk*bin_atom_len + i*8+2], bin_atom_ln[chk*bin_atom_len+i*8+3], bin_atom_ln[chk*bin_atom_len+i*8+4], bin_atom_ln[chk*bin_atom_len+i*8+5], bin_atom_ln[chk*bin_atom_len+i*8+6], bin_atom_ln[chk*bin_atom_len+i*8+7]);
  }


  int cN, cp, cxN, cyN, czN, cxsh, cysh, czsh;
  float rsh_x, rsh_y, rsh_z;
  float nl_list[Ncell][27*4];
  //float nl_list_ln[Ncell*27*4];
  float *nl_list_ln = (float *)malloc(Ncell*27*4*sizeof(float));
  const int cell_len = 27*4;
  int n_ct;

  for(c = 0; c < Ncell; c++){
     cz = c/(Lxd*Lyd);
     cp = c % (Lxd*Lyd);
     cy = cp/Lxd;
     cx = cp % Lxd;

     n_ct = 0;
     for(czN = cz-1; czN < cz+2; czN++){
       
        if(czN < 0){
          czsh = Lzd;
          rsh_z = -Lz;
        }
        else if(czN >= Lzd){
          czsh = -Lzd;
          rsh_z = Lz;
        }
        else{
          czsh = 0;
          rsh_z = 0;
        }

        for(cyN = cy-1; cyN < cy+2; cyN++){
       
           if(cyN < 0){
             cysh = Lyd;
             rsh_y = -Ly;
           }
           else if(cyN >= Lyd){
             cysh = -Lyd;
             rsh_y = Ly;
           }
           else{
             cysh = 0;
             rsh_y = 0;
           }


           for(cxN = cx-1; cxN < cx+2; cxN++){
       
              if(cxN < 0){
                cxsh = Lxd;
                rsh_x = -Lx;
              }
              else if(cxN >= Lxd){
                cxsh = -Lxd;
                rsh_x = Lx;
              }
              else{
                cxsh = 0;
                rsh_x = 0;
              }

              cN = cxN + cxsh + (cyN + cysh)*Lxd + (czN + czsh)*Lxd*Lyd;
              nl_list[c][n_ct] = cN;
              nl_list[c][n_ct+1] = rsh_x;
              nl_list[c][n_ct+2] = rsh_y;
              nl_list[c][n_ct+3] = rsh_z;
             
              nl_list_ln[c*cell_len+n_ct] = cN;
              nl_list_ln[c*cell_len+n_ct+1] = rsh_x;
              nl_list_ln[c*cell_len+n_ct+2] = rsh_y;
              nl_list_ln[c*cell_len+n_ct+3] = rsh_z;
 
              n_ct += 4;


            }

          }

       } 

  }





  int nchk = 56;
  for(int inl = 0; inl < 27; inl++){
     //printf("c = %d, cN = %f, rxsh = %f, rysh = %f, rzsh = %f\n", nchk, nl_list[nchk][inl*4], nl_list[nchk][inl*4+1], nl_list[nchk][inl*4+2], nl_list[nchk][inl*4+3]);
     //printf("c = %d, cN = %f, rxsh = %f, rysh = %f, rzsh = %f\n", nchk, nl_list_ln[nchk*cell_len+inl*4], nl_list_ln[nchk*cell_len+inl*4+1], nl_list_ln[nchk*cell_len+inl*4+2], nl_list_ln[nchk*cell_len+inl*4+3]);
  }
  printf("-----------\n");
  for(int inl = 0; inl < 27; inl++){
     //printf("c = %d, cN = %f, rxsh = %f, rysh = %f, rzsh = %f\n", nchk, nl_list[nchk][inl*4], nl_list[nchk][inl*4+1], nl_list[nchk][inl*4+2], nl_list[nchk][inl*4+3]);
     //printf("c = %d, cN = %f, rxsh = %f, rysh = %f, rzsh = %f\n", nchk, nl_list_ln[nchk*cell_len+inl*4], nl_list_ln[nchk*cell_len+inl*4+1], nl_list_ln[nchk*cell_len+inl*4+2], nl_list_ln[nchk*cell_len+inl*4+3]);
  }

  clock_gettime(CLOCK_MONOTONIC, &start_c);

  int jn;
  float Zc_ip, x_ip, y_ip, z_ip; 
  float Zc_jp, x_jp, y_jp, z_jp;
  int c_ipart, c_jpart; 
  float n_x_sh, n_y_sh, n_z_sh;
  float dx, dy, dz, r, f1, f2, f3, fr;
  int ip_cn = 0;

  for(int ic = 0; ic < Ncell; ic++){
     for(int ip = 0; ip < bin_count[ic]; ip++){
       
        c_ipart = bin_atom[ic][ip*8];
        Zc_ip = bin_atom[ic][ip*8+1];
        x_ip = bin_atom[ic][ip*8+2];
        y_ip = bin_atom[ic][ip*8+3];
        z_ip = bin_atom[ic][ip*8+4];

        for(int jc = 0; jc < 27; jc++){
           jn = (int) nl_list[ic][jc*4];
           for(int jp = 0; jp < bin_count[jn]; jp++){

              c_jpart = bin_atom[jn][jp*8];
              Zc_jp = bin_atom[jn][jp*8+1];
              x_jp = bin_atom[jn][jp*8+2];
              y_jp = bin_atom[jn][jp*8+3];
              z_jp = bin_atom[jn][jp*8+4];
              
              n_x_sh = nl_list[ic][jc*4+1];         
              n_y_sh = nl_list[ic][jc*4+2];         
              n_z_sh = nl_list[ic][jc*4+3];         
              //if(ic == 0) {printf("c = %d, c_i = %d, c_j = %d\n", ic, c_ipart, c_jpart);}
              if(c_ipart != c_jpart){
              
                //if((ic == 0) && (ip == 0)) {printf("c = %d, c_i = %d, c_j = %d\n", ic, c_ipart, c_jpart); ip_cn += 1;}
                dx = x_ip - (x_jp + n_x_sh);
                dy = y_ip - (y_jp + n_y_sh);
                dz = z_ip - (z_jp + n_z_sh);
                r = sqrt(dx*dx + dy*dy + dz*dz);

                if(r < rc){
                  f1 = (0.5/(r*r)) * exp(kappa*r) * erfc(G*r + 0.5*kappa/G) * (1.0 - kappa*r);
                  f2 = (0.5/(r*r)) * exp(-kappa*r) * erfc(G*r - 0.5*kappa/G) * (1.0 + kappa*r);
                  f3 = (G/(sqrt(pi)*r)) * (exp(-pow((G*r + 0.5*kappa/G),2)) * exp(kappa*r) + exp(-pow((G*r - 0.5*kappa/G),2)) * exp(-kappa*r) );
                  fr = Zc_ip*Zc_jp*(f1 + f2 + f3);           
               
                  bin_atom[ic][ip*8+5] = bin_atom[ic][ip*8+5] + (fr*dx/r);   
                  bin_atom[ic][ip*8+6] = bin_atom[ic][ip*8+6] + (fr*dy/r);   
                  bin_atom[ic][ip*8+7] = bin_atom[ic][ip*8+7] + (fr*dz/r);   
                  //bin_atom[jn][jp*8+5] = bin_atom[jn][jp*8+5] - (fr*dx/r);   
                  //bin_atom[jn][jp*8+6] = bin_atom[jn][jp*8+6] - (fr*dy/r);   
                  //bin_atom[jn][jp*8+7] = bin_atom[jn][jp*8+7] - (fr*dz/r);   

                }

              }    


           }


        }

     }

  }
  clock_gettime(CLOCK_MONOTONIC, &end_c);
  diff_c = (end_c.tv_sec - start_c.tv_sec)*1000000.0 + (end_c.tv_nsec - start_c.tv_nsec)/1000.0;
  printf("elapsed time = %lf micro-seconds\n", diff_c);

  printf("no. of particles for c = 0, ip = 0: %d\n", ip_cn);

  //update(kappa, G, rc, pi, Ncell, bin_count, bin_atom_ln, bin_atom_len, nl_list_ln, cell_len);

  for(int cchk = 0; cchk < Ncell; cchk++){ 
     for(i = 0; i < bin_count[cchk]; i++){
     if(bin_atom[cchk][i*8] < 10) 
       {printf("c = %d, ipart = %f, Z = %f, x = %f, y = %f, z = %f, ax = %f, ay = %f, az = %f\n",cchk, bin_atom[cchk][i*8], bin_atom[cchk][i*8+1], bin_atom[cchk][i*8+2], bin_atom[cchk][i*8+3], bin_atom[cchk][i*8+4], bin_atom[cchk][i*8+5], bin_atom[cchk][i*8+6], bin_atom[cchk][i*8+7]);}
     }
  }
  printf("-------------------\n");
  for(int cchk = 0; cchk < Ncell; cchk++){ 
     for(i = 0; i < bin_count[cchk]; i++){
     if(bin_atom[cchk][i*8] > (N-10)) 
       {printf("c = %d, ipart = %f, Z = %f, x = %f, y = %f, z = %f, ax = %f, ay = %f, az = %f\n",cchk, bin_atom[cchk][i*8], bin_atom[cchk][i*8+1], bin_atom[cchk][i*8+2], bin_atom[cchk][i*8+3], bin_atom[cchk][i*8+4], bin_atom[cchk][i*8+5], bin_atom[cchk][i*8+6], bin_atom[cchk][i*8+7]);}
     }
  }

  printf("-------------\n");

  float *d_bin_atom_ln, *d_nl_list_ln;
  int *d_bin_count;
  clock_gettime(CLOCK_MONOTONIC, &start);
  
  cudaMalloc((void **) &d_bin_count, Ncell*sizeof(int));
  cudaMalloc((void **) &d_bin_atom_ln, Ncell*bin_atom_len*sizeof(float));
  cudaMalloc((void **) &d_nl_list_ln, Ncell*27*4*sizeof(float));

  dim3 dimGrid(512, ceil(Ncell/512.0), 1);
  dim3 dimBlock(64, 1, 1);
  //clock_gettime(CLOCK_MONOTONIC, &start);
  
  cudaMemcpy(d_bin_count, bin_count, Ncell*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_bin_atom_ln, bin_atom_ln, Ncell*bin_atom_len*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_nl_list_ln, nl_list_ln, Ncell*27*4*sizeof(float), cudaMemcpyHostToDevice);

  //clock_gettime(CLOCK_MONOTONIC, &start);
  updateKernel<<<dimGrid, dimBlock>>>(kappa, G, rc, pi, Ncell, d_bin_count, d_bin_atom_ln, bin_atom_len, d_nl_list_ln, cell_len);
  //updateKernel<<<ceil(Ncell/64), 64>>>(kappa, G, rc, pi, Ncell, d_bin_count, d_bin_atom_ln, bin_atom_len, d_nl_list_ln, cell_len);
  //updateKernel<<<ceil(Ncell/128.0), 128>>>();
  cudaDeviceSynchronize(); 
  //clock_gettime(CLOCK_MONOTONIC, &end);
  //diff = (end.tv_sec - start.tv_sec)*1000000.0 + (end.tv_nsec - start.tv_nsec)/1000.0;
  //printf("elapsed time = %lf micro-seconds\n", diff);

  cudaMemcpy(bin_atom_ln, d_bin_atom_ln, Ncell*bin_atom_len*sizeof(float), cudaMemcpyDeviceToHost);

  clock_gettime(CLOCK_MONOTONIC, &end);
  diff = (end.tv_sec - start.tv_sec)*1000000.0 + (end.tv_nsec - start.tv_nsec)/1000.0;
  printf("elapsed time = %lf micro-seconds\n", diff);

  for(int cchk = 0; cchk < Ncell; cchk++){ 
     for(i = 0; i < bin_count[cchk]; i++){
     if(bin_atom_ln[cchk*bin_atom_len+i*8] < 10) 
       {printf("c = %d, ipart = %f, Z = %f, x = %f, y = %f, z = %f, ax = %f, ay = %f, az = %f\n",cchk, bin_atom_ln[cchk*bin_atom_len+i*8], bin_atom_ln[cchk*bin_atom_len+i*8+1], bin_atom_ln[cchk*bin_atom_len+i*8+2], bin_atom_ln[cchk*bin_atom_len+i*8+3], bin_atom_ln[cchk*bin_atom_len+i*8+4], bin_atom_ln[cchk*bin_atom_len+i*8+5], bin_atom_ln[cchk*bin_atom_len+i*8+6], bin_atom_ln[cchk*bin_atom_len+i*8+7]);}
     }
  }
  
  printf("----------------\n");
 
  for(int cchk = 0; cchk < Ncell; cchk++){ 
     for(i = 0; i < bin_count[cchk]; i++){
     if(bin_atom_ln[cchk*bin_atom_len+i*8] > (N-10)) 
       {printf("c = %d, ipart = %f, Z = %f, x = %f, y = %f, z = %f, ax = %f, ay = %f, az = %f\n",cchk, bin_atom_ln[cchk*bin_atom_len+i*8], bin_atom_ln[cchk*bin_atom_len+i*8+1], bin_atom_ln[cchk*bin_atom_len+i*8+2], bin_atom_ln[cchk*bin_atom_len+i*8+3], bin_atom_ln[cchk*bin_atom_len+i*8+4], bin_atom_ln[cchk*bin_atom_len+i*8+5], bin_atom_ln[cchk*bin_atom_len+i*8+6], bin_atom_ln[cchk*bin_atom_len+i*8+7]);}
     }
  }

  cudaFree(d_bin_atom_ln); cudaFree(d_nl_list_ln); cudaFree(d_bin_count);

  free(pos); free(Z);
  free(bin_count); free(bin_atom_ln); free(nl_list_ln);
  
  return 0;
} 
