// analysis sample of intensity distribution.
#include "bem3_aw_b2.h"

int main(int argc,char *argv[]) 
{
  MOBJ mo;
  FILE *fp1;
  double complex p,v[3];
  double rang,dr,r[3],*ip;
  int max,i,j,type;

  mo_dat_read(argv[1],&mo); // read datafile 
  mo_print_data(&mo);       // print data 

  max=200;
  rang=3.0*mo.md[0].aw.lambda0;
  dr=rang*2.0/(double)(max-1);
  type=1; // type setting for total_field_amsp()
  
  ip=(double *)m_alloc2(max,sizeof(double),"exampl2.c,ip");

  // x=0 plane
  if((fp1=fopen("Ip_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# y z sound_pressure_intensity");
  r[0]=0.0;
  for(i=0;i<max;i++){
    r[1]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(p,v) // omp parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      mo_pv_t(&p,v,r,type,&mo); // total field
      ip[j]=creal(p*conj(p)); // |p|^2
    } // end parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[1],r[2],ip[j]);
    }
    fprintf(fp1,"\n");
  }
  fclose(fp1);

  // y=0 plane
  if((fp1=fopen("Ip_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# x z sound_pressure_intensity");
  r[1]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(p,v) // omp parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      mo_pv_t(&p,v,r,type,&mo); // total field
      ip[j]=creal(p*conj(p));
    }// end parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[0],r[2],ip[j]);
    }
    fprintf(fp1,"\n");
  }
  fclose(fp1);

  // z=0 plane  
  if((fp1=fopen("Ip_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  } 
  fprintf(fp1,"%s\n","# x y sound_pressure_intensity");
  r[2]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(p,v) // omp parallel
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      mo_pv_t(&p,v,r,type,&mo); // total field
      ip[j]=creal(p*conj(p));
    }
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[0],r[1],ip[j]);
    }
    fprintf(fp1,"\n");
  }
  fclose(fp1);
  
  printf("Intensity plot is finished\n");
  
  free(ip);
  mo_finalize(&mo);
  return 0;
}
