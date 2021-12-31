// calculation example of far-field intensity distributions.
#include "bem3_aw_b2.h"

int main(int argc,char *argv[]) 
{
  MOBJ mo;
  FILE *fp1;
  double complex p,v[3];
  double th,ph,phd,dthd,dthr,dphr,dphd,ra,r[3],*ip,ipmax,mf;
  int i,j,sn,type;
  
  if(argc!=2 && argc!=5){
    printf("Usage : %s datafile_name [sampling_number multplier_factor type](optional)\n",argv[0]);
    printf("default sampling number 360, multiplier factor 2000 (radius = 2000*lambda0), type 0 (4 point GaussLegendre)\n");
    exit(0);
  }
  else if(argc==5){
    sn=atoi(argv[2]);
    mf=atof(argv[3]);
    type=atoi(argv[4]);
  }
  else{
    sn=360;
    mf=2000.0;
    type=0;
  }

  mo_dat_read(argv[1],&mo); // read datafile 
  mo_print_data(&mo);       // print data 

  ra=mf*mo.md[0].aw.lambda0; // radius for calculation point
  dthd=360.0/(double)sn;     // delta theta [degree]
  dthr=2.0*M_PI/(double)sn;  // delta theta [radian]
  dphd=180.0/(double)sn;     
  dphr=1.0*M_PI/(double)sn;
  
  ip=(double *)m_alloc2(sn+1,sizeof(double),"example4.c,ie");
  
  // x=0 plane, th=0 : +z-axis, th=270 : +y-axis
  if((fp1=fopen("fsIp_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","## x=0 plane, theta=0 : +z-axis, theta=270 : +y-axis ");
  fprintf(fp1,"%s %d, %s %g\n","## sampling number",sn,"multiplier factor",mf);
  fprintf(fp1,"%s\n","# theta sound_pressure_intensity normalized_intensity");
  ipmax=0.0;
  for(i=0;i<sn;i++){
    th=0.5*dthr+(double)i*dthr;
    r[0]=0.0;
    r[1]=-ra*sin(th);
    r[2]= ra*cos(th);
    mo_pv_s_dbieq(&p,v,r,type,&mo); // scattered field
    ip[i]=creal(p*conj(p));
    if(ip[i]>ipmax) ipmax=ip[i];
  }
  for(i=0;i<sn;i++){
    th=0.5*dthd+(double)i*dthd;
    fprintf(fp1,"%g %15.14e %15.14e\n",th,ip[i],ip[i]/ipmax);
  }
  fclose(fp1);

  // y=0 plane, th=0 : +z-axis, th=90 : +x-axis
  if((fp1=fopen("fsIp_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","## y=0 plane, theta=0 : +z-axis, theta=90 : +x-axis ");
  fprintf(fp1,"%s %d, %s %g\n","## sampling number",sn,"multiplier factor",mf);
  fprintf(fp1,"%s\n","# theta sound_pressure_intensity normalized_intensity");
  ipmax=0.0;
  for(i=0;i<sn;i++){
    th=0.5*dthr+(double)i*dthr;
    r[0]=ra*sin(th);
    r[1]=0.0;
    r[2]=ra*cos(th);
    mo_pv_s_dbieq(&p,v,r,type,&mo); // scattered field
    ip[i]=creal(p*conj(p));
    if(ip[i]>ipmax) ipmax=ip[i];
  }
  for(i=0;i<sn;i++){
    th=0.5*dthd+(double)i*dthd;
    fprintf(fp1,"%g %15.14e %15.14e\n",th,ip[i],ip[i]/ipmax);
  }
  fclose(fp1);
  
  // z=0 plane, th=0 : +x-axis, th=90 : +y-axis
  if((fp1=fopen("fsIp_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","## z=0 plane, theta=0 : +x-axis, theta=90 : +y-axis ");
  fprintf(fp1,"%s %d, %s %g\n","## sampling number",sn,"multiplier factor",mf);
  fprintf(fp1,"%s\n","# theta sound_pressure_intensity normalized_intensity");
  ipmax=0.0;
  for(i=0;i<sn;i++){
    th=0.5*dthr+(double)i*dthr;
    r[0]=ra*cos(th);
    r[1]=ra*sin(th);
    r[2]=0.0;
    mo_pv_s_dbieq(&p,v,r,type,&mo); // scattered field
    ip[i]=creal(p*conj(p));
    if(ip[i]>ipmax) ipmax=ip[i];
  }
  for(i=0;i<sn;i++){
    th=0.5*dthd+(double)i*dthd;
    fprintf(fp1,"%g %15.14e %15.14e\n",th,ip[i],ip[i]/ipmax);
  }
  fclose(fp1);

  // 3d plot 
  if((fp1=fopen("fsIp_3d.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","## 3d plot, x=r*sin(theta)*cos(phi), y=r*sin(theta)*sin(phi), z=r*cos(theta), r=multiplier_factor*lambda0");
  fprintf(fp1,"%s %d, %s %g\n","## sampling number",sn,"multiplier factor",mf);
  fprintf(fp1,"%s\n","# theta phi sound_pressure_intensity");
  for(i=0;i<sn;i++){
    ph =0.5*dphr+(double)i*dphr;
    phd=0.5*dphd+(double)i*dphd;
    #pragma omp parallel for schedule(dynamic) private(th,r,p,v)
    for(j=0;j<=sn;j++){
      th=(double)j*dthr;
      r[0]=ra*sin(ph)*cos(th);
      r[1]=ra*sin(ph)*sin(th);
      r[2]=ra*cos(ph);
      mo_pv_s_dbieq(&p,v,r,type,&mo); // scattered field
      ip[j]=creal(p*conj(p));
    }
    for(j=0;j<=sn;j++){
      th=(double)j*dthd;
      fprintf(fp1,"%g %g %15.14e\n",phd,th,ip[j]);
    }
    fprintf(fp1,"\n");
  }
  fclose(fp1);

  free(ip);
  mo_finalize(&mo);
  return 0;
}
