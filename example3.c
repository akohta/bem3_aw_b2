// analysis sample of instantaneous value of each field.
#include "bem3_aw_b2.h"
#include <sys/stat.h>
#include <errno.h>  
#include <png.h>

typedef struct image_data{
  char dir_name[64];      // directory name to output image
  int scale;              // number for enlarge the output image
  
  int m;                  // sampling number 
  double rang;            // range of sampling

  int ts;                 // time step per cycle
  
  int type;               // type setting for surface integral
  
  double complex *vp;     // sound pressure data
  double mp;              // maximum amplitude 
}IMD;

void directory_name(char *src,char *nn);
void make_directory(char *dir_name);

void p_field_x(IMD *id,MOBJ *ad);
void p_field_y(IMD *id,MOBJ *ad);
void p_field_z(IMD *id,MOBJ *ad);
void output_field(char *pl,IMD *id,MOBJ *ad);

// color table
png_byte ct1[9][3]={{0x00,0x00,0x90},{0x00,0x0f,0xff},{0x00,0x90,0xff},{0x0f,0xff,0xee},
                    {0xff,0xff,0xff},{0xff,0xee,0x00},{0xff,0x70,0x00},{0xee,0x00,0x00},{0x7f,0x00,0x00}};
/*                    
png_byte ct1[9][3]={{0x00,0x00,0x90},{0x00,0x0f,0xff},{0x00,0x90,0xff},{0x0f,0xff,0xee},
                    {0x90,0xff,0x70},{0xff,0xee,0x00},{0xff,0x70,0x00},{0xee,0x00,0x00},{0x7f,0x00,0x00}};  
*/

int main(int argc,char *argv[])
{
  MOBJ ad;
  IMD id;

  mo_dat_read(argv[1],&ad); // read datafile 
  mo_print_data(&ad);       // print data 
  
  directory_name(argv[1],id.dir_name); // remove file-extension from argv[1] and add "_images"
  id.scale=1;                          // number for enlarge the output image
  id.m=180;                            // sampling number 
  id.rang=3.0*ad.md[0].aw.lambda0;     // range of sampling
  id.ts=40;                            // time step per cycle
  id.type=0;                           // type setting for total_field_amsp()
  
  make_directory(id.dir_name);

  id.vp=(double complex *)m_alloc2(id.m*id.m,sizeof(double complex),"example3.c, vp");

  // x=0 plane
  p_field_x(&id,&ad);
  output_field("yz",&id,&ad);
  // y=0 plane
  p_field_y(&id,&ad);
  output_field("xz",&id,&ad);
  // z=0 plane
  p_field_z(&id,&ad);
  output_field("xy",&id,&ad);

  free(id.vp);
  
  mo_finalize(&ad);
  return 0;
}

void directory_name(char *src,char *nn)
{
  int s1,s2;
  char *sd,buf[64]="";
  
  sd=strrchr(src,'.');
  if(sd==NULL){ // no file extension
    sprintf(nn,"%s_images",src);
  }
  else {
    s1=strlen(src);
    s2=strlen(sd);
    strncpy(buf,src,s1-s2);
    sprintf(nn,"%s_images",buf);
  }
  
}

void make_directory(char *dir_name)
{
  int ret;
  
  ret=mkdir(dir_name,S_IRWXU|S_IRWXG);
  if(ret!=0 && errno!=EEXIST){
    printf("failed to make directory. Exit..");
    exit(1);
  }
}

void p_field_x(IMD *id,MOBJ *ad)
{
  double complex p,v[3];
  double x[3],dr;
  int i,j;
  
  dr=id->rang*2.0/(double)(id->m-1);
  
  id->mp=0.0;
  
  // x=0 plane  
  x[0]=0.0;
  #pragma omp parallel for schedule(dynamic) firstprivate(x) private(j,p,v) 
  for(i=0;i<id->m;i++){
    x[2]=id->rang-(double)i*dr;
    for(j=0;j<id->m;j++){
      x[1]=-id->rang+(double)j*dr;
      mo_pv_t(&p,v,x,id->type,ad); // total field
      
      #pragma omp critical
      if(cabs(p)>id->mp) id->mp=cabs(p);
      
      id->vp[i*id->m+j]=p;
    }
  }
}

void p_field_y(IMD *id,MOBJ *ad)
{
  double complex p,v[3];
  double x[3],dr;
  int i,j;
  
  dr=id->rang*2.0/(double)(id->m-1);
  
  id->mp=0.0;
  
  // y=0 plane  
  x[1]=0.0;
  #pragma omp parallel for schedule(dynamic) firstprivate(x) private(j,p,v) 
  for(i=0;i<id->m;i++){
    x[2]=id->rang-(double)i*dr;
    for(j=0;j<id->m;j++){
      x[0]=-id->rang+(double)j*dr;
      mo_pv_t(&p,v,x,id->type,ad); // total field
      
      #pragma omp critical
      if(cabs(p)>id->mp) id->mp=cabs(p);
      
      id->vp[i*id->m+j]=p;
    }
  }
}

void p_field_z(IMD *id,MOBJ *ad)
{
  double complex p,v[3];
  double x[3],dr;
  int i,j;
  
  dr=id->rang*2.0/(double)(id->m-1);
  
  id->mp=0.0;
  
  // z=0 plane  
  x[2]=0.0;
  #pragma omp parallel for schedule(dynamic) firstprivate(x) private(j,p,v) 
  for(i=0;i<id->m;i++){
    x[1]=id->rang-(double)i*dr;
    for(j=0;j<id->m;j++){
      x[0]=-id->rang+(double)j*dr;
      mo_pv_t(&p,v,x,id->type,ad); // total field
      
      #pragma omp critical
      if(cabs(p)>id->mp) id->mp=cabs(p);
      
      id->vp[i*id->m+j]=p;
    }
  }
}

void output_field(char *pl,IMD *id,MOBJ *ad)
{
  void output_png(int nt,double complex cet,char *pl,IMD *id);
  void output_color_bar(IMD *id);
  
  FILE *fp;
  char fn[128];
  double dt;
  int n;
  
  dt=1.0/(ad->md[0].aw.f*(double)id->ts);
  
  #pragma omp parallel for schedule(dynamic) 
  for(n=0;n<id->ts;n++){
    output_png(n,cexp(-I*2.0*M_PI*ad->md[0].aw.f*dt*(double)n),pl,id);
  }

  // print info
  sprintf(fn,"%s/%s_info.txt",id->dir_name,pl);
  fp=fopen(fn,"wt");
  if(fp==NULL){
    printf("Failed to open the %s file. Exit...\n",fn);
    exit(1);
  }
  fprintf(fp,"the range of color bar\n");
  fprintf(fp,"p is %8e to %8e\n",-id->mp,id->mp);
  fclose(fp);
  
  // output color bar image
  output_color_bar(id);
}

void output_png(int nt,double complex cet,char *pl,IMD *id)
{
  int color_rgb(double x,png_byte *r,png_byte *g,png_byte *b); // -1 <= x <= 1
  
  FILE *fp;
  char fname[256];
  int j,i,sj,si,m,scale;
  png_uint_32 width,height;
  png_structp png_p;
  png_infop info_p;
  png_bytepp pd_p;
  png_byte r,g,b;

  m=id->m;
  scale=id->scale;

  width =m*(scale+1);
  height=m*(scale+1);

  png_p =png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  info_p=png_create_info_struct(png_p);
  sprintf(fname,"%s/%s_p_%03d.png",id->dir_name,pl,nt);
  fp=fopen(fname,"wb");
  if(fp==NULL){
    printf("Failed to open the %s file. Exit...\n",fname);
    exit(1);
  }
   
  png_init_io(png_p,fp);
  png_set_IHDR(png_p,info_p,width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
  pd_p=(png_bytepp)png_malloc(png_p,sizeof(png_bytep)*height);
  png_set_rows(png_p,info_p,pd_p);

  for(j=0;j<height;j++){
    pd_p[j]=(png_bytep)png_malloc(png_p,sizeof(png_byte)*width*3);
  }
  
  for(i=0;i<m;i++){
    for(j=0;j<m;j++){
      
      color_rgb(creal(cet*id->vp[i*m+j])/id->mp,&r,&g,&b);
      for(si=0;si<=scale;si++){
        for(sj=0;sj<=scale;sj++){
          pd_p[i*(scale+1)+si][(j*(scale+1)+sj)*3+0]=r;
          pd_p[i*(scale+1)+si][(j*(scale+1)+sj)*3+1]=g;
          pd_p[i*(scale+1)+si][(j*(scale+1)+sj)*3+2]=b;
        }
      }
    }
  }
  
  png_write_png(png_p,info_p,PNG_TRANSFORM_IDENTITY,NULL);
    
  for(j=0;j<height;j++){
    png_free(png_p,pd_p[j]);
  }
  png_free(png_p,pd_p);
   
  fclose(fp);
}

void output_color_bar(IMD *id)
{
  int color_rgb(double x,png_byte *r,png_byte *g,png_byte *b); // -1 <= x <= 1
  
  FILE *fp;
  char fname[128];
  int j,i;
  
  png_uint_32 width,height;
  png_structp png;
  png_infop info;
  png_bytepp pdata;
  png_byte r,g,b;

  sprintf(fname,"%s/color_bar.png",id->dir_name);

  height=id->m*(id->scale+1);
  width=height/16;
  
  png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  info= png_create_info_struct(png);
  
  fp=fopen(fname,"wb");
  if(fp==NULL){
    printf("Failed to open the %s file. Exit...\n",fname);
    exit(1);
  }
  
  png_init_io(png, fp);
  png_set_IHDR(png,info,width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
  pdata=(png_bytepp)png_malloc(png, sizeof(png_bytep)*height);
  png_set_rows(png,info,pdata);

  for(j=0;j<height;j++){
    pdata[j]=(png_bytep)png_malloc(png,sizeof(png_byte)*width*3);
  }
  
  for(i=0;i<height;i++){
    color_rgb(1.0-(2.0/(double)height)*(double)i,&r,&g,&b);
    for(j=0;j<width;j++){
      pdata[i][j*3+0]=r;
      pdata[i][j*3+1]=g;
      pdata[i][j*3+2]=b;
    }
  }
  
  png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);
  
  for(j=0;j<height;j++){
    png_free(png,pdata[j]);
  }
  png_free(png,pdata);
  fclose(fp);
}

int color_rgb(double x,png_byte *r,png_byte *g,png_byte *b) // -1 <= x <= 1
{
  double i_nc,dr,dg,db;
  unsigned int i,n,nc,nd;

  if(x<-1.0 || x>1.0){
    *r=0x00;    *g=0x00;    *b=0x00;
    return -1;
  }
  
  n=(unsigned int)floor(pow(2,23)*(x+1.0));
  nc=(unsigned int)pow(2,21);
  i_nc=1.0/(double)nc;
  
  if(n<nc*1)      i=1;
  else if(n<nc*2) i=2;
  else if(n<nc*3) i=3;
  else if(n<nc*4) i=4;
  else if(n<nc*5) i=5;
  else if(n<nc*6) i=6;
  else if(n<nc*7) i=7;
  else if(n<nc*8) i=8;
  else {
    *r=ct1[8][0];    *g=ct1[8][1];    *b=ct1[8][2];
    return 0;
  }
    
  nd=n-nc*(i-1);
  dr=(double)(ct1[i][0]-ct1[i-1][0])*i_nc;
  dg=(double)(ct1[i][1]-ct1[i-1][1])*i_nc;
  db=(double)(ct1[i][2]-ct1[i-1][2])*i_nc;
  *r=(png_byte)floor((double)ct1[i-1][0]+dr*(double)nd);
  *g=(png_byte)floor((double)ct1[i-1][1]+dg*(double)nd);
  *b=(png_byte)floor((double)ct1[i-1][2]+db*(double)nd);
  
  return 0;  
}
