#include <stdio.h>
#include <stdlib.h>
#include "my_utils.h"

// using Rodrigues' rotation formula
void rotation_msh_data(char *fn,char *ofn,double *nv,double th);

int main(int argc,char **argv)
{
	double nv[3];
	int err;
	
	if(argc!=7){
		printf("This program needs command-line argument as follows\n%s msh_filename (unit_vector_rotation_axis)_x (unit_vector_rotation_axis)_y (unit_vector_rotation_axis)_z rotation_angle output_filename\nExit...\n",argv[0]);
		exit(1);
	}
	printf("mesh filename             = %s\n",argv[1]);
	printf("unit vector rotation axis = (% 8.6g, % 8.6g, % 8.6g)\n",atof(argv[2]),atof(argv[3]),atof(argv[4]));
	printf("rotation angle      [rad] = % 8.6g\n",atof(argv[5]));
	printf("output filename           = %s\n",argv[6]);
	printf("continue?(y/n) ");
	if(getchar()!='y'){ printf("Exit\n");  exit(0);}

	nv[0]=atof(argv[2]);
	nv[1]=atof(argv[3]);
	nv[2]=atof(argv[4]);
	
	err=vuni_d(nv);
	if(err<0){
		printf("absolute value of rotation axis = 0 ( zero vector ). Exit...\n");
	}
	
	rotation_msh_data(argv[1],argv[6],nv,atof(argv[5]));

	return 0;
}

void rotation_msh_data(char *fn,char *ofn,double *nv,double th)
{
	FILE *rf,*of;
	char str[256]="";
	double Rn[9],st,ct,td,vo[3];
	int N,i,j,ti;
	
	// rotation matrix
	st=sin(th);
	ct=cos(th);
	Rn[0]=ct+nv[0]*nv[0]*(1.0-ct);
	Rn[1]=nv[0]*nv[1]*(1.0-ct)-nv[2]*st;
	Rn[2]=nv[2]*nv[0]*(1.0-ct)+nv[1]*st;
	Rn[3]=nv[0]*nv[1]*(1.0-ct)+nv[2]*st;
	Rn[4]=ct+nv[1]*nv[1]*(1.0-ct);
	Rn[5]=nv[1]*nv[2]*(1.0-ct)-nv[0]*st;
	Rn[6]=nv[2]*nv[0]*(1.0-ct)-nv[1]*st;
	Rn[7]=nv[1]*nv[2]*(1.0-ct)+nv[0]*st;
	Rn[8]=ct+nv[2]*nv[2]*(1.0-ct);

	if((rf=fopen(fn,"rt"))==NULL){    printf("Can not open the %s file.\n",fn);    exit(1);  }
	if((of=fopen(ofn,"wt"))==NULL){    printf("Can not create the %s file.\n",ofn);    exit(1);  }

	for(i=0;i<4;i++){
		fgets(str,256,rf);      fprintf(of,"%s",str);
	}
	fgets(str,256,rf);  fprintf(of,"%s",str);
	N=atoi(str);

	for(i=1;i<=N;i++){
		fscanf(rf,"%d",&ti); fprintf(of,"%d",ti);
		fscanf(rf,"%lf",&td); vo[0]=td;
		fscanf(rf,"%lf",&td); vo[1]=td;
		fscanf(rf,"%lf\n",&td); vo[2]=td;

		for(j=0;j<3;j++){
			td=Rn[j*3+0]*vo[0]+Rn[j*3+1]*vo[1]+Rn[j*3+2]*vo[2];
			fprintf(of," %16.15f",td);
		}
		fprintf(of,"\n");
	}

	for(i=0;i<2;i++){
		fgets(str,256,rf);	fprintf(of,"%s",str);
	}
	fgets(str,256,rf);  fprintf(of,"%s",str);
	N=atoi(str);
	for(i=0;i<=N;i++){
		fgets(str,256,rf);	fprintf(of,"%s",str);
	}

	fclose(rf);
	fclose(of);
	printf("Done writing %s\n",ofn);
}
