#include <stdio.h>
#include <stdlib.h>
#include "my_utils.h"

void stretch_msh_data(char *fn,char *ofn,double *sv);

int main(int argc,char **argv)
{
	double sv[3];

	if(argc!=6){
		printf("This program needs command-line argument as follows\n%s msh_filename stretch_vector_x stretch_vector_y stretch_vector_z output_filename\nExit...\n",argv[0]);
		exit(1);
	}
	printf("mesh filename   = %s\n",argv[1]);
	printf("stretch vector  = (%8.6g,%8.6g,%8.6g)\n",atof(argv[2]),atof(argv[3]),atof(argv[4]));
	printf("output filename = %s\n",argv[5]);
	printf("continue?(y/n) ");
	if(getchar()!='y'){ printf("Exit\n");  exit(0);}

	sv[0]=atof(argv[2]);
	sv[1]=atof(argv[3]);
	sv[2]=atof(argv[4]);	

	stretch_msh_data(argv[1],argv[5],sv);

	return 0;
}

void stretch_msh_data(char *fn,char *ofn,double *sv)
{
	FILE *rf,*of;
	char str[256]="";
	double td,vf[3],nsv[3],tv,rv,asv,ssv;
	int N,i,ti,j;

	nsv[0]=sv[0];
	nsv[1]=sv[1];
	nsv[2]=sv[2];
	asv=vabs_d(nsv);
	if(asv==0.0){
		printf("the magnitude of stretch vector is zero. Exit...\n");
		exit(1);
	} 
	nsv[0]/=asv;
	nsv[1]/=asv;
	nsv[2]/=asv;
	//printf("stretch vector absolute value = %8.6g\n",asv);
	//printf("unit stretch vector           = (%8.6g,%8.6g,%8.6g)\n",nsv[0],nsv[1],nsv[2]);

	if((rf=fopen(fn,"rt"))==NULL){    printf("Can not open the %s file.\n",fn);    exit(1);  }
	if((of=fopen(ofn,"wt"))==NULL){    printf("Can not create the %s file.\n",ofn);    exit(1);  }

	for(i=0;i<4;i++){
		fgets(str,256,rf);	fprintf(of,"%s",str);
	}
	fgets(str,256,rf);  fprintf(of,"%s",str);
	N=atoi(str);

	for(i=1;i<=N;i++){
		fscanf(rf,"%d",&ti); fprintf(of,"%d",ti);
		fscanf(rf,"%lf",&td); vf[0]=td;
		fscanf(rf,"%lf",&td); vf[1]=td;
		fscanf(rf,"%lf\n",&td); vf[2]=td;

		td=vdot_d(vf,nsv);
		for(j=0;j<3;j++){
			tv=td*nsv[j];
			rv=vf[j]-tv;
			ssv=asv*tv+rv;
			fprintf(of," %16.15f",ssv);
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

