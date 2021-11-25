#include "bem3_aw_b2.h"

int main(int argc,char **argv)
{
  MOBJ mo;
  
  mo_initialize(argc,argv,&mo);
  mo_print_data(&mo);
  
  mo_solve_iter(&mo);
  mo_dat_write(argv[2],&mo);
  mo_output_node_particles(argv[2],&mo);
  printf("Done writing %s\n",argv[2]);
  
  mo_finalize(&mo);
  return 0;
}
  
