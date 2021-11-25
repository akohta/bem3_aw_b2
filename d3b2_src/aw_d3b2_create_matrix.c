#include "bem3_aw_b2.h"

int main(int argc,char **argv)
{
  DMDA ad;
  
  read_dmda2(argc,argv,&ad); // read datafile for single object.
  print_dmda2(&ad);          // print data.
  initialize_dmda2(&ad);     // memory allocation and initialize coefficients.
  solve_bieq_dmda2(&ad,argv[3]); // create coefficient matrix and its inverse matrix. 
  finalize_dmda2(&ad);       // memory free.
  return 0;
}
