#include <stdio.h>
#include <stdint.h>
 uint64_t s128[2];
 uint64_t s1k[16]; 
 uint64_t s[2];
main(int argc, char **argv){
  double t1, t2;
  uint32_t val;
  uint64_t val64;
  uint32_t i;
  double MPI_Wtime();
  uint32_t buffer[1000];

  MPI_Init(&argc, &argv);

  s128[0] = 123456 ; s128[1] = 456789;
  s[0] = 123456 ; s[1] = 456789;
  for(i=0 ; i<16 ; i++) s1k[i] = 2*i + 125;

  for(i=0 ; i<100000000 ; i++) val = next128r();
  t1 = MPI_Wtime();
  for(i=0 ; i<1000000000 ; i++) val = next128r();
  t2 = MPI_Wtime();
  printf("val = %d, next128r time = %f\n",val,t2-t1);

  for(i=0 ; i<100000 ; i++) Vnext128r(buffer,1000);
  t1 = MPI_Wtime();
  for(i=0 ; i<1000000000 ; i++) Vnext128r(buffer,1);
  t2 = MPI_Wtime();
  val = buffer[1];
  printf("val = %d, Vnext128r by 1 time = %f\n",val,t2-t1);

  for(i=0 ; i<100000 ; i++) Vnext128r(buffer,1000);
  t1 = MPI_Wtime();
  for(i=0 ; i<500000000 ; i++) Vnext128r(buffer,2);
  t2 = MPI_Wtime();
  val = buffer[1];
  printf("val = %d, Vnext128r by 2 time = %f\n",val,t2-t1);

  for(i=0 ; i<100000 ; i++) Vnext128r(buffer,1000);
  t1 = MPI_Wtime();
  for(i=0 ; i<200000000 ; i++) Vnext128r(buffer,5);
  t2 = MPI_Wtime();
  val = buffer[1];
  printf("val = %d, Vnext128r by 5 time = %f\n",val,t2-t1);

  for(i=0 ; i<100000 ; i++) Vnext128r(buffer,1000);
  t1 = MPI_Wtime();
  for(i=0 ; i<100000000 ; i++) Vnext128r(buffer,10);
  t2 = MPI_Wtime();
  val = buffer[9];
  printf("val = %d, Vnext128r by 10 time = %f\n",val,t2-t1);

  for(i=0 ; i<100000 ; i++) Vnext128r(buffer,1000);
  t1 = MPI_Wtime();
  for(i=0 ; i<10000000 ; i++) Vnext128r(buffer,100);
  t2 = MPI_Wtime();
  val = buffer[99];
  printf("val = %d, Vnext128r by 100 time = %f\n",val,t2-t1);

  for(i=0 ; i<100000 ; i++) Vnext128r(buffer,1000);
  t1 = MPI_Wtime();
  for(i=0 ; i<1000000 ; i++) Vnext128r(buffer,1000);
  t2 = MPI_Wtime();
  val = buffer[999];
  printf("val = %d, Vnext128r by 1000 time = %f\n",val,t2-t1);

  for(i=0 ; i<100000000 ; i++) val64 = next_orig();
  t1 = MPI_Wtime();
  for(i=0 ; i<1000000000 ; i++) val64 = next_orig();
  t2 = MPI_Wtime();
  printf("val = %d, next_orig time = %f\n",val64,t2-t1);

  for(i=0 ; i<100000000 ; i++) val64 = next128();
  t1 = MPI_Wtime();
  for(i=0 ; i<1000000000 ; i++) val64 = next128();
  t2 = MPI_Wtime();
  printf("val = %ld, next128 time = %f\n",val64,t2-t1);

  for(i=0 ; i<100000000 ; i++) val64 = next1024();
  t1 = MPI_Wtime();
  for(i=0 ; i<1000000000 ; i++) val64 = next1024();
  t2 = MPI_Wtime();
  printf("val = %ld, next1024 time = %f\n",val64,t2-t1);

  for(i=0 ; i<100000000 ; i++) val64 = next64();
  t1 = MPI_Wtime();
  for(i=0 ; i<1000000000 ; i++) val64 = next64();
  t2 = MPI_Wtime();
  printf("val = %ld, next64 time = %f\n",val64,t2-t1);
  MPI_Finalize();
}

