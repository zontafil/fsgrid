#include <iostream>
#include "fsgrid.hpp"

template<class T, int stencil> void timeit(std::array<int32_t, 3> globalSize, std::array<int32_t, 3> isPeriodic, int iterations){
   double t1,t2;   
   FsGrid<T ,stencil> testGrid(globalSize, MPI_COMM_WORLD, isPeriodic);
   int rank,size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   t1=MPI_Wtime();      
   for(int i = 0; i < iterations; i++) {
      testGrid.updateGhostCells();         
   }
   t2=MPI_Wtime();
   if(rank==0)
      printf("%g s per update: nprocs %d, grid is %d x %d x %d, stencil %d, element size %ld \n", (t2 - t1)/iterations, size, globalSize[0], globalSize[1], globalSize[2], stencil, sizeof(*testGrid.get(0,0,0)));
}


int main(int argc, char** argv) {
   
   MPI_Init(&argc,&argv);

   int rank,size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   // Create a 8×8 Testgrid
   std::array<int32_t, 3> globalSize{2000,1000,1};
   std::array<int, 3> isPeriodic{false,false,true};

   const int iterations = 200;

   timeit<std::array<double,1>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,2>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,4>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,8>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,16>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,32>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,64>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,128>, 2>(globalSize, isPeriodic, iterations);
   
      
   MPI_Finalize();
   return 0;
}
