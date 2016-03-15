#include <iostream>
#include "fsgrid.hpp"

int main(int argc, char** argv) {
   
   MPI_Init(&argc,&argv);

   // Create a 8×8 Testgrid
   std::array<uint32_t, 3> globalSize{8,8,1};
   std::array<int, 3> isPeriodic{false,false,true};
   {
      FsGrid<int,3,1> testGrid(globalSize, MPI_COMM_WORLD, isPeriodic);
   }
   MPI_Finalize();

   return 0;
}
