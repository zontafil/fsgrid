#include <iostream>
#include "fsgrid.hpp"

int main(int argc, char** argv) {
   
   MPI_Init(&argc,&argv);

   // Create a 8×8 Testgrid
   std::array<uint32_t, 2> globalSize{8,8};
   FsGrid<int,2,1> testGrid(globalSize, MPI_COMM_WORLD);

   MPI_Finalize();

   return 0;
}
