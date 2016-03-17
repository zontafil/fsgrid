#include <iostream>
#include "fsgrid.hpp"

int main(int argc, char** argv) {
   
   MPI_Init(&argc,&argv);

   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Create a 8×8 Testgrid
   std::array<int32_t, 3> globalSize{8,8,1};
   std::array<int, 3> isPeriodic{false,false,true};
   {
      FsGrid<int,1> testGrid(globalSize, MPI_COMM_WORLD, isPeriodic);

      if(rank == 0) {
         for(int i=0; i<8; i++) {
            auto taskLid = testGrid.getTaskForGlobalID(8*1*0+8*0+i);
            std::cerr << "Cell ( " << i << ", 0, 0) is located on task "
               << taskLid.first << std::endl;
            std::cerr << "   and it has LocalID " << taskLid.second << std::endl;
         }
         for(int i=0; i<8; i++) {
            auto taskLid = testGrid.getTaskForGlobalID(8*1*0+8*i+0);
            std::cerr << "Cell ( " << 0 << ", " << i << ", 0) is located on task "
               << taskLid.first << std::endl;
            std::cerr << "   and it has LocalID " << taskLid.second << std::endl;
         }
         for(int i=0; i<8; i++) {
            auto taskLid = testGrid.getTaskForGlobalID(8*1*0+8*i+i);
            std::cerr << "Cell ( " << i << ", " << i << ", 0) is located on task "
               << taskLid.first << std::endl;
         }
      }
   }
   MPI_Finalize();

   return 0;
}
