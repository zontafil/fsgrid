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
      
      std::array<int, 3>  localSize=testGrid.getLocalSize();


      
      for(int z = 0; z < localSize[2]; z++){
         for(int y = 0; y < localSize[1]; y++){
            for(int x = 0; x < localSize[0]; x++){
               testGrid.get(x, y, z) = rank;
            }
         }
      }
      
      
      if(rank==0) {
         printf("local size %d %d %d\n", localSize[0], localSize[1], localSize[2]);
         printf("----------------------------------\n");
         for(int z = -1; z < localSize[2] +1; z++){
            printf("z=%d\n", z);
            for(int y = -1; y < localSize[1] + 1; y++){
               printf("y=%d :", y);
               for(int x = -1; x < localSize[0] +1;  x++){
                  printf("%d ", testGrid.get(x, y, z));
               }
               printf("\n");
            }
         }
         printf("----------------------------------\n");
      }
      

      testGrid.updateGhostCells();

      if(rank==0) {
         printf("local size %d %d %d\n", localSize[0], localSize[1], localSize[2]);
         printf("----------------------------------\n");
         for(int z = -1; z < localSize[2] +1; z++){
            for(int y = -1; y < localSize[1] + 1; y++){
               for(int x = -1; x < localSize[0] +1;  x++){
                  printf("%d ", testGrid.get(x, y, z));
               }
               printf("\n");
            }
         }
         printf("----------------------------------\n");
      }
   }
   

   
   MPI_Finalize();

   return 0;
}
