#include <iostream>
#include "fsgrid.hpp"

int main(int argc, char** argv) {
   
   MPI_Init(&argc,&argv);

   int rank,size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   // Create a 8×8 Testgrid
   std::array<int32_t, 3> globalSize{9,2,1};
   std::array<int, 3> isPeriodic{true,true,true};
   {
      FsGrid<int,1> testGrid(globalSize, MPI_COMM_WORLD, isPeriodic);
      if(rank == 0) {
         std::cerr << " --- Test task mapping functions ---" << std::endl;
      }
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


      
      int z = 0;
      for(int y = -1; y < localSize[1] + 1; y++){
         for(int x = -1; x < localSize[0] + 1; x++){
            *(testGrid.get(x, y, z)) = rank;
         }
      }
      
      

      if(rank==1) {
         printf("local size %d %d %d\n", localSize[0], localSize[1], localSize[2]);
         printf("----------------------------------\n");
         int z = 0;
         printf("z=%d\n", z);
         for(int y = -1; y < localSize[1] + 1; y++){
            printf("y=%d :", y);
            for(int x = -1; x < localSize[0] +1;  x++){
               printf("%d ", *(testGrid.get(x, y, z)));
            }
            printf("\n");
         }
         printf("----------------------------------\n");
      }
      

      testGrid.updateGhostCells();

      if(rank==1) {
         printf("local size %d %d %d\n", localSize[0], localSize[1], localSize[2]);
         printf("----------------------------------\n");
         int z = 0;
         printf("z=%d\n", z);
         for(int y = -1; y < localSize[1] + 1; y++){
            printf("y=%d :", y);
            for(int x = -1; x < localSize[0] +1;  x++){
               printf("%d ", *(testGrid.get(x, y, z)));
            }
            printf("\n");
         }
         printf("----------------------------------\n");
      }
      

         if(rank == 0) {
            std::cerr << " --- Test data transfer into the grid ---" << std::endl;
         }

         // First, setup grid coupling to do everything from task 0
         if(rank == 0) {
            testGrid.setupForGridCoupling(8*8);
            for(int y=0; y<8; y++) {
               for(int x=0; x<8; x++) {
                  testGrid.setGridCoupling(y*8+x,0); // All cells are coupled to 0
               }
            }
         } else {
            testGrid.setupForGridCoupling(0);
         }
         testGrid.finishGridCoupling();



         // Fill in some junk data from task 0
         std::vector<int> fillData;
         if(rank == 0) {
            fillData.resize(8*8);

            // We are going to send data for all 8×8 Cells
            testGrid.setupForTransferIn(8*8);
            for(int y=0; y<8; y++) {
               for(int x=0; x<8; x++) {
                  fillData[y*8 + x] = x*y;

                  testGrid.transferDataIn(y*8+x,fillData[y*8+x]);
               }
            }
         } else {
            // The others simply recieve
            testGrid.setupForTransferIn(0);
         }
         testGrid.finishTransfersIn();


         // Now have each task output their data
         for(int i=0; i<size; i++) {
            if(i == rank) {
               std::cerr << "Contents of Task #" << rank << ": " << std::endl;
               std::array<int32_t,3> localSize = testGrid.getLocalSize();
               for(unsigned int y=0; y<localSize[1]; y++) {
                  for(unsigned int x=0; x<localSize[0]; x++) {
                     std::cerr << *testGrid.get(x,y,0) << ", ";
                  }
                  std::cerr << std::endl;
               }
            }
            MPI_Barrier(MPI_COMM_WORLD);
         }


         // Transfer it back
         std::vector<int> returnedData;
         if(rank ==0) {
            returnedData.resize(8*8);

            testGrid.setupForTransferOut(8*8);
            for(int y=0; y<8; y++) {
               for(int x=0; x<8; x++) {
                  testGrid.transferDataOut(y*8+x,returnedData[y*8+x]);
               }
            }
         } else {
            testGrid.setupForTransferOut(0);
         }
         testGrid.finishTransfersOut();

         // Validate the result
         if(rank == 0) {
            std::cerr << " --------- " << std::endl;
            std::cerr << "Returned array contents:" << std::endl;
            for(unsigned int y=0; y<8; y++) {
               for(unsigned int x=0; x<8; x++) {
                  std::cerr << returnedData[y*8+x] << ", ";
               }
               std::cerr << std::endl;
            }
         }
      }
      

      
      MPI_Finalize();

      return 0;
   }
