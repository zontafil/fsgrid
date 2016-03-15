#include <stdlib.h>
#include <stdio.h>
#include <array>
#include <math.h>
#include <algorithm>
#include <limits>


void computeDomainDecomposition(const std::array<uint32_t, 3>& GlobalSize, int nProcs, std::array<uint32_t,3>& processDomainDecomposition) {
  
   std::array<double, 3> systemDim;
   std::array<double, 3 > processBox;
   double optimValue = std::numeric_limits<double>::max();
  
   for(int i = 0; i < 3; i++) {
      systemDim[i] = (double)GlobalSize[i];
   }
   processDomainDecomposition = {1, 1, 1};
  
   for (int i = 1; i <= nProcs; i++) {
      if( i  > systemDim[0])
         continue;
      processBox[0] = std::max(systemDim[0]/i, 1.0);
    
      for (int j = 1; j <= nProcs; j++) {
         if( i * j  > nProcs || j > systemDim[1])
            continue;
       
         processBox[1] = std::max(systemDim[1]/j, 1.0);
         for (int k = 1; k <= nProcs; k++) {
            if( i * j * k > nProcs || k > systemDim[2])
               continue;
            processBox[2] = std::max(systemDim[2]/k, 1.0);
            double value = 
               1000 * processBox[0] * processBox[1] * processBox[2] + 
               (i > 1 ? processBox[1] * processBox[2]: 0) +
               (j > 1 ? processBox[0] * processBox[2]: 0) +
               (k > 1 ? processBox[0] * processBox[1]: 0);
	
            if(value < optimValue ){
               optimValue = value;
               processDomainDecomposition[0] = i;
               processDomainDecomposition[1] = j;
               processDomainDecomposition[2] = k;
            }
//            printf("%g: %d %d %d with box %g %g %g\n", value, i, j, k, processBox[0], processBox[1], processBox[2]);
         }
      }
   }
}

int main(int argc, char **argv){
  
   std::array<uint32_t,3>  sys;
   std::array<uint32_t,3> processDomainDecomposition;

   if(argc != 5) {
      printf("Usage %s size_x size_y size_z nProcesses\n", argv[0]);
      exit(1);
   }
   
      
   sys[0] = atof(argv[1]);
   sys[1] = atof(argv[2]);
   sys[2] = atof(argv[3]);
   uint nProcs = atoi(argv[4]);

   computeDomainDecomposition(sys, nProcs, processDomainDecomposition);
   printf("DD of %d %d %d for %d processes is %d %d %d \n", 
          sys[0], sys[1], sys[2], nProcs,
          processDomainDecomposition[0], processDomainDecomposition[1], processDomainDecomposition[2]);


  
}
