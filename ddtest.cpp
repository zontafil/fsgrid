/*
  Copyright (C) 2016 Finnish Meteorological Institute
  Copyright (C) 2016 CSC -IT Center for Science 

  This file is part of fsgrid
 
  fsgrid is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  lib-Slice3D is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY;
  without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with fsgrid.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <array>
#include <math.h>
#include <algorithm>
#include <limits>


void computeDomainDecomposition(const std::array<int, 3>& GlobalSize, int nProcs, std::array<int,3>& processDomainDecomposition) {
  
   std::array<double, 3> systemDim;
   std::array<double, 3 > processBox;
   double optimValue = std::numeric_limits<double>::max();
  
   for(int i = 0; i < 3; i++) {
      systemDim[i] = (double)GlobalSize[i];
   }
   processDomainDecomposition = {1, 1, 1};
   
   for (int i = 1; i <= std::min(nProcs, GlobalSize[0]); i++) {
      processBox[0] = std::max(systemDim[0]/i, 1.0);
      for (int j = 1; j <= std::min(nProcs, GlobalSize[1]) ; j++) {
         if( i * j  > nProcs )
            break;
         processBox[1] = std::max(systemDim[1]/j, 1.0);
         for (int k = 1; k <= std::min(nProcs, GlobalSize[2]); k++) {
            if( i * j * k > nProcs )
               break;
            processBox[2] = std::max(systemDim[2]/k, 1.0);
            double value = 
               10 * processBox[0] * processBox[1] * processBox[2] + 
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
  
   std::array<int,3>  sys;
   std::array<int,3> processDomainDecomposition;

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
