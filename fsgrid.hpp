#pragma once
/*
  Copyright (C) 2016 Finnish Meteorological Institute
  Copyright (C) 2016 CSC -IT Center for Science

  This file is part of fsgrid

  fsgrid is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  fsgrid is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY;
  without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with fsgrid.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <array>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <limits>
#include <stdint.h>
#include <cassert>
#ifdef __CUDACC__
   #define ARCH_HOSTDEV __host__ __device__
   #include "cuda.h"
   #include "cuda_runtime.h"
#else
   #define ARCH_HOSTDEV
#endif





struct FsGridTools{

   //! Helper function: calculate position of the local coordinate space for the given dimension
   // \param globalCells Number of cells in the global Simulation, in this dimension
   // \param ntasks Total number of tasks in this dimension
   // \param my_n This task's position in this dimension
   // \return Cell number at which this task's domains cells start (actual cells, not counting ghost cells)
   static int32_t calcLocalStart(int32_t globalCells, int ntasks, int my_n) {
      int n_per_task = globalCells / ntasks;
      int remainder = globalCells % ntasks;

      if(my_n < remainder) {
         return my_n * (n_per_task+1);
      } else {
         return my_n * n_per_task + remainder;
      }
   }
      //! Helper function: calculate size of the local coordinate space for the given dimension
      // \param globalCells Number of cells in the global Simulation, in this dimension
      // \param ntasks Total number of tasks in this dimension
      // \param my_n This task's position in this dimension
      // \return Nmuber of cells for this task's local domain (actual cells, not counting ghost cells)
      static int32_t calcLocalSize(int32_t globalCells, int ntasks, int my_n) {
         int n_per_task = globalCells/ntasks;
         int remainder = globalCells%ntasks;
         if(my_n < remainder) {
            return n_per_task+1;
         } else {
            return n_per_task;
         }
      }


      //! Helper function to optimize decomposition of this grid over the given number of tasks
      static void computeDomainDecomposition(const int GlobalSize[3], int nProcs, std::array<int,3>& processDomainDecomposition) {
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
               }
            }
         }

         if(optimValue == std::numeric_limits<double>::max() ||
               processDomainDecomposition[0] * processDomainDecomposition[1] * processDomainDecomposition[2] != nProcs) {
            std::cerr << "FSGrid domain decomposition failed, are you running on a prime number of tasks?" << std::endl;
            throw std::runtime_error("FSGrid computeDomainDecomposition failed");
         }
      }

};


struct FsGridCouplingInformation {
   std::vector<int> externalRank; //!< MPI rank that each cell is being communicated to externally

   void setCouplingSize(size_t totalStorageSize) {

      // Set up coupling information to external grid (fill with MPI_PROC_NULL to begin with)
      // Before actual grid coupling can be done, this information has to be filled in.
      if(externalRank.size() != totalStorageSize) {
         externalRank.resize(totalStorageSize);
      }
      for(uint i=0; i<externalRank.size(); i++) {
         externalRank[i] = MPI_PROC_NULL;
      }
   }
};

/*! Simple cartesian, non-loadbalancing MPI Grid for use with the fieldsolver
 *
 * \param T datastructure containing the field in each cell which this grid manages
 * \param stencil ghost cell width of this grid
 */
template <typename T, int TDim, int stencil> class FsGrid : public FsGridTools{

   public:

      typedef int64_t LocalID;
      typedef int64_t GlobalID;
      T* data;

   // Legacy constructor from coupling reference
   FsGrid(int32_t globalSize[3], MPI_Comm parent_comm, std::array<bool,3> isPeriodic, FsGridCouplingInformation& coupling) : FsGrid(globalSize, parent_comm, isPeriodic, &coupling) {}

   // Legacy constructor from coupling reference
   FsGrid(std::array<int32_t,3> globalSize, MPI_Comm parent_comm, std::array<bool,3> isPeriodic, FsGridCouplingInformation& coupling) : FsGrid(globalSize, parent_comm, isPeriodic, &coupling) {}

      /*! Constructor for this grid.
       * \param globalSize Cell size of the global simulation domain.
       * \param MPI_Comm The MPI communicator this grid should use.
       * \param isPeriodic An array specifying, for each dimension, whether it is to be treated as periodic.
       */
   FsGrid(int32_t globalSize[3], MPI_Comm parent_comm, std::array<bool,3> isPeriodic, FsGridCouplingInformation* coupling)
            : globalSize{globalSize[0], globalSize[1], globalSize[2]}, coupling(coupling) {
         int status;
         int size;

         status = MPI_Comm_size(parent_comm, &size);

         // Heuristically choose a good domain decomposition for our field size
         computeDomainDecomposition(globalSize, size, ntasks);

         //set private array
         periodic = isPeriodic;
         //set temporary int array for MPI_Cart_create
         std::array<int, 3> isPeriodicInt;
         for(unsigned int i=0; i < isPeriodic.size(); i++) {
            isPeriodicInt[i] = (int)isPeriodic[i];
         }

         // Create cartesian communicator. Note, that reorder is false so
         // ranks match the ones in parent_comm
         status = MPI_Cart_create(parent_comm, 3, ntasks.data(), isPeriodicInt.data(), 0, &comm3d);
         if(status != MPI_SUCCESS) {
            std::cerr << "Creating cartesian communicatior failed when attempting to create FsGrid!" << std::endl;
            throw std::runtime_error("FSGrid communicator setup failed");
         }

         status = MPI_Comm_rank(comm3d, &rank);
         if(status != MPI_SUCCESS) {
            std::cerr << "Getting rank failed when attempting to create FsGrid!" << std::endl;

            // Without a rank, there's really not much we can do. Just return an uninitialized grid
            // (I suppose we'll crash after this, anyway)
            return;
         }


         // Determine our position in the resulting task-grid
         status = MPI_Cart_coords(comm3d, rank, 3, taskPosition.data());
         if(status != MPI_SUCCESS) {
            std::cerr << "Rank " << rank
               << " unable to determine own position in cartesian communicator when attempting to create FsGrid!"
               << std::endl;
         }

         // Allocate the array of neighbours
         for(int i=0; i<size; i++) {
            neighbour_index.push_back(MPI_PROC_NULL);
         }
         for(int i=0; i<27; i++) {
            neighbour[i]=MPI_PROC_NULL;
         }

         // Get the IDs of the 26 direct neighbours
         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  std::array<int,3> neighPosition;

                  /*
                   * Figure out the coordinates of the neighbours in all three
                   * directions
                   */
                  neighPosition[0]=taskPosition[0]+x;
                  if(isPeriodic[0]) {
                     neighPosition[0] += ntasks[0];
                     neighPosition[0] %= ntasks[0];
                  }

                  neighPosition[1]=taskPosition[1]+y;
                  if(isPeriodic[1]) {
                     neighPosition[1] += ntasks[1];
                     neighPosition[1] %= ntasks[1];
                  }

                  neighPosition[2]=taskPosition[2]+z;
                  if(isPeriodic[2]) {
                     neighPosition[2] += ntasks[2];
                     neighPosition[2] %= ntasks[2];
                  }

                  /*
                   * If those coordinates exist, figure out the responsible CPU
                   * and store its rank
                   */
                  if(neighPosition[0]>=0 && neighPosition[0]<ntasks[0] && neighPosition[1]>=0
                        && neighPosition[1]<ntasks[1] && neighPosition[2]>=0 && neighPosition[2]<ntasks[2]) {

                     // Calculate the rank
                     int neighRank;
                     status = MPI_Cart_rank(comm3d, neighPosition.data(), &neighRank);
                     if(status != MPI_SUCCESS) {
                        std::cerr << "Rank " << rank << " can't determine neighbour rank at position [";
                        for(int i=0; i<3; i++) {
                           std::cerr << neighPosition[i] << ", ";
                        }
                        std::cerr << "]" << std::endl;
                     }

                     // Forward lookup table
                     neighbour[(x+1)*9+(y+1)*3+(z+1)]=neighRank;

                     // Reverse lookup table
                     if(neighRank >= 0 && neighRank < size) {
                        neighbour_index[neighRank]=(char) ((x+1)*9+(y+1)*3+(z+1));
                     }
                  } else {
                     neighbour[(x+1)*9+(y+1)*3+(z+1)]=MPI_PROC_NULL;
                  }
               }
            }
         }


         // Determine size of our local grid
         for(int i=0; i<3; i++) {
            localSize[i] = calcLocalSize(globalSize[i],ntasks[i], taskPosition[i]);
            localStart[i] = calcLocalStart(globalSize[i],ntasks[i], taskPosition[i]);
         }

         if(  localSize[0] == 0 || (globalSize[0] > stencil && localSize[0] < stencil)
           || localSize[1] == 0 || (globalSize[1] > stencil && localSize[1] < stencil)
           || localSize[2] == 0 || (globalSize[2] > stencil && localSize[2] < stencil)) {
            std::cerr << "FSGrid space partitioning leads to a space that is too small on Rank " << rank << "." <<std::endl;
            std::cerr << "Please run with a different number of Tasks, so that space is better divisible." <<std::endl;
            throw std::runtime_error("FSGrid too small domains");
         }

         // Allocate local storage array
         totalStorageSize=1;
         for(int i=0; i<3; i++) {
            if(globalSize[i] <= 1) {
               // Collapsed dimension => only one cell thick
               storageSize[i] = 1;
            } else {
               // Size of the local domain + 2* size for the ghost cell stencil
               storageSize[i] = (localSize[i] + stencil*2);
            }
            totalStorageSize *= storageSize[i];
         }
         #ifdef USE_GPU
         cudaMallocManaged((void**)&data, totalStorageSize * TDim * sizeof(T));
         #else
         data = (T*) malloc(totalStorageSize * TDim * sizeof(T));
         #endif
         memset(data, 0, totalStorageSize * TDim * sizeof(T));
         coupling->setCouplingSize(totalStorageSize);

         MPI_Datatype mpiTypeT;
         MPI_Type_contiguous(TDim * sizeof(T), MPI_BYTE, &mpiTypeT);
         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  neighbourSendType[(x+1) * 9 + (y + 1) * 3 + (z + 1)] = MPI_DATATYPE_NULL;
                  neighbourReceiveType[(x+1) * 9 + (y + 1) * 3 + (z + 1)] = MPI_DATATYPE_NULL;
               }
            }
         }

         // Compute send and receive datatypes

         //loop through the shifts in the different directions
         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  int subarraySize[3];
                  int subarrayStart[3];
                  const int shiftId = (x+1) * 9 + (y + 1) * 3 + (z + 1);


                  if((storageSize[0] == 1 && x!= 0 ) ||
                     (storageSize[1] == 1 && y!= 0 ) ||
                     (storageSize[2] == 1 && z!= 0 ) ||
                     (x == 0 && y == 0 && z == 0)){
                     //skip flat dimension for 2 or 1D simulations, and self
                     neighbourSendType[shiftId] = MPI_DATATYPE_NULL;
                     neighbourReceiveType[shiftId] = MPI_DATATYPE_NULL;
                     continue;
                  }

                  subarraySize[0] = (x == 0) ? localSize[0] : stencil;
                  subarraySize[1] = (y == 0) ? localSize[1] : stencil;
                  subarraySize[2] = (z == 0) ? localSize[2] : stencil;

                  if( x == 0 || x == -1 )
                     subarrayStart[0] = stencil;
                  else if (x == 1)
                     subarrayStart[0] = storageSize[0] - 2 * stencil;
                  if( y == 0 || y == -1 )
                     subarrayStart[1] = stencil;
                  else if (y == 1)
                     subarrayStart[1] = storageSize[1] - 2 * stencil;
                  if( z == 0 || z == -1 )
                     subarrayStart[2] = stencil;
                  else if (z == 1)
                     subarrayStart[2] = storageSize[2] - 2 * stencil;

                  for(int i = 0;i < 3; i++)
                     if(storageSize[i] == 1)
                        subarrayStart[i] = 0;

                  int swappedStorageSize[3] = {storageSize[0], storageSize[1], storageSize[2]};
                  swapArray(swappedStorageSize);
                  swapArray(subarraySize);
                  swapArray(subarrayStart);
                  MPI_Type_create_subarray(3,
                                           swappedStorageSize,
                                           subarraySize,
                                           subarrayStart,
                                           MPI_ORDER_C,
                                           mpiTypeT,
                                           &(neighbourSendType[shiftId]) );

                  if(x == 1 )
                     subarrayStart[0] = 0;
                  else if (x == 0)
                     subarrayStart[0] = stencil;
                  else if (x == -1)
                     subarrayStart[0] = storageSize[0] -  stencil;
                  if(y == 1 )
                     subarrayStart[1] = 0;
                  else if (y == 0)
                     subarrayStart[1] = stencil;
                  else if (y == -1)
                     subarrayStart[1] = storageSize[1] -  stencil;
                  if(z == 1 )
                     subarrayStart[2] = 0;
                  else if (z == 0)
                     subarrayStart[2] = stencil;
                  else if (z == -1)
                     subarrayStart[2] = storageSize[2] -  stencil;
                  for(int i = 0;i < 3; i++)
                     if(storageSize[i] == 1)
                        subarrayStart[i] = 0;

                  swapArray(subarrayStart);
                  MPI_Type_create_subarray(3,
                                           swappedStorageSize,
                                           subarraySize,
                                           subarrayStart,
                                           MPI_ORDER_C,
                                           mpiTypeT,
                                           &(neighbourReceiveType[shiftId]));

               }
            }
         }

         for(int i=0;i<27;i++){
            if(neighbourReceiveType[i] != MPI_DATATYPE_NULL)
               MPI_Type_commit(&(neighbourReceiveType[i]));
            if(neighbourSendType[i] != MPI_DATATYPE_NULL)
               MPI_Type_commit(&(neighbourSendType[i]));
         }


      }

      /*! Sets the data pointer to the given vector
       * \param data pointer to the data vector
       */
      void setData(T *data) {
         this->data = data;
      }

      /*! Finalize instead of destructor, as the MPI calls fail after the main program called MPI_Finalize().
       *  Cleans up the cartesian communicator
       */
      void finalize() {
         for(int i=0;i<27;i++){
            if(neighbourReceiveType[i] != MPI_DATATYPE_NULL)
               MPI_Type_free(&(neighbourReceiveType[i]));
            if(neighbourSendType[i] != MPI_DATATYPE_NULL)
               MPI_Type_free(&(neighbourSendType[i]));
         }
         MPI_Comm_free(&comm3d);
      }

      /*! Returns the task responsible, and its localID for handling the cell with the given GlobalID
       * \param id GlobalID of the cell for which task is to be determined
       * \return a task for the grid's cartesian communicator
       */
      std::pair<int,LocalID> getTaskForGlobalID(GlobalID id) {
         // Transform globalID to global cell coordinate
         std::array<int, 3> cell = globalIDtoCellCoord(id);

         // Find the index in the task grid this Cell belongs to
         std::array<int, 3> taskIndex;
         for(int i=0; i<3; i++) {
            int n_per_task = globalSize[i] / ntasks[i];
            int remainder = globalSize[i] % ntasks[i];

            if(cell[i] < remainder * (n_per_task+1)) {
               taskIndex[i] = cell[i] / (n_per_task + 1);
            } else {
               taskIndex[i] = remainder + (cell[i] - remainder*(n_per_task+1)) / n_per_task;
            }
         }

         // Get the task number from the communicator
         std::pair<int,LocalID> retVal;
         int status = MPI_Cart_rank(comm3d, taskIndex.data(), &retVal.first);
         if(status != MPI_SUCCESS) {
            std::cerr << "Unable to find FsGrid rank for global ID " << id << " (coordinates [";
            for(int i=0; i<3; i++) {
               std::cerr << cell[i] << ", ";
            }
            std::cerr << "]" << std::endl;
            return std::pair<int,LocalID>(MPI_PROC_NULL,0);
         }

         // Determine localID of that cell within the target task
         std::array<int, 3> thatTasksStart;
         std::array<int, 3> thatTaskStorageSize;
         for(int i=0; i<3; i++) {
            thatTasksStart[i] = calcLocalStart(globalSize[i], ntasks[i], taskIndex[i]);
            thatTaskStorageSize[i] = calcLocalSize(globalSize[i], ntasks[i], taskIndex[i]) + 2 * stencil;
         }

         retVal.second = 0;
         int stride = 1;
         for(int i=0; i<3; i++) {
            if(globalSize[i] <= 1) {
               // Collapsed dimension, doesn't contribute.
               retVal.second += 0;
            } else {
               retVal.second += stride*(cell[i] - thatTasksStart[i] + stencil);
               stride *= thatTaskStorageSize[i];
            }
         }

         return retVal;
      }

      /*! Transform global cell coordinates into the local domain.
       * If the coordinates are out of bounds, (-1,-1,-1) is returned.
       * \param x The cell's global x coordinate
       * \param y The cell's global y coordinate
       * \param z The cell's global z coordinate
       */
      std::array<int, 3> globalToLocal(int x, int y, int z) {
         std::array<int, 3> retval;
         retval[0] = x - localStart[0];
         retval[1] = y - localStart[1];
         retval[2] = z - localStart[2];

         if(retval[0] >= localSize[0] || retval[1] >= localSize[1] || retval[2] >= localSize[2]
               || retval[0] < 0 || retval[1] < 0 || retval[2] < 0) {
            return {-1,-1,-1};
         }

         return retval;
      }

      /*! Determine the cell's GlobalID from its local x,y,z coordinates
       * \param x The cell's task-local x coordinate
       * \param y The cell's task-local y coordinate
       * \param z The cell's task-local z coordinate
       */
      GlobalID GlobalIDForCoords(int x, int y, int z) {
         return x + localStart[0] + globalSize[0] * (y + localStart[1]) + globalSize[0] * globalSize[1] * (z + localStart[2]);
      }
      /*! Determine the cell's LocalID from its local x,y,z coordinates
       * \param x The cell's task-local x coordinate
       * \param y The cell's task-local y coordinate
       * \param z The cell's task-local z coordinate
       */
      ARCH_HOSTDEV LocalID LocalIDForCoords(int x, int y, int z) {
         LocalID index=0;
         if(globalSize[2] > 1) {
            index += storageSize[0]*storageSize[1]*(stencil+z);
         }
         if(globalSize[1] > 1) {
            index += storageSize[0]*(stencil+y);
         }
         if(globalSize[0] > 1 ) {
            index += stencil+x;
         }

         return index;
      }

      /*! Prepare for transfer of Rank information from a coupled grid.
       * Setup MPI Irecv requests to recieve data for each grid cell, and prepare
       * buffer space to hold the MPI requests of the sends.
       *
       * \param cellsToCouple How many cells are going to be coupled
       */
      void setupForGridCoupling(int cellsToCouple) {
         int status;
         // Make sure we have sufficient buffer space to store our mpi
         // requests. Here only for receives, sends are done with blocking routine.
         requests.resize(localSize[0]*localSize[1]*localSize[2]);
         numRequests=0;

         // If previous coupling information was present, remove it.
         for(uint i=0; i<coupling->externalRank.size(); i++) {
            coupling->externalRank[i] = MPI_PROC_NULL;
         }

         for(int z=0; z<localSize[2]; z++) {
            for(int y=0; y<localSize[1]; y++) {
               for(int x=0; x<localSize[0]; x++) {
                  // Calculate LocalID for this cell
                  LocalID thisCell = LocalIDForCoords(x,y,z);
                  assert(numRequests < requests.size());
                  assert(thisCell < coupling->externalRank.size());
                  status = MPI_Irecv(&coupling->externalRank[thisCell], 1, MPI_INT, MPI_ANY_SOURCE, thisCell, comm3d,
                        &requests[numRequests++]);
                  if(status != MPI_SUCCESS) {
                     std::cerr << "Error setting up MPI Irecv in FsGrid::setupForGridCoupling" << std::endl;
                  }
               }
            }
         }
      }

      /*! Set the MPI rank responsible for external communication with the given cell.
       * If, for example, a separate grid is used in another part of the code, this would be the rank
       * in that grid, responsible for the information in this cell.
       *
       * This only needs to be called if the association changes (such as: on startup, and in a load
       * balancing step)
       *
       * \param id Global cell ID to be addressed
       * \iparam cellRank Rank that owns the cell in the other external grid.
       */
      void setGridCoupling(GlobalID id, int cellRank) {
         // Determine Task and localID that this cell belongs to
         std::pair<int,LocalID> TaskLid = getTaskForGlobalID(id);
         int status;
         status = MPI_Send(&cellRank, 1, MPI_INT, TaskLid.first, TaskLid.second, comm3d);
         if(status != MPI_SUCCESS) {
            std::cerr << "Error setting up MPI Isend in FsGrid::setGridCoupling" << std::endl;
         }
      }

      /*! Called after setting up the transfers into or out of this grid.
       * Basically only does a MPI_Waitall for all requests.
       */
      void finishGridCoupling() {
         assert(numRequests == requests.size());
         MPI_Waitall(numRequests, requests.data(), MPI_STATUSES_IGNORE);
         numRequests=0;
         MPI_Barrier(comm3d);
      }

      /*! Prepare for transfer of grid cell data into this grid.
       * Setup MPI Irecv requests to recieve data for each grid cell, and prepare
       * buffer space to hold the MPI requests of the sends.
       *
       * \param cellsToSend How many cells are going to be sent by this task
       */
      void setupForTransferIn(int cellsToSend) {
         int status;
         // Make sure we have sufficient buffer space to store our mpi requests
         requests.resize(localSize[0]*localSize[1]*localSize[2] + cellsToSend);
         numRequests=0;

         for(int z=0; z<localSize[2]; z++) {
            for(int y=0; y<localSize[1]; y++) {
               for(int x=0; x<localSize[0]; x++) {
                  // Calculate LocalID for this cell
                  LocalID thisCell = LocalIDForCoords(x,y,z);
                  assert(numRequests < requests.size());
                  status = MPI_Irecv(get(thisCell), sizeof(T), MPI_BYTE, coupling->externalRank[thisCell],
                                     thisCell, comm3d, &requests[numRequests++]);
                  if(status != MPI_SUCCESS) {
                     std::cerr << "Error setting up MPI Irecv in FsGrid::setupForTransferIn" << std::endl;
                  }
               }
            }
         }
      }

      /*! Set grid cell wtih the given ID to the value specified. Note that this doesn't need to be a local cell.
       * Note that this function should be concurrently called by all tasks, to allow for point-to-point communication.
       * \param id Global cell ID to be filled
       * \iparam value New value to be assigned to it
       */
      void transferDataIn(GlobalID id, T* value) {

         // Determine Task and localID that this cell belongs to
         std::pair<int,LocalID> TaskLid = getTaskForGlobalID(id);

         // Build the MPI Isend request for this cell
         int status;
         assert(numRequests < requests.size());
         status = MPI_Isend(value, sizeof(T), MPI_BYTE, TaskLid.first, TaskLid.second, comm3d,
               &requests[numRequests++]);
         if(status != MPI_SUCCESS) {
            std::cerr << "Error setting up MPI Isend in FsGrid::transferDataIn" << std::endl;
         }
      }
      /*! Called after setting up the transfers into or out of this grid.
       * Basically only does a MPI_Waitall for all requests.
       */
      void finishTransfersIn() {
         assert(numRequests == requests.size());
         MPI_Waitall(numRequests,requests.data(),MPI_STATUSES_IGNORE);
         numRequests = 0;

      }


      /*! Prepare for transfer of grid cell data out of this grid.
       * Setup MPI Isend requests to send data for each grid cell, and prepare
       * buffer space to hold the MPI requests of the sends.
       *
       * \param cellsToSend How many cells are going to be sent by this task
       */
      void setupForTransferOut(int cellsToReceive) {
         // Make sure we have sufficient buffer space to store our mpi requests
         requests.resize(localSize[0]*localSize[1]*localSize[2] + cellsToReceive);
         numRequests=0;
      }

      /*! Get the value of the grid cell wtih the given ID. Note that this doesn't need to be a local cell.
       * \param id Global cell ID to be read
       * \iparam target Location that the result is to be stored in
       */
      void transferDataOut(GlobalID id, T* target) {

         // Determine Task and localID that this cell belongs to
         std::pair<int,LocalID> TaskLid = getTaskForGlobalID(id);

         // Build the MPI Irecv request for this cell
         int status;
         assert(numRequests < requests.size());
         status = MPI_Irecv(target, sizeof(T), MPI_BYTE, TaskLid.first, TaskLid.second, comm3d,
               &requests[numRequests++]);
         if(status != MPI_SUCCESS) {
            std::cerr << "Error setting up MPI Irecv in FsGrid::transferDataOut" << std::endl;
         }
      }

      /*! Finalize the transfer of data out of this grid, re-integrating it back
       *! into the locations where it came from.
       * This method assumes that transferDataOut() has been called for each cell beforehand,
       * so that appropriate recieve buffers exist in the target ranks.
       */
      void finishTransfersOut() {
         int status;
         for(int z=0; z<localSize[2]; z++) {
            for(int y=0; y<localSize[1]; y++) {
               for(int x=0; x<localSize[0]; x++) {
                  // Calculate LocalID for this cell
                  LocalID thisCell = LocalIDForCoords(x,y,z);
                  assert(numRequests < requests.size());
                  status = MPI_Isend(get(thisCell), sizeof(T), MPI_BYTE, coupling->externalRank[thisCell], thisCell, comm3d,
                        &requests[numRequests++]);
                  if(status != MPI_SUCCESS) {
                     std::cerr << "Error setting up MPI Isend in FsGrid::setupForTransferOut" << std::endl;
                  }
               }
            }
         }
         assert(numRequests == requests.size());
         MPI_Waitall(numRequests,requests.data(),MPI_STATUSES_IGNORE);
         numRequests = 0;
      }

      /*! Perform ghost cell communication.
       */
      void updateGhostCells() {
         //TODO, faster with simultaneous isends& ireceives?
         std::array<MPI_Request, 27> receiveRequests;
         std::array<MPI_Request, 27> sendRequests;

         for(int i = 0; i < 27; i++){
            receiveRequests[i] = MPI_REQUEST_NULL;
            sendRequests[i] = MPI_REQUEST_NULL;
         }


         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  int shiftId = (x+1) * 9 + (y + 1) * 3 + (z + 1);
                  int receiveId = (1 - x) * 9 + ( 1 - y) * 3 + ( 1 - z);
                  if(neighbour[receiveId] != MPI_PROC_NULL &&
                     neighbourSendType[shiftId] != MPI_DATATYPE_NULL) {
                     // MPI_Irecv(data.data(), 1, neighbourReceiveType[shiftId], neighbour[receiveId], shiftId, comm3d, &(receiveRequests[shiftId]));
                     MPI_Irecv(data, 1, neighbourReceiveType[shiftId], neighbour[receiveId], shiftId, comm3d, &(receiveRequests[shiftId]));
                  }
               }
            }
         }

         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  int shiftId = (x+1) * 9 + (y + 1) * 3 + (z + 1);
                  int sendId = shiftId;
                  if(neighbour[sendId] != MPI_PROC_NULL &&
                     neighbourSendType[shiftId] != MPI_DATATYPE_NULL) {
                     // MPI_Isend(data.data(), 1, neighbourSendType[shiftId], neighbour[sendId], shiftId, comm3d, &(sendRequests[shiftId]));
                     MPI_Isend(data, 1, neighbourSendType[shiftId], neighbour[sendId], shiftId, comm3d, &(sendRequests[shiftId]));
                  }
               }
            }
         }
         MPI_Waitall(27, receiveRequests.data(), MPI_STATUSES_IGNORE);
         MPI_Waitall(27, sendRequests.data(), MPI_STATUSES_IGNORE);
      }


      int32_t* getStorageSize() {
         return storageSize;
      }

      /*! Get the size of the local domain handled by this grid.
       */
      ARCH_HOSTDEV int32_t* getLocalSize() {
      // std::array<int32_t, 3>& getLocalSize() {
         return localSize;
      }

      /*! Get the sstart coordinates of the local domain handled by this grid.
       */
      ARCH_HOSTDEV int32_t* getLocalStart() {
         return localStart;
      }

      /*! Get global size of the fsgrid domain
       */
      ARCH_HOSTDEV int* getGlobalSize() {
         return globalSize;
      }

      /*! Calculate global cell position (XYZ in global cell space) from local cell coordinates.
       *
       * \param x x-Coordinate, in cells
       * \param y y-Coordinate, in cells
       * \param z z-Coordinate, in cells
       *
       * \return Global cell coordinates
       */
      ARCH_HOSTDEV void getGlobalIndices(int x, int y, int z, int32_t (&retval)[3]) {
         retval[0] = localStart[0] + x;
         retval[1] = localStart[1] + y;
         retval[2] = localStart[2] + z;
      }

      /*! Get a reference to the field data in a cell
       * \param x x-Coordinate, in cells
       * \param y y-Coordinate, in cells
       * \param z z-Coordinate, in cells
       * \return A reference to cell data in the given cell
       */
      ARCH_HOSTDEV T* get(int x, int y, int z) {

         // Keep track which neighbour this cell actually belongs to (13 = ourself)
         int isInNeighbourDomain=13;
         int coord_shift[3] = {0,0,0};
         if(x < 0) {
            isInNeighbourDomain -= 9;
            coord_shift[0] = 1;
         }
         if(x >= localSize[0]) {
            isInNeighbourDomain += 9;
            coord_shift[0] = -1;
         }
         if(y < 0) {
            isInNeighbourDomain -= 3;
            coord_shift[1] = 1;
         }
         if(y >= localSize[1]) {
            isInNeighbourDomain += 3;
            coord_shift[1] = -1;
         }
         if(z < 0) {
            isInNeighbourDomain -= 1;
            coord_shift[2] = 1;
         }
         if(z >= localSize[2]) {
            isInNeighbourDomain += 1;
            coord_shift[2] = -1;
         }

         // Santiy-Check that the requested cell is actually inside our domain
         // TODO: ugh, this is ugly.
#ifdef FSGRID_DEBUG
         bool inside=true;
         if(localSize[0] <= 1 && !periodic[0]) {
            if(x != 0) {
               std::cerr << "x != 0 despite non-periodic x-axis with only one cell." << std::endl;
               inside = false;
            }
         } else {
            if(x < -stencil || x >= localSize[0] + stencil) {
               std::cerr << "x = " << x << " is outside of [ " << -stencil <<
                  ", " << localSize[0] + stencil << "[!" << std::endl;
               inside = false;
            }
         }

         if(localSize[1] <= 1 && !periodic[1]) {
            if(y != 0) {
               std::cerr << "y != 0 despite non-periodic y-axis with only one cell." << std::endl;
               inside = false;
            }
         } else {
            if(y < -stencil || y >= localSize[1] + stencil) {
               std::cerr << "y = " << y << " is outside of [ " << -stencil <<
                  ", " << localSize[1] + stencil << "[!" << std::endl;
               inside = false;
            }
         }

         if(localSize[2] <= 1 && !periodic[2]) {
            if(z != 0) {
               std::cerr << "z != 0 despite non-periodic z-axis with only one cell." << std::endl;
               inside = false;
            }
         } else {
            if(z < -stencil || z >= localSize[2] + stencil) {
               inside = false;
               std::cerr << "z = " << z << " is outside of [ " << -stencil <<
                  ", " << localSize[2] + stencil << "[!" << std::endl;
            }
         }
         if(!inside) {
            std::cerr << "Out-of bounds access in FsGrid::get! Expect weirdness." << std::endl;
            return NULL;
         }
#endif // FSGRID_DEBUG

         if(isInNeighbourDomain != 13) {

            // Check if the corresponding neighbour exists
            if(neighbour[isInNeighbourDomain]==MPI_PROC_NULL) {
               // Neighbour doesn't exist, we must be an outer boundary cell
               // (or something is quite wrong)
               return NULL;
            } else if(neighbour[isInNeighbourDomain] == rank) {
               // For periodic boundaries, where the neighbour is actually ourself,
               // return our own actual cell instead of the ghost
               x += coord_shift[0] * localSize[0];
               y += coord_shift[1] * localSize[1];
               z += coord_shift[2] * localSize[2];
            }
            // Otherwise we return the ghost cell
         }
         LocalID index = LocalIDForCoords(x,y,z);

         return &data[index * TDim];
      }

      ARCH_HOSTDEV T* get(LocalID id) {
         if(id < 0 || (unsigned int)id > totalStorageSize) {
            #ifndef __CUDA_ARCH__
               std::cerr << "Out-of-bounds access in FsGrid::get!" << std::endl
                  << "(LocalID = " << id << ", but storage space is " << totalStorageSize
                  << ". Expect weirdness." << std::endl;
            #endif
            return NULL;
         }
         return &data[id * TDim];
      }

      ARCH_HOSTDEV T& get(int x, int y, int z, int i) {
         return (get(x,y,z)[i]);
      }

      ARCH_HOSTDEV T& getData(int i=0) const {
            return data[i];
      }

      /*! Physical grid spacing and physical coordinate space start.
       * TODO: Should this be private and have accesor-functions?
       */
      double DX,DY,DZ;
      std::array<double,3> physicalGlobalStart;

      /*! Get the physical coordinates in the global simulation space for
       * the given cell.
       *
       * \param x local x-Coordinate, in cells
       * \param y local y-Coordinate, in cells
       * \param z local z-Coordinate, in cells
       */
      std::array<double, 3> getPhysicalCoords(int x, int y, int z) {
         std::array<double, 3> coords;
         coords[0] = physicalGlobalStart[0] + (localStart[0]+x)*DX;
         coords[1] = physicalGlobalStart[1] + (localStart[1]+y)*DY;
         coords[2] = physicalGlobalStart[2] + (localStart[2]+z)*DZ;

         return coords;
      }

      /*! Debugging output helper function. Allows for nicely formatted printing
       * of grid contents. Since the grid data format is varying, the actual
       * printing should be done in a lambda passed to this function. Example usage
       * to plot |B|:
       *
       * perBGrid.debugOutput([](const std::array<Real, fsgrids::bfield::N_BFIELD>& a)->void{
       *     cerr << sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) << ", ";
       * });
       *
       * \param func Function pointer (or lambda) which is called with a cell reference,
       * in order. Use std::cerr in it to print desired value.
       */
      void debugOutput(void (func)(const T&)) {
         int xmin=0,xmax=1;
         int ymin=0,ymax=1;
         int zmin=0,zmax=1;
         if(localSize[0] > 1) {
            xmin = -stencil; xmax = localSize[0]+stencil;
         }
         if(localSize[1] > 1) {
            ymin = -stencil; ymax = localSize[1]+stencil;
         }
         if(localSize[2] > 1) {
            zmin = -stencil; zmax = localSize[2]+stencil;
         }
         for(int z=zmin; z<zmax; z++) {
            for(int y=ymin; y<ymax; y++) {
               for(int x=xmin; x<xmax; x++) {
                  func(*get(x,y,z));
               }
               std::cerr << std::endl;
            }
            std::cerr << " - - - - - - - - " << std::endl;
         }
      }

      /*! Get the rank of this CPU in the FsGrid communicator */
      int getRank() {
        return rank;
      }

      /*! Get in which directions, if any, this grid is periodic */
      std::array<bool, 3>& getPeriodic() {
        return periodic;
      }

      /*! Perform an MPI_Allreduce with this grid's internal communicator
       * Function syntax is identical to MPI_Allreduce, except the final (communicator
       * argument will not be needed) */
      int Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op) {
         return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm3d);
      }

   private:
      //! MPI Cartesian communicator used in this grid
      MPI_Comm comm3d;
      int rank; //!< This task's rank in the communicator
      std::vector<MPI_Request> requests;
      uint numRequests;

      int neighbour[27]; //!< Tasks of the 26 neighbours (plus ourselves)
      std::vector<char> neighbour_index; //!< Lookup table from rank to index in the neighbour array

      // We have, fundamentally, two different coordinate systems we're dealing with:
      // 1) Task grid in the MPI_Cartcomm
      std::array<int, 3> ntasks; //!< Number of tasks in each direction
      std::array<int, 3> taskPosition; //!< This task's position in the 3d task grid
      // 2) Cell numbers in global and local view

      std::array<bool, 3> periodic; //!< Information about whether a given direction is periodic
      int32_t totalStorageSize; //!< Total number of cells in the local storage, including ghost cells
      int32_t globalSize[3]; //!< Global size of the simulation space, in cells
      int32_t localSize[3];
      // std::array<int32_t, 3> localSize;  //!< Local size of simulation space handled by this task (without ghost cells)
      int32_t storageSize[3];  //!< Local size of simulation space handled by this task (including ghost cells)
      int32_t localStart[3]; //!< Offset of the local coordinate system against the global one

      FsGridCouplingInformation* coupling; // Information required to couple to external grids

      std::array<MPI_Datatype, 27> neighbourSendType; //!< Datatype for sending data
      std::array<MPI_Datatype, 27> neighbourReceiveType; //!< Datatype for receiving data



      //! Actual storage of field data
      // std::vector<T>* data;


      //! Helper function: given a global cellID, calculate the global cell coordinate from it.
      // This is then used do determine the task responsible for this cell, and the
      // local cell index in it.
      std::array<int, 3> globalIDtoCellCoord(GlobalID id) {

         // Transform globalID to global cell coordinate
         std::array<int, 3> cell;

         assert(id >= 0);
         assert(id < globalSize[0] * globalSize[1] * globalSize[2]);

         int stride=1;
         for(int i=0; i<3; i++) {
            cell[i] = (id / stride) % globalSize[i];
            stride *= globalSize[i];
         }

         return cell;
      }

      void swapArray(int array[3]) {
         int a = array[0];
         array[0] = array[2];
         array[2] = a;
      }
};
