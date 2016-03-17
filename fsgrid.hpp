/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * */
#include <array>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <limits>
#include <stdint.h>

static const int fsgrid_dims = 3;

/*! Simple cartesian, non-loadbalancing MPI Grid for use with the fieldsolver
 *
 * \param T datastructure containing the field in each cell which this grid manages
 * \param fsgrid_dims number of dimensions this field has
 * \param stencil ghost cell width of this grid
 */
template <typename T, int stencil> class FsGrid {

   public:

      typedef uint64_t LocalID;
      typedef uint64_t GlobalID;

      /*! Constructor for this grid.
       * \param globalSize Cell size of the global simulation domain.
       * \param MPI_Comm The MPI communicator this grid should use.
       * \param isPeriodic An array specifying, for each dimension, whether it is to be treated as periodic.
       */
      FsGrid(std::array<int32_t,fsgrid_dims> globalSize, MPI_Comm parent_comm, std::array<int,fsgrid_dims> isPeriodic)
            : globalSize(globalSize) {

         int status;
         int size;
         status = MPI_Comm_size(parent_comm, &size);

         // Tell MPI how many Tasks we want in which direction.
         // By giving a value of 0, we let MPI choose this number itself.
         //for(int i=0; i<fsgrid_dims; i++) {
         //   if(globalSize[i] <= 1) {
         //      ntasks[i] = 1;
         //   } else {
         //      ntasks[i] = 0;
         //   }
         //}
         //status = MPI_Dims_create(size,fsgrid_dims,ntasks.data());
         computeDomainDecomposition(globalSize, size, ntasks);


         // Create cartesian communicator
         status = MPI_Cart_create(parent_comm, fsgrid_dims, ntasks.data(), isPeriodic.data(), 0, &comm3d);
         if(status != MPI_SUCCESS) {
            std::cerr << "Creating cartesian communicatior failed when attempting to create FsGrid!" << std::endl;
            return;
         }

         status = MPI_Comm_rank(comm3d, &rank);
         if(status != MPI_SUCCESS) {
            std::cerr << "Getting rank failed when attempting to create FsGrid!" << std::endl;

            // Without a rank, there's really not much we can do. Just return an uninitialized grid
            // (I suppose we'll crash after this, anyway)
            return;
         }


         // Determine our position in the resulting task-grid
         status = MPI_Cart_coords(comm3d, rank, fsgrid_dims, taskPosition.data());
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
                  //TODO: This has 3 dimensions hardcoded!
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
                        for(int i=0; i<fsgrid_dims; i++) {
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
         std::array<int, 3> localStart;
         for(int i=0; i<fsgrid_dims; i++) {
            localSize[i] = calcLocalSize(globalSize[i],ntasks[i], taskPosition[i]);
            localStart[i] = calcLocalStart(globalSize[i],ntasks[i], taskPosition[i]);
         }

         //std::cerr << "Rank " << rank << " here, my space starts at [" << localStart[0] << ", " << localStart[1] << ", "
         //   << localStart[2] << "] and ends at [" << (localStart[0] + localSize[0]) << ", " <<
         //   (localStart[1]+localSize[1]) << ", " << (localStart[2]+localSize[2]) << "]" << std::endl;

         // Allocate local storage array
         size_t totalStorageSize=1;
         for(int i=0; i<fsgrid_dims; i++) {
            if(globalSize[i] <= 1) {
               // Collapsed dimension => only one cell thick
               storageSize[i] = 1;
            } else {
               // Size of the local domain + 2* size for the ghost cell stencil
               storageSize[i] = (localSize[i] + stencil*2);
            }
            totalStorageSize *= storageSize[i];
         }
         data.resize(totalStorageSize);

         MPI_Datatype mpiTypeT;
         MPI_Type_contiguous(sizeof(T), MPI_BYTE, &mpiTypeT);
         // Compute send and receive datatypes
         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  std::array<int,3> subarraySize;
                  std::array<int,3> subarrayStart;
                  
                  if((storageSize[0] == 1 && x!= 0 ) ||
                     (storageSize[1] == 1 && y!= 0 ) ||
                     (storageSize[2] == 1 && z!= 0 )){
                     //check for 2 or 1D simulations
                     neighbourSendType[(x+1) * 9 + (y + 1) * 3 + (z + 1)] = MPI_DATATYPE_NULL;
                     neighbourReceiveType[(x+1) * 9 + (y + 1) * 3 + (z + 1)] = MPI_DATATYPE_NULL;
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
                  /*
                  printf("create datatype for %d, %d, %d:\n %d %d %d\n %d %d %d\n %d %d %d\n", 
                         x, y, z,
                         storageSize[0], storageSize[1], storageSize[2], 
                         subarraySize[0], subarraySize[1], subarraySize[2], 
                         subarrayStart[0], subarrayStart[1], subarrayStart[2]);
                  */
                         
                  MPI_Type_create_subarray(3,
                                           storageSize.data(),
                                           subarraySize.data(),
                                           subarrayStart.data(),
                                           MPI_ORDER_C,
                                           mpiTypeT,
                                           &(neighbourSendType[(x+1) * 9 + (y + 1) * 3 + (z + 1)]) );
                  
                  if(x == -1 )
                     subarrayStart[0] = 0;
                  else if (x == 0)
                     subarrayStart[0] = stencil;
                  else if (x == 1)
                     subarrayStart[0] = storageSize[0] -  stencil;
                  if(y == -1 )
                     subarrayStart[1] = 0;
                  else if (y == 0)
                     subarrayStart[1] = stencil;
                  else if (y == 1)
                     subarrayStart[1] = storageSize[1] -  stencil;
                  if(z == -1 )
                     subarrayStart[2] = 0;
                  else if (z == 0)
                     subarrayStart[2] = stencil;
                  else if (z == 1)
                     subarrayStart[2] = storageSize[2] -  stencil;
                  for(int i = 0;i < 3; i++)
                     if(storageSize[i] == 1) 
                        subarrayStart[i] = 0;
                  
                  MPI_Type_create_subarray(3,
                                           storageSize.data(),
                                           subarraySize.data(),
                                           subarrayStart.data(),
                                           MPI_ORDER_C,
                                           mpiTypeT,
                                           &(neighbourReceiveType[(x+1)*9+(y+1)*3+(z+1)]) );
                  
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

      /*! Destructor, cleans up the cartesian communicator
       */
      ~FsGrid() {
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
         std::array<int, fsgrid_dims> cell = globalIDtoCellCoord(id);

         // Find the index in the task grid this Cell belongs to
         std::array<int, fsgrid_dims> taskIndex;
         for(int i=0; i<fsgrid_dims; i++) {
            int n_per_task = globalSize[i]/ntasks[i];
            int remainder = globalSize[i]%ntasks[i];

            if(cell[i] < remainder*(n_per_task+1)) {
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
            for(int i=0; i<fsgrid_dims; i++) {
               std::cerr << cell[i] << ", ";
            }
            std::cerr << "]" << std::endl;
            return std::pair<int,LocalID>(MPI_PROC_NULL,0);
         }

         // Determine localID of that cell within the target task
         std::array<int, fsgrid_dims> thatTasksStart;
         for(int i=0; i<fsgrid_dims; i++) {
            thatTasksStart[i] = calcLocalStart(globalSize[i], ntasks[i], taskIndex[i]);
         }

         retVal.second = 0;
         int stride = 1;
         for(int i=0; i<fsgrid_dims; i++) {
            if(globalSize[i] <= 1) {
               // Collapsed dimension, doesn't contribute.
               retVal.second += 0;
            } else {
               retVal.second += stride*(cell[i] - thatTasksStart[i] + stencil);
               stride *= storageSize[i];
            }
         }

         return retVal;
      }

      /*! Set grid cell wtih the given ID to the value specified. Note that this doesn't need to be a local cell.
       * Note that this function should be concurrently called by all tasks, to allow for point-to-point communication.
       * \param id Global cell ID to be filled
       * \iparam value New value to be assigned to it
       *
       * TODO: Shouldn't this maybe rather be part of the external grid-glue?
       */
      void setFieldData(GlobalID id, T& value) {

         // Determine Task and localID that this cell belongs to
         std::pair<int,LocalID> TaskLid = getTaskForGlobalID(id);

         if(TaskLid.first == rank) {
            // This is our own cell, nice!
         } else {
         }
      }

      /*! Perform ghost cell communication.
       */
      void updateGhostCells() {
         //TODO, faster with simultaneous isends& ireceives?
         for(int i = 0; i < 27; i++) {
            MPI_Sendrecv(data.data(), 1, neighbourSendType[i], neighbour[i], i, 
                         data.data(), 1, neighbourReceiveType[i], neighbour[i], i,
                         comm3d, MPI_STATUS_IGNORE);
         }
      }

      /*! Get the size of the local domain handled by this grid.
       */
      std::array<int, fsgrid_dims>& getLocalSize() {
         return localSize;
      }

      /*! Get a reference to the field data in a cell
       * \param x x-Coordinate, in cells
       * \param y y-Coordinate, in cells
       * \param z z-Coordinate, in cells
       * \return A reference to cell data in the given cell
       */
      T& get(int x, int y, int z) {
         //TODO: remove this sanity check once everything works
         if(x < -stencil || x > localSize[0] + stencil
               || y < -stencil || y > localSize[1] + stencil
               || z < -stencil || z > localSize[2] + stencil) {
            std::cerr << "Out-of bounds access in FsGrid::get! Expect weirdness." << std::endl;
         }

         int index=storageSize[0]*storageSize[1]*z
            + storageSize[0]*y
            + x;

         return data[index];
      }
      const T& get(int x, int y, int z) const;

      T& get(LocalID id) {
         if(id < 0 || id > data.get_size()) {
            std::cerr << "Out-of-bounds access in FsGrid::get! Expect weirdness." << std::endl;
         }
         return data[id];
      }

   private:
      //! MPI Cartesian communicator used in this grid
      MPI_Comm comm3d;
      int rank;

      std::array<int, 27> neighbour; //!< Tasks of the 26 neighbours (plus ourselves)
      std::vector<char> neighbour_index; //!< Lookup table from rank to index in the neighbour array

      // We have, fundamentally, two different coordinate systems we're dealing with:
      // 1) Task grid in the MPI_Cartcomm
      std::array<int,fsgrid_dims> ntasks; //!< Number of tasks in each direction
      std::array<int, fsgrid_dims> taskPosition; //!< This task's position in the 3d task grid
      // 2) Cell numbers in global and local view
      std::array<int32_t,fsgrid_dims> globalSize; //!< Global size of the simulation space, in cells
      std::array<int32_t,fsgrid_dims> localSize;  //!< Local size of simulation space handled by this task (without ghost cells)
      std::array<int32_t,fsgrid_dims> storageSize;  //!< Local size of simulation space handled by this task (including ghost cells)

      std::array<MPI_Datatype, 27> neighbourSendType; //!< Datatype for sending data
      std::array<MPI_Datatype, 27> neighbourReceiveType; //!< Datatype for receiving data


      //! Actual storage of field data
      std::vector<T> data;

      //! Helper function: given a global cellID, calculate the global cell coordinate from it.
      // This is then used do determine the task responsible for this cell, and the
      // local cell index in it.
      std::array<int, fsgrid_dims> globalIDtoCellCoord(GlobalID id) {

         // Transform globalID to global cell coordinate
         std::array<int, fsgrid_dims> cell;

         int stride=1;
         for(int i=0; i<fsgrid_dims; i++) {
            cell[i] = (id / stride) % globalSize[i];
            stride *= globalSize[i];
         }

         return cell;
      }

      //! Helper function to optimize decomposition of this grid over the given number of tasks
      void computeDomainDecomposition(const std::array<int32_t, 3>& GlobalSize, int nProcs, 
            std::array<int,3>& processDomainDecomposition) {

         std::array<double, 3> systemDim;
         std::array<double, 3 > processBox;
         double optimValue = std::numeric_limits<double>::max();

         processDomainDecomposition = {1,1,1};

         for(int i = 0; i < 3; i++) {
            systemDim[i] = (double)GlobalSize[i];
         }

         for (int i = 1; i<= nProcs; i++) {
            if( i  > systemDim[0])
               continue;
            processBox[0] = std::max(systemDim[0]/i, 1.0);

            for (int j = 1; j<= nProcs; j++) {
               if( i * j  >nProcs || j > systemDim[1])
                  continue;

               processBox[1] = std::max(systemDim[1]/j, 1.0);
               for (int k = 1; k<= nProcs; k++) {
                  if( i * j * k >nProcs || k > systemDim[2])
                     continue;
                  processBox[2] = std::max(systemDim[2]/k, 1.0);
                  double value = 
                     100 * processBox[0] * processBox[1] * processBox[2] + 
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
      }


      //! Helper function: calculate position of the local coordinate space for the given dimension
      //TODO: Inverse of these, to get task for a given position
      // \param globalCells Number of cells in the global Simulation, in this dimension
      // \param ntasks Total number of tasks in this dimension
      // \param my_n This task's position in this dimension
      // \return Cell number at which this task's domains cells start (actual cells, not counting ghost cells)
      int32_t calcLocalStart(int32_t globalCells, int ntasks, int my_n) {
         int n_per_task = globalCells / ntasks;
         int remainder = globalCells % ntasks;

         if(my_n < remainder) {
            return my_n * (n_per_task+1);
         } else {
            return my_n * n_per_task + remainder;
         }
      }
      //! Helper function: calculate size of the local coordinate space for the given dimension
      //TODO: Inverse of these, to get task for a given position
      // \param globalCells Number of cells in the global Simulation, in this dimension
      // \param ntasks Total number of tasks in this dimension
      // \param my_n This task's position in this dimension
      // \return Nmuber of cells for this task's local domain (actual cells, not counting ghost cells)
      int32_t calcLocalSize(int32_t globalCells, int ntasks, int my_n) {
         int n_per_task = globalCells/ntasks;
         int remainder = globalCells%ntasks;
         if(my_n < remainder) {
            return n_per_task+1;
         } else {
            return n_per_task;
         }
      }

};
