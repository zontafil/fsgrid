/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * */
#include <array>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <stdint.h>

/*! Simple cartesian, non-loadbalancing MPI Grid for use with the fieldsolver
 *
 * \param T datastructure containing the field in each cell which this grid manages
 * \param dims number of dimensions this field has
 * \param stencil ghost cell width of this grid
 */
template <typename T, int dims, int stencil> class FsGrid {

   public:

      typedef uint64_t LocalID;
      typedef uint64_t GlobalID;

      /*! Constructor for this grid.
       * \param globalSize Cell size of the global simulation domain.
       * \param MPI_Comm The MPI communicator this grid should use.
       * \param isPeriodic An array specifying, for each dimension, whether it is to be treated as periodic.
       */
      FsGrid(std::array<uint32_t,dims> globalSize, MPI_Comm parent_comm, std::array<int,dims> isPeriodic)
            : globalSize(globalSize) {

         int status;
         int size;
         status = MPI_Comm_size(parent_comm, &size);

         // Tell MPI how many Tasks we want in which direction.
         // By giving a value of 0, we let MPI choose this number itself.
         for(int i=0; i<dims; i++) {
            ntasks[i] = 0;
            if(ntasks.data()[i] != 0) {
               std::cerr << "Everything is horrible and all you believe is wrong!" << std::endl;
            }
         }
         status = MPI_Dims_create(size,dims,ntasks.data());

         // Create cartesian communicator
         status = MPI_Cart_create(parent_comm, dims, ntasks.data(), isPeriodic.data(), 0, &comm3d);
         if(status != MPI_SUCCESS) {
            std::cerr << "Creating cartesian communicatior failed when attempting to create FsGrid!" << std::endl;
            return;
         }

         int rank;
         status = MPI_Comm_rank(comm3d, &rank);
         if(status != MPI_SUCCESS) {
            std::cerr << "Getting rank failed when attempting to create FsGrid!" << std::endl;

            // Without a rank, there's really not much we can do. Just return an uninitialized grid
            // (I suppose we'll crash after this, anyway)
            return;
         }


         // Determine our position in the resulting task-grid
         status = MPI_Cart_coords(comm3d, rank, dims, taskPosition.data());
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
                        for(int i=0; i<dims; i++) {
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
         for(int i=0; i<dims; i++) {
            localSize[i] = calcLocalSize(globalSize[i],ntasks[i], taskPosition[i]);
         }

         // Allocate local storage array
         size_t storageSize=1;
         for(auto i : localSize) {
            // Size of the local domain + 2* size for the ghost cell stencil
            storageSize *= (i + stencil*2);
         }
         data.resize(storageSize);

      }

      /*! Destructor, cleans up the cartesian communicator
       */
      ~FsGrid() {
         MPI_Comm_free(&comm3d);
      }

      /*! Returns the task responsible for handling the cell with the given GlobalID
       * \param id GlobalID of the cell for which task is to be determined
       * \return a task for the grid's cartesian communicator
       */
      int getTaskForGlobalID(GlobalID id) {
         // Transform globalID to global cell coordinate
         std::array<int, dims> globalCell;

         int stride=1;
         for(int i=0; i<dims; i++) {
            globalCell[i] = (id / stride) % globalSize[i];
            stride *= globalSize[i];
         }

         // Find the CPU this Cell belongs to
         // TODO: continue here
         return 0;
      }
      /*! Get the localID (=offset into the data-array) for the given globalID
       * \param id the global cellID of a cell.
       */
      LocalID getLocalIdForGlobalID(GlobalID id);

      void setFieldData(GlobalID id, T& value);

      /*! Perform ghost cell communication.
       */
      void updateGhostCells() {


      }

      /*! Get the size of the local domain handled by this grid.
       */
      std::array<int, dims> getLocalSize();

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

         int index=(localSize[0]+2*stencil)*(localSize[1]+2*stencil)*z
            + (localSize[0]+2*stencil)*y
            + x;

         return data[index];
      }
      const T& get(int x, int y, int z) const;

   private:
      //! MPI Cartesian communicator used in this grid
      MPI_Comm comm3d;

      std::array<int, 27> neighbour; //!< IDs of the 26 neighbours (plus ourselves)
      std::vector<char> neighbour_index; //!< Lookup table from rank to index in the neighbour array

      // We have, fundamentally, two different coordinate systems we're dealing with:
      // 1) Task grid in the MPI_Cartcomm
      std::array<int,dims> ntasks; //!< Number of tasks in each direction
      std::array<int, dims> taskPosition; //!< This task's position in the 3d task grid
      // 2) Cell numbers in global and local view
      std::array<uint32_t,dims> globalSize; //!< Global size of the simulation space, in cells
      std::array<uint32_t,dims> localSize;  //!< Local size of simulation space handled by this task

      //! Actual storage of field data
      std::vector<T> data;

      //! Helper function: calculate position of the local coordinate space for the given dimension
      //TODO: Inverse of these, to get task for a given position
      // \param globalCells Number of cells in the global Simulation, in this dimension
      // \param ntasks Total number of tasks in this dimension
      // \param my_n This task's position in this dimension
      // \return Cell number at which this task's domains cells start (actual cells, not counting ghost cells)
      uint32_t calcLocalStart(uint32_t globalCells, int ntasks, int my_n) {
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
      uint32_t calcLocalSize(uint32_t globalCells, int ntasks, int my_n) {
         int n_per_task = globalCells/ntasks;
         int remainder = globalCells%ntasks;
         if(my_n < remainder) {
            return n_per_task+1;
         } else {
            return n_per_task;
         }
      }

};
