/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * */
#include <array>
#include <vector>
#include <mpi.h>
#include <stdint.h>

/**
 * Simple cartesian, non-loadbalancing MPI Grid for use with the fieldsolver
 *
 * @Type: datastructure containing the field in each cell which this grid manages
 * @stencil: ghost cell width of this grid
 */
template <typename T, int dims, int stencil> class FsGrid {

   public:

      typedef uint64_t LocalID;
      typedef uint64_t GlobalID;

      /**
       * Constructor for this grid, given:
       * @globalSize: Cell size of the global simulation domain
       * @MPI_Copp: The MPI communicator this grid should use
       */
      FsGrid(std::array<uint32_t,dims> globalSize, MPI_Comm parent_comm) : globalSize(globalSize) {

         int status;

         // Tell MPI how many Tasks we want in which direction.
         // By giving a value of 0, we let MPI choose this number itself.
         std::array<int, dims> periodic;
         for(int i=0; i<dims; i++) {
            ntasks[i] = 0;
            periodic[i] = 0; //TODO: Determine boundary here?
         }

         // Create cartesian communicator
         status = MPI_Cart_create(parent_comm, dims, ntasks.data(), periodic.data(), 0, &comm3d);

         if(status != MPI_SUCCESS) {
            // TODO: Error handling!
         }

         // Determine our position in the resulting task-grid
         int rank;
         status = MPI_Comm_rank(comm3d, &rank);
         if(status != MPI_SUCCESS) {
            // TODO: Error handling!
         }
         status = MPI_Cart_coords(comm3d, rank, dims, taskPosition.data());
         if(status != MPI_SUCCESS) {
            // TODO: Error handling!
         }

         // Determine size of our local grid
         for(int i=0; i<dims; i++) {
            localSize[i] = calcLocalSize(globalSize[i],ntasks[i], taskPosition[i]);
         }

         // Allocate local storage array
         size_t storageSize=1;
         for(auto i : localSize) {
            storageSize *= i;
         }
         data.resize(storageSize);
      }

      ~FsGrid() {
         MPI_Comm_free(&comm3d);
      }


      int getTaskForGlobalID(GlobalID id);
      LocalID getLocalIdForGlobalID(GlobalID id);

      void setFieldData(GlobalID id, T& value);

      /**
       * Perform ghost cell communication.
       */
      void updateGhostCells();

      /**
       * Get number of cells handeled by the local task in each spatial dimension.
       */
      std::array<int, dims> getLocalSize();

      /**
       * Get a reference to the field data in a cell
       */
      T& get(int x, int y, int z);
      const T& get(int x, int y, int z) const;

   private: 
      std::array<uint32_t,dims> globalSize;
      std::array<uint32_t,dims> localSize;

      MPI_Comm comm3d;

      // Number of tasks in each direction
      std::array<int,dims> ntasks;

      // This tasks position in the 3d task grid
      std::array<int, dims> taskPosition;

      // Actual storage of field data
      std::vector<T> data;


      // Helper functions: calculate position and size of the local coordinate space for the given dimension
      //TODO: Inverse of these, to get task for a given position
      uint32_t calcLocalStart(uint32_t globalCells, int ntasks, int my_n) {
         int n_per_task = globalCells / ntasks;
         int remainder = globalCells % ntasks;

         if(my_n < remainder) {
            return my_n * (n_per_task+1);
         } else {
            return my_n * n_per_task + remainder;
         }
      }
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
