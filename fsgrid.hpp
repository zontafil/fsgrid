/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * */
#include <array>
#include <stdint.h>

/**
 * Simple cartesian, non-loadbalancing MPI Grid for use with the fieldsolver
 *
 * @Type: datastructure containing the field in each cell which this grid manages
 * @stencil: ghost cell width of this grid
 */
template <typename T, int dims, int stencil> class FsGrid {
   
   std::array<uint32_t,dims> globalSize; // (this is identical to Dccrg's grid size anyway)
   std::array<uint32_t,dims> localSize;

   MPI_Cartcomm comm3d;

   // Number of tasks in each direction
   std::array<int,dims> nproc;


   std::vector<T> data;

public:

   typedef uint64_t LocalID;
   typedef uint64_t GlobalID;

   FsGrid(std::array<uint32_t,dims> globalSize, MPI_Comm parent_comm) : globalSize(globalSize) {
   
      int status;

      // Tell MPI how many Tasks we want in which direction.
      // By giving a value of 0, we let MPI choose this number itself.
      std::array<int, dims> periodic;
      for(int i=0; i<dims; i++) {
         nproc[i] = 0;
         periodic[i] = 0; //TODO: Determine boundary here?
      }

      // Create cartesian communicator
      status = MPI_Cart_create(dims, dims, nproc, periodic, 0, &comm3d);

      if(status != MPI_SUCCESS) {
         // TODO: Error handling!
      }



      // Allocate local storage array

   }

   ~FsGrid() {
      MPI_Comm_free(&communicator);
   }


   int getTaskForGlobalID(GlobalID id);
   LocalID getLocalIdForGlobalID(GlobalID id);

   void setFieldData(GlobalID id, T& const value);

   /**
    * Perform ghost cell communication.
    */
   void updateGhostCells();

   /**
    * Get a reference to the field data in a cell
    */
   T& get(int x, int y, int z);
   T& const get(int x, int y, int z) const;
};
