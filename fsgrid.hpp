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
       * \param globalSize Cell size of the global simulation domain
       * \param MPI_Comm The MPI communicator this grid should use
       */
      FsGrid(std::array<uint32_t,dims> globalSize, MPI_Comm parent_comm) : globalSize(globalSize) {

         int status;

         int rank;
         status = MPI_Comm_rank(comm3d, &rank);
         if(status != MPI_SUCCESS) {
            std::cerr << "Getting rank failed when attempting to create FsGrid!" << std::endl;

            // Without a rank, there's really not much we can do. Just return an uninitialized grid
            // (I suppose we'll crash after this, anyway)
            return;
         }

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
            if(rank == 0) {
               std::cerr << "Creating cartesian communicatior failed when attempting to create FsGrid!" << std::endl;
            }
            return;
         }

         // Determine our position in the resulting task-grid
         status = MPI_Cart_coords(comm3d, rank, dims, taskPosition.data());
         if(status != MPI_SUCCESS) {
            std::cerr << "Rank " << rank 
               << " unable to determine own position in cartesian communicator when attempting to create FsGrid!"
               << std::endl;
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

      /*! Destructor, cleans up the cartesian communicator
       */
      ~FsGrid() {
         MPI_Comm_free(&comm3d);
      }

      /*! Returns the task responsible for handling the cell with the given GlobalID
       * \param id GlobalID of the cell for which task is to be determined
       * \return a task for the grid's cartesian communicator
       */
      int getTaskForGlobalID(GlobalID id);
      /*! Get the localID (=offset into the data-array) for the given globalID
       * \param id the global cellID of a cell.
       */
      LocalID getLocalIdForGlobalID(GlobalID id);

      void setFieldData(GlobalID id, T& value);

      /*! Perform ghost cell communication.
       */
      void updateGhostCells();

      /*! Get the size of the local domain handled by this grid.
       */
      std::array<int, dims> getLocalSize();

      /*! Get a reference to the field data in a cell
       * \param x x-Coordinate, in cells
       * \param y y-Coordinate, in cells
       * \param z z-Coordinate, in cells
       * \return A reference to cell data in the given cell
       */
      T& get(int x, int y, int z);
      const T& get(int x, int y, int z) const;

   private: 
      std::array<uint32_t,dims> globalSize; //!< Global size of the simulation space, in cells
      std::array<uint32_t,dims> localSize;  //!< Local size of simulation space handled by this task

      //! MPI Cartesian communicator used in this grid
      MPI_Comm comm3d;

      //! Number of tasks in each direction
      std::array<int,dims> ntasks;

      //! This task's position in the 3d task grid
      std::array<int, dims> taskPosition;

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
