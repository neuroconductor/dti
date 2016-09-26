#ifndef CONVERTER_H_
#define CONVERTER_H_

#include "Voxel.h"
#include <R.h>
#include <Rinternals.h>

/**
 * @class Converter
 * @brief This class converts the data from R into the internal
 * voxel representation.
 *
 */
class Converter
{
	private:

		/// created voxels which are used in the fiber tracking algorithm
		Voxel* voxels;
		
	public:

		/**
		 * @brief convert R data to voxels (single direction per voxel)
		 *
		 * R data is provided as double* and int* .
		 * Each voxel contains one direction vector and one fractional anisotropy value.
		 * The order of all voxels is 1.
		 *
		 * @param data_dir_coords direction vector of highest anisotropy per voxel
		 * @param data_FA_values fractional anisotropy values per voxel
		 * @param mask marks where the fiber tracking may start (!=0) and where not (0)
		 * @param dim_x spatial x dimension of the data cube
		 * @param dim_y spatial y dimension of the data cube
		 * @param dim_z spatial z dimension of the data cube
		 */
		Converter(double* data_dir_coords, double* data_FA_values, int* mask, int dim_x, int dim_y, int dim_z);

		/**
		 * @brief convert R data to voxels (multiple directions per voxel)
		 *
		 * R data is provided as double* and int* .
		 * Each voxel contains multiple direction vectors and one fractional anisotropy value.
		 * The order of each voxel is extracted from data_order.
		 *
		 * @param data_dir_coords direction vectors of highest anisotropy per voxel
		 * @param data_FA_values fractional anisotropy values per voxel
		 * @param mask marks where the fiber tracking may start (1) and where not (0)
		 * @param data_order order of each voxel
		 * @param maxorder highest observed order in a voxel in this data set
		 * @param dim_x spatial x dimension of the data cube
		 * @param dim_y spatial y dimension of the data cube
		 * @param dim_z spatial z dimension of the data cube
		 */
		Converter(double* data_dir_coords, double* data_FA_values, int* mask, int* data_order, int maxorder, int dim_x, int dim_y, int dim_z);

		/**
		 * @brief return the created voxels
		 *
		 * @return voxels from R data
		 */
		Voxel& getVoxels();
};

#endif /*CONVERTER_H_*/
