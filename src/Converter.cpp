
#include "Converter.h"

/**
 * @todo re-write this to return std::shared_ptr<Voxel> or std::vector<Voxel>.
 * Returning pointers or references may result in memory leaks.
 */
Voxel& Converter::getVoxels()
{
	return *voxels;
}

Converter::Converter(double* data_dir_coords, double* data_FA_values, int* mask, int dim_x, int dim_y, int dim_z)
{
	voxels = new Voxel[dim_x*dim_y*dim_z];

	Vector* dir;

	double dir_x, dir_y, dir_z, FA;

	/*
	 * Increasing an integer was more fail proof then doing the
	 * math with the indices dim_x, dim_y, dim_z.
	 */
	int num_vector = 0;
	bool startable = false;
	
	/*
	 * iterate over all voxels and extract direction vector, FA and masking information
	 */
    for (int k = 0; k < dim_z; k++)
    {
	    for (int j = 0; j < dim_y; j++)
	    {
	    	for (int i = 0; i < dim_x; i++)
		   	{
		   		dir_x = data_dir_coords[num_vector];num_vector++;
		   		dir_y = data_dir_coords[num_vector];num_vector++;
		   		dir_z = data_dir_coords[num_vector];num_vector++;

		   		FA 	  =  data_FA_values[i + j * dim_x + k * dim_x * dim_y];
				startable = (mask[i + j * dim_x + k * dim_x * dim_y] == 0) ? false : true;
		   		
				/*
				 * try to use more modern C++11 features here like
				 * move-semantics and smart pointers
				 */
				dir = new Vector(dir_x, dir_y, dir_z);
				
				Voxel v(i, j, k, 1, *dir, FA);
				
				voxels[i + j * dim_x + k * dim_x * dim_y] = v;
				voxels[i + j * dim_x + k * dim_x * dim_y].setStartable(startable);
		   	}
	    }
    }
}


Converter::Converter(double* data_dir_coords, double* data_FA_values, int* mask, int* data_order, int maxorder, int dim_x, int dim_y, int dim_z)
{
	voxels = new Voxel[dim_x*dim_y*dim_z];

	double dir_x, dir_y, dir_z, FA;
	int num_vector = 0, order = 0;
	bool startable = false;
	

	/*
	 * iterate over all voxels and extract direction vector, FA and masking information
	 */
    for (int k = 0; k < dim_z; k++)
    {
	    for (int j = 0; j < dim_y; j++)
	    {
	    	for (int i = 0; i < dim_x; i++)
		   	{
		   		FA 	  = data_FA_values[i + j * dim_x + k * dim_x * dim_y];
		   		order = data_order[i + j * dim_x + k * dim_x * dim_y];
				startable = (mask[i + j * dim_x + k * dim_x * dim_y] == 0) ? false : true;
		   		
		   		Vector *dir = new Vector[order];

		   		/*
		   		 * Each voxel has an order which tells, how many direction vectors are stored there.
		   		 * So we have to iterate here too.
		   		 */
		   		for (int l = 0; l < order; l++)
		   		{
			   		dir_x = data_dir_coords[num_vector];num_vector++;
			   		dir_y = data_dir_coords[num_vector];num_vector++;
			   		dir_z = data_dir_coords[num_vector];num_vector++;
			   		
			   		dir[l] = Vector(dir_x, dir_y, dir_z);
		   		}
				
		   		num_vector += (maxorder - order)*3;

				/*
				 * try to use more modern C++11 features here like
				 * move-semantics and smart pointers
				 */
		   		Voxel v(i, j, k, order, *dir, FA);

				voxels[i + j * dim_x + k * dim_x * dim_y] = v;
				voxels[i + j * dim_x + k * dim_x * dim_y].setStartable(startable);
		   	}
	    }
    }
}
