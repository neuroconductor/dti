
#include "Converter.h"

Voxel& Converter::getVoxels()
{
	return *voxels;
}

Converter::Converter(double* data_dir_coords, double* data_FA_values, int x_range, int y_range, int z_range)
{
//	printf("%d, %d, %d\n", x_range, y_range, z_range);
	
	voxels = new Voxel[x_range*y_range*z_range];
	
//	printf("converting doubles to Voxels...\n");

	Vector *dir;
    
	double dir_x, dir_y, dir_z, FA;
	int i = 0, j = 0, k = 0, num_voxel = 0, counter = 0;
	
    for (k = 0; k < z_range; k++)
    {
	    for (j = 0; j < y_range; j++)
	    {
	    	for (i = 0; i < x_range; i++)
		   	{
		   		dir_x = data_dir_coords[num_voxel];num_voxel++;
		   		dir_y = data_dir_coords[num_voxel];num_voxel++;
		   		dir_z = data_dir_coords[num_voxel];
		   		FA 	  =  data_FA_values[i + j * x_range + k * x_range * y_range];
		   		
				dir = new Vector(dir_x, dir_y, dir_z);
				
				Voxel v(i, j, k, *dir, FA);
				
				voxels[i + j * x_range + k * x_range * y_range] = v;
				
				num_voxel++;
				counter++;
		   	}
	    }
    }

//    printf("Converted %d Voxels.\n", counter);
}
