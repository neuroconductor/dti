#ifndef CONVERTER_H_
#define CONVERTER_H_

#include <R.h>
#include <Rinternals.h>
#include "Voxel.h"

class Converter
{
	private:
		Voxel* voxels;
		
	public:
		Converter(double*, double*, int, int, int);
		Voxel& getVoxels();
};

#endif /*CONVERTER_H_*/
