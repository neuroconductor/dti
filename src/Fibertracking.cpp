
#include "Fibertracking.h"

using namespace std;

int n_angle = 0;
int n_visited = 0;
int n_aniso = 0;

Fibertracking::Fibertracking()
{
	this->n_e1 = *(new Vector(( 0, 0, 1)));
	this->n_e2 = *(new Vector(( 0, 1, 0)));
	this->n_e3 = *(new Vector(( 1, 0, 0)));
	this->n_e4 = *(new Vector(( 0,-1, 0)));
	this->n_e5 = *(new Vector(( 0, 0,-1)));
	this->n_e6 = *(new Vector((-1, 0, 0)));
	
	this->max_intersec_angle = 30.;
	this->change_dir = false;
}

Fibertracking::Fibertracking(Voxel& voxels, int x, int y, int z, double dim_x, double dim_y, double dim_z, double min_anisotropy, double max_angle)
{
	this->voxels = &voxels;
	
	this->n_e1 = *(new Vector( 0, 0, 1));
	this->n_e2 = *(new Vector( 0, 1, 0));
	this->n_e3 = *(new Vector( 1, 0, 0));
	this->n_e4 = *(new Vector( 0,-1, 0));
	this->n_e5 = *(new Vector( 0, 0,-1));
	this->n_e6 = *(new Vector(-1, 0, 0));
	
	x_range = x;
	y_range = y;
	z_range = z;
	
	last_start_voxel = 0;
	num_fibers = 0;
	this->dim_x = dim_x, this->dim_y = dim_y, this->dim_z = dim_z;
	
	cur_voxel_index = 0;
	
	this->min_anisotropy = min_anisotropy;
	
	intersec_angle = 0.;
	this->max_intersec_angle = max_angle;
	
	this->change_dir = false;
	
	allVectors = *new VectorList();
}

//Fibertracking::~Fibertracking()
//{
//	currentFiber.~Fiber();
//	curVectorList.~VectorList();
//	allVectors.~VectorList();
//	
//	delete voxels;
//	
//	start_o.~Vector();
//	
//	n_e1.~Vector();
//	n_e2.~Vector();
//	n_e3.~Vector();
//	n_e4.~Vector();
//	n_e5.~Vector();
//	n_e6.~Vector();
//}

int Fibertracking::getLength()
{
	int f = allVectors.getNum_Nan()+1;
	int p = (allVectors.getLength() - allVectors.getNum_Nan()) /2;
	return 12*(p - f);
}

double* Fibertracking::convertToDouble()
{

	int f = allVectors.getNum_Nan()+1;
	int p = (allVectors.getLength() - allVectors.getNum_Nan()) /2;
	int l = 2*(p - f);
	
//	printf("%d\n",f);
	
	double* vals = new double[6*l];
	
//	double* x = new double[l];
//	double* y = new double[l];
//	double* z = new double[l];
//	double* r = new double[l];
//	double* g = new double[l];
//	double* b = new double[l];
	
//	Vector *last_coords = new Vector(allVectors.getStart().getComponents()[0], allVectors.getStart().getComponents()[1], allVectors.getStart().getComponents()[2]);
	
//	printf("vals-array mit length = %d allokiert...\n", 6*l);
	
	int i = 0;
	double temp = 0.;
	
	bool fibreStart = true;
	bool fibreEnd = false;
	double interx = 0;
	double intery = 0;
	double interz = 0;
	double dirx = 0;
	double diry = 0;
	double dirz = 0;
	
	while ( (allVectors.getLength()) > 1)
	{
		temp = allVectors.getStart().getComponents()[1];
//		printf("Schleifendurchlauf: %d mit temp: %f\n", i, temp);
		
		if (isnan(temp))
		{
//			printf("\nentering IF-sequence,\ti = %d\n\n", i);
			
			i-=1;
			allVectors.del_at_start();
			
//			printf("i zurueckgesetzt auf %d\n", i);
			fibreStart = true;
		}
		else
		{
//			printf("\nentering ELSE-sequence,\ti = %d\n\n", i);
//		        printf("Punkt einmal auf %d\n", i);
			
			interx = allVectors.getStart().getComponents()[0];
			intery = allVectors.getStart().getComponents()[1];
			interz = allVectors.getStart().getComponents()[2];

			vals[    i] = interx;
//			printf("vals[%d] = %f\t(inter_x)\n", i, vals[i]);
			vals[  l+i] = intery;
//			printf("vals[%d] = %f\t(inter_y)\n", i, vals[  l+i]);
			vals[2*l+i] = interz;
//			printf("vals[%d] = %f\t(inter_z)\n", i, vals[2*l+i]);
						
			allVectors.del_at_start();
//			printf("start deleted 1\n");
			
//			allVectors.getStart().print();
			
			temp = allVectors.getStart().getComponents()[1];
//			printf("temp %d: %f\n", i, temp);

			dirx = voxels[(int)temp].getDirection().getComponents()[0];
			diry = voxels[(int)temp].getDirection().getComponents()[1];
			dirz = voxels[(int)temp].getDirection().getComponents()[2];
		
			vals[3*l+i] = dirx;
//			printf("vals[%d] = %f\t(dir_x)\n", i, vals[3*l+i]);
			vals[4*l+i] = diry;
//			printf("vals[%d] = %f\t(dir_y)\n", i, vals[4*l+i]);
			vals[5*l+i] = dirz;
//			printf("vals[%d] = %f\t(dir_z)\n", i, vals[5*l+i]);
			
			allVectors.del_at_start();

			if (!fibreStart && allVectors.getLength() > 0) // if not fibrestart re-write point and colour
			{
				i++;						
//			        printf("Punkt nochmal auf %d\n", i);

				vals[    i] = interx;
//				printf("vals[%d] = %f\t(inter_x)\n", i, vals[i]);
				vals[  l+i] = intery;
//				printf("vals[%d] = %f\t(inter_y)\n", i, vals[  l+i]);
				vals[2*l+i] = interz;
//				printf("vals[%d] = %f\t(inter_z)\n", i, vals[2*l+i]);
				vals[3*l+i] = dirx;
//				printf("vals[%d] = %f\t(dir_x)\n", i, vals[3*l+i]);
				vals[4*l+i] = diry;
//				printf("vals[%d] = %f\t(dir_y)\n", i, vals[4*l+i]);
				vals[5*l+i] = dirz;
//				printf("vals[%d] = %f\t(dir_z)\n", i, vals[5*l+i]);
			}

			i++;
			fibreStart = false;
		}
		
//		printf("i fuer naechsten Punkt auf %d\n", i);
		
	}
	
	return vals;
}

void Fibertracking::nextVoxel_forward()
{
	int cur_x = voxels[cur_voxel_index].getX();
	int cur_y = voxels[cur_voxel_index].getY();
	int cur_z = voxels[cur_voxel_index].getZ();
	
	int x = cur_x, y = cur_y, z = cur_z;
	
	int plane_dir = 0;
	
	if ( cur_x < 0 || cur_y < 0 || cur_z < 0 || cur_x > (x_range-1) || cur_y > (y_range-1) || cur_z > (z_range-1))
	{
		return;
	}
	
	Vector voxel_d = voxels[cur_voxel_index].getDirection();
	
	Vector intersection(3);
	
//	printf("\tcur_x = %d, cur_y = %d, cur_z = %d\n", cur_x, cur_y, cur_z);
	
	/**
	 *  position vectors of the plane equation
	 **/
	// vertex of the voxel which points to the point of origin
	Vector voxel_bottom(  cur_x   *dim_x,  cur_y   *dim_y,  cur_z   *dim_z );
	// opposing vertex
	Vector voxel_top   ( (cur_x+1)*dim_x, (cur_y+1)*dim_y, (cur_z+1)*dim_z );
	
	double *distances = new double[7];
	
//	printf("\tchange = %d\n", change);
	
	if (change_dir)
	{
		voxel_d = voxel_d * -1.;
	}
	
	distances[0] = HUGE_VAL;
	distances[1] = ( ( voxel_bottom - start_o ) * n_e1 ) / ( voxel_d * n_e1 );
	distances[2] = ( ( voxel_bottom - start_o ) * n_e2 ) / ( voxel_d * n_e2 );
	distances[3] = ( ( voxel_bottom - start_o ) * n_e3 ) / ( voxel_d * n_e3 );
	distances[4] = ( (    voxel_top - start_o ) * n_e4 ) / ( voxel_d * n_e4 );
	distances[5] = ( (    voxel_top - start_o ) * n_e5 ) / ( voxel_d * n_e5 );
	distances[6] = ( (    voxel_top - start_o ) * n_e6 ) / ( voxel_d * n_e6 );
	
	double dSkalar = 0.;
	
	for (int i = 1; i < 7; i++)
	{
		if (fabs(distances[i]) < fabs(distances[plane_dir]) && (distances[i] > 0.) )
		{
			dSkalar = distances[i];
			plane_dir = i;
		}
	}
	
	switch (plane_dir)
	{
		case 1:	{ z--; break; }
		case 2:	{ y--; break; }
		case 3:	{ x--; break; }
		case 4:	{ y++; break; }
		case 5:	{ z++; break; }
		case 6:	{ x++; break; }
	}
	
	if ( x < 0 || y < 0 || z < 0 || x > (x_range-1) || y > (y_range-1) || z > (z_range-1) )
	{
		return;
	}
	else
	{	
		cur_x = x; cur_y = y; cur_z = z;
		dSkalar = distances[plane_dir];
	}
	
	cur_voxel_index = cur_x + cur_y * x_range + cur_z * x_range * y_range;
	
	float zaehler = (float) (voxel_d * (voxels[x + y * x_range + z * x_range * y_range].getDirection()));
	
	intersec_angle = 180./M_PI * acos(zaehler);
	
	intersection = start_o + ( voxel_d *  dSkalar );
	
	start_o = intersection;
	
	if (zaehler < .0f)
	{
		change_dir = true;
		intersec_angle = 180-intersec_angle;
	}
	else
	{
		change_dir = false;
	}
}

void Fibertracking::nextVoxel_backward()
{
	int cur_x = voxels[cur_voxel_index].getX();
	int cur_y = voxels[cur_voxel_index].getY();
	int cur_z = voxels[cur_voxel_index].getZ();
	
	int x = cur_x, y = cur_y, z = cur_z;
	
	int plane_dir = 0;
	
	if ( cur_x < 0 || cur_y < 0 || cur_z < 0 || cur_x > (x_range-1) || cur_y > (y_range-1) || cur_z > (z_range-1))
	{
		return;
	}
	
	Vector voxel_d = voxels[cur_voxel_index].getDirection();
	
	Vector intersection(3);
	
	/**
	 *  position vectors of the plane equation
	 **/
	// vertex of the voxel which points to the point of origin
	Vector voxel_bottom(  cur_x   *dim_x,  cur_y   *dim_y,  cur_z   *dim_z );
	// opposing vertex
	Vector voxel_top   ( (cur_x+1)*dim_x, (cur_y+1)*dim_y, (cur_z+1)*dim_z );
	
	double *distances = new double[7];
		
	if (change_dir)
	{
		voxel_d = voxel_d * -1.;
	}
		
	distances[0] = HUGE_VAL;
	distances[1] = ( ( voxel_bottom - start_o ) * n_e1 ) / ( voxel_d * n_e1 );
	distances[2] = ( ( voxel_bottom - start_o ) * n_e2 ) / ( voxel_d * n_e2 );
	distances[3] = ( ( voxel_bottom - start_o ) * n_e3 ) / ( voxel_d * n_e3 );
	distances[4] = ( (    voxel_top - start_o ) * n_e4 ) / ( voxel_d * n_e4 );
	distances[5] = ( (    voxel_top - start_o ) * n_e5 ) / ( voxel_d * n_e5 );
	distances[6] = ( (    voxel_top - start_o ) * n_e6 ) / ( voxel_d * n_e6 );
	
	double dSkalar = 0.;
	
	for (int i = 1; i < 7; i++)
	{
		if (fabs(distances[i]) < fabs(distances[plane_dir]) && (distances[i] < 0.) )
		{
			dSkalar = distances[i];
			plane_dir = i;
		}
	}
	
	switch (plane_dir)
	{
		case 1:	{ z--; break; }
		case 2:	{ y--; break; }
		case 3:	{ x--; break; }
		case 4:	{ y++; break; }
		case 5:	{ z++; break; }
		case 6:	{ x++; break; }
	}
	
	if ( x < 0 || y < 0 || z < 0 || x > (x_range-1) || y > (y_range-1) || z > (z_range-1) )
	{
		return;
	}
	else
	{	
		cur_x = x; cur_y = y; cur_z = z;
		dSkalar = distances[plane_dir];
	}
	
	cur_voxel_index = cur_x + cur_y * x_range + cur_z * x_range * y_range;
	
	float zaehler = (float) (voxel_d * (voxels[x + y * x_range + z * x_range * y_range].getDirection()));
	
	intersec_angle = 180./M_PI * acos(zaehler);
	
	intersection = start_o + ( voxel_d *  dSkalar );
	
	if (intersection.getComponents()[0] < 0. ||
		intersection.getComponents()[1] < 0. ||
		intersection.getComponents()[2] < 0.)
	{
		printf("NEGATIVER SCHNITTPUNKT!!!!!\n");
		intersection.print();
	}
	
	start_o = intersection;
	
	if (zaehler < .0f)
	{
		change_dir = true;
		intersec_angle = 180-intersec_angle;
	}
	else
	{
		change_dir = false;
	}
}

void Fibertracking::trackFiber_forward()
{
//	printf("Searching forward...\n");
	
	Voxel *current = &voxels[cur_voxel_index];
	Vector *curVec;
//	current->setVisited(false);
	
	start_o = *new Vector( (current->getX()+0.5)*dim_x, (current->getY()+0.5)*dim_y, (current->getZ()+0.5)*dim_z );

	// Startpunkt einschreiben
//	curVectorList = *new VectorList(start_o);
//	curVec = new Vector(0., (double)cur_voxel_index, 0.);
//	curVectorList.add_at_end(*curVec);
	
	curVectorList = *new VectorList();
	
	while (current->getAnisotropy() >= min_anisotropy && !current->isVisited() && (fabs(intersec_angle) <= max_intersec_angle) )
	{
//		current->print();
		current->setVisited(true);
		current->setStartable(false);
		
		currentFiber.add_at_end(*current);
		
		curVectorList.add_at_end(start_o);
		
		curVec = new Vector(0., (double)cur_voxel_index, 0.);
		curVectorList.add_at_end(*curVec);

		nextVoxel_forward();

		if ( current == &voxels[cur_voxel_index]) // || voxels[cur_voxel_index].isVisited()
		{
			break;
		}
		else
		{
//			currentFiber.add_at_end(*current);
			current = &voxels[cur_voxel_index];
			
//			curVectorList.add_at_end(start_o);
			
//			curVec = new Vector(0., (double)cur_voxel_index, 0.);
//			curVectorList.add_at_end(*curVec);
		}
	}
	
//	printf("Fiber ended because of: ");
	
	if (current->isVisited())
	{
//		printf("next voxel is already visited");
		n_visited++;
	}
	
	if (current->getAnisotropy() < min_anisotropy)
	{
//		printf(" next voxel has to low anisotrpy ");
		n_aniso++;
	}
	
	if (fabs(intersec_angle) > max_intersec_angle)
	{
//		printf(" intersection angle is to large ");
		n_angle++;
	}
	
//	printf(".\n");
}

void Fibertracking::trackFiber_backward()
{
//	printf("Searching backward...\n");
	
	Voxel *current = &voxels[cur_voxel_index];
	Vector *curVec;
	current->setVisited(false);
	
	start_o = *new Vector( (current->getX()+0.5)*dim_x, (current->getY()+0.5)*dim_y, (current->getZ()+0.5)*dim_z );
	
	// Startpunkt einschreiben
	nextVoxel_backward();
	
	if ( current == &voxels[cur_voxel_index]) // || voxels[cur_voxel_index].isVisited()
	{
		return;
	}
	else
	{
		current = &voxels[cur_voxel_index];
//
//			curVec = new Vector(0., (double)cur_voxel_index, 0.);
//			curVectorList.add_at_start(*curVec);
//
//			curVectorList.add_at_start(start_o);
//			
//			currentFiber.add_at_start(*current);
	}
	
	while (current->getAnisotropy() >= min_anisotropy && !current->isVisited() && (fabs(intersec_angle) <= max_intersec_angle) )
	{
//		current->print();
		current->setVisited(true);
		current->setStartable(false);
		
		curVec = new Vector(0., (double)cur_voxel_index, 0.);
		curVectorList.add_at_start(*curVec);

		curVectorList.add_at_start(start_o);
		
		currentFiber.add_at_start(*current);
		
		nextVoxel_backward();
		
		if ( current == &voxels[cur_voxel_index]) // || voxels[cur_voxel_index].isVisited()
		{
			return;
		}
		else
		{
			current = &voxels[cur_voxel_index];
//
//			curVec = new Vector(0., (double)cur_voxel_index, 0.);
//			curVectorList.add_at_start(*curVec);
//
//			curVectorList.add_at_start(start_o);
//			
//			currentFiber.add_at_start(*current);
		}
	}
	
//	printf("Fiber ended because of: ");
	
	if (current->isVisited())
	{
//		printf("next voxel is already visited");
		n_visited++;
	}
	
	if (current->getAnisotropy() < min_anisotropy)
	{
//		printf(" next voxel has to low anisotrpy ");
		n_aniso++;
	}
	
	if (fabs(intersec_angle) > max_intersec_angle)
	{
//		printf(" intersection angle is to large ");
		n_angle++;
	}
	
//	printf(".\n");
}

void Fibertracking::findAllFibers()
{
//	printf("Searching for fibers...\n");
	
	while (last_start_voxel < x_range*y_range*z_range)
	{
		if (voxels[last_start_voxel].getAnisotropy() >= min_anisotropy && voxels[last_start_voxel].isStartable())
		{
			num_fibers++;
			
			currentFiber = *new Fiber();
			curVectorList = *new VectorList(); 
			
//			printf("Fiber found!\n");
//			printf("============\n");
			
			cur_voxel_index = voxels[last_start_voxel].getX() + voxels[last_start_voxel].getY()*x_range + voxels[last_start_voxel].getZ()*x_range*y_range;
			trackFiber_forward();
			
			// Zuruecksetzen von wichtigen Parametern
			intersec_angle = 0.;
			cur_voxel_index = voxels[last_start_voxel].getX() + voxels[last_start_voxel].getY()*x_range + voxels[last_start_voxel].getZ()*x_range*y_range;
			trackFiber_backward();
			
			allVectors.add_list(curVectorList);
			
			currentFiber.unvisit();
			
//			printf("============\n");
//			printf("Searching continued...\n");
		}
		
		last_start_voxel++;
	}
	
	allVectors.del_at_start();
	
//	printf("End of searching.\n");
}

void Fibertracking::findMarkedFibers(int* ranges)
{
	
	// Region of Interest muss festgelegt werden durch ein Voxelarray
	int length = (ranges[1]-(ranges[0]-1))*(ranges[3]-(ranges[2]-1))*(ranges[5]-(ranges[4]-1));
	
	Voxel *marked = new Voxel[length];
	
	int counter = 0;
	
	for (int cur_z = ranges[4]-1; cur_z < ranges[5]; cur_z++)
	{
		for (int cur_y = ranges[2]-1; cur_y < ranges[3]; cur_y++)
		{
			for (int cur_x = ranges[0]-1; cur_x < ranges[1]; cur_x++)
			{
				cur_voxel_index = cur_x + cur_y * x_range + cur_z * x_range * y_range;
				
				marked[counter] = voxels[cur_voxel_index];
				counter++;
			}
		}
	}
	
	counter = 0;
	
	while (last_start_voxel < length)
	{
		if (marked[last_start_voxel].getAnisotropy() >= min_anisotropy && marked[last_start_voxel].isStartable())
		{
//			num_fibers++;
			
			currentFiber = *new Fiber();
			curVectorList = *new VectorList(); 
			
//			printf("Fiber found!\n");
//			printf("============\n");
			
			cur_voxel_index = marked[last_start_voxel].getX() + marked[last_start_voxel].getY()*x_range + marked[last_start_voxel].getZ()*x_range*y_range;
			trackFiber_forward();
			
			// Zuruecksetzen von wichtigen Parametern
			intersec_angle = 0.;
			cur_voxel_index = marked[last_start_voxel].getX() + marked[last_start_voxel].getY()*x_range + marked[last_start_voxel].getZ()*x_range*y_range;
			trackFiber_backward();
			
			allVectors.add_list(curVectorList);
			
			currentFiber.unvisit();
			
//			printf("============\n");
//			printf("Searching continued...\n");
		}
		
		last_start_voxel++;
	}
	
	allVectors.del_at_start();
	
//	printf("End of searching.\n");
	
	double all_abort = n_visited+n_angle+n_aniso;

	printf("Abort fibers because of:\nvisited\t=\t%d (%f%)\naniso\t=\t%d (%f%)\nangle\t=\t%d (%f%)\n", n_visited, n_aniso, n_angle, (double)n_visited*100./all_abort, (double)n_aniso*100./all_abort, (double)n_angle*100./all_abort);
}
