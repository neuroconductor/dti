#ifndef FIBER_H_
#define FIBER_H_

#include "Voxel.h"

/*
 * The class 'Fiber' is a LinkedList of Voxel objects.
 */

class Fiber
{
	private:
		// output file
		FILE *voxels;
		
		Voxel* start;
		Voxel* end;
	
		int length;
		
		// minimum length of the fibers to be logged in the output file
		int minLength;
		
		string path_v;
		string input_v;
		
	public:
		/**  constructors & destructors  **/
		Fiber();
		Fiber(string&);
		Fiber(Voxel*);
		
//		~Fiber();
		
		/** getters  **/
		int getLength();
		
		/**  LinkedList-methods  **/
		void add_at_end(Voxel&);
		void add_at_start(Voxel&);
		void del_at_start();
		
		// sets all voxels within the fiber to visited = false
		void unvisit();
		
		//  print-method
		void print();
		
		/**  methods for file output  **/
		void addVector_forw(Vector&, Voxel&);
		void addVector_back(Vector&, Voxel&);
		void write();
};

#endif /*FIBER_H_*/
