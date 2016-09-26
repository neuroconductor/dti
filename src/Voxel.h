#ifndef VOXEL_H_
#define VOXEL_H_

#include "Vector.h"

/**
 * @class Voxel
 * @brief This class contains the direction vectors and the anisotropy
 * from the given R variables.
 */
class Voxel
{
	private:
	
		/// spatial x location
		int x;

		/// spatial y location
		int y;

		/// spatial z location
		int z;

		/// order, i.e. how many direction vectors
		int order;

		///
		int dir_index;
		
		/// position vector [x, y, z]
		// Why do we have this twice?
		Vector position;

		/// direction vectors
		Vector *directions;

		/// fractional anisotropy
		double anisotropy;
	
		/// Can a fiber start in this voxel?
		bool startable;
		

		// linked list variables, this can be deleted, when Fiber is rewritten.
		/// next voxel in linked list
		Voxel* next;

		/// previous voxel in linked list
		Voxel* prev;
		
		/// Is this voxel already part of the current fiber? (avoids circular fibers)
		bool visited;
		
		// this used to be in discussion, but not yet implemented
//		Voxel *subvoxels;
	
	public:
		/**
		 * Default constructor. Voxel::order is set to 1.
		 */
		Voxel();

		/**
		 * @brief Initialize the voxel with the given values.
		 *
		 * @param x_in spatial x location
		 * @param y_in spatial y location
		 * @param z_in spatial z location
		 * @param order_in order, i.e. how many direction vectors
		 * @param directions_in direction vectors
		 * @param anisotropy_in fractional anisotropy
		 */
		Voxel(int x_in, int y_in, int z_in, int order_in, Vector& directions_in, double anisotropy_in);
		
//		~Voxel();
		
		/**
		 *
		 * @return
		 */
		bool isVisited();

		/**
		 *
		 * @return
		 */
		bool isStartable();

		/**
		 *
		 * @return
		 */
		double getAnisotropy();

		/**
		 *
		 *
		 * @return direction
		 */
		Vector* getDirections();

		/**
		 *
		 * @return
		 */
		Vector& getPosition();
		int getX();
		int getY();
		int getZ();
		int getOrder();
		int getDir_Index();
		Voxel* getNext();
		Voxel* getPrev();
//		Voxel* getSubvoxels();
		
		/**  setters  **/
		void setDir_Index(int);
		void setVisited(bool);
		void setStartable(bool);
		void setNext(Voxel*);
		void setPrev(Voxel*);
		
		/**
		 * Not yet implemented.
		 */
		void print();
};

#endif /*VOXEL_H_*/
