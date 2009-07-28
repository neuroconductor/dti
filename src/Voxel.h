#ifndef VOXEL_H_
#define VOXEL_H_

#include "Vector.h"

using namespace std;

/*
 * The class 'Voxel' contains the direction vectors and the anisotropy
 * from the given R variables.
 */

class Voxel
{
	private:
	
		/* x: 	x-location
 		 * y: 	y-location
 		 * z: 	z-location
 		 * direction:	direction-vector
		 * anisotropy:	anisotropy 		*/
	
		int x, y, z;
		
		Vector position;
		Vector direction;
		double anisotropy;
	
		bool startable;
		
		// Color values to generate a EPS-File via GNUPLOT
		int red, green, blue;
	
		// LinkedList-Variables
		Voxel* next;
		Voxel* prev;
		
		bool visited;
	
	public:
		/**  constructors & destructor  **/
		Voxel();
		Voxel(int, int, int, Vector&, double);
		
//		~Voxel();
		
		/**  getters  **/
		bool isVisited();
		bool isStartable();
		double getAnisotropy();
		Vector& getDirection();
		Vector& getPosition();
		int getX();
		int getY();
		int getZ();
		int getRed();
		int getGreen();
		int getBlue();
		Voxel* getNext();
		Voxel* getPrev();
		
		/**  setters  **/
		void setVisited(bool);
		void setStartable(bool);
		void setNext(Voxel*);
		void setPrev(Voxel*);
		
		/**  print-methode  **/
		void print();
};

#endif /*VOXEL_H_*/
