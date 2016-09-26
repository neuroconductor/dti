#ifndef FIBER_H_
#define FIBER_H_

#include "Voxel.h"

/**
 * @class Fiber
 * @brief This class implements a double linked list of voxel objects.
 * @todo Replace this with std::deque<Voxel> or std::vector<Voxel> and add unvisit/print facilities elsewhere.
 *
 * I will not comment on this class in detail, as you can find out about double linked lists
 * at Wikipedia and this class should be replaced ASAP.
 */
class Fiber
{
	private:

		/// front of the linked list
		Voxel* start;

		/// back of the linked list
		Voxel* end;
	
		/// size of the linked list
		int length;
		
	public:
		/**
		 * @brief Default constructor, set every pointer to NULL.
		 */
		Fiber();
		
		/**
		 * @brief Initialize the fiber with one voxel.
		 *
		 * @param set_start
		 */
		Fiber(Voxel* set_start);

//		~Fiber();
		
		/**
		 * @brief Return the current length of the fiber.
		 *
		 * @return how many voxels are in the fiber
		 */
		int getLength();
		
		/**
		 * @brief Add a voxel at the end of the fiber.
		 * @see std::vector<T>::push_back()
		 *
		 * @param add voxel to add
		 */
		void add_at_end(Voxel& add);

		/**
		 * @brief Add a voxel at the start of the fiber.
		 * @see std::vector<T>::push_front()
		 *
		 * @param add voxel to add
		 */
		void add_at_start(Voxel& add);

		/**
		 * @brief Delete a voxel at the start of the fiber.
		 * @see std::vector<T>::pop_front()
		 */
		void del_at_start();
		
		/**
		 * @brief Sets all voxels within the fiber to Voxel::visited = false .
		 */
		void unvisit();
		
		/**
		 * @brief Print the current fiber from start to end.
		 *
		 * This uses the Voxel::print() function.
		 *
		 * @see Voxel::print()
		 */
		void print();
};

#endif /*FIBER_H_*/
