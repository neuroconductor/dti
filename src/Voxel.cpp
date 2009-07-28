#include "Voxel.h"

using namespace std;

Voxel::Voxel()
{
	this->position = Vector(-1, -1, -1);
	this->direction = Vector(0, 0, 0);
	this->x = -1;
	this->y = -1;
	this->z = -1;
	
	this->anisotropy = 0.;
	
	this->visited = false;
	this->startable = true;
	
	this->red = 0;
	this->green = 0;
	this->blue = 0;
	
	this->next = NULL;
	this->prev = NULL;
}

Voxel::Voxel(int x, int y, int z, Vector& direction, double anisotropy)
{
	this->direction = direction;
	this->anisotropy = anisotropy;
	this->x = x;
	this->y = y;
	this->z = z;

	this->position = Vector((double)this->x,(double)this->y,(double)this->z);
	
	visited = false;
	startable = true;
	
	red = 	(int)(255 * fabs((this->direction.getComponents()[0]) * fabs(this->anisotropy)));
	green = (int)(255 * fabs((this->direction.getComponents()[1]) * fabs(this->anisotropy)));
	blue =  (int)(255 * fabs((this->direction.getComponents()[2]) * fabs(this->anisotropy)));
	
	next = NULL;
	prev = NULL;
}

//Voxel::~Voxel()
//{
//	delete next;
//	delete prev;
//	
//	 position.~Vector();
//	direction.~Vector();
//}

/**  getters  **/

bool Voxel::isVisited()			{return visited;}
bool Voxel::isStartable()		{return startable;}
double Voxel::getAnisotropy()	{return anisotropy;}
Vector& Voxel::getDirection()	{return direction;}
Vector& Voxel::getPosition()	{return position;}
int Voxel::getX()				{return x;}
int Voxel::getY()				{return y;}
int Voxel::getZ()				{return z;}
int Voxel::getRed()				{return red;}
int Voxel::getGreen()			{return green;}
int Voxel::getBlue()			{return blue;}
Voxel* Voxel::getNext()			{return next;}
Voxel* Voxel::getPrev()			{return prev;}

/**  setters  **/
void Voxel::setVisited(bool isVisited)	{this->visited = isVisited;}
void Voxel::setStartable(bool startable){this->startable = startable;}
void Voxel::setNext(Voxel *next)		{this->next = next;}
void Voxel::setPrev(Voxel *prev)		{this->prev = prev;}

/**  print-methode  **/
void Voxel::print()
{
	printf("%d %d %d\n", x, y, z);
}
