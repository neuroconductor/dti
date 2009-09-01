#include "Fiber.h"

using namespace std;

Fiber::Fiber()
{
	start = NULL;
	end   = NULL;
	length = 0;
	minLength = 2;
}

Fiber::Fiber(string& path_v)
{
	start = NULL;
	end   = NULL;
	length = 0;
	minLength = 2;
	
	this->path_v = path_v;
	input_v = "";
	
	voxels = fopen( path_v.c_str(), "w");
}

Fiber::Fiber(Voxel* start)
{
	path_v = "";
	input_v = "";
	
	this->start = start;
	this->end   = start;
	length = 1;
	minLength = 2;
	
	voxels = fopen( path_v.c_str(), "w");
}

//Fiber::~Fiber()
//{
//	delete start;
//	delete end;
//	
//	fclose(voxels);
//	
//	delete voxels;
//}

int Fiber::getLength()		{return length;}

void Fiber::add_at_end(Voxel &add) // OLD!!!
{
	if (end == NULL)
	{
		add.setPrev(NULL);
		add.setNext(NULL);
		start = &add;
		end   = &add;
	}
	else
	{
		add.setPrev(end);
		end->setNext(&add);
		end = &add;
		end->setNext(NULL);
	}
	
	length++;
}

//void Fiber::add_at_end(Voxel &add) // NEW!!!
//{
//	if (&add == NULL)
//	{
//		return;
//	}
//	else
//	{
//		if (end == NULL)
//		{
//			start = &add;
//			end   = &add;
//			start->setNext(NULL);
//			start->setPrev(NULL);
//		}
//		else
//		{
//			add.setNext(NULL);
//			add.setPrev(end);
//			end->setNext(&add);
//			end = &add;
//		}
//		
//		length++;
//	}
//}

void Fiber::add_at_start(Voxel &add) // OLD!!!
{
	if (start == NULL)
	{
		add.setPrev(NULL);
		add.setNext(NULL);
		start = &add;
		end   = &add;
	}
	else
	{
		start->setPrev(&add);
		add.setNext(start);
		start = &add;
		start->setPrev(NULL);
	}
	
	length++;
}

//void Fiber::add_at_start(Voxel &add) //NEW!!!
//{
//	if (&add == NULL)
//	{
//		return;
//	}
//	else
//	{
//		if (start == NULL)
//		{
//			start = &add;
//			end   = &add;
//			start->setNext(NULL);
//			start->setPrev(NULL);
//		}
//		else
//		{
//			Voxel* uebergang = start;
//			start = &add;
//			start->setNext(uebergang);
//		}
//		
//		length++;
//	}
//}

void Fiber::del_at_start()
{
	start->getNext()->setPrev(NULL);
	start->setNext(NULL);
	
	delete start;
	
	length--;
}

void Fiber::unvisit()
{
	if (length == 0)
	{
		return;
	}
	
	Voxel* aktuell = start;
	
	int i;
	
	for (i = 0; i < length; i++)
	{
		aktuell->setVisited(false);
		aktuell = aktuell->getNext();
		
//		if (aktuell == NULL)
//		{
//			break;
//		}	
	}
}

void Fiber::print()
{
	if (length == 0)
	{
		return;
	}
	
	Voxel* aktuell = start;
	
	int i;
	
	for (i = 0; i < length; i++)
	{
		aktuell->print();
		aktuell = aktuell->getNext();
		
//		if (aktuell == NULL)
//		{
//			break;
//		}	
	}
	
//	end->print();
}

void Fiber::addVector_forw(Vector &v, Voxel &vo)
{
	stringstream s;
	
	s << v.getComponents()[0] << "\t" << v.getComponents()[1] << "\t" << v.getComponents()[2] << "\t" << vo.getRed() << "\t" << vo.getGreen() << "\t" << vo.getBlue() << "\n";
	
	input_v = input_v + s.str();
}

void Fiber::addVector_back(Vector &v, Voxel &vo)
{
	stringstream s;
	
	s << v.getComponents()[0] << "\t" << v.getComponents()[1] << "\t" << v.getComponents()[2] << "\t" << vo.getRed() << "\t" << vo.getGreen() << "\t" << vo.getBlue() << "\n";
	
	input_v = s.str() + input_v;
}

void Fiber::write()
{
	// If the fiber has not the minimum length, the file will be deleted here.
	
	if (length <= minLength)
	{
		input_v = "";
		fclose(voxels);
		
		printf("Fiber to short and will be deleted.\n");
		remove(path_v.c_str());
	}
	else
	{
		fprintf(voxels, input_v.c_str());
		fclose(voxels);
	}
}
