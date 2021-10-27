#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <exception>


struct PrintUnknownVolume : public std::exception
{
	const char * what () const throw ()
    {
    	return "Volume type not recognized at print";
    }
};


struct AssemblyUnknownVolume : public std::exception
{
	const char * what () const throw ()
    {
    	return "Volume type not recognized at system assembly";
    }
};


struct DistanceNotVolume : public std::exception
{
	const char * what () const throw ()
    {
    	return "Volume type is not suitable for distance measurement";
    }
};


struct EnergyBalanceUnknownVolume : public std::exception
{
	const char * what () const throw ()
    {
    	return "Volume type not recognized at energy balance";
    }
};


struct MeshUnknownVolume : public std::exception
{
	const char * what () const throw ()
    {
    	return "Volume type not recognized at mesh assembly";
    }
};


struct UnconsistemProblemDimensions : public std::exception
{
	const char * what () const throw ()
    {
    	return "The size of a matrix is not coherent with problem dimensions";
    }
};


struct UnconsistemNumberOfVolumes : public std::exception
{
	const char * what () const throw ()
    {
    	return "The size of a matrix is not coherent with the number of volumes";
    }
};


struct UnconsistemNumberOfBoundaries : public std::exception
{
	const char * what () const throw ()
    {
    	return "The size of a matrix is not coherent with the number of boundaries";
    }
};


struct BadQuantityOfAttributes : public std::exception
{
	const char * what () const throw ()
    {
    	return "The number of columns of volms_data or boundary_data is incorrect";
    }
};


#endif