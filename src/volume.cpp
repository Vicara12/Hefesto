#include "volume.h"
#include <math.h>
#include <iostream>


VType Volume::volumeType () const
{
    return vol_type_;
}


Volume::Volume (VType vol_type) :
    vol_type_(vol_type)
{
    //
}


////////////////////////////////////////////////////////////////


SolidVolume::SolidVolume (double volume, double lambda, double qv,
                          const double *surfaces, int index,
                          const double *position) :
        Volume(VType::solid), lambda_(lambda), qv_(qv), volume_(volume), index_(index)
{
    for (int i = 0; i < PROBLEM_DIM*2; i++)
        surface_[i] = surfaces[i];

    for (int i = 0; i < PROBLEM_DIM; i++)
        position_[i] = position[i];
}


void SolidVolume::setBoundaries (const Volume **boundaries)
{
    for (int i = 0; i < PROBLEM_DIM*2; i++)
        boundaries_[i] = boundaries[i];
}


void SolidVolume::getEquation (DoubleVector &coefs)
{
    int n_nodes = coefs.size()-1;

    // for each boundary (two per dimension are assumed)
    for (int i = 0; i < PROBLEM_DIM*2; i++)
    {
        VType boundary_type = boundaries_[i]->volumeType();

        // the coefficients that need to be actualized
        // depend on the type of volume or boundary
        if (boundary_type == VType::solid)
        {
            // cast volume pointer from generic volume to solid volume
            const SolidVolume *boundary = (SolidVolume*)&(boundaries_[i]);

            int boundary_i = boundary->index_;
            // this is not exactly the average on the boundary of both volumes
            // but in the midlle point between their centers, but it's close enough
            double lambda = 2/(1/this->lambda_ + 1/boundary->lambda_);
            double S = surface_[i];
            double d = distanceToVolume(boundary);

            coefs[index_] -= lambda_*S/d;
            coefs[boundary_i] += lambda_*S/d;
        }
        else if (boundary_type == VType::convection_boundary)
        {
            // cast volume pointer from generic volume to convection boundary
            const ConvectionBoundary *boundary = (ConvectionBoundary*)&(boundaries_[i]);

            double alpha = boundary->getAlpha();
            double S = surface_[i];
            double Text = boundary->getTExt();

            coefs[index_] -= alpha*S;
            coefs[n_nodes] -= alpha*S*Text;
        }
        else if (boundary_type == VType::fixed_T_boundary)
        {
            // cast volume pointer from generic volume to fixed T boundary
            const FixedTBoundary *boundary = (FixedTBoundary*)&(boundaries_[i]);
            
            double S = surface_[i];
            double d = distanceToVolume(boundary);

            coefs[index_] -= lambda_*S/d;
            coefs[n_nodes] -= lambda_*S/d*boundary->getT();
        }
        else
        {
            throw "Volume type not recognized at system assembly";
        }
    }

    // take into account internally generated heat (qv)
    coefs[n_nodes] -= qv_*volume_;
}


void SolidVolume::setLambda (double new_lambda)
{
    lambda_ = new_lambda;
}


double SolidVolume::getLambda  () const
{
    return lambda_;
}


double SolidVolume::getIndex () const
{
    return index_;
}

double SolidVolume::distanceToVolume (const Volume *other) const
{
    double sum = 0;

    if (other->volumeType() == solid)
    {
        const SolidVolume *casted_other = (SolidVolume*) other;

        for (int i = 0; i < PROBLEM_DIM; i++)
            sum += pow(casted_other->position_[i] - this->position_[i], 2);
    }
    else if (other->volumeType() == fixed_T_boundary)
    {
        const FixedTBoundary *casted_other = (FixedTBoundary*) other;

        for (int i = 0; i < PROBLEM_DIM; i++)
            sum += pow(casted_other->getCoordinate(i) - this->position_[i], 2);
    }
    else
    {
        throw "bad call to function distanceToVolume";
    }
    
    return sqrt(sum);
}


void SolidVolume::print (int index) const
{
    std::cout << " * (" << index << ")Solid volume:\n";
    std::cout << "\tVolume: " << volume_ << " m^3\n";
    std::cout << "\tInternal heat generated: " << qv_ << " W/m^3\n";
    std::cout << "\tLambda: " << lambda_ << " W/(K*m^2)\n";
    std::cout << "\tPosition: ";

    for (int i = 0; i < PROBLEM_DIM; i++)
        std::cout << (i != 0 ? ", " : "") << position_[i];
    std::cout << " m\n";

    std::cout << "\tSurfaces: ";

    for (int i = 0; i < PROBLEM_DIM*2; i++)
        std::cout << (i != 0 ? ", " : "") << surface_[i];
    std::cout << " m^2\n";

    std::cout << "\tBoundaries:\n";

    for (int i = 0; i < PROBLEM_DIM*2; i++)
    {
        VType boundary_type = boundaries_[i]->volumeType();

        if (boundary_type == solid)
        {
            SolidVolume* boundary = (SolidVolume*)boundaries_[i];
            std::cout << "\t   (" << i << ") Solid volume " << boundary->index_ << "\n";
        }
        else if (boundary_type == convection_boundary)
        {
            ConvectionBoundary* boundary = (ConvectionBoundary*)boundaries_[i];
            std::cout << "\t   (" << i << ") Convection boundary: alpha = " << boundary->getAlpha() <<
                " W/(K*m^2)    T_ext = " << boundary->getTExt() << " K\n";
        }
        else if (boundary_type == fixed_T_boundary)
        {
            FixedTBoundary* boundary = (FixedTBoundary*)boundaries_[i];
            std::cout << "\t   (" << i << ") Fixed T boundary: T = " << boundary->getT() << " K\n";
        }
        else
        {
            throw "boundary type not implemented at print";
        }
    }

    std::cout << std::endl;
}


////////////////////////////////////////////////////////////////


ConvectionBoundary::ConvectionBoundary (double T_ext, double alpha) :
        Volume(VType::convection_boundary), T_ext_(T_ext), alpha_(alpha)
{
    //
}


void ConvectionBoundary::setTExt (double new_T_ext)
{
    T_ext_ = new_T_ext;
}


double  ConvectionBoundary::getTExt () const
{
    return T_ext_;
}


void ConvectionBoundary::setAlpha (double new_alpha)
{
    alpha_ = new_alpha;
}


double ConvectionBoundary::getAlpha () const
{
    return alpha_;
}


void ConvectionBoundary::print (int index) const
{
    std::cout << " * (" << index << ") Convection boundary: alpha = " << alpha_ <<
            " W/(K*m^2)   T_ext = " << T_ext_ << " K\n";
    
    std::cout << std::endl;
}


////////////////////////////////////////////////////////////////


FixedTBoundary::FixedTBoundary (double T, const double *position) :
        Volume(VType::fixed_T_boundary), T_(T)
{
    for (int i = 0; i < PROBLEM_DIM; i++)
        position_[i] = position[i];
}

void FixedTBoundary::setT (double new_T)
{
    T_ = new_T;
}


double FixedTBoundary::getT () const
{
    return T_;
}


double FixedTBoundary::getCoordinate (int dimension) const
{
    return position_[dimension];
}


void FixedTBoundary::print (int index) const
{
    std::cout << " * (" << index << ")Fixed T boundary: T = " << T_ << " K \n";

    std::cout << std::endl;
}