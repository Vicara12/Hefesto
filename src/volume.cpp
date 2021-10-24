#include "volume.h"
#include <stdexcept>
#include <math.h>


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
                          double *surfaces, Volume **boundaries,
                          int index, double *position) :
        Volume(solid), lambda_(lambda), qv_(qv), volume_(volume), index_(index)
{
    for (int i = 0; i < PROBLEM_DIM*2; i++)
    {
        boundaries_[i] = boundaries[i];
        surface_[i] = surfaces[i];
    }

    for (int i = 0; i < PROBLEM_DIM; i++)
        position_[i] = position[i];
}


void SolidVolume::getEquation (const Volume *boundaries, double *coefs, int n_nodes)
{
    // for each boundary (two per dimension are assumed)
    for (int i = 0; i < PROBLEM_DIM*2; i++)
    {
        VType boundary_type = boundaries_[i]->volumeType();

        // the coefficients that need to be actualized
        // depend on the type of volume or boundary
        if (boundary_type == solid)
        {
            // cast volume pointer from generic volume to solid volume
            const SolidVolume *boundary = (SolidVolume*)&(boundaries[i]);

            int boundary_i = boundary->index_;
            // this is not exactly the average on the boundary of both volumes
            // but in the midlle point between their centers, but it's close enough
            double lambda = 2/(1/this->lambda_ + 1/boundary->lambda_);
            double S = surface_[i];
            double d = distanceToVolume(boundary);

            coefs[boundary_i] += lambda_*S/d;
            coefs[index_] -= lambda_*S/d;
        }
        else if (boundary_type == convection_boundary)
        {
            // cast volume pointer from generic volume to convection boundary
            const ConvectionBoundary *boundary = (ConvectionBoundary*)&(boundaries[i]);
        }
        else if (boundary_type == fixed_T_boundary)
        {
            // cast volume pointer from generic volume to fixed T boundary
            const FixedTBoundary *boundary = (FixedTBoundary*)&(boundaries[i]);
        }
        else
        {
            throw std::logic_error("boundary type not recognized");
        }
    }

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

double SolidVolume::distanceToVolume (const SolidVolume *other) const
{
    double sum = 0;

    for (int i = 0; i < PROBLEM_DIM; i++)
        sum += pow(other->position_[i] - this->position_[i], 2);
    
    return sqrt(sum);
}


////////////////////////////////////////////////////////////////


ConvectionBoundary::ConvectionBoundary (double T_ext, double alpha, double surface) :
        Volume(convection_boundary), T_ext_(T_ext), alpha_(alpha), surface_(surface)
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


void ConvectionBoundary::setSurface (double new_surface)
{
    surface_ = surface_;
}


double ConvectionBoundary::getSurface () const
{
    return surface_;
}


////////////////////////////////////////////////////////////////


FixedTBoundary::FixedTBoundary (double T) : Volume(fixed_T_boundary), T_(T)
{
    //
}

void FixedTBoundary::setT (double new_T)
{
    T_ = new_T;
}


double FixedTBoundary::getT () const
{
    return T_;
}