#include "volume.h"
#include <stdexcept>


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

SolidVolume::SolidVolume (double volume, double alpha, double qv,
                          double *surfaces, Volume **boundaries, double index) :
        Volume(solid), alpha_(alpha), qv_(qv), volume_(volume), index_(index)
{
    for (int i = 0; i < PROBLEM_DIM*2; i++)
    {
        boundaries_[i] = boundaries[i];
        surface_[i] = surfaces[i];
    }
}

void SolidVolume::getEquation (const Volume *boundaries, std::vector<double> &coefs)
{
    for (int i = 0; i < PROBLEM_DIM*2; i++)
    {
        VType boundary_type = boundaries_[i]->volumeType();

        if (boundary_type == solid)
        {
            //
        }
        else if (boundary_type == convection_boundary)
        {
            //
        }
        else
        {
            throw std::logic_error("boundary type not recognized");
        }
    }
}

double SolidVolume::setAlpha (double new_alpha)
{
    alpha_ = new_alpha;
}

double SolidVolume::getAlpha () const
{
    return alpha_;
}

double SolidVolume::getIndex () const
{
    return index_;
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