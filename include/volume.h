#ifndef VOLUME_H_
#define VOLUME_H_

#include "definitions.h"


// abstract class for a generic volume
class Volume
{
public:

    virtual void print(int index) const = 0;

    VType volumeType () const;

protected:

    Volume (VType vol_type);

private:

    VType vol_type_;
};



// class for fixed temperature boundaries
class FixedTBoundary : public Volume
{
public:

    FixedTBoundary (double T, double distance);

    void setT (double new_T);
    double getT () const;

    void print (int index) const override;

    double getDistance () const;

private:

    double T_;
    double d_; // distance to the assigned volume's center
};



// class for a solid volume (conduction heat transfer)
class SolidVolume : public Volume
{
public:

    SolidVolume (unsigned int n_dimensions, double volume, double lambda, double qv, 
                 const DoubleVector &surfaces, int index, const DoubleVector &position);
    
    void setBoundaries (const std::vector<const Volume*> &boundaries);

    // setBoundaries must be called before this method
    // coefst is the index equation with the format  sum(a_i * x_i) = b_i
    void getEquation (DoubleVector &coefs);

    void setLambda  (double new_lambda);

    void print (int index) const override;

    // Checks energy balance for the node at temperature T
    double checkEnergyBalance (const DoubleVector &T) const;

private:

    double getLambda  () const;
    double getIndex () const;
    // get distance from this volume to other solid volume or fixed T volume
    double distanceToVolume (const Volume *other) const;

    unsigned int n_dimensions_;
    std::vector<const Volume*> boundaries_;
    double volume_;
    double qv_;
    DoubleVector surface_;
    double lambda_;
    int index_;
    DoubleVector position_;
};



// class for conduction boundaries
class ConvectionBoundary : public Volume
{
public:

    ConvectionBoundary (double T_ext, double alpha);

    void setTExt (double new_T_ext);
    double  getTExt () const;

    void setAlpha (double new_alpha);
    double getAlpha () const;

    void print (int index) const override;

private:

    double T_ext_;
    double alpha_;
};


#endif