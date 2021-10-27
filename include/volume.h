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

    SolidVolume (double volume, double lambda, double qv, const double *surfaces,
                 int index, const double *position);
    
    void setBoundaries (const Volume **boundaries);

    // setBoundaries must be called before this method
    // coefst is the index equation with the format  sum(a_i * x_i) = b_i
    void getEquation (DoubleVector &coefs);

    void setLambda  (double new_lambda);

    void print (int index) const override;

private:

    double getLambda  () const;
    double getIndex () const;
    // get distance from this volume to other solid volume or fixed T volume
    double distanceToVolume (const Volume *other) const;

    const Volume *boundaries_ [PROBLEM_DIM*2];
    double volume_;
    double qv_;
    double surface_ [PROBLEM_DIM*2];
    double lambda_;
    int index_;
    double position_ [PROBLEM_DIM];
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