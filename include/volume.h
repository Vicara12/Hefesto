#ifndef VOLUME_H_
#define VOLUME_H_

#include <vector>

#define PROBLEM_DIM 2

enum VType {solid, convection_boundary, fixed_T_boundary};
enum Positions {E, W, U, D, F, B};




// abstract class for a generic volume
class Volume
{
public:

    virtual void getEquation (const Volume *boundaries,
                              double *coefs,
                              int n_nodes) = 0;

    VType volumeType () const;

protected:

    Volume (VType vol_type);

private:

    VType vol_type_;
};



// class for a solid volume (conduction heat transfer)
class SolidVolume : public Volume
{
public:

    SolidVolume (double volume, double lambda, double qv, double *surfaces,
                 Volume **boundaries, int index, double *position);

    // coefst is the index equation with the format  sum(a_i * x_i) = b_i
    void getEquation (const Volume *boundaries, double *coefs, int n_nodes);

    void setLambda  (double new_lambda);

private:

    double getLambda  () const;
    double getIndex () const;
    // get distance between two solid volume centers
    double distanceToVolume (const SolidVolume *other) const;

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

    ConvectionBoundary (double T_ext, double alpha, double surface);

    void setTExt (double new_T_ext);
    double  getTExt () const;

    void setAlpha (double new_alpha);
    double getAlpha () const;

    void setSurface (double new_surface);
    double getSurface () const;

private:
    // this method does nothing because it's a boundary, not a volume
    void getEquation (const Volume *boundaries, double *coefs, int n_nodes);

    double T_ext_;
    double alpha_;
    double surface_;
};



// class for fixed temperature boundaries
class FixedTBoundary : public Volume
{
public:

    FixedTBoundary (double T);

    void setT (double new_T);
    double getT () const;

private:
    // this method does nothing because it's a boundary, not a volume
    void getEquation (const Volume *boundaries, double *coefs, int n_nodes);

    double T_;
};


#endif