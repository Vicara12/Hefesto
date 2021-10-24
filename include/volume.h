#ifndef VOLUME_H_
#define VOLUME_H_

#include <vector>

#define PROBLEM_DIM 2

enum VType {solid, convection_boundary};
enum Positions {E, W, U, D, F, B};




// abstract class for a generic volume
class Volume
{
public:

    virtual std::vector<double> getEquation (const Volume *boundaries) = 0;

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

    SolidVolume (double volume, double alpha, double qv, double *surfaces,
                 Volume **boundaries, double index);

    // coefst is the index equation with the format  sum(a_i * x_i) = b_i
    void getEquation (const Volume *boundaries, std::vector<double> &coefs);

    double setAlpha (double new_alpha);

private:

    double getAlpha () const;
    double getIndex () const;

    const Volume *boundaries_ [PROBLEM_DIM*2];
    double volume_;
    double qv_;
    double surface_ [PROBLEM_DIM*2];
    double alpha_;
    double index_;
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
    std::vector<double> getEquation (const Volume *boundaries);

    double T_ext_;
    double alpha_;
    double surface_;
};

#endif