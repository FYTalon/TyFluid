#ifndef _FLUID_H_
#define _FLUID_H_

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/SparseCore>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>

#define _velocity fields[0]
#define _density fields[1]
#define _heat fields[2]
/*number of fields*/
#define N_f 3

using namespace std;
using namespace Eigen;


/*number of dimension*/
int K;

/*time step*/
int _dt;


class Fluid{
public:
    /*number of cells every dimension*/
    vector<int> _res, _peeled;
    /*number of cells in i dimensions*/
    vector<int> num;
    /*fields need calculate, e.g velocity, density, heat*/
    vector<VectorXf> fields;
    vector<int> dim;
    /*obstacles*/
    vector<int> obstacles;

    /*buoyancy and vorticity and dx*/
    float _buoyancy, _x, _vorticityEPS, _h;
    /* fields' diffusion matrix */
    vector<SparseMatrix<float> > diffusion_matrix;

    /* velocity's pressure matrix*/
    SparseMatrix<float> pressure_matrix;

    void initialize(vector<int>*a);

    int get_index(VectorXd c);
    VectorXd get_cell(int id);
    void calc_neighbors(int id, vector<int>&neighbors);

    void AdvectionAll(VectorXf& old_field, VectorXf& new_field, int d);
    VectorXf Advection(VectorXf &old_field, int d, VectorXf &velocity, int id);

    void BuildDiffusionMatrix(SparseMatrix<float>&Mat, int d);
    void Diffusion(VectorXf& old_field, VectorXf& new_field, int d, float coe);

    void ComputeVelocity2Divergence(SparseMatrix<float>&Mat);
    void ComputePossion(SparseMatrix<float>&Mat);
    void ComputePressure2Velocity(SparseMatrix<float>&Mat);
    void BuildPressureMatrix(SparseMatrix<float>&Mat);
    void Pressure(VectorXf& old_field, VectorXf& new_field);

    float AddBuoyancy(int id);
    void AddVorticity();
    void External();

};

#endif