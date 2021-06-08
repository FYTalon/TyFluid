#ifndef _FLUID_H_
#define _FLUID_H_ 
#include <Eigen/Eigen>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <cstring>
#include <string>

#define _velocity fields[0]
#define _density fields[1]
#define _heat fields[2]
#define DENSITY_AMOUNT 2
#define HEAT_AMOUNT 200
#define VELOCITY_AMOUNT 20

#define SCALE 128
/*number of fields*/
#define N_f 3

using namespace std;
using namespace Eigen;


class Fluid{
public:
    /*number of dimension*/
    int K;

    /*time step*/
    float _dt;
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
    SparseMatrix<float> Velocity2Divergence;
    SparseMatrix<float> Possion;
    SparseMatrix<float> Pressure2Velocity;

    void initialize(vector<int>*a);

    int get_index(VectorXd c);
    VectorXd get_cell(int id);
    void calc_neighbors(int id, vector<int>&neighbors);

    void AdvectionAll(VectorXf& old_field, VectorXf& new_field, int d);
    VectorXf Advection(VectorXf &old_field, int d, VectorXf &velocity, int id);

    void BuildDiffusionMatrix(SparseMatrix<float>&Mat, int d);
    void Diffusion(VectorXf& old_field, VectorXf& new_field, int d);

    void ComputeVelocity2Divergence(SparseMatrix<float>&Mat);
    void ComputePossion(SparseMatrix<float>&Mat);
    void ComputePressure2Velocity(SparseMatrix<float>&Mat);
    void BuildPressureMatrix(SparseMatrix<float>&Mat);
    void Pressure(VectorXf& old_field, VectorXf& new_field);
    void solvePoisson(VectorXf& x, VectorXf& b);

    float AddBuoyancy(int id);
    void AddVorticity();
    void External();

    void step();
};

void WriteMat2File(SparseMatrix<float>&Mat, string s);
void WriteFile2Mat(SparseMatrix<float>&Mat, string s);
void test_advection();
void test_diffusion();
void test_projection();
void test();

#endif