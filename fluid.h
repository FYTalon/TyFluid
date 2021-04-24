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


/*number of cells every dimension*/
vector<int> _res;
/*number of cells in i dimensions*/
vector<int> num;
/*fields need calculate, e.g velocity, density, heat*/
vector<VectorXf> fields;

/* fields' diffusion matrix */
vector<SparseMatrix<float> > diffusion_matrix;

/* velocity's pressure matrix*/
SparseMatrix<float> pressure_matrix;

void initialize(vector<int>*a);

int get_index(VectorXd c);
VectorXd get_cell(int id);


VectorXf Advection(VectorXf &old_field, int d, VectorXf &velocity, VectorXd &c);

void BuildDiffusionMatrix(sparseMatrix<float>&Mat, int d);
void Diffusion(VectorXf& old_fields, VectorXf& new_fields);

void BuildPressureMatrix(sparseMatrix<float>&Mat, int d);
void Pressure(VectorXf& old_fields, VectorXf& new_fields);

void External();

void Vorticity();