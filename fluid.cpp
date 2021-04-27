#include "fluid.h"
#include "transform.h"

/*
the i-th velocity at i+1.5 (0 v_0 v_1 ... v_n-1 0)

i-th dimension domained [0, _res[i+1]]

0,1,2,...,K-1 dimension
*/

//initialize the fields and dimension.
void Fluid::initialize(vector<int>*a){
    _res.resize(K);
    num.resize(K + 1);
    _peeled.resize(K + 1);
    fields.resize(N_f);
    
    num[0] = 1;
    _peeled[0] = 1;
    for(int i = 0; i < K; i++){
        _res[i] = (*a)[i];
        _peeled[i + 1] = _res[i] - 2;
        num[i + 1] = num[i] * _peeled[i + 1];
        _x = max(_x, (float)_res[i]);
    }
    _velocity.resize(K * num[K]);
    _velocity.setZero();
    _density.resize(num[K]);
    _density.setZero();
    _heat.resize(num[K]);
    _heat.setZero();

    _buoyancy = 0.1f;
    _vorticityEPS = 2.0f;
    _h = max(1.0f, _x / 64.0f);
    _x = 1.0 / _x;
}

//cell2index
int Fluid::get_index(VectorXd c){
    int id = 0;

    /* x + x_num*y + x_num*y_num*z */
    for(int i = 0; i < K; i++)
        id += c(i) * num[i];
    
    return id;
}

//index2cell
VectorXd Fluid::get_cell(int id){
    VectorXd c(K);
    for(int i = 0; i < K; i++){
        c(i) = id % _peeled[i + 1];
        id /= _peeled[i + 1];
    }
    return c;
}

void Fluid::AdvectionAll(VectorXf& old_field, VectorXf& new_field, int d){
    for(int i = 0; i < num[K]; i++){
        VectorXf c = Advection(old_field, d, _velocity, i);
        for(int j = 0; j < d; j++)
            new_field(i * d + j) = c(j);
    }
}
// advect a single cell
VectorXf Fluid::Advection(VectorXf &old_field, int d, VectorXf &velocity, int id){
    VectorXd c = get_cell(id);

    vector<float> Trace;
    vector<int>bound;
    Trace.resize(K);
    bound.resize(2 * K);

    VectorXf v(K);

    for(int i = 0; i < K; i++)
        v(i) = velocity(id * K + i);


    for(int i = 0; i < K; i++)
        Trace[i] = c(i) + 1.5 - _dt * v(i);

    // find the bound
    for(int i = 0; i < K; i++){
        bound[i] = max((int)(Trace[i] - 1.5), 0);
        bound[i + K] = min((int)(Trace[i] - 0.5), _res[i] - 3);
    }
    
    VectorXf ans(d);

    for(int i = 0; i < d; i++)
        ans(i) = 0;
    for(int S = 0; S < (1 << K); S++){
        double weight = 1;
        VectorXd p(K);
        // calculate S cell's weight
        for(int i = 0; i < K; i++){
            if(bound[i] != bound[i + K])
                weight *= (S >> i & 1) ? bound[i + K] + 1.5 - Trace[i] : Trace[i] - bound[i] - 1.5;
            else weight *= 0.5;
            p(i) = (S >> i & 1) ? bound[i + K] : bound[i];
        }
        int p_id = get_index(p);
        for(int i = 0; i < d; i++)
            ans(i) += weight * old_field(p_id * d + i);
    }
    return ans;
}

void Fluid::calc_neighbors(int id, vector<int>&neighbors){
    neighbors.resize(2 * K);
    VectorXd c = get_cell(id);
    //get neighbors every dimension
    for(int i = 0; i < K; i++){
        c(i) -= 1;
        neighbors[i] = get_index(c);
        c(i) += 2;
        neighbors[i + K] = get_index(c);
        c(i) -= 1;
    }
}

/* build the diffusion Matrix*/
/* d is the dimension */

void Fluid::BuildDiffusionMatrix(SparseMatrix<float>&Mat, int d){
    const float w = 0.9 / K;
    const float wp = 0.1;

    Mat.resize(d * num[K], d * num[K]);
    Mat.setZero();

    for(int i = 0; i < d * num[K]; i++)
        Mat.coeffRef(i, i) = 1;

    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);
        //filter
        for(int i = 0; i < d; i++){
            Mat.coeffRef(id * d + i, id * d + i) = wp;
            for(int j = 0; j < K; j++){
                if(c(j) != 0)
                    Mat.coeffRef(neighbors[j] * d + i, neighbors[j] * d + i) = w;
                if(c(j) != _peeled[j + 1] - 1)
                    Mat.coeffRef(neighbors[j + K] * d + i, neighbors[j + K] * d + i) = w;
            }
        }
    }
}

/* diffusion */
void Fluid::Diffusion(VectorXf& old_field, VectorXf& new_field, int d, float coe){
    new_field = coe * diffusion_matrix[d] * old_field;
}

// didn't folded away
void Fluid::ComputeVelocity2Divergence(SparseMatrix<float>&Mat){
    float coe = 0.5;
    Mat.resize(num[K], K * num[K]);
    Mat.setZero();
    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);

        //Taylor series
        for(int i = 0; i < K; i++){
            if(c(i) != _peeled[i + 1] - 1)
                Mat.coeffRef(id, K * neighbors[i + K] + i) += coe;
            else 
                Mat.coeffRef(id, K * neighbors[i] + i) += coe;
            
            if(c(i) != 0)
                Mat.coeffRef(id, K * neighbors[i] + i) -= coe;
            else 
                Mat.coeffRef(id, K * neighbors[i + K] + i) -= coe;
        }
    }
}

void Fluid::ComputePossion(SparseMatrix<float>&Mat){
    float coe = -1;
    Mat.resize(num[K], num[K]);
    Mat.setZero();
    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);
        Mat.coeffRef(id, id) = 0;
        for(int i = 0; i < K; i++){
            if(c(i) != _peeled[i + 1] - 1 && !obstacles[neighbors[i + K]])
                Mat.coeffRef(id, neighbors[i + K]) = -coe,
                Mat.coeffRef(id, id) += coe;
            if(c(i) != 0 && !obstacles[neighbors[i]])
                Mat.coeffRef(id, neighbors[i]) = -coe,
                Mat.coeffRef(id, id) += coe;
        }
    }
}

void Fluid::ComputePressure2Velocity(SparseMatrix<float>&Mat){
    float coe = -0.5;
    Mat.resize(K * num[K], num[K]);
    Mat.setZero();
    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);

        //Taylor series
        for(int i = 0; i < K; i++){
            if(c(i) != _peeled[i + 1] - 1)
                Mat.coeffRef(id + i, neighbors[i + K] + i) += coe;
            else {
                Mat.coeffRef(id + i, id) += 2 * coe;
                Mat.coeffRef(id + i, neighbors[i]) -= coe;
            }
                
            
            if(c(i) != 0)
                Mat.coeffRef(id + i, neighbors[i]) -= coe;
            else {
                Mat.coeffRef(id + i, id) -= 2 * coe;
                Mat.coeffRef(id + i, neighbors[i + K]) += coe;
            }
                
        }
    }
}

/* build the pressure matrix */
void Fluid::BuildPressureMatrix(SparseMatrix<float>&Mat){
    SparseMatrix<float> Velocity2Divergence;
    SparseMatrix<float> Possion;
    SparseMatrix<float> Pressure2Velocity;

    ComputeVelocity2Divergence(Velocity2Divergence);
    ComputePossion(Possion);
    ComputePressure2Velocity(Pressure2Velocity);

    Mat = Pressure2Velocity * (Possion.inverse()) * Velocity2Divergence;
}

/* pressure projection */
void Fluid::Pressure(VectorXf& old_field, VectorXf& new_field){
    new_field = pressure_matrix * old_field;
}

/* single cell Buoyancy */
float Fluid::AddBuoyancy(int id){
    float force = _buoyancy * _heat(id);
    return force;
}

/* vorticity */
void Fluid::AddVorticity(){
    VectorXf vor;
    vor.resize(K * num[K]);
    vor.setZero();
    float coe = _vorticityEPS * _h;

    VectorXf len;
    len.resize(num[K]);
    len.setZero();
    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c;
        c = get_cell(id);
        calc_neighbors(id, neighbors);
        VectorXf dx;
        dx.resize(K);
        dx.setZero();
        
        for(int i = 0; i < K; i++){
            int up = obstacles[neighbors[i + K]] || neighbors[i + K] == _res[i] - 3 ? id : neighbors[i + K];
            int down = obstacles[neighbors[i]] || neighbors[i] == 0 ? id : neighbors[i];
            dx(i) = (up == id || down == id) ? 1.0f : 2.0f;
            switch (K){
                case 2:
                    vor(K * id) += i ? -(_velocity(K * up) - _velocity(K * down)) / dx(i)
                                    : (_velocity(K * up + 1) - _velocity(K * down + 1)) / dx(i);
                    break;
                case 3:
                    int j = (i + 2) % 3;
                    int k = (i + 1) % 3;
                    vor(K * id + j) += (_velocity(K * up + k) - _velocity(K * down + k)) / dx(i);
                    vor(K * id + k) -= (_velocity(K * up + j) - _velocity(K * down + j)) / dx(i);
            }
        }
        for(int i = 0; i < K; i++)
            len(id) += vor(K * id + i) * vor(K * id + i);
        len(id) = sqrt(len(id));
    }
    for(int id = 0; id < num[K]; id++) if(!obstacles[id]){
        vector<int>neighbors;
        VectorXd c;
        c = get_cell(id);
        calc_neighbors(id, neighbors);
        VectorXf dx;
        dx.resize(K);
        dx.setZero();
        VectorXf N;
        N.resize(K);
        N.setZero();
        float L = 0;
        for(int i = 0; i < K; i++){
            int up = obstacles[neighbors[i + K]] || neighbors[i + K] == _res[i] - 3 ? id : neighbors[i + K];
            int down = obstacles[neighbors[i]] || neighbors[i] == 0 ? id : neighbors[i];
            dx(i) = (up == id || down == id) ? 1.0f : 2.0f;
            N(i) = (len(up) - len(down)) / dx(i);
            L += N(i) * N(i);
        }
        L = sqrt(L);
        if(L > 0.0f){
            N = L * N;
            switch (K) {
                case 2:
                    _velocity(id * K) += _dt * N(1) * vor(K * id);
                    _velocity(id * K) -= _dt * N(0) * vor(K * id);
                    break;
                case 3:
                    _velocity(id * K) += _dt * coe * (N(1) * vor(K * id + 2) - N(2) * vor(K * id + 1));
                    _velocity(id * K + 1) += _dt * coe * (N(2) * vor(K * id) - N(0) * vor(K * id + 2));
                    _velocity(id * K + 2) += _dt * coe * (N(0) * vor(K * id + 1) - N(1) * vor(K * id));
            }
        }
    }
}

void Fluid::External(){
    for(int i = 0; i < num[K]; i++)
        _velocity(i * K + 1) += _dt * AddBuoyancy(i);

    AddVorticity();
}