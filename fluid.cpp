#include <fluid.h>
#include <transform.h>

/*
the i-th grid domained [i, i+1)

i-th dimension domained [0, _res[i+1])

0,1,2,...,K-1 dimension
*/

//initialize the fields and dimension.
void initialize(vector<int>*a){
    _res.resize(K + 1);
    num.resize(K + 1); 
    fields.resize(N_f);

    _res[0] = 1;
    num.push_back(1);

    for(int i = 0; i < K; i++){
        _res[i + 1] = (*a)[i];
        num[i + 1] = num[i] * _res[i + 1];
    }
    _velocity.resize(K * num[K]);
    _density.resize(num[K]);
    _heat.resize(num[K]);
}

//cell2index
int get_index(VectorXd c){
    int id = 0;

    /* x + x_num*y + x_num*y_num*z */
    for(int i = 0; i < K; i++)
        id += c(i) * num[i];
    
    return id;
}

//index2cell
VectorXd get_cell(int id){
    VectorXd c(K);
    for(int i = 0; i < K; i++){
        c(i) = id % _res[i + 1];
        id /= _res[i + 1];
    }
    return c;
}

// advect a single cell
VectorXf Advection(VectorXf &old_field, int d, VectorXf &velocity, VectorXd &c){
    int id = get_index(c);

    vector<float> Trace;
    vector<int>bound;
    Trace.resize(K);
    bound.resize(2 * K);

    VectorXf v(K);

    for(int i = 0; i < K; i++)
        v(i) = velocity(id * K + i);


    for(int i = 0; i < K; i++)
        Trace[i] = c(i) - _dt * v(i);

    // find the bound
    for(int i = 0; i < K; i++){
        bound[i] = max((int) (Trace[i] - 0.5), 0);
        bound[i + K] = min((int) (Trace[i] + 0.5)), _res[i + 1] - 1);
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
                weight *= (S >> i & 1) ? bound[i + K] + 0.5 - Trace[i] : Trace[i] - bound[i] - 0.5;
            else weight *= 0.5;
            p(i) = (S >> i & 1) ? bound[i + K] : bound[i];
        }
        int p_id = get_index(p);
        for(int i = 0; i < d; i++)
            ans(i) += weight * old_field(p_id * d + i);
    }
    return ans;
}

void calc_neighbors(int id, vector<int>&neighbors){
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

void BuildDiffusionMatrix(sparseMatrix<float>&Mat, int d){
    const float w = 0.9 / K;
    const float wp = 0.1;

    Mat.resize(d * num[K], d * num[k]);

    for(int i = 0; i < d * num[K]; i++)
        Mat(i, i) = 1;

    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);
        //filter
        for(int i = 0; i < d; i++){
            Mat(id * d + i, id * d + i) = wp;
            for(int j = 0; j < K; j++){
                if(c(j) != 0)
                    Mat(neighbors[j] * d + i, neighbors[j] * d + i) = w;
                if(c(j) != _res[j + 1] - 1)
                    Mat(neighbors[j + K] * d + i, neighbors[j + K] * d + i) = w;
            }
        }
    }
}

/* diffusion */
void Diffusion(VectorXf& old_fields, VectorXf& new_fields){
    new_fields = diffusion_matrix * old_fields;
}

void ComputeVelocity2Divergence(SparseMatrix<float>&Mat){
    float coe = -0.5;
    Mat.resize(num[K], K * num[K]);
    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);

        //Taylor series
        for(int i = 0; i < K; i++){
            if(c(i) != res[i + 1] - 1)
                Mat(id, K * neighbors[i + K] + i) += coe;
            else {
                Mat(id, K * id + i) += 2 * coe;
                Mat(id, K * neighbors[i] + i) -= coe;
            }
                
            
            if(c(i) != 0)
                Mat(id, K * neighbors[i] + i) -= coe;
            else {
                Mat(id, K * id + i) -= 2 * coe;
                Mat(id, K * neighbors[i + K] + i) += coe;
            }
        }
    }
}

void ComputePossion(SparseMatrix<float>&Mat){
    float coe = 1;
    Mat.resize(num[K], num[K]);
    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);

        Mat(id, id) = K * 2 * coe;
        for(int i = 0; i < K; i++){
            if(c(i) != res[i + 1] - 1)
                Mat(id, neighbors[i + K]) = -coe;
            if(c(i) != 0)
                Mat(id, neighbors[i]) = -coe;
        }
    }
}

void ComputePressure2Velocity(SparseMatrix<float>&Mat){
    float coe = 0.5;
    Mat.resize(K * num[K], num[K]);
    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);

        //Taylor series
        for(int i = 0; i < K; i++){
            if(c(i) != res[i + 1] - 1)
                Mat(id + i, neighbors[i + K] + i) += coe;
            else {
                Mat(id + i, id) += 2 * coe;
                Mat(id + i, neighbors[i]) -= coe;
            }
                
            
            if(c(i) != 0)
                Mat(id + i, neighbors[i]) -= coe;
            else {
                Mat(id + i, id) -= 2 * coe;
                Mat(id + i, Neighbors[i + K]) += coe;
            }
                
        }
    }
}

/* build the pressure matrix */
void BuildPressureMatrix(SparseMatrix<float>&Mat){
    SparseMatrix<float> Velocity2Divergence;
    SparseMatrix<float> Possion;
    SparseMatrix<float> Pressure2Velocity;

    ComputeVelocity2Divergence(Velocity2Divergence);
    ComputePossion(Possion);
    ComputePressure2Velocity(Pressure2Velocity);

    Mat = Pressure2Velocity * (Possion.inverse()) * Velocity2Divergence;
}


