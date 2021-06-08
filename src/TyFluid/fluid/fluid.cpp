#include "fluid.h"
#include "transform.h"

/*
the i-th velocity at i+1.5 (0 v_0 v_1 ... v_n-1 0)

i-th dimension domained [0, _res[i+1]]

0,1,2,...,K-1 dimension
*/

void WriteMat(SparseMatrix<float>&Mat, string s){
    FILE* out = fopen(s.c_str(), "wb");
    int r = Mat.rows(), c = Mat.cols();
    float v;
    fwrite((void*)&r, sizeof(int), 1, out);
    fwrite((void*)&c, sizeof(int), 1, out);
    for(int k = 0; k < Mat.outerSize(); k++)
        for(SparseMatrix<float>::InnerIterator it(Mat, k); it ; ++it){
            //if(fabs(it.value()) < 1e-8) continue ;
            r = it.row();
            c = it.col();
            v = it.value();
            fwrite((void*)&r, sizeof(int), 1, out);
            fwrite((void*)&c, sizeof(int), 1, out);
            fwrite((void*)&v, sizeof(float), 1, out);
        }
    fclose(out);
}

void ReadMat(SparseMatrix<float>&Mat, string s){
    FILE* in = fopen(s.c_str(), "rb");
    int r, c;
    float v;
    fread((void*)&r, sizeof(int), 1, in);
    fread((void*)&c, sizeof(int), 1, in);
    Mat.resize(r, c);
    Mat.setZero();
    while(fread((void*)&r, sizeof(int), 1, in) > 0){
        fread((void*)&c, sizeof(int), 1, in);
        fread((void*)&v, sizeof(float), 1, in);
        Mat.coeffRef(r, c) = v;
    }
    fclose(in);
}

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

    printf("num: ");
    for(int i = 0; i <= K; i++)
        printf("%d ", num[i]);
    printf("\n");

    printf("_peeled: ");
    for(int i = 0; i <= K; i++)
        printf("%d ", _peeled[i]);
    printf("\n");

    printf("_res: ");
    for(int i = 0; i < K; i++)
        printf("%d ", _res[i]);
    printf("\n");

    obstacles.resize(K * num[K]);
    for(int i = 0; i < K * num[K]; i++)
        obstacles[i] = 0;
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

    diffusion_matrix.resize(2);
    double start, end;
    char buf[100];
    string suf, pre;
    sprintf(buf, "%d", SCALE);
    suf = buf;
    memset(buf, 0, sizeof(buf));
    sprintf(buf, "./matrix/");
    pre = buf;
    /*start = clock();
    BuildDiffusionMatrix(diffusion_matrix[0], K);
    BuildDiffusionMatrix(diffusion_matrix[1], 1);
    
    WriteMat(diffusion_matrix[0], pre + (string)"dm0." + suf);
    WriteMat(diffusion_matrix[1], pre + (string)"dm1." + suf);
    end = clock();
    printf("diffusion : %lf\n", end - start);
    start = clock();
    BuildPressureMatrix(pressure_matrix);
    WriteMat(Velocity2Divergence, pre + (string)"v2d." + suf);
    WriteMat(Possion, pre + (string)"p." + suf);
    WriteMat(Pressure2Velocity, pre + (string)"p2v." + suf);
    //WriteMat(pressure_matrix, pre + "pm.64" + suf);
    end = clock();
    printf("project : %lf\n", end - start);*/

    start = clock();
    ReadMat(diffusion_matrix[0], pre + (string)"dm0." + suf);
    ReadMat(diffusion_matrix[1], pre + (string)"dm1." + suf);
    end = clock();
    printf("diffusion : %lf\n", end - start);

    start = clock();
    ReadMat(Velocity2Divergence, pre + (string)"v2d." + suf);
    ReadMat(Possion, pre + (string)"p." + suf);
    ReadMat(Pressure2Velocity, pre + (string)"p2v." + suf);
    end = clock();
    printf("project : %lf\n", end - start);
    /*start = clock();
    BuildPressureMatrix(pressure_matrix);
    end = clock();
    printf("project : %lf\n", end - start);
    start = clock();
    ReadMat(pressure_matrix, "pm.64");
    end = clock();
    printf("project : %lf\n", end - start);*/
    //cout << Possion << endl;
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

void test_advection(){
    Fluid f;
    vector<int>a;
    f.K = 2;
    f._dt = 0.01;
    a.push_back(5);
    a.push_back(5);
    f.initialize(&a);

    // 0 0 0 0 0 
    // 0       0
    // 0       0
    // 0 0 0 0 0

    f._velocity[0] = 1; f._velocity[1] = 1;
    f._velocity[3] = 1; f._velocity[4] = -1;
    f._velocity[5] = 1; f._velocity[9] = 3;
    f._velocity[8] = 3;

    VectorXf ov;
    ov.resize(f.K * f.num[f.K]);
    ov.setZero();
    f.AdvectionAll(f._velocity, ov, f.K);

    printf("_velocity:\n");
    cout << f._velocity << endl;
    printf("ov:\n");
    cout << ov << endl;
}

void Fluid::AdvectionAll(VectorXf& old_field, VectorXf& new_field, int d){
    for(int i = 0; i < num[K]; i++){
        VectorXf c = Advection(old_field, d, _velocity, i);
        //cout << i << " " <<  c << endl;
        for(int j = 0; j < d; j++)
            new_field(i * d + j) = c(j);

        /*VectorXf f(d);
        for(int j = 0; j < d; j++)
            f(j) = old_field(i * d + j);
        cout << "prev:" << endl;
        cout << f << endl;
        cout << "cur:" << endl;
        cout << c << endl;*/
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

    //printf("id=%d: ", id);
    for(int i = 0; i < K; i++){
        Trace[i] = c(i) + 1 - _dt * v(i) / _x;
        Trace[i] = Trace[i] < 1.5 ? 1.5 : Trace[i];
        Trace[i] = Trace[i] > _res[i] - 2.5 ? _res[i] - 2.5 : Trace[i];
        //cout << "trace" << Trace[i] << endl;
    }
        

    // find the bound
    for(int i = 0; i < K; i++){
        bound[i] = (int) Trace[i] - 1;
        bound[i + K] = bound[i] != _peeled[i + 1] - 1 ? bound[i] + 1 : bound[i];
        //cout << "left" << bound[i] << endl;
        //cout << "right" << bound[i + K] << endl;
    }
    
    VectorXf ans(d);

    //interpolate
    for(int i = 0; i < d; i++)
        ans(i) = 0;
    for(int S = 0; S < (1 << K); S++){
        double weight = 1;
        VectorXd p(K);
        // calculate S cell's weight
        for(int i = 0; i < K; i++){
            if(bound[i] != bound[i + K])
                weight *= (S >> i & 1) ? Trace[i] - bound[i] - 1 : bound[i] + 2 - Trace[i];
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
    //printf("id = %d\n", id);
    for(int i = 0; i < K; i++){
        c(i) -= 1;
        neighbors[i] = get_index(c);
        //printf("neighbor down: %d\n", neighbors[i]);
        c(i) += 2;
        neighbors[i + K] = get_index(c);
        
        //printf("neighbor up: %d\n", neighbors[i + K]);
        c(i) -= 1;
    }
}

/* build the diffusion Matrix*/
/* d is the dimension */

void test_diffusion(){
    Fluid f;
    vector<int>a;
    f.K = 2;
    f._dt = 0.01;
    a.push_back(5);
    a.push_back(5);
    f.initialize(&a);

    SparseMatrix<float>dmatrix;
    f.BuildDiffusionMatrix(dmatrix, f.K);

    for(int k = 0; k < f.K * f.num[f.K]; k++){
        for(SparseMatrix<float>::InnerIterator it(dmatrix, k); it ; ++it)
            cout << it.value() << " ";
        cout << endl;
    }

}

void Fluid::BuildDiffusionMatrix(SparseMatrix<float>&Mat, int d){
    const float w = 0.15;
    const float wp = 1.0 - w * 2 * K;
    
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
                    Mat.coeffRef(id * d + i, neighbors[j] * d + i) = w;
                if(c(j) != _peeled[j + 1] - 1)
                    Mat.coeffRef(id * d + i, neighbors[j + K] * d + i) = w;
            }
        }
    }
}

/* diffusion */
void Fluid::Diffusion(VectorXf& old_field, VectorXf& new_field, int d){
    new_field = diffusion_matrix[d] * old_field;
}

// didn't folded away
void Fluid::ComputeVelocity2Divergence(SparseMatrix<float>&Mat){
    float coe = 0.5 * _x;
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
            else{
                if(!i)
                    Mat.coeffRef(id, K * id + i) -= coe;
                else 
                    Mat.coeffRef(id, K * id + i) += coe;
            } 
                
            
            if(c(i) != 0)
                Mat.coeffRef(id, K * neighbors[i] + i) -= coe;
            else {
                if(!i)
                    Mat.coeffRef(id, K * id + i) += coe;
                else 
                    Mat.coeffRef(id, K * id + i) -= coe;
            }
                
        }
    }
}

void Fluid::ComputePossion(SparseMatrix<float>&Mat){
    float coe = 1;
    Mat.resize(num[K], num[K]);
    Mat.setZero();
    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);
        //Mat.coeffRef(id, id) = -2 * K * coe;
        for(int i = 0; i < K; i++){
            if(c(i) != _peeled[i + 1] - 1){
                Mat.coeffRef(id, neighbors[i + K]) = coe;
                Mat.coeffRef(id, id) -= coe;
            }
            if(c(i) != 0){
                Mat.coeffRef(id, neighbors[i]) = coe;
                Mat.coeffRef(id, id) -= coe;
            }
        }
    }
    //cout << Mat << endl;
}

void Fluid::ComputePressure2Velocity(SparseMatrix<float>&Mat){
    float coe = -0.5 / _x;
    Mat.resize(K * num[K], num[K]);
    Mat.setZero();
    for(int id = 0; id < num[K]; id++){
        vector<int>neighbors;
        VectorXd c = get_cell(id);
        calc_neighbors(id, neighbors);

        //Taylor series
        for(int i = 0; i < K; i++){
            if(c(i) != _peeled[i + 1] - 1)
                Mat.coeffRef(id * K + i, neighbors[i + K]) += coe;
            else {
                Mat.coeffRef(id * K + i, id) += 2 * coe;
                Mat.coeffRef(id * K + i, neighbors[i]) += -coe;
            }
                
            
            if(c(i) != 0)
                Mat.coeffRef(id * K + i, neighbors[i]) += -coe;
            else {
                Mat.coeffRef(id * K + i, id) += - 2 * coe;
                Mat.coeffRef(id * K + i, neighbors[i + K]) += coe;
            }
                
        }
    }
}

void test_projection(){
    Fluid f;
    vector<int>a;
    f.K = 1;
    f._dt = 0.1;
    a.push_back(7);
    f.initialize(&a);

    //cout << f._velocity << endl;
    f._velocity[0] = 2;
    
    f.BuildPressureMatrix(f.pressure_matrix);
    f._velocity += f.pressure_matrix * f._velocity;
    //cout << f.pressure_matrix << endl;
    // 0 2 2 1 0 0 0
    cout << f._velocity << endl;
}

/* build the pressure matrix */
void Fluid::BuildPressureMatrix(SparseMatrix<float>&Mat){

    ComputeVelocity2Divergence(Velocity2Divergence);
    //cout << Velocity2Divergence << endl;
    ComputePossion(Possion);
    //cout << Possion << endl;
    ComputePressure2Velocity(Pressure2Velocity);
    //cout << Possion << endl;
    /*SparseMatrix<float> V;
    SimplicialCholesky<SparseMatrix<float> > Solver_sparse(Possion);
    V = Solver_sparse.solve(Velocity2Divergence);
    Mat = Pressure2Velocity * V;*/
    /*VectorXf nv = Velocity2Divergence * Pressure2Velocity * V * _velocity;
    cout << nv << endl;*/
    //cout << Pressure2Velocity * V * _velocity << endl;
}

void Fluid::solvePoisson(VectorXf& x, VectorXf& b){
    SparseMatrix<float>D, B;
    VectorXf f;
    D.resize(num[K], num[K]);
    D.setZero();
    for(int i = 0; i < num[K]; i++)
        D.coeffRef(i, i) = 1.0 / Possion.coeffRef(i, i);
    B = D * Possion;
    f = D * b;
    x.setZero();

    for(int i = 0; i < num[K]; i++)
        B.coeffRef(i, i) -= 1;
    for(int i = 0; i < 1; i++)
        x = f - B * x;
}

/* pressure projection */
void Fluid::Pressure(VectorXf& old_field, VectorXf& new_field){
    //new_field = pressure_matrix * old_field;
    SimplicialCholesky<SparseMatrix<float> > solver;
    solver.compute(Possion);
    new_field = solver.solve(Velocity2Divergence * old_field);
    /*VectorXf field = Velocity2Divergence * old_field;
    solvePoisson(new_field, field);*/
    //VectorXf check = Possion * new_field;
    
    /*cout << check.norm() << endl;
    check = Velocity2Divergence * old_field;
    cout << check.norm() << endl;
    cout << endl;*/
    new_field = old_field + Pressure2Velocity * new_field;
    
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
            int up = c(i) == _peeled[i + 1] - 1 ? id : neighbors[i + K];
            int down = c(i) == 0 ? id : neighbors[i];
            dx(i) = (up == id || down == id) ? 1.0f : 2.0f;
            switch (K){
                case 2:
                    vor(K * id) += i ? -(_velocity(K * up) - _velocity(K * down)) / dx(i) / _x
                                    : (_velocity(K * up + 1) - _velocity(K * down + 1)) / dx(i) / _x;
                    break;
                case 3:
                    int j = (i + 2) % 3;
                    int k = (i + 1) % 3;
                    vor(K * id + j) += (_velocity(K * up + k) - _velocity(K * down + k)) / dx(i) / _x;
                    vor(K * id + k) -= (_velocity(K * up + j) - _velocity(K * down + j)) / dx(i) / _x;
            }
        }
        for(int i = 0; i < K; i++)
            len(id) += vor(K * id + i) * vor(K * id + i);
        len(id) = sqrt(len(id));
    }
    //cout << "flag" << endl;
    for(int id = 0; id < num[K]; id++) {
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
            int up = c(i) == _peeled[i + 1] - 1 ? id : neighbors[i + K];
            int down = c(i) == 0 ? id : neighbors[i];
            dx(i) = (up == id || down == id) ? 1.0f : 2.0f;
            N(i) = (len(up) - len(down)) / dx(i) / _x;
            L += N(i) * N(i);
        }
        L = sqrt(L);
        if(L > 0.0f){
            N.normalize();
            switch (K) {
                case 2:
                    _velocity(id * K) += _dt * N(1) * vor(K * id) * _x * coe;
                    _velocity(id * K + 1) -= _dt * N(0) * vor(K * id) * _x * coe;
                    break;
                case 3:
                    _velocity(id * K) += _dt * coe * (N(1) * vor(K * id + 2) - N(2) * vor(K * id + 1)) * _x;
                    _velocity(id * K + 1) += _dt * coe * (N(2) * vor(K * id) - N(0) * vor(K * id + 2)) * _x;
                    _velocity(id * K + 2) += _dt * coe * (N(0) * vor(K * id + 1) - N(1) * vor(K * id)) * _x;
            }
        }
    }
}

void Fluid::External(){
    for(int i = 0; i < num[K]; i++)
        _velocity(i * K + 1) += _dt * AddBuoyancy(i);

    AddVorticity();
}

void test(){
    
}

void Fluid::step(){
    VectorXf ov, oh, od;
    ov.resize(K * num[K]);
    oh.resize(num[K]);
    od.resize(num[K]);
    ov.setZero();
    oh.setZero();
    od.setZero();
    /*for(int id = 0; id < num[K]; id++)
        _velocity(id * K + 1) += _dt * AddBuoyancy(id);*/
    //printf("step\n");
    /*External();
    //printf("E\n");
    
    AdvectionAll(_velocity, ov, K);
    //printf("A\n");
    //od = _density;
    Diffusion(ov, _velocity, 0);
    //printf("D\n");
    //Pressure(_velocity, ov);
    //printf("P\n");
    //_velocity = ov;
    Diffusion(_heat, oh, 1);
    Diffusion(_density, od, 1);
    AdvectionAll(oh, _heat, 1);
    //oh = _heat;
    //ov = _velocity;
    AdvectionAll(od, _density, 1);*/
    
    
    External();
    
    AdvectionAll(_velocity, ov, K);
    _velocity = ov;
    Diffusion(_velocity, ov, 0);
    
    Diffusion(_heat, oh, 1);
    Diffusion(_density, od, 1);
    AdvectionAll(oh, _heat, 1);
    AdvectionAll(od, _density, 1);
    Pressure(ov, _velocity);
    //_velocity = ov;
    
    

    
    
    //cout << _velocity << endl;
    
    /*VectorXd c(2);
    c(0) = _peeled[1] / 2; c(1) = _peeled[2] / 2 + 1;
    int id = get_index(c);
    cout << _velocity(id * K) << " " << _velocity(id * K + 1) << endl;*/
    
    //_density /= 1.02;
}

