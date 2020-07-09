#include <iostream>
#include <Eigen/Dense>
#include <cstdio>
#include <vector>
#include <cmath>
#include <random>
#include "particle_filter.h"

#define PI 3.14159265359
#define MAX_RANGE 20.0
#define NP 100
#define NTh NP/2

using namespace std;
using namespace Eigen;

Particle_Filter::Particle_Filter() 
{
    // Time
    ts = 0.0f;
    dt = 0.1f;
    tf = 50.0f;
    // State Vector [x y yaw]'
    xDR << 0, 0, 0, 0;
    xEst << 0, 0, 0, 0;
    xTrue << 0, 0, 0, 0;
    // Covariance Matrix for motion
    Q = 0.5f;
    // Covariance Matrix for observation 
    R << pow(10.0, 2.0), 0,
         0, pow(toRadian(30), 2.0);
    // Simulation parameter
    Qsigma = 0.05f;
    Rsigma << pow(1.5, 2.0), 0,
              0, pow(toRadian(5), 2.0);
    PEst = Matrix4f::Identity();
    // Landmark
    RFID << 10, 0,
            10, 10,
            0, 15,
            -5, 20;
    // particle stor
    px = Matrix<float, 4, NP>::Zero();
    pw = Matrix<float, NP, 1>::Ones() * 1.0 / NP;
}

Particle_Filter::~Particle_Filter()
{
    cout << "Finish" << endl;
}

float Particle_Filter::toRadian(float degree)
{
    // degree to radian
    float radian = degree / 180.0 * PI;
    return radian;
}

Vector2f Particle_Filter::doControl(float time) 
{
    // Calc Input Parameter

    float T = 10.0f;
    float V = 1.0f; 
    float yawrate = 5.0f; 
    // Input
    u << V * (1 - exp(-time / T)), toRadian(yawrate) * (1 - exp(-time / T));
    
    return u;
}

Vector2f Particle_Filter::noise_distribution(Vector2f u) 
{
    // gaussian distribution
    random_device rd{};
    mt19937 gen{rd()};
    normal_distribution<> gaussian_d{0, 1};

    u(0, 0) = u(0, 0) + gaussian_d(gen) * Rsigma(0, 0);
    u(1, 0) = u(1, 0) + gaussian_d(gen) * Rsigma(1, 1);

    return u;
}

Vector4f Particle_Filter::motion_model(Vector4f x, Vector2f u) {
    F <<  1.0, 0, 0, 0,
          0, 1.0, 0, 0,
          0, 0, 1.0, 0,
          0, 0, 0, 1.0;
            
    B << dt* cos(x(2, 0)), 0,
         dt * sin(x(2, 0)), 0,
         0.0, dt,
         1.0, 0.0;

  return F * x + B * u;
}

float Particle_Filter::gauss_likelihood(float x, float sigma){
  float p = 1.0 / sqrt(2.0 * PI * sigma * sigma) * exp(-x * x / (2 * sigma * sigma));
  return p;
}

Matrix4f Particle_Filter::calc_covariance(Vector4f xEst, Matrix<float, 4, NP> px, Matrix<float, NP, 1> pw) {
    Matrix4f PEst_ = Matrix4f::Zero();

    for(int i = 0; i < px.cols(); i++){
        Vector4f dx = px.col(i) - xEst;
        PEst_ += pw(i) * dx * dx.transpose();
    }

    return PEst_;
}

void Particle_Filter::pf_localization(Matrix<float, 4, NP>& px, Matrix<float, NP, 1>& pw, Vector4f& xEst, Matrix4f& PEst, vector<RowVector3f> z, Vector2f u, Matrix2f Rsim, float Q) {

    for(int ip = 0; ip < NP; ip++){
        Vector4f x = px.col(ip);
        float w    = pw(ip);
        // noise distribution
        u = noise_distribution(u);
        x = motion_model(x, u);

        for(unsigned int i = 0; i < z.size(); i++){
            RowVector3f temp = z[i];
            float dx = x(0) - temp(1);
            float dy = x(1) - temp(2);
            float prez = sqrt(dx*dx + dy*dy);
            float dz = prez - temp(0);
            w = w * gauss_likelihood(dz, Q);
        }
        px.col(ip) = x;
        pw(ip)     = w;
  }

  pw   = pw / pw.sum();
  xEst = px * pw;
  PEst = calc_covariance(xEst, px, pw);
}

Matrix<float, NP, 1> Particle_Filter::cumsum(Matrix<float, NP, 1> pw) {
    Matrix<float, NP, 1> cum;
    cum(0) = pw(0);
    for(int i = 1; i < pw.rows(); i++){
        cum(i) = cum(i - 1) + pw(i);
    }

    return cum;
}

void Particle_Filter::resampling(Matrix<float, 4, NP>& px, Matrix<float, NP, 1>& pw){
    // noise
    random_device rd2{};
    mt19937 gen{rd2()};
    uniform_real_distribution<> uni_d{1.0, 2.0};

    float Neff = 1.0 / (pw.transpose() * pw);

    if (Neff < NTh){
        Matrix<float, NP, 1> wcum = cumsum(pw);
        Matrix<float, NP, 1> base = cumsum(pw * 0.0 +  Matrix<float, NP, 1>::Ones() * 1.0 / NP) - Matrix<float, NP, 1>::Ones() * 1.0 / NP;
        Matrix<float, NP, 1> resampleid;
        Matrix<float, 4, NP> output;

        for(int j = 0; j < pw.rows(); j++){
            resampleid(j) = base(j) + uni_d(gen) / NP;
        }

        int ind = 0;
        for(int i = 0; i < NP; i++){
            while (resampleid(i) > wcum(ind) && ind < NP - 1) {
            ind += 1;
            }
            output.col(i) = px.col(ind);
        }
        px = output;
        pw = Matrix<float, NP, 1>::Ones() * 1.0 / NP;
    }
}

void Particle_Filter::simulation()
{
    // save data
    FILE *fp;
    if ((fp = fopen("data.txt", "w")) == NULL) {
        printf("Error\n");
        exit(1);
    }
    fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", ts, xTrue(0, 0), xTrue(1, 0), xEst(0, 0), xEst(1, 0));
    // main loop
    for (ts; ts <= tf; ts += dt) {
        // Input
        u     = doControl(ts);
        xTrue = motion_model(xTrue, u);
        // noise distribution
        u     = noise_distribution(u);
        xDR   = motion_model(xDR, u);
        // gaussian distribution
        random_device rd{};
        mt19937 gen{rd()};
        normal_distribution<> gaussian_d{0, 1};

        z.clear();
        for(int i = 0; i < RFID.rows(); i++){
            float dx = xTrue(0) - RFID(i, 0);
            float dy = xTrue(1) - RFID(i, 1);
            float d  = sqrt(dx*dx + dy*dy);
            if (d <= MAX_RANGE){
                float dn = d + gaussian_d(gen) * Qsigma;
                RowVector3f zi;
                zi << dn, RFID(i, 0), RFID(i, 1);
                z.push_back(zi);
            }
        }
        // ------Particle Filter --------
        pf_localization(px, pw, xEst, PEst, z, u, Rsigma, Q);
        resampling(px, pw);
        
        fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", ts, xTrue(0, 0), xTrue(1, 0), xEst(0, 0), xEst(1, 0));
    }

    fclose(fp);
}
