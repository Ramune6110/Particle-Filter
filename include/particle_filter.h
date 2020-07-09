#ifndef PARTICLE_FILTER
#define PARTICLE_FILTER

#include <vector>
#include <Eigen/Dense>

#define NP 100
#define NTh NP/2

using namespace std;
using namespace Eigen;

class Particle_Filter 
{
    public:
        float ts;
        float dt;
        float tf;
    private:
        // system matrix
        Matrix<float, 4, 4> F;
        Matrix<float, 4, 2> B;
        Matrix<float, 2, 4> H;
        // control input
        Vector2f u;
        // observation z
        vector<RowVector3f> z;
        // RFID remarks
        Matrix<float, 4, 2> RFID;
        // dead reckoning
        Vector4f xDR;
        // ground truth reading
        Vector4f xTrue;
        // Estimation
        Vector4f xEst;
        Matrix<float, 4, 4> PEst;
        // Motional model covariance
        float Q;
        // Observation model covariance
        Matrix<float, 2, 2> R;
        // Motion model simulation error
        float Qsigma;
        // Observation model simulation error
        Matrix<float, 2, 2> Rsigma;
        // particle stor
        Matrix<float, 4, NP> px;
        Matrix<float, NP, 1> pw;
    private:
        float toRadian(float degree);
        Vector2f doControl(float time);
        Vector2f observation_model(Vector4f x);
        Vector2f noise_distribution(Vector2f u);
        float gauss_likelihood(float x, float sigma);
        Vector4f motion_model(Vector4f x, Vector2f u);
        Matrix<float, NP, 1> cumsum(Matrix<float, NP, 1> pw);
        void resampling(Matrix<float, 4, NP>& px, Matrix<float, NP, 1>& pw);
        Matrix4f calc_covariance(Vector4f xEst, Matrix<float, 4, NP> px, Matrix<float, NP, 1> pw);
        void pf_localization(Matrix<float, 4, NP>& px, Matrix<float, NP, 1>& pw, Vector4f& xEst, Matrix4f& PEst, vector<RowVector3f> z, Vector2f u, Matrix2f Rsigma, float Q);
    public:
        Particle_Filter();
        ~Particle_Filter();
        void simulation();
};

#endif // PARTICLE_FILTER

