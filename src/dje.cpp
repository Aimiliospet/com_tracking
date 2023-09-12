#include <DJE.h>

DJE::DJE()
{
    // Gravity Vector
    g -> setZero();
    (*g)(2) = -9.80;
}   

void DJE::init()
{
    firstrun = true;

    If->setIdentity();
    
    Qf->setZero();
    R->setZero();
    Hf->setIdentity();
    Af->setIdentity();
    x->setZero();
    z->setZero();
    dxf->setZero();
    temp->setZero();
    Kf->setZero();

    s->setZero();
    sv->setZero();
    
    P = DARE(Af, Hf, Qf, R, maxiterations, tolerance);


    /*
    bgyr->setZero();
    bacc->setZero();
    gyro->setZero();
    acc->setZero();
    angle->setZero();

    fhat->setZero();
    omegahat->setZero();
    v -> setZero();
    r-> setZero();
    angleX = 0.0;
    angleY = 0.0;
    angleZ = 0.0;
    gyroX = 0.0;
    gyroY = 0.0;
    gyroZ = 0.0;
    accX = 0.0;
    accY = 0.0;
    accZ = 0.0;
    rX = 0.0;
    rY = 0.0;
    rZ = 0.0;
    velX = 0.0;
    velY = 0.0;
    velZ = 0.0;
    Tib->setIdentity();
    */
    std::cout << "joint state estimator Initialized Successfully" << std::endl;
}

std::unique_ptr<MatrixXd> DJE::DARE(std::unique_ptr<MatrixXd>& F, std::unique_ptr<MatrixXd>& H, std::unique_ptr<MatrixXd>& Q, std::unique_ptr<MatrixXd>& R, int maxiterations, double tolerance)
{
    *P = *Q;
    for (int i = 0; i < maxiterations; ++i) {
        MatrixXd P_next = F->transpose() * (*P) * (*F) - F->transpose() * (*P) * (*H) *
                          (H->transpose() * (*P) * (*H) + (*R)).inverse() * H->transpose() * (*P) * (*F) + (*Q);

        if ((P_next - (*P)).norm() < tolerance) {
            (*P) = P_next;
            break;
        }

        (*P) = P_next;
    }
    (*Kf) = (*P) * H->transpose() * ((*R) + (*H) * (*P) * H->transpose()).inverse();
    return P;
}

std::unique_ptr<MatrixXd> DJE::computeDiscreteDyn(const std::unique_ptr<MatrixXd>& x_, const std::unique_ptr<MatrixXd>& jointvelocities)     
{
    std::unique_ptr<MatrixXd> res = x_;
    *res = *x_ + (*jointvelocities)*dt;
}         


void DJE::predict(const std::unique_ptr<MatrixXd>& jointvelocities)
{
   x = computeDiscreteDyn(x, jointvelocities);
    
}

void DJE::update(const std::unique_ptr<MatrixXd>& jointpositions)
{

dxf = (*Kf) * jointpositions;
x->noalias() += dxf;


}

