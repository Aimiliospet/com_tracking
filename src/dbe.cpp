#include <dbe.h>

DBE::DBE()
{
    // Gravity Vector
    g -> setZero();
    (*g)(2) = -9.80;
}   

void DBE::init()
{
    firstrun = true;


    If->setIdentity();
    P->setZero();
    //set velocity, rotational and positional errors and gyro, acc biases
    P->block<3, 3>(0, 0).setConstant(1e-3);
    P->block<3, 3>(3, 3).setConstant(1e-3);
    P->block<3, 3>(6, 6).setConstant(1e-5);
    P->block<3, 3>(9, 9).setConstant(1e-3);
    P->block<3, 3>(12, 12).setConstant(1e-3);

    Hf->setZero();
    Hf->block<3, 3>(0, 6).setIdentity();
    Hf->block<3, 3>(3, 3).setIdentity();
    Hvf->setZero();
    Hvf->block<3, 3>(3, 3).setIdentity();
    Hv->setZero();

    Rib->setIdentity();
    x->setZero();

    z->setZero();
    zv->setZero();

    v->setZero();
    dxf->setZero();
    temp->setZero();
    Kf->setZero();
    Kv->setZero();

    s->setZero();
    sv->setZero();

    R->setZero();
    Rv->setZero();

    Acf->setZero();
    //Qff->setZero();
    Qf->setZero();
    Af->setZero();

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

    std::cout << "Base EKF Initialized Successfully" << std::endl;
}


std::unique_ptr<Eigen::Matrix<double, 15, 15>> DBE::computeTrans(const std::unique_ptr<Eigen::Matrix<double, 15, 1>>& x_,
                                                                    const std::unique_ptr<Eigen::Matrix<double, 3, 3>>& Rib_,
                                                                    const std::unique_ptr<Quaterniond>& qib, 
                                                                    const std::unique_ptr<Eigen::Vector3d>& omega_,
                                                                    const std::unique_ptr<Eigen::Vector3d>& f_)
{

std::unique_ptr<Eigen::Matrix<double, 15, 15>> res = std::make_unique<Eigen::Matrix<double, 15, 15>>();
(*res).setZero();
(*res).block<3, 3>(0, 0)->setIdentity();
(*res).block<3, 3>(6, 6)->setIdentity();
(*res).block<3, 3>(9, 9)->setIdentity();
(*res).block<3, 3>(12, 12)->setIdentity();
std::unique_ptr<Eigen::Matrix<double, 3, 3>> eye = std::make_unique<Eigen::Matrix<double, 3, 3>>();
(*eye).block<3, 3>(0, 0)->setIdentity();
(*res).block<3, 3>(0, 6) = (*eye)*dt;
(*res).block<3, 3>(3, 3) = expmap(dt*(x->segment<3>(9)));
(*res).block<3, 3>(6, 3) = dt*(*wedge(Rib_->transpose()*(*f_ - x->segment<3>(12))));  //have to check if here we use Rib or R(pq) (R: fct from quaternion to rotation matrix)
(*res).block<3, 3>(6, 12) = -dt*Rib_->transpose();
return res;
}

void DBE::euler(const std::unique_ptr<Vector3d>& omega_,const std::unique_ptr<Vector3d>& f_)
{
    Acf = computeTrans(x, Rib, qib, omega_, f_);
    // Euler Discretization - First order Truncation
    (*Af) = (*If);
    Af->noalias() += (*Acf) * dt;
    x = computeDiscreteDyn(x, Rib, qib, omega_, f_);
    // x.noalias() += computeContinuousDyn(x,Rib,omega_,f_)*dt;
}

//maybe have to change Rib with qib !!!! and update it separately !!!!
std::unique_ptr<Matrix<double, 15, 1>> DBE::computeDiscreteDyn(const std::unique_ptr<Matrix<double, 15, 1>>& x_,
                                                          const std::unique_ptr<Quaterniond>& qib, 
														  const std::unique_ptr<Eigen::Matrix3d>& Rib_,
														  const std::unique_ptr<Vector3d>& omega_,
														  const std::unique_ptr<Vector3d>& f_)
{
    std::unique_ptr<Matrix<double, 15, 1>> res = std::make_unique<Eigen::Matrix<double, 15, 1>>();
    res->segment<3>(0).noalias() = x_->segment<3>(0) + x_->segment<3>(6)*dt;
    //here we have to find a solution with the quaternion!!! segment<3> not gonna fit the quaternion !!!!
     //here x_->segment<3>(9) must change to the RATE of the angular velocity! (in the paper it sais quaternion check this !!!)
    res->segment<3>(6).noalias() = x_->segment<3>(6) + ((quaternionToRotationMatrix(*qib)->transpose())*((*f_) - x_->segment<3>(12))-(*g))*dt;
    res->segment<3>(9).noalias() = x_->segment<3>(9) + ((*omega_) - (*omega_old)) ;
    res->segment<3>(12).noalias() = x_->segment<3>(12);
    return res; 
    
}                                                          

void DBE::predict(const std::unique_ptr<Vector3d>& omega_, const std::unique_ptr<Vector3d>& f_)
{
    (*omega_old) = (*omega);
    (*omega) = (*omega_);
    (*f) = (*f_);
    //omegahat->noalias() = (*omega) - x->segment<3>(9);
    //(*v) = x->segment<3>(0);

    // Update the Input-noise Jacobian
    //Lcf->block<3, 3>(0, 0).noalias() = -(*wedge(v));

    euler(omega_, f_);

    // Covariance Q with full state + biases
    (*Qf)(0, 0) = gyr_qx * gyr_qx;
    (*Qf)(1, 1) = gyr_qy * gyr_qy;
    (*Qf)(2, 2) = gyr_qz * gyr_qz;
    (*Qf)(3, 3) = acc_qx * acc_qx;
    (*Qf)(4, 4) = acc_qy * acc_qy;
    (*Qf)(5, 5) = acc_qz * acc_qz;
    (*Qf)(6, 6) = gyrb_qx * gyrb_qx;
    (*Qf)(7, 7) = gyrb_qy * gyrb_qy;
    (*Qf)(8, 8) = gyrb_qz * gyrb_qz;
    (*Qf)(9, 9) = accb_qx * accb_qx;
    (*Qf)(10, 10) = accb_qy * accb_qy;
    (*Qf)(11, 11) = accb_qz * accb_qz;

    //Qff->noalias() = (*Af) * (*Lcf) * (*Qf) * Lcf->transpose() * Af->transpose() * dt;
    /** Predict Step: Propagate the Error Covariance  **/
    P->noalias() = (*Af) * (*P) * Af->transpose() + (*Qf); //have to check how do we exactly compute Q !!!
    // update separately Rib !
    // Propagate only if non-zero input

   
    if (!omegahat->isZero())
    {   
        std::unique_ptr<Vector3d> temp;
        (*temp) = (*omegahat) * dt;
        //update Rib
        (*Rib) *= (*expMap(temp));
        //update qib
         *qib +=  0.5*vectorRate(omega, omega_)*(*qib)*dt;
    }

    x->segment<3>(3).setZero();
    updateVars();
}

void DBE::updateWithTwistRotation(const std::unique_ptr<Vector3d>& y,const std::unique_ptr<Eigen::Quaterniond>& qy)
{






}

void DBE::updateVars()
{

    (*pos) = x->segment<3>(6);
    rX = (*pos)(0);
    rY = (*pos)(1);
    rZ = (*pos)(2);
    Tib->linear() = (*Rib);
    Tib->translation() = (*pos);
    (*qib) = Quaterniond(Tib->linear());

    // Update the biases
    (*bgyr) = x->segment<3>(9);
    (*bacc) = x->segment<3>(12);
    bias_gx = (*x)(9);
    bias_gy = (*x)(10);
    bias_gz = (*x)(11);
    bias_ax = (*x)(12);
    bias_ay = (*x)(13);
    bias_az = (*x)(14);

    (*omegahat) = (*omega) - (*bgyr);
    (*fhat) = (*f) - (*bacc);

    (*gyro) = (*Rib) * (*omegahat);
    gyroX = (*gyro)(0);
    gyroY = (*gyro)(1);
    gyroZ = (*gyro)(2);

    (*acc) = (*Rib) * (*fhat);
    accX = (*acc)(0);
    accY = (*acc)(1);
    accZ = (*acc)(2);

    (*v) = x->segment<3>(0);
    (*vel) = (*Rib) * (*v);
    velX = (*vel)(0);
    velY = (*vel)(1);
    velZ = (*vel)(2);

    // ROLL - PITCH - YAW
    angle = getEulerAngles(Rib);
    angleX = (*angle)(0);
    angleY = (*angle)(1);
    angleZ = (*angle)(2);
}