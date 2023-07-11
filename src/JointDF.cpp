#include <JointDF.h>

void JointDF::init(string JointName_, double fsampling, double fcutoff)
{
    JointName = JointName_;
    JointPosition = 0.000;
    JointVelocity = 0.000;

    bw.init(JointName, fsampling, fcutoff);
    df.init(JointName, 1.00 / fsampling);
    std::cout << JointName << " velocity filter initialized successfully"
        << std::endl;
}

void JointDF::reset()
{
    JointPosition = 0.000;
    JointVelocity = 0.000;
    std::cout << JointName << " velocity filter reset" << std::endl;
}

double JointDF::filter(double JointPosMeasurement)
{
    JointPosition = JointPosMeasurement;
    JointVelocity = bw.filter(df.diff(JointPosMeasurement));

    return JointVelocity;
}
