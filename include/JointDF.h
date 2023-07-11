#pragma once
#include <differentiator.h>
#include <butterworthLPF.h>
#include <iostream>

using namespace std;
class JointDF
{

private:
    /// 2nd order butterworth filter to smooth the angular velocity
    butterworthLPF bw;
    /// linear differentiator filter to compute the angular velocity
    Differentiator df;

public:
    /// Joint angular position
    double JointPosition;
    /// Joint angular velocity
    double JointVelocity;
    /// Joint name
    string JointName;

    /** @fn double filter(double JointPosMeasurement);
     *  @brief estimates the Joint Velocity using the Joint Position measurement
     *  by the encoders
     */
    double filter(double JointPosMeasurement);
    void reset();
    /** @fn void init(string JointName_,double fsampling, double fcutoff);
     *  @brief initializes the differentiator filter
     *  @param JointName_ the name of the filter e.g. "LHipPitch"
     *  @param fsampling the sampling frequency of the sensor e.g. 100hz
     *  @param fcutoff the cut-off frequency of the  2nd order Low Pass
     *  Butterworth Filter filter e.g. 10hz
     */
    void init(string JointName_, double fsampling, double fcutoff);
};
