#ifndef  __DIFFERENTIATOR_H__
#define  __DIFFERENTIATOR_H__


#include <iostream>
#include <string>
using namespace std;

class Differentiator
{

private:
    double x_, dt;
    bool firstrun;
    string name;
public:
    double x;
    double xdot;

    /** @fn void setParams(double dt_)
     *  @brief differentiates the measurement with finite differences
     *  @param dt_ Sampling time in seconds e.g. 0.01s
    */  
    void setParams(double dt_)
    {
        dt=dt_;
    }
    /** @fn void diff(double x)
     *  @brief differentiates the measurement with finite differences
     *  @param x signal to be differentiatied
    */
    double diff(double x);

    /** @fn void init(string name_,double dt_);
     *  @brief initializes the numerical differentiator
     *  @param name_ name of the signal e.g LHipYawPitch
     *  @param dt_ Sampling time in seconds e.g. 0.01s
    */  
    void init(string name_,double dt_);

    /** @fn void reset();
     *  @brief  resets the the numerical differentiator's state
    */    
    void reset();
};
#endif
