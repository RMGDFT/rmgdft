#include <string>
#include <iostream>
#include "CheckValue.h"
#include "RmgException.h"

// These are used to check input values. CheckValueAndReturn should be used
// when the caller is able to correct bad input values while the second form
// CheckValueAndTerminate is used when correction is not possible.

template bool CheckAndReturn(int, int, int);
template bool CheckAndReturn(float, float, float);
template bool CheckAndReturn(double, double, double);

template void CheckAndTerminate(int, int, int, std::string);
template void CheckAndTerminate(float, float, float, std::string);
template void CheckAndTerminate(double, double, double, std::string);

template void CheckAndFix(int *, int, int, int, std::string);
template void CheckAndFix(float *, float, float, float, std::string);
template void CheckAndFix(double *, double, double, double, std::string);

template <typename DataType>
bool CheckAndReturn(DataType val, DataType min, DataType max)
{
    if((val < min) || (val > max)) {
        return false;
    }
    return true;
}


template <typename DataType>
void CheckAndTerminate(DataType val, DataType min, DataType max, std::string msg)
{
    if((val < min) || (val > max)) {
        throw RmgFatalException() << msg << " in " << __FILE__ << " at line " << __LINE__ << "\n";
    }
}


template <typename DataType>
void CheckAndFix(DataType *val, DataType min, DataType max, DataType fix, std::string msg)
{
    if((*val < min) || (*val > max)) {
       *val = fix;
       std::cout << msg;
    }
}

