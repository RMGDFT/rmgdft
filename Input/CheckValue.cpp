#include <string>
#include "CheckValue.h"
#include "RmgException.h"

// These are used to check input values. CheckValueAndReturn should be used
// when the caller is able to correct bad input values while the second form
// CheckValueAndTerminate is used when correction is not possible.

template bool CheckValueAndReturn(int, int, int);
template bool CheckValueAndReturn(float, float, float);
template bool CheckValueAndReturn(double, double, double);

template void CheckValueAndTerminate(int, int, int, std::string);
template void CheckValueAndTerminate(float, float, float, std::string);
template void CheckValueAndTerminate(double, double, double, std::string);


template <typename DataType>
bool CheckValueAndReturn(DataType val, DataType min, DataType max)
{
    if((val < min) || (val > max)) {
        return false;
    }
    return true;
}


template <typename DataType>
void CheckValueAndTerminate(DataType val, DataType min, DataType max, std::string msg)
{
    if((val < min) || (val > max)) {
        throw RmgFatalException() << msg;
    }
}

