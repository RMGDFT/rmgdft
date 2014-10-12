#ifndef RMG_CheckValue_H
#define RMG_CheckValue_H

#include <string>

template <typename DataType>
bool CheckAndReturn(DataType val, DataType min, DataType max);
template <typename DataType>
void CheckAndTerminate(DataType val, DataType min, DataType max, std::string msg);
template <typename DataType>
void CheckAndFix(DataType *val, DataType min, DataType max, DataType fix, std::string msg);

#endif
