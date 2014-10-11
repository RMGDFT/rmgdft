#ifndef RMG_CheckValue_H
#define RMG_CheckValue_H

#include <string>

template <typename DataType>
bool CheckValueAndReturn(DataType val, DataType min, DataType max);
template <typename DataType>
void CheckValueAndTerminate(DataType val, DataType min, DataType max, std::string msg);

#endif
