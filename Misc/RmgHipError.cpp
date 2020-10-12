
#include "rmg_error.h"
#include "ErrorFuncs.h"


#if HIP_ENABLED

void RmgHipError(const char *file, int line, const hipError_t hipStatus, const char * errorMessage)
{
    if (hipStatus != hipSuccess) {
        rmg_error_handler(file, line, errorMessage);
    }
}


void RmgHipError(const char *file, int line, const hipblasStatus_t status, const char * errorMessage)
{
    if (status != HIPBLAS_STATUS_SUCCESS) {
        rmg_error_handler(file, line, errorMessage);
    }
}

#endif

