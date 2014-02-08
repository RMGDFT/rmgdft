#ifndef RMG_ERROR_H
#define RMG_ERROR_H 1

#if __cplusplus
extern "C" {
    void rmg_error_handler(const char *message);
}
#else
void rmg_error_handler(char *message);
#endif


#endif
