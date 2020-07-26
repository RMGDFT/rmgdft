#ifndef RMG_rmgfiles_h
#define RMG_rmgfiles_h


int FileOpenAndCreate(std::string &pathname, int flags, mode_t mode);
void *CreateMmapArray(int &fd, size_t length);
void DeleteNvmeArrays(void);

#endif
