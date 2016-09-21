#include "transition.h"
#include "common_prototypes.h"

// C-binding used for ON which needs to reinitialize the q functions
// in the fast relax routine which is still in C and several layers deep
extern std::unordered_map<std::string, InputKey *> ControlMap;

extern "C" void init_qfunct(void)
{
    InitQfunct(ControlMap);
}

