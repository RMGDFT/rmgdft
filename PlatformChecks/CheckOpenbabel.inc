# Default to failure then run compile checks
set(CHECK_OPENBABEL_INTERNAL 0)

set(CMAKE_REQUIRED_INCLUDES "${OPENBABEL_INCLUDES}")
set(CMAKE_REQUIRED_LIBRARIES "${OPENBABEL_LIBRARIES}")
CHECK_CXX_SOURCE_COMPILES("
#include <iostream>
#include <openbabel/obconversion.h>
using namespace std;
int main(int argc,char **argv)
{
  int m;
  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);
  OpenBabel::OBConversion conv(&ifs, &ofs);
  conv.SetInAndOutFormats(\"CML\",\"MOL\");
  int n = conv.Convert();
  return 0;
}
" CHECK_OPENBABEL_INTERNAL1)

if(CHECK_OPENBABEL_INTERNAL1)
    set(CHECK_OPENBABEL_INTERNAL 1)
else(NOT CHECK_OPENBABEL_INTERNAL1)

    # fix screwy openbabel path naming setup
    string(LENGTH ${OPENBABEL_INCLUDES} OPENBABEL_INCLUDES_STRLEN)
    math(EXPR NEW_OPENBABEL_INCLUDES_STRLEN "${OPENBABEL_INCLUDES_STRLEN} - 9")
    string(SUBSTRING ${OPENBABEL_INCLUDES} 0 ${NEW_OPENBABEL_INCLUDES_STRLEN} OPENBABEL_INCLUDES)
    set(CMAKE_REQUIRED_INCLUDES "${OPENBABEL_INCLUDES}")
    set(CMAKE_REQUIRED_LIBRARIES "${OPENBABEL_LIBRARIES}")
    CHECK_CXX_SOURCE_COMPILES("
    #include <iostream>
    #include <openbabel/obconversion.h>
    using namespace std;
    int main(int argc,char **argv)
    {
      ifstream ifs(argv[1]);
      ofstream ofs(argv[2]);
      OpenBabel::OBConversion conv(&ifs, &ofs);
      conv.SetInAndOutFormats(\"CML\",\"MOL\");
      int n = conv.Convert();
      return 0;
    }
    " CHECK_OPENBABEL_INTERNAL2)
    if(CHECK_OPENBABEL_INTERNAL2)
        set(CHECK_OPENBABEL_INTERNAL 1)
    endif(CHECK_OPENBABEL_INTERNAL2)

endif(CHECK_OPENBABEL_INTERNAL1)
