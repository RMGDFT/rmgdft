#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <cfloat>
#include <climits>
#include <unordered_map>
#include <typeinfo>
#include "InputKey.h"

namespace po = boost::program_options;

#define CHECK_AND_FIX 0
#define CHECK_AND_TERMINATE 1
#define REQUIRED true
#define OPTIONAL false


class RmgInputFile {
    public:
        RmgInputFile(char *inputfile);
        ~RmgInputFile(void);

        template <typename T>
        void RegisterInputKey(std::string KeyName, T *Readval, T Minval, T Maxval, T Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg);
        void RegisterInputKey(std::string KeyName, std::string *Readstr, const char *defstr, bool Fix, bool Required, const char *helpmsg, const char *errmsg);
        void RegisterInputKey(std::string KeyName, RmgInput::ReadVector<int> *V, int count, bool Required, const char *helpmsg, const char *errmsg);
        void RegisterInputKey(std::string KeyName, RmgInput::ReadVector<double> *V, int count, bool Required, const char *helpmsg, const char *errmsg);
        void RegisterInputKey(std::string KeyName, bool *ReadVal, bool Defval, const char *helpmsg);


        void LoadInputKeys(void);

    private:
        std::string PreprocessInputFile(char *cfile);
        std::string sfile;
        po::options_description control;
        std::unordered_map<std::string, InputKey *> InputMap;
        po::variables_map vm;


};


