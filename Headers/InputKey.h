#ifndef RMG_InputKey_H
#define RMG_InputKey_H 1

#include <string>
#include <vector>
#include <unordered_map>

namespace RmgInput {

    template <typename VectorType>
    class ReadVector
    {
        public:
            std::vector<VectorType> vals;
    };

}


class InputKey {
    public:
        // Scalar types
        InputKey(std::string& KeyName, int *ReadVal, int Minval, int Maxval, int Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg);
        InputKey(std::string& KeyName, double *ReadVal, double Minval, double Maxval, double Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg);
        InputKey(std::string& KeyName, bool *ReadVal, bool Defval, const char *helpmsg);

        // Strings
        InputKey(std::string& KeyName, std::string *ReadStr, const char *defstr, bool Fix, bool Required, const char *helpmsg, const char *errmsg);

        // Enumerated strings
        InputKey(std::string& KeyName, std::string *ReadStr, const char *Defstr, bool Fix, bool Required, const std::unordered_map<std::string, int> Allowed, const char *helpmsg, const char *errmsg);


        // Vectors
        InputKey(std::string& KeyName, RmgInput::ReadVector<int> *V , size_t count, bool Required, const char* helpmsg, const char *errmsg);
        InputKey(std::string& KeyName, RmgInput::ReadVector<double> *V, size_t count, bool Required, const char* helpmsg, const char *errmsg);
        

        std::string& KeyName;
        std::unordered_map<std::string, int> Range;
        bool MapPresent;
        size_t KeyType;
        bool Fix;
        bool Required;
        int *Readintval;
        int Minintval;
        int Maxintval;
        int Defintval;
        double *Readdoubleval;
        double Mindoubleval;
        double Maxdoubleval;
        double Defdoubleval;
        bool *Readboolval;
        bool Defboolval;
        std::string* Readstr;
        RmgInput::ReadVector<int> *Vint;
        RmgInput::ReadVector<double> *Vdouble;
        size_t count;
        const char *Defstr;
        const char* helpmsg;
        const char *errmsg;
};

#endif
