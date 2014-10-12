#ifndef RMG_InputKey_H
#define RMG_InputKey_H 1

#include <string>


class InputKey {
    public:
        // Scalar types
        InputKey(std::string KeyName, int *ReadVal, int Minval, int Maxval, int Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg);
        InputKey(std::string KeyName, float *ReadVal, float Minval, float Maxval, float Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg);
        InputKey(std::string KeyName, double *ReadVal, double Minval, double Maxval, double Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg);
        InputKey(std::string KeyName, bool *ReadVal, bool Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg);

        // Strings
        InputKey(std::string KeyName, std::string *ReadStr, const char *defstr, bool Fix, bool Required, const char *helpmsg, const char *errmsg);

        // Vectors
//        template <typename T>
//        InputKey(std::string KeyName, Ri::ReadVector<T>, bool Required, const char *helpmsg, const char *errmsg);
        

        std::string KeyName;
        size_t KeyType;
        bool Fix;
        bool Required;
        int *Readintval;
        int Minintval;
        int Maxintval;
        int Defintval;
        float *Readfloatval;
        float Minfloatval;
        float Maxfloatval;
        float Deffloatval;
        double *Readdoubleval;
        double Mindoubleval;
        double Maxdoubleval;
        double Defdoubleval;
        bool *Readboolval;
        bool Defboolval;
        std::string *Readstr;
        const char *Defstr;
        const char *helpmsg;
        const char *errmsg;
};

#endif
