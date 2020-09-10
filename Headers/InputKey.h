#ifndef RMG_InputKey_H
#define RMG_InputKey_H 1

#include <string>
#include <vector>
#include <unordered_map>

#define CONTROL_OPTIONS         0
#define CELL_OPTIONS            1
#define PSEUDO_OPTIONS          2
#define KS_SOLVER_OPTIONS       3
#define XC_OPTIONS              4
#define OCCUPATION_OPTIONS      5
#define MIXING_OPTIONS          6
#define MD_OPTIONS              7
#define DIAG_OPTIONS            8
#define PERF_OPTIONS            9
#define LDAU_OPTIONS           10
#define TDDFT_OPTIONS          11
#define POISSON_OPTIONS        12
#define OUTPUT_OPTIONS         13
#define TESTING_OPTIONS        14
#define MISC_OPTIONS           15

// These can be or'd with the above to attach a modifier
#define EXPERT_OPTION 16777216
#define EXPERIMENTAL_OPTION 33554432

// This is used to mask off the modifiers
#define MASKOFF_MODIFIERS 16777215



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
        InputKey(std::string& NewKeyName, int *ReadVal, int Minval, int Maxval, int Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg, size_t group);
        InputKey(std::string& NewKeyName, double *ReadVal, double Minval, double Maxval, double Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg, size_t group);
        InputKey(std::string& NewKeyName, bool *ReadVal, bool Defval, const char *helpmsg, size_t group);

        // Strings
        InputKey(std::string& NewKeyName, std::string *ReadStr, const char *defstr, bool Fix, bool Required, const char *helpmsg, const char *errmsg, size_t group);

        // Enumerated strings
        InputKey(std::string& NewKeyName, std::string *ReadStr, int *ReadVal, const char *Defstr, bool Fix, bool Required, const std::unordered_map<std::string, int> Allowed, const char *helpmsg, const char *errmsg, size_t group);


        // Vectors
        InputKey(std::string& NewKeyName, RmgInput::ReadVector<int> *V , RmgInput::ReadVector<int> *Defintvec, size_t count, bool Required, const char* helpmsg, const char *errmsg, size_t group);
        InputKey(std::string& NewKeyName, RmgInput::ReadVector<double> *V, RmgInput::ReadVector<double> *Defdblvec, size_t count, bool Required, const char* helpmsg, const char *errmsg, size_t group);
        
        // Reads stored value
        std::string Print(void);

        std::string KeyName;
        std::unordered_map<std::string, int> Range;
        bool MapPresent;  // For enumerated strings
        bool allocated=false;
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
        std::string Readstr;
        std::string *Readstrorig;
        RmgInput::ReadVector<int> *Vintorig;
        RmgInput::ReadVector<double> *Vdoubleorig;
        RmgInput::ReadVector<int> Vint;
        RmgInput::ReadVector<double> Vdouble;
        RmgInput::ReadVector<int> Defintvec;
        RmgInput::ReadVector<double> Defdblvec;
        size_t count;
        const char *Defstr;
        const char* helpmsg;
        const char *errmsg;
        size_t grouping;
};

// Used to sort keys by group and then alphabetically
struct keycompare
{
    bool operator()(InputKey *ik_lhs, InputKey *ik_rhs)
    {
        // First sort by group
        size_t left = ik_lhs->grouping & MASKOFF_MODIFIERS;
        size_t right = ik_rhs->grouping & MASKOFF_MODIFIERS;

        if(left < right) return true;
        if(left > right) return false;

        // Now alphabetically by key name within group
        if(ik_lhs->KeyName.compare(ik_rhs->KeyName) < 0) return true;
        return false;
    }
};

#endif
