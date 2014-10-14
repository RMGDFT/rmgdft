#include <typeinfo>
#include <iostream>
#include <unordered_map>
#include "InputKey.h"

// This is a rather kludgy way of doing this but I didn't have much luck getting
// the boost::program_options class working with run time type info so I went with this for now
InputKey::InputKey(std::string& KeyName, int *ReadVal, int Minval, int Maxval, int Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg) : KeyName(KeyName) {
    InputKey::Readintval = ReadVal;
    InputKey::Minintval = Minval;
    InputKey::Maxintval = Maxval;
    InputKey::Defintval = Defval;
    InputKey::Fix = Fix;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::KeyType = typeid(int).hash_code();
}

InputKey::InputKey(std::string& KeyName, double *ReadVal, double Minval, double Maxval, double Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg) : KeyName(KeyName) {
    InputKey::Readdoubleval = ReadVal;
    InputKey::Mindoubleval = Minval;
    InputKey::Maxdoubleval = Maxval;
    InputKey::Defdoubleval = Defval;
    InputKey::Fix = Fix;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::KeyType = typeid(double).hash_code();
}

InputKey::InputKey(std::string& KeyName, bool *ReadVal, bool Defval, const char *helpmsg) : KeyName(KeyName) {
    InputKey::Readboolval = ReadVal;
    InputKey::Defboolval = ReadVal;
    InputKey::helpmsg = helpmsg;
    InputKey::KeyType = typeid(bool).hash_code();
}

// Regular string
InputKey::InputKey(std::string& KeyName, std::string *ReadStr, const char *Defstr, bool Fix, bool Required, const char *helpmsg, const char *errmsg) : KeyName(KeyName) {
    InputKey::Readstr = ReadStr;
    InputKey::Defstr = Defstr;
    InputKey::Fix = Fix;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::MapPresent = false;
    InputKey::KeyType = typeid(std::string).hash_code();
}

// Enumerated string
InputKey::InputKey(std::string& KeyName, std::string *ReadStr, const char *Defstr, bool Fix, bool Required, const std::unordered_map<std::string, int> Allowed, const char *helpmsg, const char *errmsg) : KeyName(KeyName) {
    InputKey::Readstr = ReadStr;
    InputKey::Range = Allowed;
    InputKey::Defstr = Defstr;
    InputKey::Fix = Fix;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::MapPresent = true;
    InputKey::KeyType = typeid(std::string).hash_code();
}

InputKey::InputKey(std::string& KeyName, RmgInput::ReadVector<int> *V , size_t count, bool Required, const char* helpmsg, const char *errmsg) : KeyName(KeyName) {

    InputKey::Vint = V;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::count = count;
    InputKey::KeyType = typeid(RmgInput::ReadVector<int>).hash_code();
}

InputKey::InputKey(std::string& KeyName, RmgInput::ReadVector<double> *V, size_t count, bool Required, const char* helpmsg, const char *errmsg) : KeyName(KeyName) {

    InputKey::Vdouble = V;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::count = count;
    InputKey::KeyType = typeid(RmgInput::ReadVector<double>).hash_code();
}

