#include <typeinfo>
#include <iostream>
#include <unordered_map>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <iomanip>
#include "InputKey.h"

// This is a rather kludgy way of doing this but I didn't have much luck getting
// the boost::program_options class working with run time type info so I went with this for now
// 
// If the passed value of ReadVal is NULL then we allocate memory to store the item inside
// the constructor. Right now this is not freed inside the constructor since the typical
// program usage is such that we don't care. (Only called once and items live for the
// duration of the program). But we may want to change this in the future in which case
// some means of tracking whether or not the memory needs to be freed in the des

// ints
InputKey::InputKey(std::string& NewKeyName, int *ReadVal, int Minval, int Maxval, int Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg) : KeyName(NewKeyName) {

    if(!ReadVal) {
        allocated = true;
        ReadVal = new int[1];
    }
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

// doubles
InputKey::InputKey(std::string& NewKeyName, double *ReadVal, double Minval, double Maxval, double Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg) : KeyName(NewKeyName) {

    if(!ReadVal) {
        allocated = true;
        ReadVal = new double[1];
    }
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

// booleans
InputKey::InputKey(std::string& NewKeyName, bool *ReadVal, bool Defval, const char *helpmsg) : KeyName(NewKeyName) {

    if(!ReadVal) {
        allocated = true;
        ReadVal = new bool[1];
        *ReadVal = Defval;
    }
    InputKey::Required = false;
    InputKey::Readboolval = ReadVal;
    InputKey::Defboolval = Defval;
    InputKey::helpmsg = helpmsg;
    InputKey::KeyType = typeid(bool).hash_code();

}

// Regular string
InputKey::InputKey(std::string& NewKeyName, std::string *ReadStr, const char *Defstr, bool Fix, bool Required, const char *helpmsg, const char *errmsg) : KeyName(NewKeyName) {

    if(ReadStr) {
        InputKey::Readstr = *ReadStr;
    }
    InputKey::Readstrorig = ReadStr;
    InputKey::Defstr = Defstr;
    InputKey::Fix = Fix;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::MapPresent = false;
    InputKey::KeyType = typeid(std::string).hash_code();

}

// Enumerated string
InputKey::InputKey(std::string& NewKeyName, std::string *ReadStr, int *ReadVal, const char *Defstr, bool Fix, bool Required, const std::unordered_map<std::string, int> Allowed, const char *helpmsg, const char *errmsg) : KeyName(NewKeyName) {

    if(ReadStr) {
        InputKey::Readstr = *ReadStr;
    }
    InputKey::Readstrorig = ReadStr;
    InputKey::Readintval = ReadVal;
    InputKey::Range = Allowed;
    InputKey::Defstr = Defstr;
    InputKey::Fix = Fix;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::MapPresent = true;
    InputKey::KeyType = typeid(std::string).hash_code();

}

// Integer vector
InputKey::InputKey(std::string& NewKeyName, RmgInput::ReadVector<int> *V , RmgInput::ReadVector<int> *Defintvec, size_t count, bool Required, const char* helpmsg, const char *errmsg) : KeyName(NewKeyName) {

    if(!V) {
        allocated = true;
        V = new RmgInput::ReadVector<int>[1];
    }
    InputKey::Vintorig = V;
    InputKey::Vint = *V;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::count = count;
    InputKey::Defintvec = *Defintvec;
    InputKey::KeyType = typeid(RmgInput::ReadVector<int>).hash_code();
}

// Double vector
InputKey::InputKey(std::string& NewKeyName, RmgInput::ReadVector<double> *V,  RmgInput::ReadVector<double> *Defdblvec, size_t count, bool Required, const char* helpmsg, const char *errmsg) : KeyName(NewKeyName) {

    if(!V) {
        allocated = true;
        V = new RmgInput::ReadVector<double>[1];
    }
    InputKey::Vdoubleorig = V;
    InputKey::Vdouble = *V;
    InputKey::Required = Required;
    InputKey::helpmsg = helpmsg;
    InputKey::errmsg = errmsg;
    InputKey::count = count;
    InputKey::Defdblvec = *Defdblvec;
    InputKey::KeyType = typeid(RmgInput::ReadVector<double>).hash_code();

}

std::string InputKey::Print(void) {

    if(this->KeyType == typeid(int).hash_code())
        return boost::lexical_cast<std::string>(*this->Readintval);

    if(this->KeyType == typeid(double).hash_code()) {
        std::ostringstream ss;
        if(fabs(*this->Readdoubleval) > 0.0001) {
            ss << std::fixed << std::setprecision(8);
            ss << *this->Readdoubleval;
        }
        else {
            ss << std::scientific << std::setprecision(8);
            ss << *this->Readdoubleval;
        }
        return ss.str();
        //return boost::lexical_cast<std::string>(*this->Readdoubleval);
    }

    if(this->KeyType == typeid(bool).hash_code()) {
        if(*this->Readboolval) return "true";
        return "false";
    }

    if(this->KeyType == typeid(std::string).hash_code()) {
        return this->Readstr;
    }

    if(this->KeyType == typeid(RmgInput::ReadVector<int>).hash_code()) {
        std::string rstr;
        for(size_t i = 0;i < this->count;i++)
            rstr = rstr + boost::lexical_cast<std::string>(this->Vint.vals.at(i)) + " ";

        return rstr;
    }

    if(this->KeyType == typeid(RmgInput::ReadVector<double>).hash_code()) {
        std::string rstr;
        for(size_t i = 0;i < this->count;i++)
            rstr = rstr + boost::lexical_cast<std::string>(this->Vdouble.vals.at(i)) + " ";

    }

    return "Not done yet";
}
