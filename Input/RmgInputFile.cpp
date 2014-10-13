#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
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
#include "RmgInputFile.h"
#include "RmgException.h"
#include "CheckValue.h"

namespace po = boost::program_options;


template void RmgInputFile::RegisterInputKey<int>(std::string, int *, int, int, int, bool, bool, const char *, const char *);
template void RmgInputFile::RegisterInputKey(std::string, double *, double, double, double, bool, bool, const char *, const char *);

// Custom validators for boost program options. 
namespace RmgInput {

    // handles integer vectors
    void validate(boost::any& v, const std::vector<std::string>& values, ReadVector<int>*, int)
    {
        ReadVector<int> A;
        po::validators::check_first_occurrence(v);
        const std::string& s = po::validators::get_single_string(values);
        std::string t1 = s;
        boost::trim_if(t1, boost::algorithm::is_any_of("\" \t"));

        while(t1.size() > 0) {
            A.vals.push_back( std::atoi(t1.c_str()));
            size_t f1 = t1.find_first_of(" \t");
            if(f1 == std::string::npos) break;
            t1 = t1.substr(f1, std::string::npos);
            boost::trim(t1);
        }
        v = A;
    }

    // handles double vectors
    void validate(boost::any& v, const std::vector<std::string>& values, ReadVector<double>*, double)
    {
        ReadVector<double> A;
        po::validators::check_first_occurrence(v);
        const std::string& s = po::validators::get_single_string(values);
        std::string t1 = s;
        boost::trim_if(t1, boost::algorithm::is_any_of("\" \t"));

        while(t1.size() > 0) {
            A.vals.push_back( std::atof(t1.c_str()));
            size_t f1 = t1.find_first_of(" \t");
            if(f1 == std::string::npos) break;
            t1 = t1.substr(f1, std::string::npos);
            boost::trim(t1);
        }
        v = A;
    }

}


RmgInputFile::RmgInputFile(char *inputfile) {
    sfile = PreprocessInputFile(inputfile);

}


template <typename T>
void RmgInputFile::RegisterInputKey(std::string KeyName, T *Readval, T Minval, T Maxval, T Defval, bool Fix, bool Required, const char *helpmsg, const char *errmsg) {
    InputKey *NewKey = new InputKey(KeyName, Readval, Minval, Maxval, Defval, Fix, Required, helpmsg, errmsg);
    std::pair <std::string, InputKey *> NewEntry(KeyName, NewKey);
    InputMap.insert(NewEntry);
}

void RmgInputFile::RegisterInputKey(std::string KeyName, bool *ReadVal, bool Defval, const char *helpmsg) {
    InputKey *NewKey = new InputKey(KeyName, ReadVal, Defval, helpmsg);
    std::pair <std::string, InputKey *> NewEntry(KeyName, NewKey);
    InputMap.insert(NewEntry);
}

void RmgInputFile::RegisterInputKey(std::string KeyName, std::string *Readstr, const char *defstr, bool Fix, bool Required, const char *helpmsg, const char *errmsg) {
    InputKey *NewKey = new InputKey(KeyName, Readstr, defstr, Fix, Required, helpmsg, errmsg);
    std::pair <std::string, InputKey *> NewEntry(KeyName, NewKey);
    InputMap.insert(NewEntry);
}

void RmgInputFile::RegisterInputKey(std::string KeyName, RmgInput::ReadVector<int> *V, int count, bool Required, const char *helpmsg, const char *errmsg) {
    InputKey *NewKey = new InputKey(KeyName, V, count, Required, helpmsg, errmsg);
    std::pair <std::string, InputKey *> NewEntry(KeyName, NewKey);
    InputMap.insert(NewEntry);
}

void RmgInputFile::RegisterInputKey(std::string KeyName, RmgInput::ReadVector<double> *V, int count, bool Required, const char *helpmsg, const char *errmsg) {
    InputKey *NewKey = new InputKey(KeyName, V, count, Required, helpmsg, errmsg);
    std::pair <std::string, InputKey *> NewEntry(KeyName, NewKey);
    InputMap.insert(NewEntry);
}

void RmgInputFile::LoadInputKeys(void) {

    // Load the keys into the map
    for (auto item = InputMap.begin();item != InputMap.end();item++) {

        std::string KeyName = item->first;
        InputKey *Ik = item->second;
        if(Ik->Required) {
            if(Ik->KeyType == typeid(int).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Readintval)->required(), Ik->helpmsg);
            if(Ik->KeyType == typeid(double).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Readdoubleval)->required(), Ik->helpmsg);
            if(Ik->KeyType == typeid(std::string).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Readstr)->required(), Ik->helpmsg);
            if(Ik->KeyType == typeid(RmgInput::ReadVector<int>).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Vint)->required(), Ik->helpmsg);
            if(Ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Vdouble)->required(), Ik->helpmsg);
        }
        else {
            if(Ik->KeyType == typeid(int).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Readintval)->default_value(Ik->Defintval), Ik->helpmsg);
            if(Ik->KeyType == typeid(double).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Readdoubleval)->default_value(Ik->Defdoubleval), Ik->helpmsg);
            if(Ik->KeyType == typeid(bool).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Readboolval)->default_value(Ik->Defboolval), Ik->helpmsg);
            if(Ik->KeyType == typeid(std::string).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Readstr)->default_value(Ik->Defstr), Ik->helpmsg);
            if(Ik->KeyType == typeid(RmgInput::ReadVector<int>).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Vint), Ik->helpmsg);
            if(Ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code())
                control.add_options() (KeyName.c_str(), po::value(Ik->Vdouble), Ik->helpmsg);
        }
        

    }

    std::stringstream ss;
    ss << sfile;

    // Parse the input file
    auto parsedOptions = po::parse_config_file(ss, control, true);
    po::store(parsedOptions, vm);
    po::notify(vm);
    auto unregistered = po::collect_unrecognized(parsedOptions.options, po::include_positional);


    // Process the results
    for (auto item = InputMap.begin();item != InputMap.end();item++) {

        std::string KeyName = item->first;
        InputKey *Ik = item->second;

        if((Ik->Fix == CHECK_AND_TERMINATE) && (Ik->KeyType==typeid(int).hash_code())) CheckAndTerminate<int>(*Ik->Readintval, Ik->Minintval, Ik->Maxintval, Ik->errmsg);
        if((Ik->Fix == CHECK_AND_TERMINATE) && (Ik->KeyType==typeid(double).hash_code())) CheckAndTerminate<double>(*Ik->Readdoubleval, Ik->Mindoubleval, Ik->Maxdoubleval, Ik->errmsg);

        if((Ik->Fix == CHECK_AND_FIX) && (Ik->KeyType==typeid(int).hash_code())) CheckAndFix<int>(Ik->Readintval, Ik->Minintval, Ik->Maxintval, Ik->Defintval, Ik->errmsg);
        if((Ik->Fix == CHECK_AND_FIX) && (Ik->KeyType==typeid(double).hash_code())) CheckAndFix<double>(Ik->Readdoubleval, Ik->Mindoubleval, Ik->Maxdoubleval, Ik->Defdoubleval, Ik->errmsg);
        
        if(Ik->KeyType == typeid(RmgInput::ReadVector<int>).hash_code()) {
            if(Ik->Vint->vals.size() != Ik->count) {
                throw RmgFatalException() << Ik->KeyName << Ik->errmsg;
            }
        }

        if(Ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code()) {
            if(Ik->Vint->vals.size() != Ik->count) {
                throw RmgFatalException() << Ik->KeyName << Ik->errmsg;
            }
        }


    }

}


std::string RmgInputFile::PreprocessInputFile(char *cfile)
{
    std::string config_file(cfile);
    std::string outbuf;
    std::string tbuf;

    std::ifstream ifs(cfile);

    // First pass to clean it up
    for(std::string line; std::getline(ifs, line); ) {

        // Strip leading and trailing whitespace
        boost::trim(line);

        // Filter lines and only include non-empty and non-comment lines
        std::size_t len = line.size();
        if(len > 0) {
            if( line.compare(0, 1, "#") ) {
                // Chop trailing comments. Might need to fix this for # characters embedded in strings
                std::size_t found;
                while(std::string::npos != (found = line.find_last_of("#"))) line.resize(found);
                boost::trim(line);
                tbuf = tbuf + line + "\n";
            }
        }

    }
    // Split into individual lines again
    std::string delims = "\n";
    std::vector<std::string> lines;
    std::vector<std::string> lines1;
    std::vector<std::string>::iterator it;
    boost::algorithm::split( lines, tbuf, boost::is_any_of(delims), boost::token_compress_on );
    tbuf.clear();
    int idx = -1;
    bool join = false;
    for (it = lines.begin(); it != lines.end(); ++it) {
        std::string line = *it;

        if(join) {
            lines1.at(idx) = lines1.at(idx) + line;
            join = false;
        }

        // If line contains "=" then it's a new key-value pair
        if(std::string::npos != line.find("=")) {
            idx++;
            // new key-value pair
            lines1.push_back(line);
            join = false;
        }

        // Forward line
        std::vector<std::string>::iterator it1;
        it1 = it + 1;

        // Any quotes in this line?
        std::size_t f1 = line.find("\"", 1);
        if(f1 != std::string::npos) {
            std::size_t f2 = line.find("\"", f1 + 1);
            // If two quotes then we are done, if one quote then join with next line unless it contains an =
            if(f2 != std::string::npos) {
                join = false;
            }
            else {
                join = true;
                if(it1 != lines.end()) {
                    std::string fline = *it1;
                    f1 = fline.find("=", 1);
                    if(f1 != std::string::npos) join = false;
                }
            }
        }
        else {
            // no quotes so join unless next line contains an equal sign
            join = true;
            if(it1 != lines.end()) {
                std::string fline = *it1;
                f1 = fline.find("=", 1);
                if(f1 != std::string::npos) join = false;
            }
        }
    }

    for (it = lines1.begin(); it != lines1.end(); ++it) {
        // if *it does not contain quotes but does contain true or false then it's a boolean field so
        // we turn the actual value into a 0 or 1 for the next level parsing routines
        std::string& fline = *it;
        std::size_t f1 = fline.find("\"", 1);
        if(f1 == std::string::npos) {
            // no quotes
            std::size_t f2 = fline.find("true");
            std::size_t f3 = fline.find("false");
            if((f2 != std::string::npos) && (f3 != std::string::npos)) {
                // Both true and false so some sort of error            
                throw RmgFatalException() << "Syntax error in " << cfile << " near " << fline;
            }
            if(f2 != std::string::npos) {
               fline.replace(f2, 4, "1");
            }
            if(f3 != std::string::npos) {
               fline.replace(f3, 5, "0");
            }
        }

        outbuf = outbuf + *it + "\n";
    }

//    std::cout << outbuf << std::endl;
    return outbuf;

}
