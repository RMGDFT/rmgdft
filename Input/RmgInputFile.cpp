#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
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
#include "portability.h"
#include "InputKey.h"
#include "RmgInputFile.h"
#include "RmgException.h"
#include "CheckValue.h"

#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "transition.h"

namespace po = boost::program_options;


template void RmgInputFile::RegisterInputKey<int>(std::string, int *, int, int, int, bool, bool, const char *, const char *);
template void RmgInputFile::RegisterInputKey<double>(std::string, double *, double, double, double, bool, bool, const char *, const char *);


// Custom validators for boost program options. 
namespace RmgInput {

    // handles integer vectors
    void validate(boost::any& v, const std::vector<std::string>& values, ReadVector<int>*, int)
    {
        ReadVector<int> *A = new ReadVector<int>;
        po::validators::check_first_occurrence(v);
        const std::string& s = po::validators::get_single_string(values);
        std::string t1 = s;
        boost::trim_if(t1, boost::algorithm::is_any_of("\" \t"));

        while(t1.size() > 0) {
            A->vals.push_back( std::atoi(t1.c_str()));
            size_t f1 = t1.find_first_of(" \t");
            if(f1 == std::string::npos) break;
            t1 = t1.substr(f1, std::string::npos);
            boost::trim(t1);
        }
        v = *A;
    }

    // handles double vectors
    void validate(boost::any& v, const std::vector<std::string>& values, ReadVector<double>*, double)
    {
        ReadVector<double> *A = new ReadVector<double>;
        po::validators::check_first_occurrence(v);
        const std::string& s = po::validators::get_single_string(values);
        std::string t1 = s;
        boost::trim_if(t1, boost::algorithm::is_any_of("\" \t"));

        while(t1.size() > 0) {
            A->vals.push_back( std::atof(t1.c_str()));
            size_t f1 = t1.find_first_of(" \t");
            if(f1 == std::string::npos) break;
            t1 = t1.substr(f1, std::string::npos);
            boost::trim(t1);
        }
        v = *A;
    }

}


RmgInputFile::RmgInputFile(char *inputfile, std::unordered_map<std::string, InputKey *>& Map, MPI_Comm comm) : InputMap(Map)  {
    PreprocessInputFile(inputfile, comm);
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

void RmgInputFile::RegisterInputKey(std::string KeyName, std::string *Readstr, int *ReadVal, const char *defstr, bool Fix, bool Required, const std::unordered_map<std::string, int>& Allowed, const char *helpmsg, const char *errmsg) {
    InputKey *NewKey = new InputKey(KeyName, Readstr, ReadVal, defstr, Fix, Required, Allowed, helpmsg, errmsg);
    std::pair <std::string, InputKey *> NewEntry(KeyName, NewKey);
    InputMap.insert(NewEntry);
}

void RmgInputFile::RegisterInputKey(std::string KeyName, RmgInput::ReadVector<int> *V, RmgInput::ReadVector<int> *Defintvec, int count, bool Required, const char *helpmsg, const char *errmsg) {
    InputKey *NewKey = new InputKey(KeyName, V, Defintvec, count, Required, helpmsg, errmsg);
    if(!Required) {
//        for(int i = 0;i < count;i++) NewKey->Vint.vals = Defintvec->vals;
NewKey->Vint.vals = Defintvec->vals;
    }
    std::pair <std::string, InputKey *> NewEntry(KeyName, NewKey);
    InputMap.insert(NewEntry);
}

void RmgInputFile::RegisterInputKey(std::string KeyName, RmgInput::ReadVector<double> *V, RmgInput::ReadVector<double> *Defdblvec, int count, bool Required, const char *helpmsg, const char *errmsg) {
    InputKey *NewKey = new InputKey(KeyName, V, Defdblvec, count, Required, helpmsg, errmsg);
    if(!Required) {
//        for(int i = 0;i < count;i++) NewKey->Vdouble.vals = Defdblvec->vals;
NewKey->Vdouble.vals = Defdblvec->vals;
    }
    std::pair <std::string, InputKey *> NewEntry(KeyName, NewKey);
    InputMap.insert(NewEntry);
}

void RmgInputFile::LoadInputKeys(void) {

    std::string sfile;
    InputKey *Key;
    for (auto item = InputPairs.begin();item != InputPairs.end();item++) {

        std::string PairName = item->first;
        std::string PairValue = item->second;

        try {
            Key = InputMap.at(PairName);
            std::string val(PairValue);
            boost::trim_if(val, boost::algorithm::is_any_of("\" "));
            if((Key->KeyType == typeid(int).hash_code()) || (Key->KeyType == typeid(double).hash_code()) || (Key->KeyType == typeid(bool).hash_code())) {
                sfile = sfile + PairName + "=" + val + "\n";
            }
            else {
                sfile = sfile + PairName + "=" + "\"" + val + "\"\n";
            }
        }
        catch (const std::out_of_range& oor) {
//            std::cout << "Warning!! Unknown input tag: " << PairName << std::endl;
        }

    }
    //std::cout << sfile << std::endl;exit(0);

    // Load the keys into the map
    for (auto item = InputMap.begin();item != InputMap.end();item++) {

        std::string KeyName = item->first;
        InputKey *Ik = item->second;
        if(Ik->Required) {

            if(Ik->KeyType == typeid(int).hash_code())
                control.add_options() (KeyName.c_str(), po::value<int>(Ik->Readintval)->required(), Ik->helpmsg);

            if(Ik->KeyType == typeid(double).hash_code())
                control.add_options() (KeyName.c_str(), po::value<double>(Ik->Readdoubleval)->required(), Ik->helpmsg);

            if(Ik->KeyType == typeid(std::string).hash_code()) {
                boost::trim_if(Ik->Readstr, boost::algorithm::is_any_of("\"^"));
                control.add_options() (KeyName.c_str(), po::value<std::string>(&Ik->Readstr)->required(), Ik->helpmsg);
            }

            if(Ik->KeyType == typeid(RmgInput::ReadVector<int>).hash_code())
                control.add_options() (KeyName.c_str(), po::value(&Ik->Vint)->required(), Ik->helpmsg);

            if(Ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code())
                control.add_options() (KeyName.c_str(), po::value(&Ik->Vdouble)->required(), Ik->helpmsg);

        }
        else {

            if(Ik->KeyType == typeid(int).hash_code())
                control.add_options() (KeyName.c_str(), po::value<int>(Ik->Readintval)->default_value(Ik->Defintval), Ik->helpmsg);

            if(Ik->KeyType == typeid(double).hash_code())
                control.add_options() (KeyName.c_str(), po::value<double>(Ik->Readdoubleval)->default_value(Ik->Defdoubleval), Ik->helpmsg);

            if(Ik->KeyType == typeid(bool).hash_code())
                control.add_options() (KeyName.c_str(), po::value<bool>(Ik->Readboolval)->default_value(Ik->Defboolval), Ik->helpmsg);

            if(Ik->KeyType == typeid(std::string).hash_code()) {
                boost::trim_if(Ik->Readstr, boost::algorithm::is_any_of("\"^"));
                control.add_options() (KeyName.c_str(), po::value<std::string>(&Ik->Readstr)->default_value(Ik->Defstr), Ik->helpmsg);
            }

            if(Ik->KeyType == typeid(RmgInput::ReadVector<int>).hash_code())
                control.add_options() (KeyName.c_str(), po::value(&Ik->Vint), Ik->helpmsg);

            if(Ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code())
                control.add_options() (KeyName.c_str(), po::value(&Ik->Vdouble), Ik->helpmsg);

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
            if(Ik->Vint.vals.size() != Ik->count) {
                throw RmgFatalException() << Ik->KeyName << " " << Ik->errmsg
<< " in " << __FILE__ << " at line " << __LINE__ << "\n";
            }
            *Ik->Vintorig = Ik->Vint;
        }

        if(Ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code()) {
            if(Ik->Vdouble.vals.size() != Ik->count) {
                throw RmgFatalException() << Ik->KeyName << " " << Ik->errmsg
<< " in " << __FILE__ << " at line " << __LINE__ << "\n";
            }
            *Ik->Vdoubleorig = Ik->Vdouble;
        }

        if(Ik->KeyType == typeid(std::string).hash_code()) {

            // Trim quotes
            if(Ik->MapPresent) {
                // enumerated string so check if the value is allowed
                boost::trim_if(Ik->Readstr, boost::algorithm::is_any_of("\"^"));
                if(Ik->Range.count(Ik->Readstr) == 0) {
                    throw RmgFatalException() << Ik->KeyName << " " << Ik->errmsg;
                }
                // If integer val associated with the enum is desired set it
                if(Ik->Readintval) *Ik->Readintval = Ik->Range[Ik->Readstr];
                if(Ik->Readstrorig) *Ik->Readstrorig = Ik->Readstr;
                
            }
            else {
                // regular string
                boost::trim_if(Ik->Readstr, boost::algorithm::is_any_of("\"^"));
                if(Ik->Readstrorig) *Ik->Readstrorig = Ik->Readstr;
                //std::cout << *Ik->Readstr << std::endl;

            }

        }


    }

}


void RmgInputFile::PreprocessInputFile(char *cfile, MPI_Comm comm)
{
    std::string config_file(cfile);
    std::string outbuf;
    std::string tbuf;
    std::string Msg;
    std::stringstream ifs;
    char *input_buffer;
    int input_buffer_len;

    int rank = pct.imgpe;
    if(comm == MPI_COMM_WORLD) rank = pct.worldrank;

    // Open on one pe and read entire file into a character buffer
    if(rank == 0) {

        // Check for file existence
        boost::filesystem::path input_filepath(cfile);
        if( !boost::filesystem::exists(input_filepath) ) {

            Msg = "Input file " + boost::lexical_cast<std::string>(cfile) + " does not exist.\n";

        }
        else {

            try {
                std::ifstream input_file;
                input_file.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
                input_file.open(cfile);
                input_file.seekg(0, std::ios::end);
                input_buffer_len = input_file.tellg();
                input_file.seekg(0, std::ios::beg);
                input_buffer = new char[input_buffer_len + 1]();
                input_file.read(input_buffer, input_buffer_len);       // read the whole file into the buffer
                input_file.close();
            }
            // Catch any file io errors and rethrow later with our own error message
            catch (std::exception &e) {
                Msg = "Unable to read from input file " + boost::lexical_cast<std::string>(cfile) + "\n";
            }
        }

    }
    int openfail = Msg.length();
    MPI_Bcast(&openfail, 1, MPI_INT, 0, comm);
    if(openfail)
        throw RmgFatalException() << Msg << " in " << __FILE__ << " at line " << __LINE__ << "\n";

    // Send it to everyone else
    MPI_Bcast (&input_buffer_len, 1, MPI_INT, 0, comm);
    if(rank != 0) {
        input_buffer = new char[input_buffer_len + 1]();
    }
    MPI_Bcast (input_buffer, input_buffer_len, MPI_CHAR, 0, comm);
    std::string input_string(input_buffer);
    ifs << input_string;

    // First pass to remove leading and trailing comments
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

    // First tokenizer pass to get rid of empty tokens
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep("", "\"=");
    tokenizer tokens(tbuf, sep);
    bool breakflag = false;
    for (tokenizer::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter) {
        std::string line(*tok_iter);
        boost::trim(line);

        while(!line.size()) {
            tok_iter++;
            if(tok_iter == tokens.end()) {
                breakflag = true;
                break;
            }
            line.erase();
            line = *tok_iter;
            boost::trim(line);
        }
        if(breakflag) break;

        // Should be an input tag so throw if it's an = or a "
        if(std::string::npos != line.find("=")) throw RmgFatalException() << "Malformed input file near " << line << "\n";
        if(std::string::npos != line.find("\"")) throw RmgFatalException() << "Malformed input file near " << line << "\n";
        outbuf = outbuf + line; 

        // Advance to next non empty token
        line.erase();
        while(!line.size()) {
            tok_iter++;
            if(tok_iter == tokens.end()) {
                throw RmgFatalException() << "Malformed input file near " << line << "\n";
            }
            line.erase();
            line = *tok_iter;
            boost::trim(line);
        }

        // Should be an = so throw if not
        if(std::string::npos == line.find("=")) throw RmgFatalException() << "Malformed input file near " << line << "\n";
        outbuf = outbuf + line; 

        // Advance to next non empty token
        line.erase();
        while(!line.size()) {
            tok_iter++;
            if(tok_iter == tokens.end()) {
                throw RmgFatalException() << "Malformed input file near " << line << "\n";
            }
            line.erase();
            line = *tok_iter;
            boost::trim(line);
        }

        // Should be a " so throw if not

        outbuf = outbuf + line;
        if(std::string::npos == line.find("\"")) 
            throw RmgFatalException() << "Malformed input file near " << outbuf << "\n";
        tok_iter++;
        if(tok_iter == tokens.end()) throw RmgFatalException() << "Malformed input file near " << line << "\n";
        line.erase();
        line = *tok_iter;
        boost::trim(line);
        while(std::string::npos == line.find("\"")) {
             outbuf = outbuf + line;
             if(std::string::npos != line.find("=")) 
                throw RmgFatalException() << "Malformed input file near " << outbuf << "\n";
             tok_iter++;
             if(tok_iter == tokens.end()) break;
             line.erase();
             line = *tok_iter;
             boost::trim(line);
        }
        outbuf = outbuf + "\"\n";

    }

    // Second tokenizer pass to make the InputPairs map
    boost::char_separator<char> pairsep("\"=");
    tokenizer pairtokens(outbuf, pairsep);
    for (tokenizer::iterator tok_iter = pairtokens.begin();tok_iter != tokens.end(); ++tok_iter) {
        std::string line1(*tok_iter);
        tok_iter++;
        if(tok_iter == tokens.end()) break;
        std::string line2(*tok_iter);
        boost::trim_if(line1, boost::algorithm::is_any_of("\" \t\n\r"));
        boost::trim_if(line2, boost::algorithm::is_any_of("\" \t\n\r"));
        std::replace( line2.begin(), line2.end(), '\n', '^');
        std::pair <std::string, std::string> NewEntry(line1, line2);
        InputPairs.insert(NewEntry);
    }

//    if(pct.imgpe == 0) {

        // Write out options file
//        std::string OptionsFile(ct.basename);
//        OptionsFile = OptionsFile + ".options";

//        FILE *fhand = fopen(OptionsFile.c_str(), "w");
//        if (!fhand)
//            throw RmgFatalException() <<  "Unable to write file in " << __FILE__ << " at line " << __LINE__ << "\n";
//        fprintf(fhand, "%s", outbuf.c_str());
//        fclose(fhand);
//    }
    //std::cout << outbuf << std::endl;exit(0);

}

// Destructor frees up the map
RmgInputFile::~RmgInputFile(void) {

//    for (auto item = InputMap.begin();item != InputMap.end();item++) {
//        InputKey *Ik = item->second;
//        delete Ik;        
//    }

}
