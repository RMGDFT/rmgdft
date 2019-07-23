/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <cfloat>
#include <climits>
#include <regex>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include "transition.h"

#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "CheckValue.h"
#include "RmgException.h"
#include "RmgInputFile.h"
#include "InputOpts.h"


void wordwrap(std::string const &input, size_t width, size_t indent)
{ 

    std::vector<std::string> lines;
    boost::split(lines, input, boost::is_any_of("&"), boost::token_compress_on);

//    for(auto line = lines.begin();line != lines.end(); ++line)
    size_t last = lines.size() - 1;
    for(size_t i=0;i < lines.size();i++)
    {
        std::istringstream in(lines[i]);
        size_t current = indent;
        std::string tstr;
        while (in >> tstr)
        {
            if (current + tstr.size() > width) {
                std::cout << "\n" << std::string(indent, ' ');
                current = indent;
            }
            std::cout << tstr << ' ';
            current += tstr.size() + 1;
        }
        std::cout << "\n";
        if(i < last) std::cout << std::string(indent, ' ');
    }
}

// Writes a key to stdout
void WriteKeyStdout(InputKey *ik)
{
    const std::string& whitespace("\t\n\v\f\r ");
    const std::string yesno[2] = {"no", "yes"};
    const std::string truefalse[2] = {"false", "true"};
    std::string KeyName = ik->KeyName;
    std::string KeyType;
    char *pre, *post;
    char *pre_terminal = "    ";
    char *post_terminal = "";
    char *pre_markdown = "<b>";
    char *post_markdown = "</b>\n";

    pre = pre_terminal;
    post = post_terminal;

    if(ik->KeyType == typeid(int).hash_code()) KeyType = "integer";
    if(ik->KeyType == typeid(double).hash_code()) KeyType = "double";
    if(ik->KeyType == typeid(bool).hash_code()) KeyType = "boolean";
    if(ik->KeyType == typeid(std::string).hash_code()) KeyType = "string";
    if(ik->KeyType == typeid(RmgInput::ReadVector<int>).hash_code()) KeyType = "integer array";
    if(ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code()) KeyType = "double array";

    
    printf("%sKey name:%s     %s\n", pre, post, KeyName.c_str());
    printf("%sRequired:%s     %s\n", pre, post, yesno[ik->Required].c_str());
    printf("%sKey type:%s     %s\n", pre, post, KeyType.c_str());
    if(ik->KeyType == typeid(int).hash_code())
    {
        printf("%sMin value:%s    %d\n", pre, post, ik->Minintval);
        printf("%sMax value:%s    %d\n", pre, post, ik->Maxintval);
        printf("%sDefault:%s      %d\n", pre, post, ik->Defintval);
    }
    if(ik->KeyType == typeid(double).hash_code())
    {
        if(ik->Mindoubleval == -DBL_MAX)
            printf("%sMin value:%s    %s\n", pre, post, "-unlimited");
        else
        {
            if(fabs(ik->Mindoubleval) > 0.01)
               printf("%sMin value:%s    %f\n", pre, post, ik->Mindoubleval);
            else
                printf("%sMin value:%s    %e\n", pre, post, ik->Mindoubleval);
        }

        if(ik->Maxdoubleval == DBL_MAX)
            printf("%sMax value:%s    %s\n", pre, post, "unlimited");
        else
        {
            if(fabs(ik->Maxdoubleval) > 0.01)
                printf("%sMax value:%s    %f\n", pre, post, ik->Maxdoubleval);
            else
                printf("%sMax value:%s    %e\n", pre, post, ik->Maxdoubleval);
        }
        if(fabs(ik->Defdoubleval) > 0.01)
            printf("%sDefault:%s      %f\n", pre, post, ik->Defdoubleval);
        else
            printf("%sDefault:%s      %e\n", pre, post, ik->Defdoubleval);
    }
    if(ik->KeyType == typeid(bool).hash_code())
    {
        printf("%sDefault:%s      \"%s\"\n", pre, post, truefalse[ik->Defboolval].c_str());
    }
    if(ik->KeyType == typeid(std::string).hash_code())
    {
        printf("%sDefault:%s      \"%s\"\n", pre, post, ik->Defstr);
        printf("%sAllowed:%s      ", pre, post);
        std::string tstr;
        for(auto it = ik->Range.begin();it != ik->Range.end(); ++it) tstr = tstr + "\"" + it->first + "\"   ";
        boost::trim_right(tstr);
        wordwrap(tstr, 84, 18);
        //printf("\n");
    }
    if(ik->KeyType == typeid(RmgInput::ReadVector<int>).hash_code()) 
    {
        std::string str = ik->Print();
        printf("%sDefault:%s      \"%s\"\n", pre, post, str.c_str());
    }
    if(ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code()) 
    {
        std::string str = ik->Print();
        printf("%sDefault:%s      \"%s\"\n", pre, post, str.c_str());
    }

    // Strip trailing whitespace
    std::string Description = std::string(ik->helpmsg);
    boost::trim_right(Description);
    printf("%sDescription:%s  ", pre, post);
    wordwrap(Description, 84, 18);
    printf("\n");
}

// Writes out input options for command line help and documentation
void WriteInputOptions(std::unordered_map<std::string, InputKey *>& InputMap)
{

    printf ("\n");
    printf ("               * * * * * * * * * *\n");
    printf ("               *    R   M   G    *\n");
    printf ("               * * * * * * * * * *\n");
    printf ("\n");
    printf (" -- A Real Space Multigrid Electronic structure code --\n");
    printf (" --      More information at www.rmgdft.org          --\n");
    printf("\n\n");
    printf(""
"    The RMG input file consists of a set of key-value pairs of the form.\n"
"    \n"
"        name = \"scalar\"\n\n"
"    \n"
"    where scalar can be an integer, double or boolean value.\n\n"
"        period_of_diagonalization = \"1\"\n"
"        charge_density_mixing = \"0.5\"\n"
"        initial_diagonalization = \"true\"\n"
"    \n"
"    There are also strings and arrays which are delineated by double quotes so an\n"
"    integer array with three elements would be.\n"
"    \n"
"        processor_grid = \"2 2 2\"\n"
"    \n"
"    while a string example would be\n"
"    \n"
"        description = \"64 atom diamond cell test run at the gamma point\"\n"
"    \n"
"    strings can span multiple lines so the following would be valid as well.\n"
"    \n"
"        description = \"64 atom diamond cell test run at gamma point\n"
"        using a Vanderbilt ultrasoft pseudopotential\"\n");

    printf("\n\n");

    std::map<std::string, InputKey *> SortedMap;
    for(auto it = InputMap.begin();it != InputMap.end(); ++it)
    {
        std::string KeyName = it->first;
        InputKey *ik = it->second;
        SortedMap.insert(std::make_pair(KeyName, ik));
    }

    for(auto it = SortedMap.begin();it != SortedMap.end(); ++it)
    {
        InputKey *ik = it->second;
        WriteKeyStdout(ik);
    }
}

