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


// Writes a key to stdout
void WriteKeyStdout(InputKey *ik)
{
    const std::string& whitespace("\t\n\v\f\r ");
    const std::string yesno[2] = {"no", "yes"};
    const std::string truefalse[2] = {"false", "true"};
    std::string KeyName = ik->KeyName;
    std::string KeyType;

    if(ik->KeyType == typeid(int).hash_code()) KeyType = "integer";
    if(ik->KeyType == typeid(double).hash_code()) KeyType = "double";
    if(ik->KeyType == typeid(bool).hash_code()) KeyType = "boolean";
    if(ik->KeyType == typeid(std::string).hash_code()) KeyType = "string";
    if(ik->KeyType == typeid(RmgInput::ReadVector<int>).hash_code()) KeyType = "integer array";
    if(ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code()) KeyType = "double array";


    // Strip trailing whitespace and expand embedded new lines.
    std::string Description = std::string(ik->helpmsg);
    boost::trim_right(Description);
    boost::replace_all(Description, "\n", "\n              ");
    
    printf("Key name:     %s\n", KeyName.c_str());
    printf("Required:     %s\n", yesno[ik->Required].c_str());
    printf("Key type:     %s\n", KeyType.c_str());
    printf("Description:  %s\n", Description.c_str());
    if(ik->KeyType == typeid(int).hash_code())
    {
        printf("Min value:    %d\n", ik->Minintval);
        printf("Max value:    %d\n", ik->Maxintval);
        printf("Default:      %d\n", ik->Defintval);
    }
    if(ik->KeyType == typeid(double).hash_code())
    {
        if(ik->Mindoubleval == -DBL_MAX)
            printf("Min value:    %s\n", "-unlimited");
        else
        {
            if(fabs(ik->Mindoubleval) > 0.01)
               printf("Min value:    %f\n", ik->Mindoubleval);
            else
                printf("Min value:    %e\n", ik->Mindoubleval);
        }

        if(ik->Maxdoubleval == DBL_MAX)
            printf("Max value:    %s\n", "unlimited");
        else
        {
            if(fabs(ik->Maxdoubleval) > 0.01)
                printf("Max value:    %f\n", ik->Maxdoubleval);
            else
                printf("Max value:    %e\n", ik->Maxdoubleval);
        }
        if(fabs(ik->Defdoubleval) > 0.01)
            printf("Default:      %f\n", ik->Defdoubleval);
        else
            printf("Default:      %e\n", ik->Defdoubleval);
    }
    if(ik->KeyType == typeid(bool).hash_code())
    {
        printf("Default:      \"%s\"\n", truefalse[ik->Defboolval].c_str());
    }
    if(ik->KeyType == typeid(std::string).hash_code())
    {
        printf("Default:      \"%s\"\n", ik->Defstr);
        printf("Allowed:      ");
        int counter=0;
        for(auto it = ik->Range.begin();it != ik->Range.end(); ++it)
        {
            printf("\"%s\"  ",it->first.c_str()); 
            counter++;
            if(!(counter % 4)) printf("\n              ");
        }
        printf("\n");
    }
    if(ik->KeyType == typeid(RmgInput::ReadVector<int>).hash_code()) 
    {
        std::string str = ik->Print();
        printf("Default:      \"%s\"\n", str.c_str());
    }
    if(ik->KeyType == typeid(RmgInput::ReadVector<double>).hash_code()) 
    {
        std::string str = ik->Print();
        printf("Default:      \"%s\"\n", str.c_str());
    }
    printf("\n");
}

// Writes out input options for command line help and documentation
void WriteInputOptions(std::unordered_map<std::string, InputKey *>& InputMap)
{

    printf("\n\n");
    printf("The RMG input file consists of a set of key-value pairs of the form.\n");
    printf("name = \"scalar\"    where scalar can be an integer, double or boolean value.\n");
    printf("period_of_diagonalization = \"1\"\n");
    printf("charge_density_mixing = \"0.5\"\n");
    printf("initial_diagonalization = \"true\"\n");
    printf("\n\n");
    printf("There are also strings and arrays which are delineated by double quotes so an\n");
    printf("integer array with three elements would be.\n");
    printf("    processor_grid = \"2 2 2\"\n\n");
    printf("while a string example would be\n\n");
    printf("description = \"64 atom diamond cell test run at the gamma point\"\n\n");
    printf("strings can span multiple lines so the following would be valid as well.\n\n");
    printf("description = \"64 atom diamond cell test run at gamma point\n");
    printf("using a Vanderbilt ultrasoft pseudopotential\"\n");
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

