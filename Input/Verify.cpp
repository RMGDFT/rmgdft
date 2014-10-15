#include <typeinfo>
#include <iostream>
#include <unordered_map>
#include "InputKey.h"

// Verify function for std::string
bool Verify(const std::string& KeyName, const std::string& KeyVal, const std::unordered_map<std::string, InputKey *>& Map)
{

    // Look for the key first
    std::unordered_map<std::string, InputKey *>::const_iterator found = Map.find (KeyName);

    // not found so return false
    if(found == Map.end()) return false;

    // found so need to check if types match
    InputKey *Key = found->second;
    if(Key->KeyType != typeid(std::string).hash_code()) return false;

    // Types match so check if vals match
    if(KeyVal.compare(*Key->Readstr)) return false;
    
    return true;
}

// Verify function for char *
bool Verify(const std::string& KeyName, const char *keyval, const std::unordered_map<std::string, InputKey *>& Map)
{
    std::string KeyVal(keyval);    
    return Verify(KeyName, KeyVal, Map);

}

// Verify function for booleans
bool Verify(const std::string& KeyName, const bool& KeyVal, const std::unordered_map<std::string, InputKey *>& Map)
{

    // Look for the key first
    std::unordered_map<std::string, InputKey *>::const_iterator found = Map.find (KeyName);

    // not found so return false
    if(found == Map.end()) return false;

    // found so need to check if types match
    InputKey *Key = found->second;
    if(Key->KeyType != typeid(bool).hash_code()) return false;

    // Types match so check if vals match
    if(KeyVal == *Key->Readboolval) return true;
    
    return false;
}

