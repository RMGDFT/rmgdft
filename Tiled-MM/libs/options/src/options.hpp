#pragma once
#include <iostream>
#include <string>
#include <regex>
#include <stdexcept>

namespace options {
    static std::string command_line;
    void initialize(int argc, char** argv);

    // finds the position after the defined flag or throws an exception 
    // if flag is not found in the line.
    int find_flag(const std::string& short_flag, const std::string& long_flag, 
                     const std::string& message, const std::string& line, 
                     bool throw_exception=true);

    // looks for the defined flag in the line
    // if found return true, otherwise returns false
    bool flag_exists(const std::string& short_flag, const std::string& long_flag, 
            const std::string& line);

    // finds the next int after start in the line
    int next_int(int start, const std::string& line);
    long long next_long_long(int start, const std::string& line);
    double next_double(int start, const std::string& line);

    int next_int(const std::string& short_flag, const std::string& long_flag,
                 const std::string& message, int default_value);

    long long next_long_long(const std::string& short_flag, const std::string& long_flag,
                 const std::string& message, long long default_value);

    double next_double(const std::string& short_flag, const std::string& long_flag,
                 const std::string& message, double default_value);

    bool flag_exists(const std::string& short_flag, const std::string& long_flag);

    std::pair<int, int> next_int_pair(const std::string& short_flag,
                                      const std::string& long_flag,
                                      const std::string& message,
                                      int default_value);
};
