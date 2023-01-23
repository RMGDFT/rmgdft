#include "options.hpp"

void options::initialize(int argc, char** argv) {
    for (auto i = 1; i < argc; ++i) {
        command_line += std::string(argv[i]) + " ";
    }
}

// finds the position after the defined flag or throws an exception 
// if flag is not found in the line.
int options::find_flag(const std::string& short_flag, const std::string& long_flag, 
        const std::string& message, const std::string& line, bool compulsory) {
    auto position = line.find(short_flag);
    auto length = short_flag.length();
    if (position == std::string::npos
            // check if it is not some other long flag
            // with same starting letter
            || (position > 0 && 
                line[position] == '-' &&
                line[position-1] == '-')) {
        position = line.find(long_flag);
        length = long_flag.length();
        if (position == std::string::npos) {
            if (compulsory)
                throw std::runtime_error(message);
            else
                return -1;
        }
    }
    while (line[position + length] == ' ') length++;
    return position + length;
}

// looks for the defined flag in the line
// if found return true, otherwise returns false
bool options::flag_exists(const std::string& short_flag, const std::string& long_flag, 
        const std::string& line) {
    auto position = line.find(long_flag);
    if (position == std::string::npos) {
        position = line.find(short_flag);
        if (position == std::string::npos) {
            return false;
        } else {
            auto after_position = position + short_flag.length();
            if (after_position < line.length()) {
                if (line[after_position] != ' ') {
                    return false;
                }
            }
        }
    } else {
        auto after_position = position + long_flag.length();
        if (after_position < line.length()) {
            if (line[after_position] != ' ') {
                return false;
            }
        }
    }
    return true;
}

// finds the next int after start in the line
int options::next_int(int start, const std::string& line) {
    if (start < 0) 
        return -1;
    std::regex int_expr("([0-9]+)");
    auto it = std::sregex_iterator(line.begin() + start, line.end(), int_expr);
    int result = std::stoi(it->str());
    return result;
}

// finds the next long long after start in the line
long long options::next_long_long(int start, const std::string& line) {
    if (start < 0) 
        return -1;
    std::regex int_expr("([0-9]+)");
    auto it = std::sregex_iterator(line.begin() + start, line.end(), int_expr);
    long long result = std::stoll(it->str());
    return result;
}

// finds the next double after start in the line
double options::next_double(int start, const std::string& line) {
    if (start < 0)
        return -1;
    std::regex double_expr("/^[0-9]+(\\.[0-9]+)?$");
    auto it = std::sregex_iterator(line.begin() + start, line.end(), double_expr);
    long long result = std::stod(it->str());
    return result;
}

int options::next_int(const std::string& short_flag, const std::string& long_flag,
         const std::string& message, int default_value) {
    try {
        auto it = options::find_flag(short_flag, long_flag, message, command_line);
        return options::next_int(it, command_line);
    } catch (const std::runtime_error& e) {
        return default_value;
    }
}

long long options::next_long_long(const std::string& short_flag, const std::string& long_flag,
         const std::string& message, long long default_value) {
    try {
        auto it = options::find_flag(short_flag, long_flag, message, command_line);
        return options::next_long_long(it, command_line);
    } catch (const std::runtime_error& e) {
        return default_value;
    }
}

double options::next_double(const std::string& short_flag, const std::string& long_flag,
         const std::string& message, double default_value) {
    try {
        auto it = options::find_flag(short_flag, long_flag, message, command_line);
        return options::next_double(it, command_line);
    } catch (const std::runtime_error& e) {
        return default_value;
    }
}

bool options::flag_exists(const std::string& short_flag, const std::string& long_flag) {
    try {
        auto it = options::find_flag(short_flag, long_flag, "", command_line);
        return options::flag_exists(short_flag, long_flag, command_line);
    } catch (const std::runtime_error& e) {
        return false;
    }
}

std::pair<int, int> options::next_int_pair(
        const std::string& short_flag,
        const std::string& long_flag,
        const std::string& message,
        int default_value) {
    try {
        auto it = options::find_flag(short_flag, long_flag, message, command_line);
        std::size_t delimiter_found = command_line.find('x', it);
        if (delimiter_found != std::string::npos) {
            int first = next_int(it, command_line);
            int second = next_int(delimiter_found+1, command_line);
            return {first, second};
        }
        return {default_value, default_value};
    } catch (const std::runtime_error& e) {
        return {default_value, default_value};
    }
}

