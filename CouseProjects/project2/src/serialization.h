#pragma once
#include <fstream>

constexpr u_int32_t FILE_VERSION = 0x53595356; // "SYSV" magic number

void writeString(std::ofstream &ofs, const std::string &str);
std::string readString(std::ifstream &ifs);
