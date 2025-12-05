
#include "serialization.h"

// Writes a std::string to an ofstream, prefixing it with its length (uint32_t)
void writeString(std::ofstream &ofs, const std::string &str)
{
    u_int32_t len = static_cast<u_int32_t>(str.length());
    ofs.write((char *)&len, sizeof(len));
    if (len > 0)
    {
        ofs.write(str.data(), len);
    }
}

// Reads a std::string from an ifstream
std::string readString(std::ifstream &ifs)
{
    u_int32_t len = 0;
    ifs.read((char *)&len, sizeof(len));
    if (len > 0)
    {
        std::string str(len, '\0');
        ifs.read(str.data(), len);
        return str;
    }
    return "";
}