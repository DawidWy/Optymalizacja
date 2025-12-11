#include "csv.h"

CSVStream::CSVStream(const std::string &filename, const std::vector<std::string> &headers, char sep) : separator(sep)
{
    file.open(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }
    for (const auto &header : headers)
    {
        file << header;
        file << separator;
    }
    file << std::endl; 
    headers_num = headers.size();
}
CSVStream::~CSVStream()
{
    if (file.is_open())
    {
        file.close();
    }
}
CSVStream &CSVStream::operator<<(std::ostream &(*manip)(std::ostream &))
{
    file << manip;
    return *this;
}
CSVStream &CSVStream::operator<<(CSVStream &(*manip)(CSVStream &))
{
    return manip(*this);
}
CSVStream &CSVStream::sep()
{
    file << separator;
    return *this;
}

void CSVStream::newline()
{
    file << std::endl;
}