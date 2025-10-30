#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

class CSVStream
{

public:
    std::ofstream file;
    char separator;

    CSVStream(const std::string &filename, char sep = ',', const std::vector<std::string> &headers = {});
    ~CSVStream();

    template <typename T>
    CSVStream &operator<<(const T &value)
    {
        file << value;
        file << separator;
        return *this;
    }

    CSVStream &operator<<(std::ostream &(*manip)(std::ostream &));

    CSVStream &operator<<(CSVStream &(*manip)(CSVStream &));

    // Custom manipulator for separator
    CSVStream &sep();

    void newline();
};