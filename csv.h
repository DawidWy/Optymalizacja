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
    int headers_num;

    int current_column = 0;


    CSVStream(const std::string &filename, char sep = ',', const std::vector<std::string> &headers = {});
    ~CSVStream();

    template <typename T>
    CSVStream &operator<<(const T &value)
    {
        current_column++;
        file << value;
        if (current_column < headers_num){
        file << separator;
        } else {
            current_column = 0;
            file << std::endl;
        }
        return *this;
    }

    CSVStream &operator<<(std::ostream &(*manip)(std::ostream &));

    CSVStream &operator<<(CSVStream &(*manip)(CSVStream &));

    // Custom manipulator for separator
    CSVStream &sep();

    void newline();
};