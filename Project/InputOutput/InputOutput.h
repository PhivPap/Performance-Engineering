#ifndef INPUT_OUTPUT
#define INPUT_OUTPUT

#include <fstream>
#include <string>
#include <list>

namespace IO {
    class Body {
    public:
        std::string id;
        double mass, x, y, vel_x, vel_y;
        Body(void);
        Body(const std::string &id, double mass, double x, double y, double vel_x, double vel_y);
        void copy_from(const Body &b);
    };


    class Parser {  
    public:
        Parser(const std::string &path);
        ~Parser(void);
        int next_body_info(Body &io_body);

    private:
        bool failed = false;
        uint32_t parsed_bodies = 0;
        std::ifstream input_file;
    };


    class Writer {
    public:
        Writer(const std::string &path);
        ~Writer(void);
        void write_body(const Body& io_body);

    private:
        bool failed = false;
        std::ofstream output_file;
    };
}

#endif