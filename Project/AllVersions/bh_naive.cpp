#include <iostream>
#include <chrono>
#include <vector>
#include <math.h>
#include "../InputOutput/InputOutput.h"

const std::string DEF_IN = "BodyFiles/in/in0.tsv";
const std::string DEF_OUT = "BodyFiles/out/out0.tsv";
const double G = 6.67e-11;                              // Gravitational constant
const uint32_t total_time_steps = 50;                   // 50 hours
const double time_step_length = 3600;                   // 1 hour


void parse_input(const std::string& input_path /*, ... */){
    try {
        auto parser = IO::Parser(input_path);
        auto io_body = IO::Body();
        while (parser.next_body_info(io_body) == 0){
            // TODO: Parse values from io_body to correct the structure (dont use IO::Body for simulation)
        }
    }
    catch (const std::string& e) {
        std::cout << "Parse error: " << e << std::endl;
        exit(1);
    }
}

void write_output(const std::string& output_path /*, ... */){
    try {
        auto writer = IO::Writer(output_path);
        uint32_t id = 0;
        // TODO: for (body in bodies)
        //     writer.write_body(IO::Body(id++, body.mass, body.x, body.y, body.vel_x, body.vel_y));
    }
    catch (const std::string& e) {
        std::cout << "Write error: " << e << std::endl;
        exit(1);
    }
}

void simulate(/* ..., */ double time_step, uint32_t iterations){
    
}

int main(int argc, const char** argv){
    parse_input(DEF_IN /*, ... */);

    const auto start = std::chrono::system_clock::now();
    simulate(/* ..., */ time_step_length, total_time_steps);
    const auto end = std::chrono::system_clock::now();
    
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e9;
    std::cout << "Seconds: " << elapsed << std::endl;

    write_output(DEF_OUT /*, ... */);
}
