#include <iostream>
#include <chrono>
#include "../InputOutput/InputOutput.h"

const std::string DEF_IN = "BodyFiles/in/in0.tsv";
const std::string DEF_OUT = "BodyFiles/out/out0.tsv";

void parse_input(const std::string &input_path /*, ??? */){
    try {
        auto parser = IO::Parser(input_path);
        auto io_body = IO::Body();
        while (parser.next_body_info(io_body) == 0){
            /* ---->: create body copying values from 'io_body'. 
            Do not use IO::Body for the simulation. */
        }
    }
    catch (const std::string &e) {
        std::cout << "Parse error: " << e << std::endl;
        exit(1);
    }
}

void write_output(const std::string &output_path /*, ??? */){
    try {
        auto writer = IO::Writer(output_path);
        auto io_body = IO::Body();
        for (/* const auto &body : bodies */;0;) {
            /* ---->: copy values from body to io_body */
            writer.write_body(io_body);
        }
    }
    catch (const std::string &e) {
        std::cout << "Write error: " << e << std::endl;
        exit(1);
    }
}

void simulate(/* ??? */){
    // ???
}

int main(int argc, const char **argv){   
    try {
        parse_input(DEF_IN /*, ??? */);
        const auto start = std::chrono::system_clock::now();
        simulate(/* ??? */);
        const auto end = std::chrono::system_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e9;
        std::cout << "Seconds: "  << elapsed << std::endl;
        write_output(DEF_OUT /*, ??? */);
    }
    catch (const std::string &e) {
        std::cout << e << std::endl;
    }
}
