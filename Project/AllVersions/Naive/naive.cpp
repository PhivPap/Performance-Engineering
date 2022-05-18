#include <iostream>
#include <chrono>
#include <vector>
#include <math.h>
#include "../../InputOutput/InputOutput.h"


// default values can be altered with the main args
struct CFG {
    std::string input_file = "BodyFiles/in/def_in.tsv";
    std::string output_file = "BodyFiles/out/naive_out_050.tsv";
    uint32_t iterations = 50;
    double iter_len = 3600;

    void print(){
        std::cout << "Configuration:" << std::endl;
        std::cout   << "Input: " << input_file << "\nOutput: " << output_file 
                    << "\nIterations: " << iterations << "\nIteration legth: "
                    << iter_len << "s" << std::endl << std::endl;
    }
};

const double G = 6.67e-11;        // Gravitational constant
CFG config;

struct Body {
    double mass, x, y, vel_x, vel_y;
};



void parse_input(const std::string& input_path, std::vector<Body>& bodies){
    try {
        IO::Parser parser(input_path);
        auto io_body = IO::Body();
        while (parser.next_body_info(io_body) == 0)
            bodies.push_back(Body{ io_body.mass, io_body.x, io_body.y, io_body.vel_x, io_body.vel_y });
    }
    catch (const std::string& e) {
        std::cout << "Parse error: " << e << std::endl;
        exit(1);
    }
}

void write_output(const std::string& output_path, std::vector<Body>& bodies){
    try {
        IO::Writer writer(output_path);
        uint32_t id = 0;
        for (const auto& body : bodies)
            writer.write_body(IO::Body(id++, body.mass, body.x, body.y, body.vel_x, body.vel_y));
    }
    catch (const std::string& e) {
        std::cout << "Write error: " << e << std::endl;
        exit(1);
    }
}

void update_body_positions(Body* bodies, uint32_t body_count, double time_step){
    for (uint32_t i = 0; i < body_count; i++){
        Body& body = bodies[i];
        auto delta_x = body.vel_x * time_step;
        auto delta_y = body.vel_y * time_step;
        body.x += delta_x;
        body.y += delta_y;
    }
}

void update_body_velocities(Body* bodies, uint32_t body_count, double time_step){
    for (uint32_t i = 0; i < body_count; i++){
        Body& body = bodies[i];
        double Fx = 0.0, Fy = 0.0, ax, ay;
        for (uint32_t j = 0; j < body_count; j++){
            if (i == j) continue;
            Body& other_body = bodies[j];
            double distance = sqrt((body.x - other_body.x) * (body.x - other_body.x) +
                (body.y - other_body.y) * (body.y - other_body.y));
            double F = (G * body.mass * other_body.mass) / (distance * distance);
            Fx += F * (other_body.x - body.x) / distance;
            Fy += F * (other_body.y - body.y) / distance;
        }
        ax = Fx / body.mass;
        ay = Fy / body.mass;
        body.vel_x += ax * time_step;
        body.vel_y += ay * time_step;
    }
}

void simulate(Body* bodies, uint32_t body_count, double time_step, uint32_t iterations){
    double elapsed_p0 = 0.0, elapsed_p1 = 0.0;
    for (uint32_t i = 0; i < iterations; i++){
        const auto cp0 = std::chrono::system_clock::now();
        update_body_positions(bodies, body_count, time_step);
        const auto cp1 = std::chrono::system_clock::now();
        update_body_velocities(bodies, body_count, time_step);
        const auto cp2 = std::chrono::system_clock::now();

        elapsed_p0 += std::chrono::duration_cast<std::chrono::nanoseconds>(cp1 - cp0).count();
        elapsed_p1 += std::chrono::duration_cast<std::chrono::nanoseconds>(cp2 - cp1).count();
    }
    std::cout << "Update positions (per iteration): " << elapsed_p0 / (1e9 * iterations) << std::endl;
    std::cout << "Velocity computation (per iteration): " << elapsed_p1 / (1e9 * iterations) << std::endl;
}

void parse_args(int argc, const char** argv, CFG& config){
    int32_t idx;
    IO::ArgParser arg_parser(argc, argv);
    try {
        if ((idx = arg_parser.get_next_idx("-in")) > 0)
            config.input_file = arg_parser.get(idx);

        if ((idx = arg_parser.get_next_idx("-out")) > 0)
            config.output_file = arg_parser.get(idx);

        if ((idx = arg_parser.get_next_idx("-it")) > 0)
            config.iterations = std::stoi(arg_parser.get(idx));

        if ((idx = arg_parser.get_next_idx("-it_len")) > 0)
            config.iter_len = std::stod(arg_parser.get(idx));
    }
    catch (const std::string& ex){
        std::cout << "Argument parsing exception: " << ex << std::endl;
        exit(1);
    }
    catch (const std::invalid_argument& e){
        std::cout << "Argument parsing exception: invalid argument." << std::endl;
        exit(1);
    }
    catch (const std::out_of_range& e){
        std::cout << "Argument parsing exception: out of range." << std::endl;
        exit(1);
    }
}

int main(int argc, const char** argv){
    std::vector<Body> bodies;
    parse_args(argc, argv, config);
    config.print();
    parse_input(config.input_file, bodies);

    const auto start = std::chrono::system_clock::now();
    simulate(bodies.data(), bodies.size(), config.iter_len, config.iterations);
    const auto end = std::chrono::system_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e9;
    std::cout << "Seconds: " << elapsed << std::endl;
    std::cout << "Seconds per iteration: " << elapsed / config.iterations << std::endl;

    write_output(config.output_file, bodies);
}
