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

struct Body {
    double mass, x, y, vel_x, vel_y;
};



void parse_input(const std::string& input_path, std::vector<Body>& bodies){
    try {
        auto parser = IO::Parser(input_path);
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
        auto writer = IO::Writer(output_path);
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
    for (uint32_t i = 0; i < iterations; i++){
        update_body_positions(bodies, body_count, time_step);
        update_body_velocities(bodies, body_count, time_step);
    }
}

int main(int argc, const char** argv){
    std::vector<Body> bodies;
    parse_input(DEF_IN, bodies);

    const auto start = std::chrono::system_clock::now();
    simulate(bodies.data(), bodies.size(), time_step_length, total_time_steps);
    const auto end = std::chrono::system_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e9;
    std::cout << "Seconds: " << elapsed << std::endl;

    write_output(DEF_OUT, bodies);
}
