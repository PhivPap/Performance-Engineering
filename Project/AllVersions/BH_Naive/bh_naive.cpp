#include <iostream>
#include <chrono>
#include <vector>
#include <math.h>
#include "../../InputOutput/InputOutput.h"
#include "body.h"
#include "quad.h"
#include "util.h"
#include <assert.h>
#include <algorithm>
#include "float.h"


// default values can be altered with the main args
struct CFG {
    std::string input_file = "BodyFiles/in/def_in.tsv";
    std::string output_file = "BodyFiles/out/bh_naive_out_050.tsv";
    uint32_t iterations = 50;
    double iter_len = 3600;
    double theta = 0.50;

    void print(){
        std::cout << "Configuration:" << std::endl;
        std::cout   << "Input: " << input_file << "\nOutput: " << output_file 
                    << "\nIterations: " << iterations << "\nIteration legth: "
                    << iter_len << "s\nTheta: " << theta << std::endl << std::endl;
    }
};

const double G = 6.67e-11;        // Gravitational constant
CFG config;
//uint64_t approximations, total_visits;


void parse_input(const std::string& input_path, std::vector<Body>& bodies){
    try {
        IO::Parser parser(input_path);
        auto io_body = IO::Body();
        while (parser.next_body_info(io_body) == 0)
            bodies.push_back(Body(io_body.mass, { io_body.x, io_body.y }, io_body.vel_x, io_body.vel_y));
    }
    catch (const std::string& e) {
        std::cout << "Parse error: " << e << std::endl;
        exit(1);
    }
}

void write_output(const std::string& output_path, const std::vector<Body>& bodies){
    try {
        IO::Writer writer(output_path);
        uint32_t id = 0;
        for (const auto body : bodies)
            writer.write_body(IO::Body(id++, body.mass, body.coords.x, body.coords.y, body.vel_x, body.vel_y));
    }
    catch (const std::string& e) {
        std::cout << "Write error: " << e << std::endl;
        exit(1);
    }
}

Area update_body_positions_and_get_area(Body* bodies, uint32_t body_count, double time_step){
    Area area(DBL_MAX, DBL_MIN, DBL_MAX, DBL_MIN);

    for (uint32_t i = 0; i < body_count; i++){
        Body& body = bodies[i];
        auto delta_x = body.vel_x * time_step;
        auto delta_y = body.vel_y * time_step;
        body.coords.x += delta_x;
        body.coords.y += delta_y;
        
        area.x1 = std::min(area.x1, body.coords.x);
        area.x2 = std::max(area.x2, body.coords.x);
        area.y1 = std::min(area.y1, body.coords.y);
        area.y2 = std::max(area.y2, body.coords.y);
    }
    return area;
}

void compute_body2body_attraction(const Body* body1, const Body* body2, double& Fx, double& Fy) {
    if (body1 == body2) // assume that it's almost always false and it's correctly predicted by the cpu
        return;
    const double distance = body1->coords.distance_to(body2->coords);
    const double F = (G * body1->mass * body2->mass) / (distance * distance);
    Fx += F * (body2->coords.x - body1->coords.x) / distance;
    Fy += F * (body2->coords.y - body1->coords.y) / distance;
}

void compute_body2quad_attraction(const Body* body, const Quad* quad, double& Fx, double& Fy) {
    const double distance = body->coords.distance_to(quad->center_of_mass);
    const double F = (G * body->mass * quad->mass) / (distance * distance);
    Fx += F * (quad->center_of_mass.x - body->coords.x) / distance;
    Fy += F * (quad->center_of_mass.y - body->coords.y) / distance;
}

void compute_body_forces(Quad* quad, Body* body, double& Fx, double& Fy){
    //total_visits++;
    const auto quad_body_count = quad->body_count;
    if (quad_body_count == 0)
        return;
    if (quad_body_count == 1) {
        //approximations += 1;
        compute_body2body_attraction(body, quad->contained_bodies.front(), Fx, Fy);
        return;
    }

    // magic formula check https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation
    const auto quad_diag = quad->diag_len;
    const auto distance = body->coords.distance_to(quad->center_of_mass);

    if (quad_diag / distance < config.theta){
        //approximations += 1;
        compute_body2quad_attraction(body, quad, Fx, Fy);
    }
    else {
        compute_body_forces(quad->top_left_quad, body, Fx, Fy);
        compute_body_forces(quad->top_right_quad, body, Fx, Fy);
        compute_body_forces(quad->bot_left_quad, body, Fx, Fy);
        compute_body_forces(quad->bot_right_quad, body, Fx, Fy);
    }
}

void update_body_velocities(Quad* root, Body* bodies, uint32_t body_count, double time_step){
    for (uint32_t i = 0; i < body_count; i++){
        Body& body = bodies[i];
        double Fx = 0, Fy = 0;
        compute_body_forces(root, &body, Fx, Fy);
        body.vel_x += (Fx / body.mass) * time_step;
        body.vel_y += (Fy / body.mass) * time_step;
    }
}


void simulate(Body* bodies, uint32_t body_count, double time_step, uint32_t iterations){
    double elapsed_p0 = 0.0, elapsed_p1 = 0.0, elapsed_p2 = 0.0;
    for (int i = 0; i < iterations; i++) {
        //approximations = total_visits = 0;

        //const auto cp0 = std::chrono::system_clock::now();
        Area area = update_body_positions_and_get_area(bodies, body_count, time_step);
        //const auto cp1 = std::chrono::system_clock::now();
        Quad root(bodies, body_count, area);
        //const auto cp2 = std::chrono::system_clock::now();
        update_body_velocities(&root, bodies, body_count, time_step);
        //const auto cp3 = std::chrono::system_clock::now();

        // elapsed_p0 += std::chrono::duration_cast<std::chrono::nanoseconds>(cp1 - cp0).count();
        // elapsed_p1 += std::chrono::duration_cast<std::chrono::nanoseconds>(cp2 - cp1).count();
        // elapsed_p2 += std::chrono::duration_cast<std::chrono::nanoseconds>(cp3 - cp2).count();

        // std::cout << "Iteration: " << i << std::endl;
        // std::cout << "\tDiagonal Length: " << area.diagonal_length() << std::endl;
        // std::cout << "\tApprox factor: " << (double)approximations / (body_count * body_count) << std::endl;
        // std::cout << "\tAvg visits (per body): " << (double)total_visits / body_count << std::endl;
        
    }
    // std::cout << "Update positions (per iteration): " << elapsed_p0 / (1e9 * iterations) << "s" << std::endl;
    // std::cout << "QuadTree generation (per iteration): " << elapsed_p1 / (1e9 * iterations) << "s" << std::endl;
    // std::cout << "Velocity computation (per iteration): " << elapsed_p2 / (1e9 * iterations) << "s" << std::endl;
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

        if ((idx = arg_parser.get_next_idx("-theta")) > 0)
            config.theta = std::stod(arg_parser.get(idx));
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
