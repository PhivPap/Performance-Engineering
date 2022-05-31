#include "quad.h"
#include <iostream>
#include <chrono>

// double Quad::elapsed0 = 0.0;
// double Quad::elapsed1 = 0.0;
// double Quad::counter = 0.0;


Quad::Quad(const Area& area) : area(area), mass(0), body_count(0), center_of_mass({0,0}) {
    top_left_quad = top_right_quad = bot_left_quad = bot_right_quad = nullptr;
    diag_len_2 = area.diagonal_length_2();
}

// this constructor is used to generate the root of the quad tree
Quad::Quad(Body* bodies, uint32_t body_count, const Area& area) : center_of_mass({0,0}), body_count(0),  
            mass(0), area(area) {
    diag_len_2 = area.diagonal_length_2();
    for (uint32_t i = 0; i < body_count; i++)
        insert_body(&(bodies[i]));
    compute_bhtree_recursive();
    recursive_center_of_mass_computation();
}

Quad::~Quad() {
    if (body_count > 1){
        delete top_left_quad;
        delete top_right_quad;
        delete bot_left_quad;
        delete bot_right_quad;
    }
}

void Quad::insert_body(Body* body){
    contained_bodies.push_back(body);
    mass += body->mass;
    body_count++;
}

void Quad::compute_bhtree_recursive(void){
    if (body_count <= 1)
        return;

    const Point& center = area.get_center();

    //const auto p0 = std::chrono::high_resolution_clock::now();

    top_left_quad = new Quad({ area.x1, center.x, area.y1, center.y });
    top_right_quad = new Quad({ center.x, area.x2, area.y1, center.y });
    bot_left_quad = new Quad({ area.x1, center.x, center.y, area.y2 });
    bot_right_quad = new Quad({ center.x, area.x2, center.y, area.y2 });

    //const auto p1 = std::chrono::high_resolution_clock::now();

    for (Body* body : contained_bodies){
        const Point& coords = body->coords;

        if (coords.x > center.x) { // right
            if (coords.y > center.y) // bottom
                bot_right_quad->insert_body(body);
            else // top
                top_right_quad->insert_body(body);
        }
        else { // left
            if (coords.y > center.y) // bottom
                bot_left_quad->insert_body(body);
            else  // top
                top_left_quad->insert_body(body);
        }
    }

    //const auto p2 = std::chrono::high_resolution_clock::now();

    //elapsed0 += std::chrono::duration_cast<std::chrono::nanoseconds>(p1 - p0).count();
    //elapsed1 += std::chrono::duration_cast<std::chrono::nanoseconds>(p2 - p1).count() / contained_bodies.size();
    //counter += 1;

        // std::cout << "Subtree generation: " << elapsed0 << "ns" << std::endl;
        // std::cout << "Body insertion: " << elapsed1 / contained_bodies.size() << "ns" << std::endl;


    //level++;
    top_left_quad->compute_bhtree_recursive();
    top_right_quad->compute_bhtree_recursive();
    bot_left_quad->compute_bhtree_recursive();
    bot_right_quad->compute_bhtree_recursive();
}

// must be called on the root of the tree once the tree is generated.
void Quad::recursive_center_of_mass_computation(void){
    if (body_count == 0)
        return;
    else if (body_count == 1)
        center_of_mass = contained_bodies.front()->coords;
    else {
        top_left_quad->recursive_center_of_mass_computation();
        top_right_quad->recursive_center_of_mass_computation();
        bot_left_quad->recursive_center_of_mass_computation();
        bot_right_quad->recursive_center_of_mass_computation();

        // computes the center of mass by using the center of mass of it's four children
        center_of_mass = Point::get_center_of_mass(
            Point::get_center_of_mass(
                top_left_quad->center_of_mass, top_left_quad->mass,
                top_right_quad->center_of_mass, top_right_quad->mass
            ), top_left_quad->mass + top_right_quad->mass,
            Point::get_center_of_mass(
                bot_left_quad->center_of_mass, bot_left_quad->mass,
                bot_right_quad->center_of_mass, bot_right_quad->mass
            ), bot_left_quad->mass + bot_right_quad->mass
        );
    }
}