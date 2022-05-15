#include "quad.h"
#include <iostream>

Quad::Quad(const Area& area) : area(area), mass(0), body_count(0), center_of_mass({0,0}), diag_len(0) {
    top_left_quad = top_right_quad = bot_left_quad = bot_right_quad = nullptr;
    diag_len = area.diagonal_length();
}

// this constructor is used to generate the root of the quad tree
Quad::Quad(Body* bodies, uint32_t body_count, const Area& area) : center_of_mass({0,0}), body_count(0),  
            mass(0), area(area) {
    diag_len = area.diagonal_length();
    for (uint32_t i = 0; i < body_count; i++)
        insert_body(&(bodies[i]));
    compute_bhtree_recursive(0);
    recursive_center_of_mass_computation(0);
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

void Quad::compute_bhtree_recursive(int level){
    // std::cout << "body_count: " << body_count << ", level: " << level << std::endl;

    if (body_count <= 1)
        return;

    const Point& center = area.get_center();
    top_left_quad = new Quad({ area.x1, center.x, area.y1, center.y });
    top_right_quad = new Quad({ center.x, area.x2, area.y1, center.y });
    bot_left_quad = new Quad({ area.x1, center.x, center.y, area.y2 });
    bot_right_quad = new Quad({ center.x, area.x2, center.y, area.y2 });

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

    top_left_quad->compute_bhtree_recursive(level + 1);
    top_right_quad->compute_bhtree_recursive(level + 1);
    bot_left_quad->compute_bhtree_recursive(level + 1);
    bot_right_quad->compute_bhtree_recursive(level + 1);
}

// must be called on the root of the tree once the tree is generated.
void Quad::recursive_center_of_mass_computation(int level){
    if (body_count == 0)
        return;
    else if (body_count == 1)
        center_of_mass = contained_bodies.front()->coords;
    else {
        top_left_quad->recursive_center_of_mass_computation(level + 1);
        top_right_quad->recursive_center_of_mass_computation(level + 1);
        bot_left_quad->recursive_center_of_mass_computation(level + 1);
        bot_right_quad->recursive_center_of_mass_computation(level + 1);

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


    // std::string prelude = "";
    // for (int i = 0; i < level; i++)
    //     prelude += "\t";

    // std::cout << prelude << "x: " << area.x1 << " <= " <<center_of_mass.x << " <= " << area.x2 << std::endl;
    // std::cout << prelude <<  "y: " << area.y1 << " <= " <<center_of_mass.y << " <= " << area.y2 << std::endl;
    
}