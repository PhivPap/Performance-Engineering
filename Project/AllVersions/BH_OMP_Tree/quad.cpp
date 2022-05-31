#include "quad.h"
#include <iostream>
#include <assert.h>

uint8_t Quad::MAX_TASK_GEN_DEPTH = 0;
Quad *Quad::pool;
uint32_t Quad::pool_size;
uint32_t Quad::pool_idx;

Quad::Quad(void) : mass(0), body_count(0), center_of_mass({0, 0})
{
    contained_bodies.reserve(100);
    top_left_quad = top_right_quad = bot_left_quad = bot_right_quad = nullptr;
}

// this constructor is used to generate the root of the quad tree
Quad::Quad(Body *bodies, uint32_t body_count, const Area &area) : center_of_mass({0, 0}), body_count(0),
                                                                 mass(0), area(area)
{
    diag_len_2 = area.diagonal_length_2();
    contained_bodies.reserve(body_count);
    for (uint32_t i = 0; i < body_count; i++)
        insert_body(&(bodies[i]));

    compute_bhtree_recursive();
}

Quad::~Quad()
{
    // if (body_count > 1)
    //     delete[] top_left_quad;
}

void Quad::set_area(const Area &area)
{
    this->area = area;
    diag_len_2 = area.diagonal_length_2();
}

void Quad::insert_body(Body *body)
{
    contained_bodies.push_back(body);
    mass += body->mass;
    body_count++;
}

void Quad::compute_bhtree_recursive(void)
{
    if (body_count == 0)
        return;
    else if (body_count == 1)
    {
        center_of_mass = contained_bodies[0]->coords;
        return;
    }

    const Point &center = area.get_center();

    auto my_pool_index = Quad::pool_get_idx();
    auto subtree = &(Quad::pool[my_pool_index]);

    top_left_quad = subtree;
    top_right_quad = &(subtree[1]);
    bot_left_quad = &(subtree[2]);
    bot_right_quad = &(subtree[3]);

    top_left_quad->set_area({area.x1, center.x, area.y1, center.y});
    top_right_quad->set_area({center.x, area.x2, area.y1, center.y});
    bot_left_quad->set_area({area.x1, center.x, center.y, area.y2});
    bot_right_quad->set_area({center.x, area.x2, center.y, area.y2});

    const uint32_t n = contained_bodies.size();
    auto cb =  contained_bodies.data();
    for (uint32_t i = 0; i < n; i++)
    {
        Body* body = cb[i];
        const Point &coords = body->coords;

        if (coords.x > center.x)
        {                            // right
            if (coords.y > center.y) // bottom
                bot_right_quad->insert_body(body);
            else // top
                top_right_quad->insert_body(body);
        }
        else
        {                            // left
            if (coords.y > center.y) // bottom
                bot_left_quad->insert_body(body);
            else // top
                top_left_quad->insert_body(body);
        }
    }

    top_left_quad->compute_bhtree_recursive();
    top_right_quad->compute_bhtree_recursive();
    bot_left_quad->compute_bhtree_recursive();
    bot_right_quad->compute_bhtree_recursive();

    center_of_mass = Point::get_center_of_mass(
        Point::get_center_of_mass(
            top_left_quad->center_of_mass, top_left_quad->mass,
            top_right_quad->center_of_mass, top_right_quad->mass),
        top_left_quad->mass + top_right_quad->mass,
        Point::get_center_of_mass(
            bot_left_quad->center_of_mass, bot_left_quad->mass,
            bot_right_quad->center_of_mass, bot_right_quad->mass),
        bot_left_quad->mass + bot_right_quad->mass);
}

void Quad::set_max_task_generation_depth(uint8_t max_task_gen_depth)
{
    MAX_TASK_GEN_DEPTH = max_task_gen_depth;
}

void Quad::set_pool(uint32_t init_size)
{
    assert(init_size % 4 == 0);
    pool_size = init_size;
    pool = new Quad[init_size];
}

void Quad::reset_pool(void)
{
    delete[] pool;
    pool_idx = 0;
}

uint32_t Quad::pool_get_idx(void)
{
    double my_idx = pool_idx;
    pool_idx += 4;
    return my_idx;
}