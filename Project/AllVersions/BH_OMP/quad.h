#ifndef QUAD
#define QUAD

#include "util.h"
#include "body.h"
#include <list>
#include <inttypes.h>

// The main quadtree class
class Quad{
private:
    void compute_bhtree_recursive(uint8_t depth);
    void compute_bhtree_recursive_seq(void);
    //void recursive_center_of_mass_computation(void);
    static uint8_t MAX_TASK_GEN_DEPTH;

public:
    Area area; 
    uint32_t body_count;
    double mass, diag_len;
    Point center_of_mass;
    std::list<Body*> contained_bodies;
    
    // Children of this tree
    Quad* top_left_quad;
    Quad* top_right_quad;
    Quad* bot_left_quad;
    Quad* bot_right_quad;

    Quad(void) = delete;
    Quad(const Area& area);
    Quad(Body* bodies, uint32_t body_count, const Area& area);
    ~Quad();
    void insert_body(Body* body);
    static void set_max_task_generation_depth(uint8_t max_task_gen_depth);
};


#endif