#ifndef QUAD
#define QUAD

#include "util.h"
#include "body.h"
#include <list>
#include <inttypes.h>

// The main quadtree class
class Quad{
private:
    void compute_bhtree_recursive(int level);

public:
    Area area; 
    uint32_t body_count;
    double mass;
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
    void recursive_center_of_mass_computation(void);
};


#endif