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
    Quad *quad;
};

// https://www.geeksforgeeks.org/quad-tree/

// Used to hold details of a point
struct Point
{
    double x;
    double y;
    Point(double _x, double _y)
    {
        x = _x;
        y = _y;
    }
    Point()
    {
        x = 0;
        y = 0;
    }
};
  
// The objects that we want stored in the quadtree
// struct Node
// {
//     int id;
//     Point pos;
//     double mass;
//     Point vel;
//     Node(int _id, Point _pos, double _mass, Point _vel)
//     {
//         id = _id;
//         pos = _pos;
//         mass = _mass;
//         vel = _vel;
//     }
// };
  
// The main quadtree class
class Quad
{
    // Hold details of the boundary of this node
    Point topLeft;
    Point botRight;

    Quad *parent;
  
    // Children of this tree
    Quad *topLeftTree;
    Quad *topRightTree;
    Quad *botLeftTree;
    Quad *botRightTree;
  
public:
  
    // Contains details of node
    Body *body;
    Point *center;
    double mass;

    Quad()
    {
        topLeft = Point(0, 0);
        botRight = Point(0, 0);
        mass = 0;
        center = NULL;
        parent = NULL;
        body = NULL;
        topLeftTree  = NULL;
        topRightTree = NULL;
        botLeftTree  = NULL;
        botRightTree = NULL;
    }
    Quad(Point topL, Point botR)
    {
        mass = 0;
        center = NULL;
        parent = NULL;
        body = NULL;
        topLeftTree  = NULL;
        topRightTree = NULL;
        botLeftTree  = NULL;
        botRightTree = NULL;
        topLeft = topL;
        botRight = botR;
    }
    Quad(Point topL, Point botR, Quad _parent)
    {
        mass = 0;
        center = NULL;
        parent = _parent;
        body = NULL;
        topLeftTree  = NULL;
        topRightTree = NULL;
        botLeftTree  = NULL;
        botRightTree = NULL;
        topLeft = topL;
        botRight = botR;
    }
    void insert(Node*);
    Node* search(Point);
    bool inBoundary(Point);
};
  
// Insert a node into the quadtree
void Quad::insert(Body *b)
{
    if (b == NULL)
        return;
  
    // Current quad cannot contain it
    if (!inBoundary(b->pos))
        return;

    mass += b.mass;

    if (pos == NULL)
        pos = Point(b.x, b.y);
    else {
        weight = b.mass / mass;

        pos.x = weight * b.x + (1 - weight) * pos.x
    }
  
    // We are at a quad of unit area
    // We cannot subdivide this quad further
    if (abs(topLeft.x - botRight.x) <= 1 &&
        abs(topLeft.y - botRight.y) <= 1)
    {
        if (body == NULL)
        {
            body = b;
            b.quad = this;
        }
        return;
    }
  
    if ((topLeft.x + botRight.x) / 2 >= b->pos.x)
    {
        // Indicates topLeftTree
        if ((topLeft.y + botRight.y) / 2 >= b->pos.y)
        {
            if (topLeftTree == NULL)
                topLeftTree = new Quad(
                    Point(topLeft.x, topLeft.y),
                    Point((topLeft.x + botRight.x) / 2,
                        (topLeft.y + botRight.y) / 2),
                    this);
            topLeftTree->insert(b);
        }
  
        // Indicates botLeftTree
        else
        {
            if (botLeftTree == NULL)
                botLeftTree = new Quad(
                    Point(topLeft.x,
                        (topLeft.y + botRight.y) / 2),
                    Point((topLeft.x + botRight.x) / 2,
                        botRight.y),
                    this);
            botLeftTree->insert(b);
        }
    }
    else
    {
        
        // Indicates topRightTree
        if ((topLeft.y + botRight.y) / 2 >= b->pos.y)
        {
            if (topRightTree == NULL)
                topRightTree = new Quad(
                    Point((topLeft.x + botRight.x) / 2,
                        topLeft.y),
                    Point(botRight.x,
                        (topLeft.y + botRight.y) / 2),
                    this);
            topRightTree->insert(b);
        }
  
        // Indicates botRightTree
        else
        {
            if (botRightTree == NULL)
                botRightTree = new Quad(
                    Point((topLeft.x + botRight.x) / 2,
                        (topLeft.y + botRight.y) / 2),
                    Point(botRight.x, botRight.y),
                    this);
            botRightTree->insert(b);
        }
    }
}
  
// Check if current quadtree contains the point
bool Quad::inBoundary(Point p)
{
    return (p.x >= topLeft.x &&
        p.x <= botRight.x &&
        p.y >= topLeft.y &&
        p.y <= botRight.y);
}

void parse_input(const std::string& input_path, Quad *root, std::vector<Body>& bodies){
    Quad quadRoot = Quad(Point(-1e13, 1e13), Point(1e13, -1e13));
    
    try {
        auto parser = IO::Parser(input_path);
        auto io_body = IO::Body();
        while (parser.next_body_info(io_body) == 0){
            Body body = Body{ io_body.mass, io_body.x, io_body.y, io_body.vel_x, io_body.vel_y }
            quadRoot.insert(&body);
            bodies.push_back(body);
            // TODO: Parse values from io_body to correct the structure (dont use IO::Body for simulation)
        }
    }
    catch (const std::string& e) {
        std::cout << "Parse error: " << e << std::endl;
        exit(1);
    }

    root = &quadRoot;
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

void update_body_positions(Body* bodies, Quad quad, double time_step){
    for (uint32_t i = 0; i < body_count; i++){
        Body& body = bodies[i];
        auto delta_x = body.vel_x * time_step;
        auto delta_y = body.vel_y * time_step;
        body.x += delta_x;
        body.y += delta_y;
    }
}

void simulate(Body* bodies, Quad root, uint32_t body_count, double time_step, uint32_t iterations){
    for (int i = 0; i < iterations; i++) {
        update_body_positions(bodies, root, body_count, time_step);
        update_body_velocities(bodies, body_count, time_step);
    }
}

int main(int argc, const char** argv){
    std::vector<Body> bodies;
    Quad *quad = NULL;
    parse_input(DEF_IN, quad);

    const auto start = std::chrono::system_clock::now();
    simulate(/* ..., */ time_step_length, total_time_steps);
    const auto end = std::chrono::system_clock::now();
    
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e9;
    std::cout << "Seconds: " << elapsed << std::endl;

    write_output(DEF_OUT /*, ... */);
}
