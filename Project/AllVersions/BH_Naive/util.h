#ifndef UTIL
#define UTIL

class Point{
public:
    double x, y;
    Point() = delete;
    Point(double x, double y);
    Point(const Point& p);
    double distance_to(const Point& p2) const;
    static Point get_center_of_mass(const Point& p1, double mass1, const Point& p2, double mass2);
};

class Area{
public:
    double x1, x2, y1, y2;
    Area() = delete;
    Area(double x1, double x2, double y1, double y2);
    bool contains(const Point& p) const;
    double side_length() const;
    double diagonal_length() const;
    Point get_center() const;
};


#endif