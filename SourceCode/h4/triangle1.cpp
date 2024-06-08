#include <iostream>
#include <cmath>

class Vector {
    float x, y, z;
public:
    Vector(float x, float y, float z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    float length() const {
        return sqrt(x * x + y * y + z * z);
    }
    void normalize() {
        float cachedLength = length();
        x = x / cachedLength;
        y = y / cachedLength;
        z = z / cachedLength;
    }

    float getX() const {
        return x;
    }
    float getY() const {
        return y;
    }
    float getZ() const {
        return z;
    }

    void setX(float x) {
        this->x = x;
    }
    void setY(float y) {
        this->y = y;
    }
    void setZ(float z) {
        this->z = z;
    }

    void print() const {
        std::cout << x
        << " " << y
        << " " << z
        << std::endl;
    }
};

namespace Math {
    Vector add(Vector a, Vector b) {
        return Vector(a.getX() + b.getX(), 
                        a.getY() + b.getY(), 
                        a.getZ() + b.getZ());
    }

    Vector subtract(Vector a, Vector b) {
        return Vector(a.getX() - b.getX(), 
                        a.getY() - b.getY(), 
                        a.getZ() - b.getZ());
    }

    Vector crossProduct(Vector a, Vector b) {
        float x = a.getY() * b.getZ() - a.getZ() * b.getY();
        float y = a.getZ() * b.getX() - a.getX() * b.getZ();
        float z = a.getX() * b.getY() - a.getY() * b.getX();
        return Vector(x, y, z);
    }
};
// Task 1
class Triangle {
    Vector vertex0, vertex1, vertex2;
public:
    Triangle(Vector vertex0, Vector vertex1, Vector vertex2) 
        : vertex0(vertex0), 
          vertex1(vertex1), 
          vertex2(vertex2)
        {}
    
    Vector getNormal() const {
        Vector first = Math::subtract(vertex1, vertex0);
        Vector second = Math::subtract(vertex2, vertex0);
        return Math::crossProduct(first, second);
    }
    
    Vector getVertex0() const {
        return vertex0;
    }
    Vector getVertex1() const {
        return vertex1;
    }
    Vector getVertex2() const {
        return vertex2;
    }

    void setVertex0(Vector v) {
        vertex0 = v;
    }
    void setVertex1(Vector v) {
        vertex1 = v;
    }
    void setVertex2(Vector v) {
        vertex2 = v;
    }
};

struct Ray {
    Vector direction {0, 0, -1};
    Vector origin {0, 0, 0};
};

int main() {
    // Task 2
    std::cout << "Task 2:" << std::endl;
    Vector result(0, 0, 0);
    result = Math::crossProduct(Vector(3.5, 0, 0), Vector(1.75, 3.5, 0));
    result.print();

    result = Math::crossProduct(Vector(3, -3, 1), Vector(4, 9, 3));
    result.print();
    std::cout << result.length() << std::endl;

    result = Math::crossProduct(Vector(3, -3, 1), Vector(-12, 12, -4));
    std::cout << result.length() << std::endl;


    // Task 3
    std::cout << std::endl << "Task 3:" << std::endl;
    Triangle tr(Vector(-1.75, -1.75, -3), Vector(1.75, -1.75, -3), Vector(0, 1.75, -3));
    result = tr.getNormal();
    result.print();
    std::cout << "Area: " << result.length() / 2 << std::endl;

    tr.setVertex0(Vector(0, 0, -1));
    tr.setVertex1(Vector(1, 0, 1));
    tr.setVertex2(Vector(-1, 0, 1));
    result = tr.getNormal();
    result.print();
    std::cout << "Area: " << result.length() / 2 << std::endl;

    tr.setVertex0(Vector(0.56, 1.11, 1.23));
    tr.setVertex1(Vector(0.44, -2.368, -0.54));
    tr.setVertex2(Vector(-1.56, 0.15, -1.92));
    result = tr.getNormal();
    result.print();
    std::cout << "Area: " << result.length() / 2 << std::endl;
    return 0;
}