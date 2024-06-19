#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

class Vector {
    float x, y, z;
public:
    Vector(const float x, const float y, const float z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    float length() const {
        return sqrt(x * x + y * y + z * z);
    }
    void normalize() {
        float cachedLength = length();
        cachedLength = 1 / cachedLength;
        x = x * cachedLength;
        y = y * cachedLength;
        z = z * cachedLength;
    }
    void scale(const float num) {
        x *= num;
        y *= num;
        z *= num;
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
    Vector add(const Vector a, const Vector b) {
        return Vector(a.getX() + b.getX(), 
                        a.getY() + b.getY(), 
                        a.getZ() + b.getZ());
    }

    Vector subtract(const Vector a, const Vector b) {
        return Vector(a.getX() - b.getX(), 
                        a.getY() - b.getY(), 
                        a.getZ() - b.getZ());
    }

    Vector crossProduct(const Vector a, const Vector b) {
        float x = a.getY() * b.getZ() - a.getZ() * b.getY();
        float y = a.getZ() * b.getX() - a.getX() * b.getZ();
        float z = a.getX() * b.getY() - a.getY() * b.getX();
        return Vector(x, y, z);
    }

    float dotProduct(const Vector a, const Vector b) {
        return a.getX() * b.getX() 
                + a.getY() * b.getY() 
                + a.getZ() * b.getZ();
    }
};

#define NUM_OF_TRIANGLE_VERTICES 3

class Triangle {
    Vector vertex[NUM_OF_TRIANGLE_VERTICES];
    Vector normal {0, 0, 0};
public:
    Triangle(Vector vertex0, Vector vertex1, Vector vertex2) 
        : vertex {vertex0, vertex1, vertex2}
        {}
    
    Vector calculateNormal() const {
        Vector first = Math::subtract(vertex[1], vertex[0]);
        Vector second = Math::subtract(vertex[2], vertex[0]);
        return Math::crossProduct(first, second);
    }

    void cacheNormal() {
        normal = calculateNormal();
        normal.normalize();
    }
    void cacheNormal(const Vector v) {
        normal = v;
    }
    
    Vector getVertex0() const {
        return vertex[0];
    }
    Vector getVertex1() const {
        return vertex[1];
    }
    Vector getVertex2() const {
        return vertex[2];
    }
    Vector getVertex(const unsigned int index) const {
        return vertex[index];
    }
    Vector getNormal() const {
        return normal;
    }

    void setVertex0(const Vector v) {
        vertex[0] = v;
    }
    void setVertex1(const Vector v) {
        vertex[1] = v;
    }
    void setVertex2(const Vector v) {
        vertex[2] = v;
    }
    void setVertex(const unsigned int index, const Vector v) {
        vertex[index] = v;
    }
};

static const int imageWidth = 1920;
static const int imageHeight = 1080;
static const int maxColorComponent = 255;

#define SET_COLORS(COLORS, TR_INDEX)    for(unsigned int i = 0; i < 3; i++) { \
                                            COLORS[i] = ((TR_INDEX + 1) * (i + 16) * 64 \
                                                            + (i+1) * 64) % 256; \
                                        }

#define CLEAR_COLORS(COLORS)    for(unsigned int i = 0; i < 4; i++) { \
                                    COLORS[i] = 0; \
                                }

int main() {
    const unsigned int trCount = 6;
    Triangle triangles[trCount] {
        Triangle(Vector(-1, -1.75, -3.3), Vector(1, -1.75, -3.3), Vector(0, 0, -3)),
        Triangle(Vector(0, 0, -3), Vector(1, -1.75, -3.3), Vector(1, 1.75, -3.3)),
        Triangle(Vector(-1, -1.75, -3.3), Vector(0, 0, -3), Vector(-1, 1.75, -3.3)),
        Triangle(Vector(-1, 1.75, -3.3), Vector(0, 0, -3), Vector(1, 1.75, -3.3)),
        Triangle(Vector(-1, -1.75, -3.3), Vector(1, -1.75, -3.3), Vector(-1, 1.75, -3.3)),
        Triangle(Vector(-1, 1.75, -3.3), Vector(1, -1.75, -3.3), Vector(1, 1.75, -3.3))
    };
    for(unsigned int i = 0; i < trCount; i++) {
        triangles[i].cacheNormal();
    }

    std::ofstream ppmFileStream("3D_figure_1920_1080.ppm", 
                    std::ios::out | std::ios::binary);
    ppmFileStream << "P3\n";
    ppmFileStream << imageWidth << " " << imageHeight << "\n";
    ppmFileStream << maxColorComponent << "\n";
    for (int rowIdx = 0; rowIdx < imageHeight; ++rowIdx) {
        for (int colIdx = 0; colIdx < imageWidth; ++colIdx) {
            float x = colIdx + 0.5, y = rowIdx + 0.5;
            x /= imageWidth;
            y /= imageHeight;
            x = (2 * x) - 1;
            y = 1 - (2 * y);
            x *= imageWidth / imageHeight;
            float colors[4] = {0, 0, 0, 0};
            for(unsigned int trIndex = 0; trIndex < trCount; trIndex++) {
                Vector cameraDirection(x, y, -1);
                cameraDirection.normalize();
                if(Math::dotProduct(cameraDirection, triangles[trIndex].getNormal()) != 0 
                    && Math::dotProduct(triangles[trIndex].getNormal(), 
                                            triangles[trIndex].getVertex0()) < 0)
                {
                    float scaleFactor = fabs(Math::dotProduct(
                        triangles[trIndex].getVertex0(), 
                        triangles[trIndex].getNormal()
                    ) / Math::dotProduct(cameraDirection, triangles[trIndex].getNormal()));
                    cameraDirection.scale(scaleFactor);
                    bool isIn = true;
                    for(unsigned int i = 0; i < NUM_OF_TRIANGLE_VERTICES; i++) {
                        unsigned int tail = i;
                        unsigned int head = (i + 1) % NUM_OF_TRIANGLE_VERTICES;
                        Vector temp = Math::crossProduct(
                            Math::subtract(triangles[trIndex].getVertex(head), 
                                            triangles[trIndex].getVertex(tail)), 
                            Math::subtract(cameraDirection, 
                                            triangles[trIndex].getVertex(tail))
                        );
                        float checkSum = Math::dotProduct(triangles[trIndex].getNormal(), 
                                                            temp);
                        if(checkSum <= 0) {
                            //if(colors[3] == scaleFactor) CLEAR_COLORS(colors);
                            isIn = false;
                            break;
                        }
                    }
                    if(isIn && (colors[3] == 0 || colors[3] > scaleFactor)) {
                        SET_COLORS(colors, trIndex);
                        colors[3] = scaleFactor;
                    }
                }
            }
            ppmFileStream << colors[0] 
            << " " << colors[1] 
            << " " << colors[2] 
            << "\t";
        }
        ppmFileStream << "\n";
    }
    ppmFileStream.close();
    return 0;
}