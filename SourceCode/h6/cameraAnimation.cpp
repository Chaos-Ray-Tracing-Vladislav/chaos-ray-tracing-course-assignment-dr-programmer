#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iterator>

#define MAX_MATRIX_ROWS 4
#define MAX_MATRIX_COLUMNS 4

class Matrix {
    float fields[MAX_MATRIX_ROWS][MAX_MATRIX_COLUMNS] = {0};
    unsigned int rows, columns;
public:
    Matrix(const float * const matrix, const unsigned int rows, const unsigned int columns) {
        if(!matrix) return;
        this->rows = rows;
        this->columns = columns;
        for(unsigned int i = 0; i < rows; i++) {
            for(unsigned int j = 0; j < columns; j++) {
                fields[i][j] = matrix[j + columns * i];
            }
        }
    }
    Matrix(const unsigned int rows, const unsigned int columns) {
        this->rows = rows;
        this->columns = columns;
    }

    unsigned int getRows() const {
        return rows;
    }
    unsigned int getColumns() const {
        return columns;
    }
    float getField(const unsigned int i, const unsigned int j) const {
        return fields[i][j];
    }
    float (&getFields())[MAX_MATRIX_ROWS][MAX_MATRIX_COLUMNS] {
        return fields;
    }

    void setRows(const unsigned int rows) {
        this->rows = rows;
    }
    void setColumns(const unsigned int columns) {
        this->columns = columns;
    }
    void setField(const unsigned int i, const unsigned int j, const float value) {
        fields[i][j] = value;
    }
    void setFields(const float * const matrix) {
        for(unsigned int i = 0; i < rows; i++) {
            for(unsigned int j = 0; j < columns; j++) {
                fields[i][j] = matrix[j + columns * i];
            }
        }
    }

    void print() const {
        for(unsigned int i = 0; i < rows; i++) {
            for(unsigned int j = 0; j < columns; j++) {
                std::cout << fields[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
};

class Vector : public Matrix {
public:
    Vector(const float x, const float y, const float z) 
        : Matrix(std::begin({x, y, z}), 1, 3) 
        {}
    Vector(const Matrix &m) 
        : Matrix(
            std::begin({m.getField(0, 0), m.getField(0, 1), m.getField(0, 2)}), 
            1, 
            3
        )
        {}

    float length() const {
        return sqrt(getField(0, 0) * getField(0, 0) 
                        + getField(0, 1) * getField(0, 1)
                        + getField(0, 2) * getField(0, 2));
    }
    void normalize() {
        float cachedLength = length();
        cachedLength = 1 / cachedLength;
        setField(0, 0, getField(0, 0) * cachedLength);
        setField(0, 1, getField(0, 1) * cachedLength);
        setField(0, 2, getField(0, 2) * cachedLength);
    }
    void scale(const float num) {
        setField(0, 0, getField(0, 0) * num);
        setField(0, 1, getField(0, 1) * num);
        setField(0, 2, getField(0, 2) * num);
    }

    float getX() const {
        return getField(0, 0);
    }
    float getY() const {
        return getField(0, 1);
    }
    float getZ() const {
        return getField(0, 2);
    }
    float get(const unsigned int index) const {
        return getField(0, index);
    }

    void setX(const float x) {
        setField(0, 0, x);
    }
    void setY(const float y) {
        setField(0, 1, y);
    }
    void setZ(const float z) {
        setField(0, 2, z);
    }
    void set(const unsigned int index, const float value) {
        setField(0, index, value);
    }

    void print() const {
        Matrix::print();
    }

    Vector &operator=(const Matrix& m) {
        for(unsigned int i = 0; i < 3; i++) {
            setField(0, i, m.getField(0, i));
        }
        return *this;
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

    Matrix matrixMultiply(Matrix a, Matrix b) {
        Matrix result(a.getRows(), b.getColumns());
        auto getValue = [&a, &b](const unsigned int i, const unsigned int j) -> float {
            float result = 0;
            for(unsigned int add = 0; add < a.getColumns(); add++) {
                result += a.getField(i, add) * b.getField(add, j);
            }
            return result;
        };
        for(unsigned int i = 0; i < a.getRows(); i++) {
            for(unsigned int j = 0; j < b.getColumns(); j++) {
                result.setField(
                    i, 
                    j, 
                    getValue(i, j)
                );
            }
        }
        return result;
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
    
    Vector& getVertex0() {
        return vertex[0];
    }
    Vector& getVertex1() {
        return vertex[1];
    }
    Vector& getVertex2() {
        return vertex[2];
    }
    Vector& getVertex(const unsigned int index) {
        return vertex[index];
    }
    Vector& getNormal() {
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

class Camera {
    Vector position, imagePlane, direction, up;

    void applyRotation(Matrix& transformation) {
        imagePlane = Math::matrixMultiply(imagePlane, transformation);
        direction = Math::matrixMultiply(direction, transformation);
        up = Math::matrixMultiply(up, transformation);
    }

public:
    Camera(const Vector position, const Vector imagePlane) 
        : position(position), 
          imagePlane(imagePlane), 
          direction(Vector(0, 0, imagePlane.getZ())), 
          up(Vector(0, -imagePlane.getZ(), 0))
        {
            this->imagePlane.normalize();
            this->direction.normalize();
            this->up.normalize();
        }
    
    void pan(const float angle) {
        const float rAngle = angle * (M_PI / 180);
        Matrix transformation(
            std::begin({
                cosf(rAngle), 0.0f, -sinf(rAngle), 
                0.0f, 1.0f, 0.0f, 
                sinf(rAngle), 0.0f, cosf(rAngle)
            }),
            3,
            3
        );
        applyRotation(transformation);
    }
    void tilt(const float angle) {
        const float rAngle = angle * (M_PI / 180);
        Matrix transformation(
            std::begin({
                1.0f, 0.0f, 0.0f, 
                0.0f, cosf(rAngle), sinf(rAngle), 
                0.0f, -sinf(rAngle), cosf(rAngle)
            }),
            3, 
            3
        );
        applyRotation(transformation);
    }
    void roll(const float angle) {
        const float rAngle = angle * (M_PI / 180);
        Matrix transformation(
            std::begin({
                cosf(rAngle), sinf(rAngle), 0.0f, 
                -sinf(rAngle), cosf(rAngle), 0.0f, 
                0.0f, 0.0f, 1.0f
            }),
            3, 
            3
        );
        applyRotation(transformation);
    }
    void dolly(const float distance) {
        Vector translation(direction.getX(), direction.getY(), direction.getZ());
        translation.scale(-distance);
        position = Math::add(position, translation);
    }
    void truck(const float distance) {
        const Vector left = Math::crossProduct(direction, up);
        Vector translation(left.getX(), left.getY(), left.getZ());
        translation.scale(distance);
        position = Math::add(position, translation);
    }
    void pedestal(const float distance) {
        Vector translation(up.getX(), up.getY(), up.getZ());
        translation.scale(distance);
        position = Math::add(position, translation);
    }

    Vector& getPosition() {
        return position;
    }
    Vector& getImagePlane() {
        return imagePlane;
    }

    void setPosition(const Vector position) {
        this->position = position;
    }
    void setImagePlane(const Vector imagePlane) {
        this->imagePlane = imagePlane;
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
        //Triangle(Vector(-1.75, -1.75, -3), Vector(1.75, -1.75, -3), Vector(0, 1.75, -3))
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

    std::ofstream ppmFileStream("task3_roll_after_1920_1080.ppm", 
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
                Camera camera(Vector(0, 0, 0), Vector(x, y, -1));

                // camera.dolly(3);
                // camera.truck(3);
                // camera.pan(16);
                // camera.pedestal(6);
                // camera.tilt(-16);
                // camera.truck(3);
                // camera.pan(30);
                camera.roll(16);

                if(Math::dotProduct(camera.getImagePlane(), 
                                        triangles[trIndex].getNormal()) != 0 
                    && Math::dotProduct(triangles[trIndex].getNormal(), 
                                            triangles[trIndex].getVertex0()) < 0)
                {
                    float scaleFactor = Math::dotProduct(
                        Math::subtract(
                            triangles[trIndex].getVertex0(), 
                            camera.getPosition()
                        ), 
                        triangles[trIndex].getNormal()
                    ) / Math::dotProduct(camera.getImagePlane(), 
                                            triangles[trIndex].getNormal());
                    camera.getImagePlane().scale(scaleFactor);
                    camera.setImagePlane(
                        Math::add(
                            camera.getPosition(), 
                            camera.getImagePlane()
                        )
                    );
                    bool isIn = true;
                    for(unsigned int i = 0; i < NUM_OF_TRIANGLE_VERTICES; i++) {
                        unsigned int tail = i;
                        unsigned int head = (i + 1) % NUM_OF_TRIANGLE_VERTICES;
                        Vector temp = Math::crossProduct(
                            Math::subtract(triangles[trIndex].getVertex(head), 
                                            triangles[trIndex].getVertex(tail)), 
                            Math::subtract(camera.getImagePlane(), 
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

    Matrix a(std::begin({1.0f, 2.0f, 3.0f, 2.0f, 3.0f, 4.0f, 3.0f, 4.0f, 5.0f}), 3, 3);
    a.print();
    Matrix b(std::begin({1.0f, 2.0f, 1.0f, 1.0f, 1.0f, 4.0f, 1.0f, 1.0f, 5.0f}), 3, 3);
    b.print();
    Matrix result = Math::matrixMultiply(a, b);
    result.print();

    Vector test(1.0, 2.0, 3.0);
    Vector test1 = Math::matrixMultiply(test, a);
    test1.print();
    return 0;
}