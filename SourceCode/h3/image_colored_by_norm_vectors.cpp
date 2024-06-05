#include <fstream>
#include <cmath>

// Output image resolution
static const int imageWidth = 6000;
static const int imageHeight = 2000;
static const int maxColorComponent = 255;

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
        x = x / length();
        y = y / length();
        z = z / length();
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
};

struct Ray {
    Vector direction {0, 0, -1};
    Vector origin {0, 0, 0};
};

int main() {
    std::ofstream ppmFileStream("norm_vectors_image_6000_2000.ppm", 
                    std::ios::out | std::ios::binary);
    ppmFileStream << "P3\n";
    ppmFileStream << imageWidth << " " << imageHeight << "\n";
    ppmFileStream << maxColorComponent << "\n";
    for (int rowIdx = 0; rowIdx < imageHeight; ++rowIdx) {
        for (int colIdx = 0; colIdx < imageWidth; ++colIdx) {
            struct Ray ray;
            ray.direction.setX((colIdx + 0.5) / imageWidth);
            ray.direction.setY((rowIdx + 0.5) / imageHeight);

            ray.direction.setX(ray.direction.getX() * 2 - 1);
            ray.direction.setY(1 - ray.direction.getY() * 2);

            unsigned int aspectRatio = imageWidth / imageHeight;
            ray.direction.setX(ray.direction.getX() * aspectRatio);
            ray.direction.setY(ray.direction.getY() * aspectRatio);

            ray.direction.normalize();

            ppmFileStream << (int)abs(ray.direction.getX() * 1000) % maxColorComponent
            << " " << (int)abs(ray.direction.getY() * 1000) % maxColorComponent
            << " " << (int)abs(ray.direction.getZ() * 1000) % maxColorComponent
            << "\t";
        }
        ppmFileStream << "\n";
    }
    ppmFileStream.close();
    return 0;
}