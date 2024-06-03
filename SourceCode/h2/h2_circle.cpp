#include <fstream>
#include <cmath>

// Output image resolution
static const int imageWidth = 1920;
static const int imageHeight = 1080;
static const int maxColorComponent = 255;

#define CIRCLE_RADIUS_PX 333
#define COLOR_DEFORMATION 255

int main() {
    std::ofstream ppmFileStream("crt_output_image_circle.ppm", 
                    std::ios::out | std::ios::binary);
    ppmFileStream << "P3\n";
    ppmFileStream << imageWidth << " " << imageHeight << "\n";
    ppmFileStream << maxColorComponent << "\n";
    for (int rowIdx = 0; rowIdx < imageHeight; ++rowIdx) {
        for (int colIdx = 0; colIdx < imageWidth; ++colIdx) {
            ppmFileStream << "0 " 
            << ((int)sqrt(pow(imageHeight/2 - rowIdx, 2)
            + pow(imageWidth/2 - colIdx, 2)) <= CIRCLE_RADIUS_PX)
            * COLOR_DEFORMATION
            << " 255\t";
        }
        ppmFileStream << "\n";
    }
    ppmFileStream.close();
    return 0;
}