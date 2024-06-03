#include <fstream>

// Output image resolution
static const int imageWidth = 1920;
static const int imageHeight = 1080;
static const int maxColorComponent = 255;

#define NUMBER_OF_HORIZONTAL_DIVISIONS 6
#define NUMBER_OF_VERTICAL_DIVISIONS 6

int main() {
    std::ofstream ppmFileStream("crt_output_image_squares.ppm", 
                    std::ios::out | std::ios::binary);
    ppmFileStream << "P3\n";
    ppmFileStream << imageWidth << " " << imageHeight << "\n";
    ppmFileStream << maxColorComponent << "\n";
    for (int rowIdx = 0; rowIdx < imageHeight; ++rowIdx) {
        for (int colIdx = 0; colIdx < imageWidth; ++colIdx) {
            ppmFileStream << (rowIdx / (imageHeight / NUMBER_OF_HORIZONTAL_DIVISIONS) + 1) 
                                * (255 / NUMBER_OF_HORIZONTAL_DIVISIONS) 
            << " " << (colIdx / (imageWidth / NUMBER_OF_VERTICAL_DIVISIONS) + 1) 
                                * (255 / NUMBER_OF_VERTICAL_DIVISIONS)  
            << " 255\t";
        }
        ppmFileStream << "\n";
    }
    ppmFileStream.close();
    return 0;
}