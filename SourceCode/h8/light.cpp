#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iterator>
#include <vector>

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

    void scale(const float num) {
        for(unsigned int i = 0; i < rows; i++) {
            for(unsigned int j = 0; j < columns; j++) {
                fields[i][j] *= num;
            }
        }
    }
    void bound(const float num) {
        for(unsigned int i = 0; i < rows; i++) {
            for(unsigned int j = 0; j < columns; j++) {
                if(fields[i][j] > num) fields[i][j] = num;
            }
        }
    }
    void transpose() {
        for(unsigned int i = 0; i < rows-1; i++) {
            for(unsigned int j = i+1; j < columns; j++) {
                const float temp = fields[i][j];
                fields[i][j] = fields[j][i];
                fields[j][i] = temp;
            }
        }
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
    Vector(const float n) 
        : Matrix(std::begin({n, n, n}), 1, 3) 
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
    Vector &operator=(const float n) {
        for(unsigned int i = 0; i < 3; i++) {
            set(i, n);
        }
        return *this;
    }

    Vector operator*(const Vector& v) {
        Vector result(0, 0, 0);
        for(unsigned int i = 0; i < 3; i++) {
            result.set(i, get(i) * v.get(i));
        }
        return result;
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

    void applyRotation(Matrix& transformation) {
        imagePlane = Math::matrixMultiply(imagePlane, transformation);
        direction = Math::matrixMultiply(direction, transformation);
        up = Math::matrixMultiply(up, transformation);
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

#include "../include/rapidjson/document.h"
#include "../include/rapidjson/istreamwrapper.h"
#include "../include/rapidjson/rapidjson.h"

class Mesh {
    std::vector<Vector> vertices;
    std::vector<Triangle> triangles;
public:
    Mesh() = default;

    std::vector<Vector>& getVertices() {
        return vertices;
    }
    std::vector<Triangle>& getTriangles() {
        return triangles;
    }
};

class Light {
    Vector position;
    unsigned int intensity;
public:
    Light(const Vector position, const unsigned int intensity) 
        : position(position), intensity(intensity) 
        {}
    
    Vector getPosition() const {
        return position;
    }
    unsigned int getIntensity() const {
        return intensity;
    }
};

class Scene {
    std::string fileName;
    struct {
        Vector bgColor {0, 0, 0};
        struct {
            unsigned int width;
            unsigned int height;
        } imageSettings;
    } settings;
    Vector cameraPosition {0, 0, 0};
    Matrix initialTransformation {3, 3};
    std::vector<Mesh> meshes;
    std::vector<Light> lights;

    rapidjson::Document getJsonDocument() const {
        std::ifstream ifile(fileName);
        assert(ifile.is_open());

        rapidjson::IStreamWrapper iwrap(ifile);
        rapidjson::Document document;
        document.ParseStream(iwrap);

        if(document.HasParseError()) {
            std::cout << "Error: " << document.GetParseError() << std::endl;
            std::cout << "Offset: " << document.GetErrorOffset() << std::endl;
            assert(false);
        }

        assert(document.IsObject());
        return document;
    }
    Vector loadVector(const rapidjson::Value::ConstArray& arr) const {
        assert(arr.Size() == 3);
        return Vector(
            static_cast<float>(arr[0].GetDouble()), 
            static_cast<float>(arr[1].GetDouble()), 
            static_cast<float>(arr[2].GetDouble())
        );
    }
    Matrix loadMatrix(const rapidjson::Value::ConstArray& arr) const {
        assert(arr.Size() == 9);
        Matrix result(3, 3);
        for(unsigned int i = 0; i < result.getRows(); i++) {
            for(unsigned int j = 0; j < result.getColumns(); j++) {
                result.setField(
                    i, 
                    j, 
                    static_cast<float>(arr[j + result.getColumns() * i].GetDouble())
                );
            }
        }
        return result;
    }
    void loadVertices(const rapidjson::Value::ConstArray& arr) {
        assert(arr.Size() % 3 == 0);
        meshes.emplace_back(Mesh());
        for(unsigned int i = 0; i < arr.Size(); i+=3) {
            meshes.back().getVertices().emplace_back(
                Vector(
                    static_cast<float>(arr[i].GetDouble()), 
                    static_cast<float>(arr[i+1].GetDouble()), 
                    static_cast<float>(arr[i+2].GetDouble())
                )
            );
        }
    }
    void loadTriangles(const rapidjson::Value::ConstArray& arr) {
        assert(arr.Size() % 3 == 0);
        for(unsigned int i = 0; i < arr.Size(); i+=3) {
            meshes.back().getTriangles().emplace_back(
                Triangle(
                    meshes.back().getVertices().at(arr[i].GetInt()), 
                    meshes.back().getVertices().at(arr[i+1].GetInt()), 
                    meshes.back().getVertices().at(arr[i+2].GetInt())
                )
            );
            meshes.back().getTriangles().back().cacheNormal();
        }
    }

public:
    Scene(const std::string fileName) : fileName(fileName) {}

    void parseSceneFile() {
        rapidjson::Document document = getJsonDocument();

        const rapidjson::Value& settingsValue = document.FindMember("settings")->value;
        if(!settingsValue.IsNull() && settingsValue.IsObject()) {
            const rapidjson::Value& bgColorValue 
                        = settingsValue.FindMember("background_color")->value;
            assert(!bgColorValue.IsNull() && bgColorValue.IsArray());
            settings.bgColor = loadVector(bgColorValue.GetArray());

            const rapidjson::Value& imageSettingsValue 
                        = settingsValue.FindMember("image_settings")->value;
            if(!imageSettingsValue.IsNull() && imageSettingsValue.IsObject()) {
                const rapidjson::Value& imageWidthValue 
                            = imageSettingsValue.FindMember("width")->value;
                const rapidjson::Value& imageHeightValue 
                            = imageSettingsValue.FindMember("height")->value;
                assert(!imageWidthValue.IsNull() && imageWidthValue.IsInt() 
                            && !imageHeightValue.IsNull() && imageHeightValue.IsInt());
                settings.imageSettings.width = imageWidthValue.GetInt();
                settings.imageSettings.height = imageHeightValue.GetInt();
            }
        }

        const rapidjson::Value& cameraValue = document.FindMember("camera")->value;
        if(!cameraValue.IsNull() && cameraValue.IsObject()) {
            const rapidjson::Value& matrixValue = cameraValue.FindMember("matrix")->value;
            assert(!matrixValue.IsNull() && matrixValue.IsArray());
            initialTransformation = loadMatrix(matrixValue.GetArray());

            const rapidjson::Value& positionValue = cameraValue.FindMember("position")->value;
            assert(!positionValue.IsNull() && positionValue.IsArray());
            cameraPosition = loadVector(positionValue.GetArray());
        }

        const rapidjson::Value& ligthsValue = document.FindMember("lights")->value;
        if(!ligthsValue.IsNull() && ligthsValue.IsArray()) {
            for(rapidjson::Value::ConstValueIterator itr = ligthsValue.Begin(); 
                    itr != ligthsValue.End(); 
                    itr++)
            {
                const rapidjson::Value& intensityValue 
                        = itr->FindMember("intensity")->value;
                assert(!intensityValue.IsNull() && intensityValue.IsInt());

                const rapidjson::Value& positionValue 
                        = itr->FindMember("position")->value;
                assert(!positionValue.IsNull() && positionValue.IsArray());
                lights.emplace_back(
                    loadVector(positionValue.GetArray()), 
                    intensityValue.GetInt()
                );
            }
        }

        const rapidjson::Value& objectsValue = document.FindMember("objects")->value;
        if(!objectsValue.IsNull() && objectsValue.IsArray()) {
            for(rapidjson::Value::ConstValueIterator itr = objectsValue.Begin(); 
                    itr != objectsValue.End(); 
                    itr++)
            {
                const rapidjson::Value& verticesValue 
                            = itr->FindMember("vertices")->value;
                assert(!verticesValue.IsNull() && verticesValue.IsArray());
                loadVertices(verticesValue.GetArray());

                const rapidjson::Value& trianglesValue 
                            = itr->FindMember("triangles")->value;
                assert(!trianglesValue.IsNull() && trianglesValue.IsArray());
                loadTriangles(trianglesValue.GetArray());
            }
        }
    }

    std::string getFileName() const {
        return fileName;
    }
    unsigned int getImageWidth() const {
        return settings.imageSettings.width;
    }
    unsigned int getImageHeight() const {
        return settings.imageSettings.height;
    }
    Vector getBgColor() const {
        return settings.bgColor;
    }
    std::vector<Mesh> getMeshes() const {
        return meshes;
    }
    Vector getCameraPosition() const {
        return cameraPosition;
    }
    Matrix getInitialTransformation() const {
        return initialTransformation;
    }
    std::vector<Light> getLights() const {
        return lights;
    }

    void setFileName(const std::string fileName) {
        this->fileName = fileName;
    }
};

class Renderer {
    struct {
        Vector bgColor;
        struct {
            unsigned int width;
            unsigned int height;
            unsigned int maxColorComponent;
        } imageSettings;
    } settings;
    std::vector<Mesh> meshes;
    Vector cameraPosition;
    Matrix initialTransformation;
    std::vector<Light> lights;

    bool checkForIntersections(const Vector start, const Vector end) {
        //ray.normalize();
        for(auto& mesh : meshes) {
            for(auto& triangle : mesh.getTriangles()) {
                Vector ray = Math::subtract(end, start);
                if(Math::dotProduct(ray, triangle.getNormal()) != 0 
                        && Math::dotProduct(Math::subtract(triangle.getVertex0(), 
                                                                start), 
                                                            triangle.getNormal()) >= 0) {
                    float scaleFactor = Math::dotProduct(
                        Math::subtract(
                            triangle.getVertex0(), 
                            start
                        ), 
                        triangle.getNormal()
                    ) / Math::dotProduct(ray, triangle.getNormal());
                    ray.scale(scaleFactor);
                    ray = Math::add(start, ray);
                    bool isIn = true;
                    for(unsigned int i = 0; i < NUM_OF_TRIANGLE_VERTICES; i++) {
                        unsigned int tail = i;
                        unsigned int head = (i + 1) % NUM_OF_TRIANGLE_VERTICES;
                        Vector temp = Math::crossProduct(
                            Math::subtract(triangle.getVertex(head), 
                                            triangle.getVertex(tail)), 
                            Math::subtract(ray, 
                                            triangle.getVertex(tail))
                        );
                        float checkSum = Math::dotProduct(
                                            triangle.getNormal(), 
                                            temp
                                         );
                        if(checkSum <= 0) {
                            isIn = false;
                            break;
                        }
                    }
                    if(isIn) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
    void shade(Triangle& triangle, Camera& camera, Vector& colors) {
        colors = 0;
        Vector albedo(1, 1, 1);
        for(auto& light : lights) {
            bool hasIntersections = checkForIntersections(
                                        Math::add(
                                            camera.getImagePlane(), 
                                            triangle.getNormal() * Vector(0.36)
                                        ), 
                                        light.getPosition()
                                    );
            if(hasIntersections == true) continue;
            Vector lightDirection = Math::subtract(
                                        light.getPosition(), 
                                        camera.getImagePlane()
                                    );
            float lightRadius = lightDirection.length();
            lightDirection.normalize();
            float attenuation = 4 * M_PI * lightRadius * lightRadius;
            Vector temp = std::fmax(
                0, 
                Math::dotProduct(
                    triangle.getNormal(), 
                    lightDirection
                ) * 255 * (light.getIntensity() / attenuation)
            );
            temp = temp * albedo;
            colors = Math::add(colors, temp);
        }
        colors.bound(255);
    }
public:
    Renderer(const Scene scene) 
        : settings  {
                        scene.getBgColor(), 
                        {
                            scene.getImageWidth(),
                            scene.getImageHeight(),
                            255
                        }
                    }, 
        meshes(scene.getMeshes()), 
        cameraPosition(scene.getCameraPosition()), 
        initialTransformation(scene.getInitialTransformation()), 
        lights(scene.getLights()) 
        {}

    void render(const std::string fileName) {
        std::ofstream ppmFileStream(fileName, 
                        std::ios::out | std::ios::binary);
        ppmFileStream << "P3\n";
        ppmFileStream << settings.imageSettings.width 
        << " " << settings.imageSettings.height << "\n";
        ppmFileStream << settings.imageSettings.maxColorComponent << "\n";
        for (int rowIdx = 0; rowIdx < settings.imageSettings.height; ++rowIdx) {
            for (int colIdx = 0; colIdx < settings.imageSettings.width; ++colIdx) {
                float x = colIdx + 0.5, y = rowIdx + 0.5;
                x /= settings.imageSettings.width;
                y /= settings.imageSettings.height;
                x = (2 * x) - 1;
                y = 1 - (2 * y);
                x *= (float) settings.imageSettings.width / settings.imageSettings.height;
                struct {
                    Vector colors {0, 0, 0};
                    unsigned int zValue = 0;
                } pixel;
                pixel.colors = {
                    (float)(int)
                    (settings.bgColor.getX() * settings.imageSettings.maxColorComponent), 
                    (float)(int)
                    (settings.bgColor.getY() * settings.imageSettings.maxColorComponent), 
                    (float)(int)
                    (settings.bgColor.getZ() * settings.imageSettings.maxColorComponent)
                };
                for(auto &mesh : meshes) {
                    std::vector<Triangle>& triangles = mesh.getTriangles();
                    const unsigned int trCount = triangles.size();
                    for(unsigned int trIndex = 0; trIndex < trCount; trIndex++) {
                        Camera camera(cameraPosition, Vector(x, y, -1));
                        camera.applyRotation(initialTransformation);

                        if(Math::dotProduct(camera.getImagePlane(), 
                                                triangles[trIndex].getNormal()) != 0 
                            && Math::dotProduct(triangles[trIndex].getNormal(), 
                                                    Math::subtract(
                                                        triangles[trIndex].getVertex0(), 
                                                        camera.getPosition()
                                                    )) < 0)
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
                                float checkSum = Math::dotProduct(
                                                    triangles[trIndex].getNormal(), 
                                                    temp
                                                 );
                                if(checkSum <= 0) {
                                    isIn = false;
                                    break;
                                }
                            }
                            if(isIn && (pixel.zValue == 0 
                                            || pixel.zValue > scaleFactor))
                            {
                                shade(triangles[trIndex], camera, pixel.colors);
                                pixel.zValue = scaleFactor;
                            }
                        }
                    }
                }
                ppmFileStream << (int)pixel.colors.getX() 
                << " " << (int)pixel.colors.getY() 
                << " " << (int)pixel.colors.getZ() 
                << "\t";
            }
            ppmFileStream << "\n";
        }
        ppmFileStream.close();
    }
};

int main() {
    Scene scene("scene3.crtscene");
    scene.parseSceneFile();
    Renderer renderer(scene);
    renderer.render("scene_3.ppm");
    return 0;
}