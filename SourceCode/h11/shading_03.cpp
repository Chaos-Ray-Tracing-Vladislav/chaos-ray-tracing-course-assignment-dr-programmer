#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iterator>
#include <vector>

constexpr unsigned int MAX_MATRIX_ROWS = 4;
constexpr unsigned int MAX_MATRIX_COLUMNS = 4;

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
    Vector* normal = nullptr;
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

    Vector* getNormal() const {
        return normal;
    }
    void setNormal(Vector * const normal) {
        this->normal = normal;
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

constexpr unsigned int NUM_OF_TRIANGLE_VERTICES = 3;

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
    unsigned int materialIndex;
public:
    Mesh() = default;
    Mesh(const unsigned int materialIndex) : materialIndex(materialIndex) {}

    std::vector<Vector>& getVertices() {
        return vertices;
    }
    std::vector<Triangle>& getTriangles() {
        return triangles;
    }
    unsigned int getMaterialIndex() const {
        return materialIndex;
    }
};

class Light {
    Vector position;
    unsigned int intensity;
public:
    Light(const Vector position, const unsigned int intensity) 
        : position(position), intensity(intensity) 
        {}
    
    Vector& getPosition() {
        return position;
    }
    unsigned int getIntensity() const {
        return intensity;
    }
};

class Material {
    std::string type;
    Vector albedo {1};
    float ior = 1;
    bool smoothShading;
public:
    Material(const std::string& type, const Vector& albedo, const bool smoothShading) 
        : type(type), albedo(albedo), smoothShading(smoothShading) 
        {}
    Material(const std::string& type, const float ior, const bool smoothShading) 
        : type(type), ior(ior), smoothShading(smoothShading) 
        {}
    
    std::string getType() const {
        return type;
    }
    Vector& getAlbedo() {
        return albedo;
    }
    float getIOR() const {
        return ior;
    }
    bool hasSmoothShading() const {
        return smoothShading;
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
    std::vector<Material> materials;

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
        for(unsigned int i = 0; i < arr.Size(); i+=3) {
            meshes.back().getVertices().emplace_back(
                Vector(
                    static_cast<float>(arr[i].GetDouble()), 
                    static_cast<float>(arr[i+1].GetDouble()), 
                    static_cast<float>(arr[i+2].GetDouble())
                )
            );
            meshes.back().getVertices().back().setNormal(new Vector(0));
        }
    }
    void loadTriangles(const rapidjson::Value::ConstArray& arr) {
        assert(arr.Size() % 3 == 0);
        for(unsigned int i = 0; i < arr.Size(); i+=3) {
            Vector& one = meshes.back().getVertices().at(arr[i].GetInt());
            Vector& two = meshes.back().getVertices().at(arr[i+1].GetInt());
            Vector& three = meshes.back().getVertices().at(arr[i+2].GetInt());

            Vector trNormal = Math::crossProduct(
                Math::subtract(two, one), 
                Math::subtract(three, one)
            );

            *one.getNormal() = Math::add(*one.getNormal(), trNormal);
            *two.getNormal() = Math::add(*two.getNormal(), trNormal);
            *three.getNormal() = Math::add(*three.getNormal(), trNormal);

            trNormal.normalize();

            meshes.back().getTriangles().emplace_back(
                Triangle(
                    one, 
                    two, 
                    three
                )
            );
            meshes.back().getTriangles().back().cacheNormal(trNormal);
        }
    }
    void normalizeVertexNormals() {
        for(auto& vertex : meshes.back().getVertices()) {
            vertex.getNormal()->normalize();
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

        const rapidjson::Value& materialsValue = document.FindMember("materials")->value;
        if(!materialsValue.IsNull() && materialsValue.IsArray()) {
            for(rapidjson::Value::ConstValueIterator itr = materialsValue.Begin(); 
                    itr != materialsValue.End(); 
                    itr++)
            {
                const rapidjson::Value& typeValue = itr->FindMember("type")->value;
                assert(!typeValue.IsNull() && typeValue.IsString());

                Vector albedo(1);
                float ior = 0;
                const std::string type = typeValue.GetString();
                if(type == "refractive") {
                    const rapidjson::Value& iorValue = itr->FindMember("ior")->value;
                    assert(!iorValue.IsNull() && iorValue.IsDouble());
                    ior = static_cast<float>(iorValue.GetDouble());
                }
                else {
                    const rapidjson::Value& albedoValue = itr->FindMember("albedo")->value;
                    assert(!albedoValue.IsNull() && albedoValue.IsArray());
                    albedo = loadVector(albedoValue.GetArray());
                }

                const rapidjson::Value& smoothValue 
                        = itr->FindMember("smooth_shading")->value;
                assert(!smoothValue.IsNull() && smoothValue.IsBool());

                if(type == "refractive") {
                    materials.emplace_back(
                        type, 
                        ior, 
                        smoothValue.GetBool()
                    );
                }
                else {
                    materials.emplace_back(
                        type, 
                        albedo, 
                        smoothValue.GetBool()
                    );
                }
            }
        }

        const rapidjson::Value& objectsValue = document.FindMember("objects")->value;
        if(!objectsValue.IsNull() && objectsValue.IsArray()) {
            for(rapidjson::Value::ConstValueIterator itr = objectsValue.Begin(); 
                    itr != objectsValue.End(); 
                    itr++)
            {
                const rapidjson::Value& materialIndexValue 
                        = itr->FindMember("material_index")->value;
                assert(!materialIndexValue.IsNull() && materialIndexValue.IsInt());

                meshes.emplace_back(Mesh(materialIndexValue.GetInt()));

                const rapidjson::Value& verticesValue 
                            = itr->FindMember("vertices")->value;
                assert(!verticesValue.IsNull() && verticesValue.IsArray());
                loadVertices(verticesValue.GetArray());

                const rapidjson::Value& trianglesValue 
                            = itr->FindMember("triangles")->value;
                assert(!trianglesValue.IsNull() && trianglesValue.IsArray());
                loadTriangles(trianglesValue.GetArray());

                normalizeVertexNormals();
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
    Vector& getBgColor() {
        return settings.bgColor;
    }
    std::vector<Mesh>& getMeshes() {
        return meshes;
    }
    Vector& getCameraPosition() {
        return cameraPosition;
    }
    Matrix& getInitialTransformation() {
        return initialTransformation;
    }
    std::vector<Light>& getLights() {
        return lights;
    }
    std::vector<Material>& getMaterials() {
        return materials;
    }

    void setFileName(const std::string fileName) {
        this->fileName = fileName;
    }
};

constexpr float rayStartOffset = 0.3;
constexpr float rayStartOffsetReflective = 0.3;
constexpr float rayStartOffsetRefractive = 0.3;
constexpr unsigned int maxRayDepth = 6;

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
    std::vector<Material> materials;

    Vector getBarycentricCoords(Triangle& triangle, Vector& intersection) {
        float triangleArea = Math::crossProduct(
                                Math::subtract(
                                    triangle.getVertex1(), 
                                    triangle.getVertex0()
                                ), 
                                Math::subtract(
                                    triangle.getVertex2(), 
                                    triangle.getVertex0()
                                )
                             ).length();
        float u = (
            Math::crossProduct(
                Math::subtract(
                    intersection, 
                    triangle.getVertex0()
                ), 
                Math::subtract(
                    triangle.getVertex2(), 
                    triangle.getVertex0()
                )
            ).length() 
            / 
            triangleArea
        );
        float v = (
            Math::crossProduct(
                Math::subtract(
                    triangle.getVertex1(), 
                    triangle.getVertex0()
                ), 
                Math::subtract(
                    intersection, 
                    triangle.getVertex0()
                )
            ).length() 
            / 
            triangleArea
        );
        float w = 1 - u - v;
        return Vector(u, v, w);
    }
    bool checkRayForIntersections(const Vector start, 
                                    const Vector originalRay, 
                                    const bool applyShade = false, 
                                    Vector* colors = nullptr, 
                                    const unsigned int rayDepth = 0) 
    {
        float zValue = 0;
        unsigned int nearestTriangleIndex = 0;
        unsigned int nearestMeshIndex = 0;

        unsigned int meshIndex = 0;
        for(auto& mesh : meshes) {
            unsigned int triangleIndex = 0;
            for(auto& triangle : mesh.getTriangles()) {
                Vector ray = originalRay;
                if(Math::dotProduct(ray, triangle.getNormal()) != 0) {
                    float scaleFactor = Math::dotProduct(
                        Math::subtract(
                            triangle.getVertex0(), 
                            start
                        ), 
                        triangle.getNormal()
                    ) / Math::dotProduct(ray, triangle.getNormal());
                    if(scaleFactor <= 0 || (scaleFactor > 1 && applyShade == false)) {
                        triangleIndex++;
                        continue;
                    }
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
                        if(checkSum < 0) {
                            isIn = false;
                            break;
                        }
                    }
                    if(isIn) {
                        if(applyShade == true && (zValue == 0 || zValue > scaleFactor)) {
                            zValue = scaleFactor;
                            nearestTriangleIndex = triangleIndex;
                            nearestMeshIndex = meshIndex;
                        }
                        else if(applyShade == false 
                                && materials.at(mesh.getMaterialIndex()).getType() 
                                    != "refractive") 
                            return true;
                    }
                }
                triangleIndex++;
            }
            meshIndex++;
        }
        if(applyShade == true && zValue != 0) {
            Vector intersection = originalRay;
            intersection.scale(zValue);
            intersection = Math::add(start, intersection);
            *colors = shade(meshes[nearestMeshIndex].getTriangles()[nearestTriangleIndex], 
                                meshes[nearestMeshIndex], 
                                intersection, 
                                rayDepth);
            return true;
        }
        return false;
    }
    bool checkForIntersections(const Vector start, const Vector end) {
        return checkRayForIntersections(start, Math::subtract(end, start));
    }

    Vector reflect(Vector& intersection, Vector& normal, Vector& albedo, 
                        const unsigned int rayDepth) 
    {
        Vector colors(0);

        float rayProjectionOnNormal = Math::dotProduct(intersection, 
                                                            normal);
        Vector newRay = Math::subtract(
            intersection, 
            Vector(2) * normal * Vector(rayProjectionOnNormal)
        );
        bool hasIntersection = checkRayForIntersections(
            Math::add(
                intersection, 
                normal * Vector(rayStartOffsetReflective)
            ), 
            newRay, 
            true, 
            &colors, 
            rayDepth + 1
        );
        if(hasIntersection == false) return settings.bgColor * albedo;
        return colors * albedo;
    }
    Vector shade(Triangle& triangle, Mesh& mesh, Vector& intersection, 
                    const unsigned int rayDepth = 0) 
    {
        Vector colors(0);
        Vector albedo = materials.at(mesh.getMaterialIndex()).getAlbedo();
        Vector normal(0);

        if(materials.at(mesh.getMaterialIndex()).hasSmoothShading() == true) {
            Vector barCoords = getBarycentricCoords(triangle, intersection);

            Vector v0Normal = *triangle.getVertex0().getNormal();
            Vector v1Normal = *triangle.getVertex1().getNormal();
            Vector v2Normal = *triangle.getVertex2().getNormal();
            v1Normal.scale(barCoords.getX());
            v2Normal.scale(barCoords.getY());
            v0Normal.scale(barCoords.getZ());

            normal = Math::add(
                Math::add(
                    v1Normal, 
                    v2Normal
                ), 
                v0Normal
            );
        }
        else normal = triangle.getNormal();

        if(materials.at(mesh.getMaterialIndex()).getType() == "reflective") {
            if(rayDepth > maxRayDepth) return Vector(0);

            colors = reflect(intersection, normal, albedo, rayDepth);
            return colors;
        }
        if(materials.at(mesh.getMaterialIndex()).getType() == "refractive") {
            if(rayDepth > maxRayDepth) return Vector(0);

            Vector reflectionColor(0);
            Vector refractionColor(0);
            Vector incident = intersection;
            incident.normalize();
            float n1 = 1, n2 = 1;
            if(Math::dotProduct(incident, normal) > 0) {
                n1 = materials.at(mesh.getMaterialIndex()).getIOR();
                normal = normal * Vector(-1);
            }
            else n2 = materials.at(mesh.getMaterialIndex()).getIOR();
            std::cout << n1 << " " << n2;

            float cosAlpha = -Math::dotProduct(incident, normal);
            float sinAlpha = sqrt(1 - (cosAlpha * cosAlpha));
            if(sinAlpha < n2 / n1) {
                float sinBeta = (sinAlpha * n1) / n2;
                Vector A = Vector(sqrt(1 - (sinBeta * sinBeta))) * normal * Vector(-1);
                Vector C = Math::add(incident, Vector(cosAlpha) * normal);
                C.normalize();
                Vector B = C * sinBeta;
                Vector refracted = Math::add(A, B);
                refracted = Math::add(
                    intersection, 
                    refracted
                );

                bool hit = checkRayForIntersections(
                    Math::add(
                        intersection, 
                        Vector(-1) * normal * Vector(rayStartOffsetRefractive)
                    ), 
                    refracted, 
                    true, 
                    &refractionColor, 
                    rayDepth + 1
                );
                if(hit == false) refractionColor = settings.bgColor;
            }
            else return reflect(intersection, normal, albedo, rayDepth);

            reflectionColor = reflect(intersection, normal, albedo, rayDepth);

            const float fresnel = 0.5 * pow((1.0 + -cosAlpha), 5);
            colors = Math::add(
                Vector(fresnel) * reflectionColor, 
                Vector(1 - fresnel) * refractionColor
            );

            return colors;
        }

        for(auto& light : lights) {
            bool hasIntersections = checkForIntersections(
                                        Math::add(
                                            intersection, 
                                            normal * Vector(rayStartOffset)
                                        ), 
                                        light.getPosition()
                                    );
            if(hasIntersections == true) continue;
            Vector lightDirection = Math::subtract(
                                        light.getPosition(), 
                                        intersection
                                    );
            float lightRadius = lightDirection.length();
            lightDirection.normalize();
            float attenuation = 4 * M_PI * lightRadius * lightRadius;
            Vector temp = std::fmax(
                0, 
                Math::dotProduct(
                    normal, 
                    lightDirection
                ) * 255 * (light.getIntensity() / attenuation)
            );
            temp = temp * albedo;
            colors = Math::add(colors, temp);
        }
        colors.bound(255);
        return colors;
    }
public:
    Renderer(Scene scene) 
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
        lights(scene.getLights()), 
        materials(scene.getMaterials()) 
        {}

    void render(const std::string fileName) {

        settings.bgColor = {
            (float)(int)
            (settings.bgColor.getX() * settings.imageSettings.maxColorComponent), 
            (float)(int)
            (settings.bgColor.getY() * settings.imageSettings.maxColorComponent), 
            (float)(int)
            (settings.bgColor.getZ() * settings.imageSettings.maxColorComponent)
        };

        std::ofstream ppmFileStream(fileName, 
                        std::ios::out | std::ios::binary);
        ppmFileStream << "P3\n";
        ppmFileStream << settings.imageSettings.width 
        << " " << settings.imageSettings.height << "\n";
        ppmFileStream << settings.imageSettings.maxColorComponent << "\n";
        for (int rowIdx = 0; rowIdx < settings.imageSettings.height; ++rowIdx) {
            for (int colIdx = 0; colIdx < settings.imageSettings.width; ++colIdx) {
                //system("clear");
                std::cout << "row: " << rowIdx << " column: " << colIdx << std::endl;
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
                pixel.colors = settings.bgColor;

                Camera camera(cameraPosition, Vector(x, y, -1));
                camera.applyRotation(initialTransformation);
                checkRayForIntersections(camera.getPosition(), 
                                            camera.getImagePlane(), 
                                            true, 
                                            &pixel.colors);

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

int main() {/*
    for(unsigned int i = 0; i < 9; i++) {
        std::string ifile = "scene" + std::to_string(i) + ".crtscene";
        std::string ofile = "scene_" + std::to_string(i) + ".ppm";

        Scene scene(ifile);
        scene.parseSceneFile();
        Renderer renderer(scene);
        renderer.render(ofile);
    }*/

    Scene scene("scene8.crtscene");
    scene.parseSceneFile();
    Renderer renderer(scene);
    renderer.render("scene_8.ppm");

    return 0;
}