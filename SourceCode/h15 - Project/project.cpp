#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iterator>
#include <vector>
#include <map>
#include <memory>
#include <chrono>
#include <thread>
#include <mutex>
#include <list>
#include <limits>
#include <stack>

std::mutex mtx;

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
    std::shared_ptr<Vector> normal = nullptr;
    std::shared_ptr<unsigned int> index = nullptr;
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

    std::shared_ptr<Vector> getNormal() const {
        return normal;
    }
    void setNormal(std::shared_ptr<Vector> normal) {
        this->normal = normal;
    }
    std::shared_ptr<unsigned int> getIndex() const {
        return index;
    }
    void setIndex(std::shared_ptr<unsigned int> index) {
        this->index = index;
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
    std::vector<Vector> uvs;
    std::vector<Triangle> triangles;
    unsigned int materialIndex;
public:
    Mesh() = default;
    Mesh(const unsigned int materialIndex) : materialIndex(materialIndex) {}

    std::vector<Vector>& getUVs() {
        return uvs;
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
    std::string albedoTexture {1};
    bool smoothShading;
public:
    Material(const std::string& type, const Vector& albedo, const bool smoothShading) 
        : type(type), albedo(albedo), smoothShading(smoothShading) 
        {}
    Material(const std::string& type, const float ior, const bool smoothShading) 
        : type(type), ior(ior), smoothShading(smoothShading) 
        {}
    Material(const std::string& type, 
                const std::string albedoTexture, 
                const bool smoothShading) 
        : type(type), albedoTexture(albedoTexture), smoothShading(smoothShading) 
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
    std::string getAlbedoTexture() const {
        return albedoTexture;
    }
    bool hasSmoothShading() const {
        return smoothShading;
    }
};

#define STB_IMAGE_IMPLEMENTATION
#include "../stb_image.h"

class Texture {
    std::string name;
    std::string type;
public:
    Texture() = default;
    Texture(const std::string name, const std::string type) 
        : name(name), type(type) 
        {}
    
    virtual Vector getColor(const Vector& barCoords, Triangle& triangle, Mesh& mesh) const {
        return Vector(1);
    }

    std::string getName() const {
        return name;
    }
    std::string getType() const {
        return type;
    }
};

class TextureColor : public Texture {
    Vector albedo;
public:
    TextureColor(const std::string name, const std::string type, const Vector& albedo) 
        : Texture(name, type), 
          albedo(albedo) 
        {}
    
    Vector getColor(const Vector& barCoords, Triangle& triangle, Mesh& mesh) const {
        return albedo;
    }
};

class TextureEdge : public Texture {
    Vector edgeColor;
    Vector innerColor;
    float edgeWidth;
public:
    TextureEdge(const std::string name, 
                    const std::string type, 
                    const Vector& edgeColor, 
                    const Vector& innerColor, 
                    const float edgeWidth) 
        : Texture(name, type), 
          edgeColor(edgeColor), 
          innerColor(innerColor), 
          edgeWidth(edgeWidth) 
        {}
    
    Vector getColor(const Vector& barCoords, Triangle& triangle, Mesh& mesh) const {
        if(barCoords.getX() < edgeWidth 
            || barCoords.getY() < edgeWidth 
            || barCoords.getZ() < edgeWidth)
            return edgeColor;
        else return innerColor;
    }
};

class TextureChecker : public Texture {
    Vector colorA;
    Vector colorB;
    float squareSize;
public:
    TextureChecker(const std::string name, 
                    const std::string type, 
                    const Vector& colorA, 
                    const Vector& colorB, 
                    const float squareSize) 
        : Texture(name, type), 
          colorA(colorA), 
          colorB(colorB), 
          squareSize(squareSize) 
        {}

    Vector getColor(const Vector& barCoords, Triangle& triangle, Mesh& mesh) const {
        Vector interpolatedCoords = Math::add(
            Math::add(
                Vector(barCoords.getX()) 
                    * mesh.getUVs().at(*triangle.getVertex1().getIndex()), 
                Vector(barCoords.getY()) 
                    * mesh.getUVs().at(*triangle.getVertex2().getIndex())
            ), 
            Vector(barCoords.getZ()) 
                    * mesh.getUVs().at(*triangle.getVertex0().getIndex())
        );

        int numOfSquares = 1 / squareSize;
        int currentSqaureX = int(numOfSquares * interpolatedCoords.getX());
        int currentSqaureY = int(numOfSquares * interpolatedCoords.getY());

        if(currentSqaureX % 2 == 0 
            && currentSqaureY % 2 == 0) 
            return Vector(0);
        else if(currentSqaureX % 2 == 0 
            && currentSqaureY % 2 == 1) 
            return Vector(1);
        else if(currentSqaureX % 2 == 1 
            && currentSqaureY % 2 == 0) 
            return Vector(1);
        else if(currentSqaureX % 2 == 1 
            && currentSqaureY % 2 == 1) 
            return Vector(0);
        
        return Vector(0);
    }
};

constexpr bool hasTextures = false;

class TextureBitmap : public Texture {
    std::string filePath;
    int width, height, channels;
    unsigned char *image = nullptr;

    Vector getBitmapColor(unsigned char* image, 
                            unsigned int rowIdx, 
                            unsigned int colIdx) const 
    {
        unsigned int i = rowIdx * 3;
        unsigned int j = colIdx * 3;
        Vector result = Vector(image[j + i * width], 
                                    image[j+1 + i * width], 
                                    image[j+2 + i * width]);
        result.scale((float)1/255);
        return result;
    }
public:
    TextureBitmap(const std::string name, 
                    const std::string type, 
                    const std::string filePath) 
        : Texture(name, type), 
          filePath(filePath) 
    {
        image = stbi_load(filePath.c_str(), &width, &height, &channels, 0);
        if(stbi_failure_reason()) {
            std::cout << stbi_failure_reason() << std::endl;
        }
    }

    Vector getColor(const Vector& barCoords, Triangle& triangle, Mesh& mesh) const {
        Vector interpolatedCoords = Math::add(
            Math::add(
                Vector(barCoords.getX()) 
                    * mesh.getUVs().at(*triangle.getVertex1().getIndex()), 
                Vector(barCoords.getY()) 
                    * mesh.getUVs().at(*triangle.getVertex2().getIndex())
            ), 
            Vector(barCoords.getZ()) 
                    * mesh.getUVs().at(*triangle.getVertex0().getIndex())
        );
        unsigned int rowIdx = (1 - interpolatedCoords.getY()) * height;
        unsigned int colIdx = interpolatedCoords.getX() * width;
        return getBitmapColor(image, rowIdx, colIdx);
    }
};

class Scene {
    std::string fileName;
    struct {
        Vector bgColor {0, 0, 0};
        struct {
            unsigned int width;
            unsigned int height;
            unsigned int bucketSize;
        } imageSettings;
    } settings;
    Vector cameraPosition {0, 0, 0};
    Matrix initialTransformation {3, 3};
    std::vector<Vector> vertices;
    std::vector<Mesh> meshes;
    std::vector<Light> lights;
    std::vector<Material> materials;
    std::map<std::string, std::shared_ptr<Texture>> textures;

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
            vertices.emplace_back(
                Vector(
                    static_cast<float>(arr[i].GetDouble()), 
                    static_cast<float>(arr[i+1].GetDouble()), 
                    static_cast<float>(arr[i+2].GetDouble())
                )
            );
            vertices.back().setNormal(std::shared_ptr<Vector>(new Vector(0)));
            vertices.back().setIndex(std::shared_ptr<unsigned int>(new unsigned int(i/3)));
        }
    }
    void loadUVs(const rapidjson::Value::ConstArray& arr) {
        assert(arr.Size() % 3 == 0);
        for(unsigned int i = 0; i < arr.Size(); i+=3) {
            meshes.back().getUVs().emplace_back(
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
            Vector& one = vertices.at(arr[i].GetInt());
            Vector& two = vertices.at(arr[i+1].GetInt());
            Vector& three = vertices.at(arr[i+2].GetInt());

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
        for(auto& vertex : vertices) {
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

                const rapidjson::Value& bucketSizeValue 
                            = imageSettingsValue.FindMember("bucket_size")->value;
                assert(!bucketSizeValue.IsNull() && bucketSizeValue.IsInt());

                settings.imageSettings.width = imageWidthValue.GetInt();
                settings.imageSettings.height = imageHeightValue.GetInt();
                settings.imageSettings.bucketSize = bucketSizeValue.GetInt();
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

        const rapidjson::Value& texturesValue = document.FindMember("textures")->value;
        if(!texturesValue.IsNull() && texturesValue.IsArray()) {
            for(rapidjson::Value::ConstValueIterator itr = texturesValue.Begin(); 
                    itr != texturesValue.End(); 
                    itr++) 
            {
                const rapidjson::Value& nameValue = itr->FindMember("name")->value;
                assert(!nameValue.IsNull() && nameValue.IsString());

                const rapidjson::Value& typeValue = itr->FindMember("type")->value;
                assert(!typeValue.IsNull() && typeValue.IsString());

                std::string name = nameValue.GetString();
                std::string type = typeValue.GetString();
                if(type == "albedo") {
                    const rapidjson::Value& albedoValue = itr->FindMember("albedo")->value;
                    assert(!albedoValue.IsNull() && albedoValue.IsArray());

                    textures[name] = std::unique_ptr<Texture>(new TextureColor(
                        name, 
                        type, 
                        loadVector(albedoValue.GetArray())
                    ));
                }
                else if(type == "edges") {
                    const rapidjson::Value& v1 
                            = itr->FindMember("edge_color")->value;
                    assert(!v1.IsNull() && v1.IsArray());

                    const rapidjson::Value& v2 
                            = itr->FindMember("inner_color")->value;
                    assert(!v2.IsNull() && v2.IsArray());

                    const rapidjson::Value& f 
                            = itr->FindMember("edge_width")->value;
                    assert(!f.IsNull() && f.IsDouble());

                    textures[name] = std::unique_ptr<Texture>(new TextureEdge(
                        name, 
                        type, 
                        loadVector(v1.GetArray()), 
                        loadVector(v2.GetArray()), 
                        static_cast<float>(f.GetDouble())
                    ));
                }
                else if(type == "checker") {
                    const rapidjson::Value& v1 
                            = itr->FindMember("color_A")->value;
                    assert(!v1.IsNull() && v1.IsArray());

                    const rapidjson::Value& v2 
                            = itr->FindMember("color_B")->value;
                    assert(!v2.IsNull() && v2.IsArray());

                    const rapidjson::Value& f 
                            = itr->FindMember("square_size")->value;
                    assert(!f.IsNull() && f.IsDouble());

                    textures[name] = std::unique_ptr<Texture>(new TextureChecker(
                        name, 
                        type, 
                        loadVector(v1.GetArray()), 
                        loadVector(v2.GetArray()), 
                        static_cast<float>(f.GetDouble())
                    ));
                }
                else if(type == "bitmap") {
                    const rapidjson::Value& fileValue = itr->FindMember("file_path")->value;
                    assert(!fileValue.IsNull() && fileValue.IsString());

                    textures[name] = std::unique_ptr<Texture>(new TextureBitmap(
                        name, 
                        type, 
                        fileValue.GetString()
                    ));
                }
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
                std::string texture;
                const std::string type = typeValue.GetString();
                if(type == "refractive") {
                    const rapidjson::Value& iorValue = itr->FindMember("ior")->value;
                    assert(!iorValue.IsNull() && iorValue.IsDouble());
                    ior = static_cast<float>(iorValue.GetDouble());
                }
                else if(type == "diffuse" && hasTextures) {
                    const rapidjson::Value& textureValue = itr->FindMember("albedo")->value;
                    assert(!textureValue.IsNull() && textureValue.IsString());
                    texture = textureValue.GetString();
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
                else if(type == "diffuse" && hasTextures) {
                    materials.emplace_back(
                        type, 
                        texture, 
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
                vertices.clear();

                const rapidjson::Value& materialIndexValue 
                        = itr->FindMember("material_index")->value;
                assert(!materialIndexValue.IsNull() && materialIndexValue.IsInt());

                meshes.emplace_back(Mesh(materialIndexValue.GetInt()));

                const rapidjson::Value& verticesValue 
                            = itr->FindMember("vertices")->value;
                assert(!verticesValue.IsNull() && verticesValue.IsArray());
                loadVertices(verticesValue.GetArray());

                if(hasTextures) {
                    const rapidjson::Value& uvValues = itr->FindMember("uvs")->value;
                    assert(!uvValues.IsNull() && uvValues.IsArray());
                    loadUVs(uvValues.GetArray());
                }

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
    unsigned int getBucketSize() const {
        return settings.imageSettings.bucketSize;
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
    std::map<std::string, std::shared_ptr<Texture>>& getTextures() {
        return textures;
    }

    void setFileName(const std::string fileName) {
        this->fileName = fileName;
    }
};

constexpr float rayStartOffset = 0.3;
constexpr float rayStartOffsetReflective = 0.3;
constexpr float rayStartOffsetRefractive = 0.3;
constexpr unsigned int maxRayDepth = 6;

struct AABB {
    Vector min {std::numeric_limits<float>::max()};
    Vector max {std::numeric_limits<float>::lowest()};

    Vector& get(unsigned int index) {
        if(index == 0) return min;
        else return max;
    }
};

struct MeshToTriangle {
    std::shared_ptr<Mesh> mesh;
    std::vector<std::shared_ptr<Triangle>> triangles;
};

class KDNode {
    int children[2] = {-1, -1};
    struct AABB aabb;
    std::vector<struct MeshToTriangle> intersected;
public:
    KDNode() = default;

    void addValues(AABB& aabb, std::vector<struct MeshToTriangle>& intersected) {
        this->aabb = aabb;
        this->intersected = intersected;
    }

    void setChildren(const int index0, const int index1) {
        children[0] = index0;
        children[1] = index1;
    }
    void setAABB(struct AABB& aabb) {
        this->aabb = aabb;
    }
    void setIntersected(const std::vector<struct MeshToTriangle>& intersected) {
        this->intersected = intersected;
    }

    int* getChildren() {
        return children;
    }
    struct AABB& getAABB() {
        return aabb;
    }
    std::vector<struct MeshToTriangle>& getIntersected() {
        return intersected;
    }
};

class KDTree {
    const int maxDepth = 16;

    std::vector<Mesh>& meshes;
    std::vector<struct MeshToTriangle> allTriangles;
    KDNode *nodes = new KDNode[1000000];

    void setAABB(struct AABB& aabb, Triangle& triangle) {
        for(unsigned int i = 0; i < NUM_OF_TRIANGLE_VERTICES; i++) {
            for(unsigned int j = 0; j < 3; j++) {
                float currentValue = triangle.getVertex(i).get(j);
                if(aabb.max.get(j) < currentValue) {
                    aabb.max.set(j, currentValue);
                }
                if(aabb.min.get(j) > currentValue) {
                    aabb.min.set(j, currentValue);
                }
            }
        }
    }
    void setupRoot() {
        struct AABB aabb;
        for(auto &mesh : meshes) {
            struct MeshToTriangle temp;
            temp.mesh = std::shared_ptr<Mesh>(new Mesh(mesh));
            for(auto &triangle : mesh.getTriangles()) {
                temp.triangles.emplace_back(
                    std::shared_ptr<Triangle>(new Triangle(triangle))
                );
                setAABB(aabb, triangle);
            }
            allTriangles.emplace_back(temp);
        }
        nodes[0].addValues(aabb, allTriangles);
    }
    void splitTriangles(const AABB& aabb, 
                            const std::vector<MeshToTriangle>& triangles, 
                            std::vector<MeshToTriangle>& left, 
                            std::vector<MeshToTriangle>& right, 
                            const int axis) 
    {
        float midPoint = (aabb.min.get(axis) + aabb.max.get(axis)) / 2;

        for (auto& meshToTriangle : triangles) {
            MeshToTriangle leftPart = {meshToTriangle.mesh};
            MeshToTriangle rightPart = {meshToTriangle.mesh};
            for (auto& triangle : meshToTriangle.triangles) {
                bool toLeft = false, toRight = false;
                for (int i = 0; i < NUM_OF_TRIANGLE_VERTICES; i++) {
                    float vertexCoord = triangle->getVertex(i).get(axis);
                    if (vertexCoord <= midPoint) {
                        toLeft = true;
                    }
                    if (vertexCoord >= midPoint) {
                        toRight = true;
                    }
                }
                if (toLeft) {
                    leftPart.triangles.emplace_back(triangle);
                }
                if (toRight) {
                    rightPart.triangles.emplace_back(triangle);
                }
            }
            if (!leftPart.triangles.empty()) {
                left.emplace_back(leftPart);
            }
            if (!rightPart.triangles.empty()) {
                right.emplace_back(rightPart);
            }
        }
    }
public:
    KDTree(std::vector<Mesh>& meshes) : meshes(meshes) {
        setupRoot();
    }

    void build(int parentIndex = 0, int depth = 0) {
        if (depth >= maxDepth) {
            return;
        }

        AABB parentAABB = nodes[parentIndex].getAABB();
        std::vector<MeshToTriangle>& intersected = nodes[parentIndex].getIntersected();

        int axis = depth % 3;

        std::vector<MeshToTriangle> left, right;
        splitTriangles(parentAABB, intersected, left, right, axis);

        if (left.empty() || right.empty()) {
            return;
        }

        AABB leftAABB = parentAABB;
        AABB rightAABB = parentAABB;
        float midPoint = (parentAABB.min.get(axis) + parentAABB.max.get(axis)) / 2;
        leftAABB.max.set(axis, midPoint);
        rightAABB.min.set(axis, midPoint);

        int leftChildIndex = 2 * parentIndex + 1;
        int rightChildIndex = 2 * parentIndex + 2;

        nodes[leftChildIndex].addValues(leftAABB, left);
        nodes[rightChildIndex].addValues(rightAABB, right);
        nodes[parentIndex].setChildren(leftChildIndex, rightChildIndex);

        build(leftChildIndex, depth + 1);
        build(rightChildIndex, depth + 1);
    }

    KDNode *getNodes() {
        return nodes;
    }

    void print() const {
        for(unsigned int i = 0; i < 16; i++) {
            KDNode node = nodes[i];
            std::cout << i << "\n max = ";
            node.getAABB().max.print();
            std::cout << "min = ";
            node.getAABB().min.print();

            int *children = node.getChildren();
            std::cout << "children = " << children[0] << " | " << children[1] << "\n";
            /*
            for(int mT = 0; mT < node.getIntersected().size(); mT++) {
                const MeshToTriangle &current = node.getIntersected().at(mT);
                for(auto& triangle : current.triangles) {
                    for(unsigned int k = 0; k < 3; k++) {
                        triangle->getVertex(k).print();
                    }
                }
            }
            */
        }
    }

    ~KDTree() {
        delete [] nodes;
    }
};

class Renderer {
    struct {
        Vector bgColor;
        struct {
            unsigned int width;
            unsigned int height;
            unsigned int maxColorComponent;
            unsigned int bucketSize;
        } imageSettings;
    } settings;
    std::vector<Mesh> meshes;
    Vector cameraPosition;
    Matrix initialTransformation;
    std::vector<Light> lights;
    std::vector<Material> materials;
    std::map<std::string, std::shared_ptr<Texture>> textures;

    std::list<unsigned int> buckets;
    KDTree kdTree {meshes};

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
    bool checkAABBcollision(const Vector& start, const Vector& ray, const AABB& aabb) {
        Vector tMin(0), tMax(0);

        for (int i = 0; i < 3; ++i) {
            if (ray.get(i) != 0) {
                float invD = 1.0f / ray.get(i);
                float t0 = (aabb.min.get(i) - start.get(i)) * invD;
                float t1 = (aabb.max.get(i) - start.get(i)) * invD;
                if (invD < 0.0f) std::swap(t0, t1);
                tMin.set(i, t0);
                tMax.set(i, t1);
            } else {
                if (start.get(i) < aabb.min.get(i) || start.get(i) > aabb.max.get(i)) {
                    return false;
                } else {
                    tMin.set(i, std::numeric_limits<float>::lowest());
                    tMax.set(i, std::numeric_limits<float>::max());
                }
            }
        }

        float tClose = std::max(std::max(tMin.getX(), tMin.getY()), tMin.getZ());
        float tFar = std::min(std::min(tMax.getX(), tMax.getY()), tMax.getZ());

        return tClose <= tFar && tFar >= 0.0f;
    }
    std::vector<KDNode> findIntersectingNodes(const Vector& start, const Vector& ray) {
        std::stack<int> stack;
        stack.push(0);
        std::vector<KDNode> intersectingNodes;

        while (!stack.empty()) {
            int nodeIndex = stack.top();
            stack.pop();

            KDNode& currentNode = kdTree.getNodes()[nodeIndex];
            if (!checkAABBcollision(start, ray, currentNode.getAABB())) {
                continue;
            }

            int* children = currentNode.getChildren();

            if(children[0] == -1 && children[1] == -1) {
                intersectingNodes.emplace_back(currentNode);
            }
            else {
                if (children[0] != -1) {
                    stack.push(children[0]);
                }
                if (children[1] != -1) {
                    stack.push(children[1]);
                }
            }
        }

        return intersectingNodes;
    }
    bool checkRayForIntersections(const Vector& start, 
                                    const Vector& originalRay, 
                                    const bool applyShade = false, 
                                    Vector* colors = nullptr, 
                                    const unsigned int rayDepth = 0) 
    {
        float zValue = 0;
        unsigned int nearestTriangleIndex = 0;
        unsigned int nearestMeshIndex = 0;
        unsigned int nearestNodeIndex = 0;

        std::vector<KDNode> nodes = findIntersectingNodes(start, originalRay);
        if (nodes.empty()) {
            return false;
        }
        for(unsigned int nodeIndex = 0; nodeIndex < nodes.size(); nodeIndex++) {
        const std::vector<MeshToTriangle>& intersected = nodes[nodeIndex].getIntersected();

        for (unsigned int meshIndex = 0; meshIndex < intersected.size(); meshIndex++) {
            const MeshToTriangle& meshToTriangles = intersected[meshIndex];
            for (unsigned int triangleIndex = 0; 
                    triangleIndex < meshToTriangles.triangles.size(); 
                    triangleIndex++) 
            {
                Triangle& triangle = *meshToTriangles.triangles[triangleIndex];
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
                            nearestNodeIndex = nodeIndex;
                        }
                        else if(applyShade == false 
                                && materials.at(meshToTriangles.mesh->getMaterialIndex())
                                    .getType() 
                                        != "refractive") 
                            return true;
                    }
                }
            }
        }
        }
        if(applyShade == true && zValue != 0) {
            Vector intersection = originalRay;
            intersection.scale(zValue);
            intersection = Math::add(start, intersection);
            *colors = shade(*nodes[nearestNodeIndex].getIntersected()[nearestMeshIndex]
                                    .triangles[nearestTriangleIndex], 
                                *nodes[nearestNodeIndex].getIntersected()[nearestMeshIndex]
                                    .mesh, 
                                intersection, 
                                rayDepth);
            return true;
        }
        return false;
    }
    bool checkForIntersections(const Vector& start, const Vector& end) {
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
            if(rayDepth > maxRayDepth) return settings.bgColor * albedo;

            colors = reflect(intersection, normal, albedo, rayDepth);
            return colors;
        }
        if(materials.at(mesh.getMaterialIndex()).getType() == "refractive") {
            if(rayDepth > maxRayDepth) return settings.bgColor;

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

        if(hasTextures == true) {
            albedo = textures[materials.at(mesh.getMaterialIndex()).getAlbedoTexture()]
                                ->getColor(
                                    getBarycentricCoords(triangle, intersection), 
                                    triangle, 
                                    mesh
                                );
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
                            255, 
                            scene.getBucketSize()
                        }
                    }, 
        meshes(scene.getMeshes()), 
        cameraPosition(scene.getCameraPosition()), 
        initialTransformation(scene.getInitialTransformation()), 
        lights(scene.getLights()), 
        materials(scene.getMaterials()), 
        textures(scene.getTextures()) 
    {
        unsigned int numOfBuckets = settings.imageSettings.width 
                                        / settings.imageSettings.bucketSize
                                        * settings.imageSettings.height 
                                        / settings.imageSettings.bucketSize;
        for(unsigned int i = 0; i < numOfBuckets; i++) {
            buckets.emplace_back(i);
            //std::cout << buckets.back() << std::endl;
        }

        kdTree.build();
        //kdTree.print();
    }

    void render(const std::string fileName, const unsigned int frame = 0) {

        settings.bgColor = {
            (float)(int)
            (settings.bgColor.getX() * settings.imageSettings.maxColorComponent), 
            (float)(int)
            (settings.bgColor.getY() * settings.imageSettings.maxColorComponent), 
            (float)(int)
            (settings.bgColor.getZ() * settings.imageSettings.maxColorComponent)
        };

        unsigned int imageBufferSize = settings.imageSettings.width 
                                            * settings.imageSettings.height 
                                            * 3;
        unsigned int *imageBuffer = new unsigned int[imageBufferSize];
        for(unsigned int i = 0; i < imageBufferSize; i++) {
            imageBuffer[i] = 0;
        }

        std::vector<std::thread> threads;
        unsigned int numOfThreads = std::thread::hardware_concurrency();
        if(numOfThreads == 0) numOfThreads = 1;

        for(unsigned int t = 0; t < numOfThreads; t++) {
            threads.emplace_back([this, imageBuffer, frame]{
            while(1) {
            mtx.lock();
            if(buckets.size() == 0) {
                mtx.unlock();
                return;
            }
            int bucket = buckets.front();
            buckets.pop_front();
            mtx.unlock();
            int bucketX = settings.imageSettings.bucketSize * (bucket 
                            / (settings.imageSettings.width 
                                        / settings.imageSettings.bucketSize));
            int bucketY = settings.imageSettings.bucketSize * (bucket
                            % (settings.imageSettings.width 
                                        / settings.imageSettings.bucketSize));
            //std::cout << bucket << " " << bucketX << " " << bucketY << std::endl;
            for (int rowIdx = bucketX; 
                        rowIdx < bucketX + settings.imageSettings.bucketSize; ++rowIdx) {
                for (int colIdx = bucketY; 
                            colIdx < bucketY + settings.imageSettings.bucketSize; ++colIdx) {
                    //system("clear");
                    //std::cout << "row: " << rowIdx << " column: " << colIdx << std::endl;
                    float x = colIdx + 0.5, y = rowIdx + 0.5;
                    x /= settings.imageSettings.width;
                    y /= settings.imageSettings.height;
                    x = (2 * x) - 1;
                    y = 1 - (2 * y);
                    x *= (float) settings.imageSettings.width 
                            / settings.imageSettings.height;

                    Vector colors(0);
                    colors = settings.bgColor;

                    Camera camera(cameraPosition, Vector(x, y, -1));
                    camera.applyRotation(initialTransformation);
                    camera.pan(((float)16 / 72) * frame);
                    camera.dolly(((float)10 / 72) * frame);
                    camera.truck(((float)10 / 72) * frame);
                    camera.dolly(((float)-28 / 72) * frame);
                    camera.pedestal(((float)-6 / 72) * frame);
                    checkRayForIntersections(camera.getPosition(), 
                                                camera.getImagePlane(), 
                                                true, 
                                                &colors);
                    
                    imageBuffer[colIdx*3 + settings.imageSettings.width * rowIdx * 3] 
                            = colors.getX();
                    imageBuffer[colIdx*3+1 + settings.imageSettings.width * rowIdx * 3] 
                            = colors.getY();
                    imageBuffer[colIdx*3+2 + settings.imageSettings.width * rowIdx * 3] 
                            = colors.getZ();
                }
            }
            }});
        }

        for(auto &th : threads) th.join();

        std::ofstream ppmFileStream(fileName, 
                        std::ios::out | std::ios::binary);
        ppmFileStream << "P3\n";
        ppmFileStream << settings.imageSettings.width 
        << " " << settings.imageSettings.height << "\n";
        ppmFileStream << settings.imageSettings.maxColorComponent << "\n";

        for (int rowIdx = 0; rowIdx < settings.imageSettings.height; ++rowIdx) {
            for (int colIdx = 0; colIdx < settings.imageSettings.width; ++colIdx) {
                ppmFileStream 
                << imageBuffer[colIdx*3 + settings.imageSettings.width * rowIdx * 3] 
                << " " 
                << imageBuffer[colIdx*3+1 + settings.imageSettings.width * rowIdx * 3] 
                << " " 
                << imageBuffer[colIdx*3+2 + settings.imageSettings.width * rowIdx * 3] 
                << "\t";
            }
            ppmFileStream << "\n";
        }

        ppmFileStream.close();

        delete [] imageBuffer;
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

    Scene scene("scene1.crtscene");
    scene.parseSceneFile();

    for(unsigned int frame = 0; frame < 72; frame++) {
    Renderer renderer(scene);
    auto start = std::chrono::high_resolution_clock::now();
    renderer.render("scene_1.ppm", frame);
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Execution time: " << duration.count() << " us (microseconds)" << std::endl;
    }
    return 0;
}