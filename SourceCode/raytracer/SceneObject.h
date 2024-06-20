#pragma once
#include <memory>

#include "CRTTypes.h"
#include "Animation.h"

struct SceneObject {
    SceneObject() = default;
    SceneObject(const Vec3& pos, const Matrix3x3& mat) : pos(pos), mat(mat) {}
    Vec3 pos{ 0.f, 0.f, 0.f };
    Matrix3x3 mat {{ 1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f }};
    std::vector<std::unique_ptr<Animation>> animations {};

    template <typename T, typename... Args>
    void addAnimation(Args&&... args) {
        animations.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    }
};