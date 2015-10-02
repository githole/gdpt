#ifndef _RADIANCE_H_
#define _RADIANCE_H_

#include <algorithm>

#include "ray.h"
#include "scene.h"
#include "sphere.h"
#include "hitpoint.h"
#include "random.h"

namespace gemspt {

// ray方向からの放射輝度を求める
Color radiance(const Ray &ray, Random &random, const int depth) {
    const Color kBackgroundColor = Color(0.0f, 0.0f, 0.0f);
    const int kDepthLimit = 10;
    // 打ち切りチェック
    if (depth >= kDepthLimit)
        return Color();
    
    // シーンと交差判定
    Hitpoint hitpoint;
    const SceneSphere *now_object = intersect_scene(ray, &hitpoint);
    // 交差チェック
    if (now_object == NULL)
        return kBackgroundColor;

    // マテリアル取得
    const Material *now_material = now_object->get_material();
    const Color emission = now_material->emission();
    if (emission.x > 0.0 || emission.y > 0.0 || emission.z > 0.0) {

        // 光源にヒットしたら放射項だけ返して終わる。
        // （今回、光源は反射率0と仮定しているため）
        return emission;
    }
    
    // 次の方向をサンプリング。
    double pdf = -1;
    Color brdf_eval;
    const Vec dir_out = now_material->sample(random, ray.dir, hitpoint.normal, &pdf, &brdf_eval);
    
    if (now_object->get_sphere()->position().y < -1000) {
        /*
        double f = 0.2 * (sin(hitpoint.position.x * 2.0) + sin(hitpoint.position.z * 2.0)) + 0.4;
        brdf_eval = f  * brdf_eval;
        */
        
        int X = std::floor(hitpoint.position.x * 4.0);
        int Y = std::floor(hitpoint.position.z * 4.0);

        if ((X + Y) % 2 == 0) {
            brdf_eval = 0.2 * brdf_eval;
        }
        
    }

    // cos項。
    const double cost = dot(hitpoint.normal, dir_out);

    // レンダリング方程式をモンテカルロ積分によって再帰的に解く。
    const Color L = multiply(
        brdf_eval,
        radiance(Ray(hitpoint.position, dir_out),random, depth + 1))
        * cost / pdf;
    return L;
}

};

#endif
