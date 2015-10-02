#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "radiance.h"
#include "ppm.h"
#include "hdr.h"
#include "random.h"
#include "timer.h"

namespace gemspt {


void save_binary_file(const char *filename, Color *image, int width, int height) {
    FILE *fp = fopen(filename, "wb");
    if (fp != NULL) {
        unsigned char *buf = (unsigned char*)image;
        fwrite(&width, sizeof(int), 1, fp);
        fwrite(&height, sizeof(int), 1, fp);
        fwrite(buf, sizeof(Color), width * height, fp);
        fclose(fp);
    }
}

int render(const char *filename, const int width, const int height, const int num_sample_per_subpixel, const int num_thread, int dx, int dy, unsigned long long seed) {
#ifdef _OPENMP
    omp_set_num_threads(num_thread);
#endif // _OPENMP
    
    // カメラ位置。
    const Vec camera_position = Vec(7.0, 3.0, 7.0);
    const Vec camera_lookat   = Vec(0.0, 1.5, 0.0);
    const Vec camera_dir      = normalize(camera_lookat - camera_position);
    const Vec camera_up       = Vec(0.0, 1.0, 0.0);

    // ワールド座標系でのイメージセンサーの大きさ。
    const double sensor_width = 30.0 * width / height; // アスペクト比調整。
    const double sensor_height= 30.0;
    // イメージセンサーまでの距離。
    const double sensor_dist  = 45.0;
    // イメージセンサーを張るベクトル。
    const Vec sensor_x_vec = normalize(cross(camera_dir, camera_up)) * sensor_width;
    const Vec sensor_y_vec = normalize(cross(sensor_x_vec, camera_dir)) * sensor_height;
    const Vec sensor_center = camera_position + camera_dir * sensor_dist;

//    Color *image = new Color[width * height];
    std::vector<Color> image(width * height);
    std::cout << width << "x" << height << " " << num_sample_per_subpixel << " spp" << std::endl;
    
    Vec offset[2] = {
        Vec(0, 0, 0), // base
        Vec(dx, dy, 0), // offset
    };

    for (int y = 0; y < height; ++y) {
        std::cerr << "Rendering (y = " << y << ", " << (100.0 * y / (height - 1)) << " %)          \r";
#pragma omp parallel for schedule(static) // OpenMP
        for (int x = 0; x < width; ++x) {
            Random random(y * width + x + 1 + seed);

            const int image_index = (height - y - 1) * width + x;
            
            Color accumulated_radiance = Color();
            for (int sample = 0; sample < num_sample_per_subpixel; ++sample) {
                double sx = random.next01();
                double sy = random.next01();
                
                // イメージセンサー上の位置。
                const Vec position_on_sensor = 
                    sensor_center + 
                    sensor_x_vec * ((sx + x) / width - 0.5) +
                    sensor_y_vec * ((sy + y) / height- 0.5);
                
                // レイを飛ばす方向。
                const Vec dir = normalize(position_on_sensor - camera_position);

                accumulated_radiance +=
                    radiance(Ray(camera_position, dir), random, 0) / (double)num_sample_per_subpixel;
            }
            image[image_index] = image[image_index] + accumulated_radiance;
        }
    }
    std::cout << std::endl;
    
    // 出力
//    save_ppm_file(filename, image, width, height);

    char buf[256];

    sprintf(buf, "%s.hdr", filename);
    save_hdr_file(buf, &image[0], width, height);
    sprintf(buf, "%s.bin", filename);
    save_binary_file(buf, &image[0], width, height);
    
    return 0;
}

enum ReflectionType {
    Reflection,
    Refraction,
};

struct Vertex {
    const Vec position;
    const Vec normal;
    const Material *material;
    const Vec half;
    const Vec out_dir;

    const ReflectionType reflection_type;
    const SceneSphere *object;

    Vertex(const Vec &aPosition, const Vec &aNormal, const Material *aMaterial, const Vec &aHalf, const SceneSphere *aObject, const Vec &aOutDir = Vec(), const ReflectionType &aReflectionType = Reflection):
    position(aPosition), normal(aNormal), material(aMaterial), half(aHalf), out_dir(aOutDir), reflection_type(aReflectionType), object(aObject) {}
};

Color radiance_gpt(const Ray &ray, Random &random, const int depth, double *path_pdf, std::vector<Vertex> *verts) {
    const Color kBackgroundColor = Color(0.0f, 0.0f, 0.0f);
    const int kDepthLimit = 4;
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
        // 光源で終わりなのでhalf-vector定義できない
        verts->push_back(Vertex(hitpoint.position, hitpoint.normal, now_material, Vec(0, 0, 0), now_object));

        // 光源にヒットしたら放射項だけ返して終わる。
        // （今回、光源は反射率0と仮定しているため）
        return emission;
    }
    
    // 次の方向をサンプリング。
    // BRDFの値。
    Color brdf_eval;
    double pdf = -1;
    const Vec dir_out = now_material->sample(random, ray.dir, hitpoint.normal, &pdf, &brdf_eval);

    //brdf_eval=now_material->eval(ray.dir, hitpoint.normal, dir_out);

    const Vec now_normal = dot(hitpoint.normal, ray.dir) < 0.0 ? hitpoint.normal: -hitpoint.normal; // 交差位置の法線（物体からのレイの入出を考慮。
    const ReflectionType reflection_type = dot(now_normal, dir_out) > 0.0 ? Reflection : Refraction;
    
    if (now_object->get_sphere()->position().y < -1000) {
        
        int X = std::floor(hitpoint.position.x * 4.0);
        int Y = std::floor(hitpoint.position.z * 4.0);

        if ((X + Y) % 2 == 0)
            brdf_eval = 0.2 * brdf_eval;
            
        /*
        double f = 0.2 * (sin(hitpoint.position.x * 2.0) + sin(hitpoint.position.z * 2.0)) + 0.4;
        brdf_eval = f  * brdf_eval;
       */
    }

    // 頂点列構築
    const Vec half = normalize(normalize(-ray.dir) + normalize(dir_out));
    verts->push_back(Vertex(hitpoint.position, hitpoint.normal, now_material, half, now_object, dir_out, reflection_type));

    // cos項。
    const double cost = dot(hitpoint.normal, dir_out);

    *path_pdf *= pdf;

    // レンダリング方程式をモンテカルロ積分によって再帰的に解く。
    Color L = multiply(brdf_eval, radiance_gpt(Ray(hitpoint.position, dir_out),random, depth + 1, path_pdf, verts));
    L *= cost;
    return L;
}

bool is_diffuse(const Material *a) {
    bool diffuse = (typeid(*a) == typeid(LambertianMaterial)) || (typeid(*a) == typeid(Lightsource));
    return diffuse;
}
bool is_lightsource(const Material *a) {
    return typeid(*a) == typeid(Lightsource);
}
bool check_same_material(const Material *a, const Material *b) {
    bool a_diffuse = is_diffuse(a);
    bool b_diffuse = is_diffuse(b);

    return a_diffuse == b_diffuse;
}
bool is_specular(const Material *a) {
    return !is_diffuse(a);
}

enum Result {
    Invertible,
    NotInveritble,
    NotSymmetric,
};

Result shift(const std::vector<Vertex> &baseVertices, const Ray &ray, std::vector<Vertex> *offsetVertices, Color *f, double *J, double *sub_pdf) {
    // いきなり光源にヒット or いきなり外へ飛んで行ったらnot invertibleということにしてnaiveにgradient計算する
    if (baseVertices.size() <= 2) 
        return NotInveritble;

    // (i) Initial Offset Segment
    {
        Hitpoint hitpoint;
        const SceneSphere *now_object = intersect_scene(ray, &hitpoint);
        if (now_object == NULL)
            return NotInveritble;
        offsetVertices->push_back(Vertex(hitpoint.position, hitpoint.normal, now_object->get_material(), Vec(), now_object)); // offset-pathではhalf-vector使わないので0を入れる
    }
    
    *J = 1;
    *sub_pdf = 1;
    *f = Color(1, 1, 1);

    int index = 1; // index = 0 はカメラ位置
    // (ii) Tracing Additional Segments
    for (;;++index) {
        if (!check_same_material(baseVertices[index].material, (*offsetVertices)[index].material)){
            return NotInveritble;
        }
        
        if (index + 1 == baseVertices.size() && !is_lightsource(baseVertices[index].material)) {
            // 光源に到達しないで途中で打ち切られるパスであった
            return NotInveritble;
        }

        if (is_specular(baseVertices[index].material) || (index + 1 < baseVertices.size() && is_specular(baseVertices[index + 1].material))) {
            // half-vectorコピー

            // refraction -> internal reflectionチェック

            // めんどいのでちょっと簡易的な実装にするよ
            Vec dir_out;

            ReflectionType reflection_type;

            if (is_diffuse(baseVertices[index].material)) {
                // half-vectorコピー
                const Vec now_half = baseVertices[index].half;
                const Vec in_dir = normalize((*offsetVertices)[index].position - (*offsetVertices)[index - 1].position);
                const double A = dot(now_half, -in_dir);
                const double t = 2 * A;
                dir_out = normalize(t * now_half + in_dir);

                reflection_type = Reflection;
            } else {
                // スペキュラであったので、普通にスペキュラ処理する
                // 完全スペキュラだけ仮定する
                const GlassMaterial *glass = dynamic_cast<const GlassMaterial*>(baseVertices[index].material);

                Vec reflection_dir;
                Vec refraction_dir;
                double Fr, Ft;
                const Vec in_dir = normalize((*offsetVertices)[index].position - (*offsetVertices)[index - 1].position);
                bool internal_reflection = !glass->calc_dir(in_dir, (*offsetVertices)[index].normal, &reflection_dir, &refraction_dir, &Fr, &Ft);

                if (baseVertices[index].reflection_type == Refraction && internal_reflection) {
                    // not invertible
                    return NotInveritble;
                }
                // baseパスの反射・屈折にそろえる
                if (baseVertices[index].reflection_type == Reflection) {
                    dir_out = reflection_dir;
                    reflection_type = Reflection;
                } else if (baseVertices[index].reflection_type == Refraction) {
                    dir_out = refraction_dir;
                    reflection_type = Refraction;
                }
            }

            // 次のoffset頂点へ向かってレイ飛ばす
            {
                Hitpoint hitpoint;
                const SceneSphere *now_object = intersect_scene(Ray((*offsetVertices)[index].position, dir_out), &hitpoint);
                if (now_object == NULL)
                    return NotInveritble;
                
                offsetVertices->push_back(Vertex(hitpoint.position, hitpoint.normal, now_object->get_material(), Vec(), now_object)); // offset-pathではhalf-vector使わないので0を入れる

                // ヤコビアン計算
                const Vec omega_o_x = normalize(baseVertices[index - 1].position - baseVertices[index].position);
                const Vec omega_i_x = normalize(baseVertices[index + 1].position - baseVertices[index].position);
                const Vec omega_o_y = normalize((*offsetVertices)[index - 1].position - (*offsetVertices)[index].position);
                const Vec omega_i_y = normalize((*offsetVertices)[index + 1].position - (*offsetVertices)[index].position);
                const Vec normal_x = baseVertices[index].normal;
                const Vec normal_y = (*offsetVertices)[index].normal;

                if (reflection_type == Reflection) {
                    const Vec h_x = normalize(omega_o_x + omega_i_x);
                    const Vec h_y = normalize(omega_o_y + omega_i_y);
                    *J *= dot(omega_o_y, h_y) / dot(omega_o_x, h_x);
                } else if (reflection_type == Refraction) {
                    const GlassMaterial *glass_x = dynamic_cast<const GlassMaterial*>(baseVertices[index].material);
                    const GlassMaterial *glass_y = dynamic_cast<const GlassMaterial*>((*offsetVertices)[index].material);
                    const Vec now_normal_x = dot(normal_x, -omega_o_x) < 0.0 ? normal_x: -normal_x; 
                    const Vec now_normal_y = dot(normal_y, -omega_o_y) < 0.0 ? normal_y: -normal_y; 

                    double n1_x, n2_x, n1_y, n2_y;
                    if (dot(normal_x, -omega_o_x) < 0.0) {
                        // 外から中に入ってくる
                        n1_x = 1;
                        n2_x = glass_x->ior();
                    } else {
                        // 中から外に入ってくる
                        n2_x = 1;
                        n1_x = glass_x->ior();
                    }
                    
                    if (dot(normal_y, -omega_o_y) < 0.0) {
                        // 外から中に入ってくる
                        n1_y = 1;
                        n2_y = glass_y->ior();
                    } else {
                        // 中から外に入ってくる
                        n2_y = 1;
                        n1_y = glass_y->ior();
                    }

                    const Vec h_x = normalize(omega_o_x * n1_x + omega_i_x * n2_x);
                    const Vec h_y = normalize(omega_o_y * n1_y + omega_i_y * n2_y);

                    const double rJ = 
                        (omega_i_y + (n2_y / n1_y) * omega_o_y).length_squared() / (omega_i_x + (n2_x / n1_x) * omega_o_x).length_squared()
                         * dot(omega_i_x, h_x) / dot(omega_i_y, h_y);

                    *J *= rJ;
                }
            }
        } else {

            if (index + 1 < baseVertices.size()) {
                // 現在と次のbase頂点マテリアルがdiffuse
                // このとき、現在のoffset頂点もdiffuseになってるはず
                // Reconnect
                const Vec next_offset_position = baseVertices[index + 1].position;
                const Vec now_offset_position = (*offsetVertices)[index].position;

                // occlude判定
                Hitpoint th;
                const SceneSphere *s = intersect_scene(Ray(now_offset_position, normalize(next_offset_position - now_offset_position)), &th);
                const Vec real_next_offset_position = th.position;

                if ((real_next_offset_position - next_offset_position).length_squared() > 1e-3) {
                    // occludeされた!
                    return NotSymmetric;
                }

                // 無事offsetがbaseパスに合流
                // 残りの頂点全部コピー
                for (int i = index + 1; i < baseVertices.size(); ++i) {
                    offsetVertices->push_back(baseVertices[i]);
                }
        
                // 最後のヤコビアンのdet、J計算する
                {
                    const Vec X_1y = now_offset_position;
                    const Vec X_2y = next_offset_position;
                    const Vec X_1x = baseVertices[index].position;
                    const Vec X_2x = next_offset_position;
                    
                    const Vec normal_x = baseVertices[index + 1].normal;
                    const Vec normal_y = baseVertices[index + 1].normal;
                    
                    *J *= dot(normal_y, normalize(X_1y - X_2y)) / dot(normal_x, normalize(X_1x - X_2x)) * (X_1x - X_2x).length_squared() / (X_1y - X_2y).length_squared();
                    
                }
            } else {
                // offset-pathが最後までbaseに合流しなかったときにくる。
            }

            // パスがReconnectされたはずなのでスループットf(y)を計算して終了する
            // 一緒にpdfも計算する
            Color f_subpath(1, 1, 1);
            double pdf = 1;
            bool debug = false;
            for (int i = 1; i < offsetVertices->size(); ++i) {
                Vertex prev_v = (*offsetVertices)[i - 1];
                Vertex v = (*offsetVertices)[i];
            
                const Color emission = v.material->emission();
                if (emission.x > 0.0 || emission.y > 0.0 || emission.z > 0.0) {
                    if (i + 1 >= offsetVertices->size())
                        f_subpath = multiply(f_subpath, emission);
                    else
                        f_subpath = Color(); // 光源の反射率は0を仮定しているため、途中で光源にヒットした場合スループットは0
                    break;
                }

                if (i + 1 >= offsetVertices->size()) {
                    f_subpath = multiply(f_subpath, Color());
                    break;
                }
            
                Vertex next_v = (*offsetVertices)[i + 1];
                const Vec in_dir = normalize(v.position - prev_v.position);
                const Vec normal = v.normal;
                const Vec dir_out = normalize(next_v.position - v.position);

                if (typeid(*v.material) == typeid(GlassMaterial)) {
                    if (offsetVertices->size() <= 5) {
                        // debug
//                        std::cout << "*";
                    }

                    Color brdf_eval = v.material->eval(in_dir, normal, dir_out); // Fr, Ftは抜かれている
                    const GlassMaterial *glass = dynamic_cast<const GlassMaterial*>(v.material);

                    Vec reflection_dir;
                    Vec refraction_dir;
                    double Fr, Ft;
                    glass->calc_dir(in_dir, (*offsetVertices)[i].normal, &reflection_dir, &refraction_dir, &Fr, &Ft);
                    
                    const Vec now_normal = dot(normal, in_dir) < 0.0 ? normal: -normal; 
                    const double cost = dot(normal, dir_out);
                    
                    const double probability  = Fr; // GlassMaterialからのマジックナンバー
                    if (dot(now_normal, dir_out) > 0.0) {
                        // reflection
                        f_subpath = multiply(f_subpath, cost * brdf_eval * Fr);
                        pdf *= probability;
                    } else {
                        // refraction
                        f_subpath = multiply(f_subpath, cost * brdf_eval * Ft);
                        pdf *= (1 - probability);
                    }
                } else {
                    Color brdf_eval = v.material->eval(in_dir, normal, dir_out);
            
                    if (v.object->get_sphere()->position().y < -1000) {
                        
                        int X = std::floor(v.position.x * 4.0);
                        int Y = std::floor(v.position.z * 4.0);

                        if ((X + Y) % 2 == 0)
                            brdf_eval = 0.2 * brdf_eval;
                        /*
                        double f = 0.2 * (sin(v.position.x * 2.0) + sin(v.position.z * 2.0)) + 0.4;
                        brdf_eval = f  * brdf_eval;
                       */
                    }

                    const double now_pdf = v.material->calc_pdf(in_dir, normal, dir_out);

                    const double cost = dot(normal, dir_out);
                    f_subpath = multiply(f_subpath, cost * brdf_eval);
                    pdf *= now_pdf;
                }
            }
            *f = f_subpath;
            *sub_pdf = pdf;

            if (debug && f->x > 0) {
                // debug
//                std::cout << "*";
            }

            return Invertible;
        }
    }
}

void render_gpt(const char *filename, const int width, const int height, const int num_sample_per_subpixel, const int num_thread, unsigned long long seed) {
#ifdef _OPENMP
    omp_set_num_threads(num_thread);
#endif

    // カメラ位置。
    const Vec camera_position = Vec(7.0, 3.0, 7.0);
    const Vec camera_lookat   = Vec(0.0, 1.5, 0.0);
    const Vec camera_dir      = normalize(camera_lookat - camera_position);
    const Vec camera_up       = Vec(0.0, 1.0, 0.0);

    // ワールド座標系でのイメージセンサーの大きさ。
    const double sensor_width = 30.0 * width / height; // アスペクト比調整。
    const double sensor_height= 30.0;
    // イメージセンサーまでの距離。
    const double sensor_dist  = 45.0;
    // イメージセンサーを張るベクトル。
    const Vec sensor_x_vec = normalize(cross(camera_dir, camera_up)) * sensor_width;
    const Vec sensor_y_vec = normalize(cross(sensor_x_vec, camera_dir)) * sensor_height;
    const Vec sensor_center = camera_position + camera_dir * sensor_dist;

    std::cout << width << "x" << height << " " << num_sample_per_subpixel << " spp" << std::endl;
    
    std::vector<Color> coarse_image(width * height);
    std::vector<Color> debug_image(width * height);
    std::vector<Color> diff_image[4];
    for (int i = 0; i < 4; ++i)
        diff_image[i].resize(width * height);

    for (int y = 0; y < height; ++y) {
        std::cerr << "Rendering (y = " << y << ", " << (100.0 * y / (height - 1)) << " %)          \r";
        
#pragma omp parallel for schedule(static) // OpenMP
        for (int x = 0; x < width; ++x) {
            Random random(y * width + x + 1 + seed);

            const int image_index = (height - y - 1) * width + x;

            for (int sample = 0; sample < num_sample_per_subpixel; ++sample) {
                double sx = random.next01();
                double sy = random.next01();

                // base path
                Color I = Color();
                std::vector<Vertex> baseVertices;
                baseVertices.reserve(16);
                double basePathPDF = 1;
                Color f_base;
                double cx = sx + x;
                double cy = sy + y;
                {
                    const Vec position_on_sensor = 
                        sensor_center + 
                        sensor_x_vec * ((sx + x) / width - 0.5) +
                        sensor_y_vec * ((sy + y) / height- 0.5);
                    const Vec dir = normalize(position_on_sensor - camera_position);
                    
                    baseVertices.push_back(Vertex(camera_position, Vec(), NULL, Vec(), NULL));
                    f_base = radiance_gpt(Ray(camera_position, dir), random, 0, &basePathPDF, &baseVertices);
                    I += f_base / basePathPDF;
                }

                coarse_image[image_index] += I;

                // 4 sub-path
                const Vec offset[4] = {
                    Vec( 1, 0, 0),
                    Vec(-1, 0, 0),
                    Vec( 0, 1, 0),
                    Vec( 0,-1, 0)
                };
                std::vector<Vertex> offsetVertices[4];
                for (int offsetDir = 0; offsetDir < 4; ++offsetDir) {
                    const Vec position_on_sensor = 
                        sensor_center + 
                        sensor_x_vec * ((sx + x + offset[offsetDir].x) / width - 0.5) +
                        sensor_y_vec * ((sy + y + offset[offsetDir].y) / height- 0.5);
                    const Vec dir = normalize(position_on_sensor - camera_position);
                    offsetVertices[offsetDir].reserve(16);
                    offsetVertices[offsetDir].push_back(Vertex(camera_position, Vec(), NULL, Vec(), NULL));

                    Color f_subpath;
                    double J;
                    double subPathPDF;
                    Result result = shift(baseVertices, Ray(camera_position, dir), &offsetVertices[offsetDir], &f_subpath, &J, &subPathPDF);
                    Color accum;
                    if (result == Result::NotInveritble) {
                        // naiveにgradient求める！
                        std::vector<Vertex> tmp;
                        subPathPDF = 1; // 初期化
                        f_subpath = radiance_gpt(Ray(camera_position, dir), random, 0, &subPathPDF, &tmp);
                        accum = f_base / basePathPDF - f_subpath / subPathPDF;
                    } else if (result == Result::Invertible) {
                        // calc MIS Weight
                        double w_ij = basePathPDF / (basePathPDF + subPathPDF * J);
                        accum= w_ij * (f_base - f_subpath * J) / basePathPDF;
                    } else if (result == Result::NotSymmetric) {
                        const double w_ij = 1;
                        accum = w_ij * f_base / basePathPDF;
                    }
                    
                    // 前方差分 or 後方差分
                    if (offsetDir == 1 || offsetDir == 2)
                        diff_image[offsetDir][image_index] += accum;
                    else
                        diff_image[offsetDir][image_index] += -accum;

                    debug_image[image_index] += J * Color(1, 1, 1) / num_sample_per_subpixel;
                }
            }
            coarse_image[image_index] = coarse_image[image_index] / num_sample_per_subpixel;
            for (int i = 0; i < 4; ++i)
                diff_image[i][image_index] = diff_image[i][image_index] / num_sample_per_subpixel;

        }
    }
    std::cout << std::endl;

    // 以下、各種バッファを画像として出力
    char buf[256];
    char name[4][256] = {
        "__1_0",
        "_-1_0",
        "__0_1",
        "__0-1"
    };
    sprintf(buf, "%s_J.hdr", filename);
    save_hdr_file(buf, &debug_image[0], width, height);

    sprintf(buf, "%s_coarse.hdr", filename);
    save_hdr_file(buf, &coarse_image[0], width, height);
    sprintf(buf, "%s_coarse.bin", filename);
    save_binary_file(buf, &coarse_image[0], width, height);

    for (int i = 0; i < 4; ++i) {
        sprintf(buf, "%s_%s.hdr", filename, name[i]);
        save_hdr_file(buf, &diff_image[i][0], width, height);
        sprintf(buf, "%s_%s.bin", filename, name[i]);
        save_binary_file(buf, &diff_image[i][0], width, height);
        
#if 0
        sprintf(buf, "%s_%s_p.hdr", filename, name[i]);
        save_hdr_file_positive(buf, &diff_image[i][0], width, height);
        sprintf(buf, "%s_%s_n.hdr", filename, name[i]);
        save_hdr_file_negative(buf, &diff_image[i][0], width, height);
#endif
    }
}

}

int main() {
    std::cout << "gdpt"<< std::endl;

    const int width = 640;
    const int height = 480;
    // const int sample = 8192;
    const int nthread = 8;
    //gemspt::render("output/reference", width, height, sample, nthread,  0, 0, 0);
#if 0
    {
        const int sample = 8192;
        Timer timer;
        timer.begin();
        
        char buf[256];
        sprintf(buf, "lfreq/reference_%03ds", sample);
        gemspt::render(buf, width, height, sample, nthread,  0, 0, 0);

        const float sec = timer.end() / 1000.0f;
        sprintf(buf, "lfreq/reference_%03ds_log.txt", sample);
        FILE *fp = fopen(buf, "wt");
        fprintf(fp, "%f sec", sec);
        fclose(fp);
    }
#endif

#if 0
    Timer timer;
    timer.begin();

    const int sample = 32;

    char buf[256];
    sprintf(buf, "output2/result2_%03ds", sample);
    gemspt::render_gpt(buf, width, height, sample, nthread, 0);


    const float sec = timer.end() / 1000.0f;
    sprintf(buf, "output2/result2_%03ds_log.txt", sample);
    FILE *fp = fopen(buf, "wt");
    fprintf(fp, "%f sec", sec);
    fclose(fp);
#endif

#if 1
    int sample_table[] = {
        4,
        16,
        64,
        256,
    };

    for (int i = 0; i < 4; ++i) {
        Timer timer;
        timer.begin();

        const int sample = sample_table[i];

        char buf[256];
        sprintf(buf, "output/result2_%03ds", sample);
        gemspt::render_gpt(buf, width, height, sample, nthread, 0);


        const float sec = timer.end() / 1000.0f;
        sprintf(buf, "output/result2_%03ds_log.txt", sample);
        FILE *fp = fopen(buf, "wt");
        fprintf(fp, "%f sec", sec);
        fclose(fp);
    }
#endif

    return 0;
}
