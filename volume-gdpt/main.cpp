#include <iostream>
#include <cmath>

#include "hlib\image.h"
#include "hlib\sampling.h"
#include "hlib\random.h"
#include "hlib\rt\ray.h"
#include "hlib\rt\ibl.h"

#include "hlib\memfile.h"
#include "hlib\hmath.h"

#include "hlib\rt\sphere.h"

#include "hlib\timer.h"
#include "hlib\rt\camera.h"
#include "hlib\rt\bbox.h"
#include "hlib\rt\objMesh.h"

#include "hlib\ext\bmpexporter.h"

using namespace hstd;

void save_binary_file(const char *filename, Color *image, int width, int height) {
    FILE *fp = fopen(filename, "wb");
    if (fp != NULL) {
        struct DColor {
            double x, y, z;
        };
        std::vector<DColor> col(width * height);
        for (int iy = 0; iy < height; ++iy)
            for (int ix = 0; ix < width; ++ix) {
                col[iy * width + ix].x = image[iy * width + ix].x;
                col[iy * width + ix].y = image[iy * width + ix].y;
                col[iy * width + ix].z = image[iy * width + ix].z;
            }

        unsigned char *buf = (unsigned char*)&col[0];
        fwrite(&width, sizeof(int), 1, fp);
        fwrite(&height, sizeof(int), 1, fp);
        fwrite(buf, sizeof(DColor), width * height, fp);
        fclose(fp);
    }
}

void save_hdr_file(const std::string &filename, const Color *image, const int width, const int height) {
    hstd::Image nimage(width, height);

    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            nimage.at(ix, iy) = image[iy * width + ix];
        }
    }

    hstd::HDROperator::save(filename.c_str(), &nimage, true);
}
class SimpleVolume {
private:
    std::vector<float> density_;

    int Nx_, Ny_, Nz_;

    float max_density_;
    float max_scattering_;
    float max_transmittance_;
    float average_transmittance_;
    float albedo_;

    Float3 org_; // cell(0, 0, 0)の端っこの点の世界座標系における座標
    Float3 scale_; // m / cell

    Float3 factor_;
    
    inline unsigned int I(unsigned int ix, unsigned int iy, unsigned int iz) const {
        return iz * (Nx_ * Ny_) + iy * Nx_ + ix;
    }
    
    inline Float3 is2ws(const Int3& is_pos) const {
        return (is_pos + Float3(0.5f, 0.5f, 0.5f));
    }

    inline Float3 uvw(const Float3& ws_pos) const {
        Float3 pos;
        
        pos = ws_pos - org_;

        Float3 u = times(factor_, pos);

        // texture wraping
        // repeat
        u.x = u.x - floor(u.x);
        u.y = u.y - floor(u.y);
        u.z = u.z - floor(u.z);

        return u;
    }


    rt::BBox bbox_;
public:
    inline float sample_transmittance(const Float3& ws_pos) const {
        const Float3 uv_pos = uvw(ws_pos);
        const Float3 is_pos = times(uv_pos, Float3(Nx_, Ny_, Nz_));
        
        float ix = is_pos.x;
        float iy = is_pos.y;
        float iz = is_pos.z;

        if (ix < 0 || iy < 0 || iz < 0 ||
            Nx_ <= ix || Ny_ <= iy || Nz_ <= iz)
            return 0;

        const float density =  density_[I(ix, iy, iz)];

        return density / albedo_;
    }

    float average_transmittance() const {
        return average_transmittance_;
    }

    SimpleVolume(std::string filename, float density_scale, float albedo, Float3 &org, Float3 &scale) {
        FileManager fm;

        fm.load(filename.c_str());

        std::vector<unsigned char>& buf = fm.buffer();
        unsigned int *u_buf = (unsigned int*)&buf[0];


        Nx_ = u_buf[1];
        Ny_ = u_buf[2];
        Nz_ = u_buf[3];

        printf("Nx, Ny, Nz: %d %d %d\n", Nx_, Ny_, Nz_);

        density_.resize(Nx_ * Ny_ * Nz_);

        // float版
        float *f_buf = (float*)&u_buf[4];
        max_density_ = -1;

        double average_density = 0;

        for (int i = 0; i < Nx_ * Ny_ * Nz_; ++i) {
            density_[i] = density_scale * f_buf[i];
            average_density += density_[i];

            if (max_density_ < density_[i])
                max_density_ = density_[i];
        }
        average_density /= (Nx_ * Ny_ * Nz_);


        albedo_ = albedo;
        max_scattering_ = max_density_;
        max_transmittance_ = max_scattering_ / albedo_;
        average_transmittance_ = average_density / albedo_;
        printf("Average Density: %f\n", average_density);
        printf("Average Transmittance: %f\n", average_transmittance_);

        org_ = org;
        scale_ = scale;
        bbox_ = rt::BBox(org_, org_ + times(scale_, Float3(Nx_, Ny_, Nz_)));

        factor_ = Float3(1 / scale_.x / Nx_, 1 / scale_.y / Ny_, 1 / scale_.z / Nz_);

        printf("Max Scattering: %f\n", max_scattering_);
    }

    bool inside(const Float3& pt) const {
        return bbox_.inside(pt);
    }
    bool check_intersect(const rt::Ray& ray, float* hitt0, float* hitt1) const {
        return bbox_.check_intersect(ray, hitt0, hitt1);
    }
    
    float next_event_simple(const rt::Ray& ray, Random& random, float *pdf) const {
        /*
        float d_max = 0;
        float hitt0 = -1, hitt1 = -1;
        if (bbox_.check_intersect(ray, &hitt0, &hitt1)) {
            d_max = std::max(hitt0, hitt1);
        }

        const float d = -log(random.nexto01()) / average_transmittance_;
        *pdf = average_transmittance_ * exp(-average_transmittance_ * d);

        if (d >= d_max) {
            *pdf = -1.0f;
            return -1.0f;
        }

        return d;
        */
        
        float u = 0;
        Float3 next_pos;

//		float d_max = hitt0 > 0 ? hitt0 : 1e32f;
        float d_max = 0;
        
        float hitt0 = -1, hitt1 = -1;

        if (bbox_.check_intersect(ray, &hitt0, &hitt1)) {
            d_max = std::max(hitt0, hitt1);
        }

        for (;;) {
            u += -log(random.nexto01()) / max_transmittance_;
            next_pos = ray.org + u * ray.dir;
            if (u >= d_max) {
                *pdf = -1.0f;
                return -1.0f;
            }
            if (sample_transmittance(next_pos) / max_transmittance_ >= random.next01())
                break;
        }

        // pdfをMC積分する
        const int kNumSample = 16;
        float sum = 0;

        float S[kNumSample] = {0};
        float B[kNumSample] = {0};

        const float cTransmittance = sample_transmittance(ray.org);

        for (int i = 0; i < kNumSample; ++i) {
            const float X = -log(random.nexto01()) / cTransmittance;
            if (hitt1 < X)
                continue;

            const float cpdf = cTransmittance * exp(-cTransmittance * X);

            sum += sample_transmittance(ray.org + X * ray.dir) / cpdf;

            S[i] = cpdf;
            B[i] = X;
        }
        *pdf = sample_transmittance(ray.org + u * ray.dir) * exp(-sum / kNumSample);

        /*
        if (*pdf <= 0) {
            printf("%f %f %f %d\n",  sample_transmittance(ray.org + u * ray.dir), sum, exp(-sum / kNumSample), kNumSample);
            for (int i = 0; i < kNumSample;++i)
                printf("> %f %f\n", S[i], B[i]);
        }*/

        return u;
    }
    
    float next_event(const rt::Ray& ray, Random& random) const {
        float u = 0;
        Float3 next_pos;

//		float d_max = hitt0 > 0 ? hitt0 : 1e32f;
        float d_max = 0;
        
        float hitt0 = -1, hitt1 = -1;

        // hitt0は0を仮定（ボリューム内だとして）
        if (bbox_.check_intersect(ray, &hitt0, &hitt1)) {
            d_max = std::max(hitt0, hitt1);
        }

        for (;;) {
            u += -log(random.nexto01()) / max_transmittance_;
            next_pos = ray.org + u * ray.dir;
            if (u >= d_max) {
                return -1.0f;
            }
            if (sample_transmittance(next_pos) / max_transmittance_ >= random.next01())
                break;

        }

        return u;
    }

};

const float kAlbedo = 0.999f;

Color outside(const rt::Ray &ray) {
    Color L;
    const float factor = dot(ray.dir, normalize(Float3(0.2, 0.7, 0.1)));
    if (factor >= 0.7) {
        L += Float3(2, 1.5, 1.3); // sun
    }

    L += Float3(0.5, 0.5, 0.5); // outside
    return L;
}

Float3 radiance(Random& random, const SimpleVolume& volume, const rt::Ray& initial_ray, float *path_pdf) {
    const int kMaxScattering = 256;

    rt::Ray current_ray = initial_ray;
    
    Color througput(1, 1, 1);
    Color L(0, 0, 0);
    for (int scattering = 0; scattering < kMaxScattering; ++scattering) {
        // VolumeのBBoxにヒットするか否か
        float hitt0 = -1, hitt1 = -1;
        bool is_volume_hit = false;
        if (volume.check_intersect(current_ray, &hitt0, &hitt1)) {
            is_volume_hit = true;
        }
        
        rt::Hitpoint current_hp;
        if (!volume.inside(current_ray.org)) {
            if (is_volume_hit) {
                current_hp.distance = hitt0;
            }
        }

        if (!is_volume_hit) {
            L += outside(current_ray);
            break;
        }

        const Float3 new_org0 = current_ray.org + (hitt0 + 1e-3f) * current_ray.dir;
        const Float3 new_org1 = current_ray.org + (hitt1 + 1e-3f) * current_ray.dir;
        if (!volume.inside(current_ray.org)) {
            // Volumeの外だったら中からスタート
            current_ray.org = new_org0;
        }

        const float tu = volume.next_event(current_ray, random);
        if (tu < 0) {
            // ボリューム外
            current_ray.org = new_org1;
            continue;
        }

        // pdf_transmittance
        float transmittance = 0;
        const int kNumSample = 16;
        for (int i = 0; i < kNumSample; ++i) {
            const float X = volume.next_event(current_ray, random);
            if (!(0 <= X && X < tu)) {
                transmittance += 1.0f / kNumSample;
            }
        }

        const Float3 next_pos = current_ray.org + tu * current_ray.dir;
        const bool cond2 = random.next01() < kAlbedo;
        if (!cond2) {
            througput = Float3(0, 0, 0); // absorp
            break;
        }


        // 位相関数でインポータンスサンプリング
        // 今回はとりあえず等方散乱としておく
        const Float3 next_dir = Sampling::uniformSphereSurface(random);
        const float phase = 1.0f / (4.0f * kPI);
        const float pdf_phase = 1.0f / (4.0f * kPI);
          
        /*
        througput *= 0.99f;
        *path_pdf *= 0.99f;
        */
        
        througput *= transmittance * phase * kAlbedo;
        *path_pdf *= transmittance * pdf_phase * kAlbedo;

        

        // throughtput = (phase / pdf) * throughtput;	
        current_ray = rt::Ray(next_pos, next_dir);
    }

    return times(L, througput);
}


struct Vertex {
    const Float3 position;
    const Float3 out_dir;

    const float transmittance; // 一個前と今の頂点の間の減衰率


    Vertex(const Float3 &aPosition, const Float3 &aOutDir, const float aTransmittance):
    position(aPosition), out_dir(aOutDir), transmittance(aTransmittance) {}
};


Float3 radiance_gdpt(Random& random, Random& random2, const SimpleVolume& volume, const rt::Ray& initial_ray, std::vector<Vertex>* vertices, float *path_pdf, Color *finalL) {
    const int kMaxScattering = 256;

    rt::Ray current_ray = initial_ray;
    
    Color througput(1, 1, 1);
    Color L(0, 0, 0);


    for (int scattering = 0; scattering < kMaxScattering; ++scattering) {
        // VolumeのBBoxにヒットするか否か
        float hitt0 = -1, hitt1 = -1;
        bool is_volume_hit = false;
        if (volume.check_intersect(current_ray, &hitt0, &hitt1)) {
            is_volume_hit = true;
        }
        
        rt::Hitpoint current_hp;
        if (!volume.inside(current_ray.org)) {
            if (is_volume_hit) {
                current_hp.distance = hitt0;
            }
        }

        if (!is_volume_hit) {
            L += outside(current_ray);
            break;
        }

        const Float3 new_org0 = current_ray.org + (hitt0 + 1e-3f) * current_ray.dir;
        const Float3 new_org1 = current_ray.org + (hitt1 + 1e-3f) * current_ray.dir;
        if (!volume.inside(current_ray.org)) {
            // Volumeの外だったら中からスタート
            current_ray.org = new_org0;
        }
        if (scattering == 0) {
            // 一個目の頂点を入れる
            vertices->push_back(Vertex(current_ray.org, current_ray.dir, 1));
        }
        const float tu = volume.next_event(current_ray, random);
        if (tu < 0) {
            // ボリューム外
            current_ray.org = new_org1;
            continue;
        }

        // 散乱した！

        // pdf_transmittance
        float transmittance = 0;
        const int kNumSample = 16;
        /*
        {
            float sum = 0;
            for (int i = 0; i < kNumSample; ++i) {
                const float length = tu;

                const float current = length / kNumSample * (i + random.next01());

                sum += volume.sample_transmittance(current_ray.org + current * current_ray.dir) / kNumSample;
            }

            transmittance = exp(-sum);
        }*/
        for (int i = 0; i < kNumSample; ++i) {
            const float X = volume.next_event(current_ray, random2);
            if (!(0 <= X && X < tu)) {
                transmittance += 1.0f / kNumSample;
            }
        }

        const Float3 next_pos = current_ray.org + tu * current_ray.dir;
        const bool cond2 = random.next01() < kAlbedo;
        if (!cond2) {
            *path_pdf = 0;
            througput = Float3(0, 0, 0); // absorp
            break;
        }

        // 位相関数でインポータンスサンプリング
        // 今回はとりあえず等方散乱としておく
        const Float3 next_dir = Sampling::uniformSphereSurface(random);
        const float phase = 1.0f / (4.0f * kPI);
        const float pdf_phase = 1.0f / (4.0f * kPI);
          
        /*
        througput *= 0.99f;
        *path_pdf *= 0.99f;
        */
        
        througput *= transmittance * phase * kAlbedo;
        *path_pdf *= transmittance * pdf_phase * kAlbedo;
        

        // 頂点追加    
        vertices->push_back(Vertex(next_pos, next_dir, transmittance));
        

        // throughtput = (phase / pdf) * throughtput;	
        current_ray = rt::Ray(next_pos, next_dir);
    }
    *finalL = L;

    return times(L, througput);
}


void render(const char *filename, const int width, const int height, const int num_sample_per_subpixel, const int num_thread, unsigned long long seed,
    const SimpleVolume& ovolume, 
    const rt::PerspectiveCamera& camera) {
    Image image(width, height);
    Timer timer;
    timer.begin();

    for (int iy = 0; iy < height; ++iy) {
        std::cout << "y: " << iy << std::endl;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(8)
        for (int ix = 0; ix < width; ++ix) {
            Random random((unsigned long long)(iy * width + ix));

            for (int ss = 0; ss < num_sample_per_subpixel; ++ss) {
                const float u = (ix + random.next01()) / width;
                const float v = (iy + random.next01()) / height;
                rt::Ray ray = camera.get_ray(u, v);


                
                /*
                float pdf = 1;
                Color L;
                std::vector<Vertex> vertices;
                Color contribution = radiance_gdpt(random, ovolume, ray, &vertices, &pdf, &L);
                if (pdf > 0)
                    image.at(ix, iy) += contribution / (float)pdf / num_sample_per_subpixel;
                */

                float pdf = 1;
                const Float3 contribution = radiance(random, ovolume, ray, &pdf);
                if (pdf > 0)
                    image.at(ix, iy) += contribution / (float)pdf / num_sample_per_subpixel;
            }
        }
    }

    std::cout << "Time: " << (timer.end() / 1000.0f) << " sec" << std::endl;

    HDROperator::save(filename, &image);
}

void render_gdpt(const char *filename, const int width, const int height, const int num_sample_per_subpixel, const int num_thread, unsigned long long seed,
    const SimpleVolume& ovolume, 
    const rt::PerspectiveCamera& camera) {
        
    std::vector<Color> coarse_image(width * height);
    std::vector<Color> debug_image(width * height);
    std::vector<Color> diff_image[4];
    for (int i = 0; i < 4; ++i)
        diff_image[i].resize(width * height);

    Timer timer;
    timer.begin();

    for (int iy = 0; iy < height; ++iy) {
        std::cout << "y: " << iy << std::endl;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(8)
        for (int ix = 0; ix < width; ++ix) {
            Random random((unsigned long long)(iy * width + ix));
            
            const int image_index = iy * width + ix;
            for (int ss = 0; ss < num_sample_per_subpixel; ++ss) {
                float sx = random.next01();
                float sy = random.next01();
                
                std::vector<Vertex> baseVertices;
                baseVertices.reserve(16);
                float basePathPDF = 1;
                Color baseContritbuion;
                Color finalL;

                float baseX = ix + sx;
                float baseY = iy + sy;
                const unsigned long long random2_seed = (unsigned long long)(iy * width + ix) * num_sample_per_subpixel + ss;
                rt::Ray base_ray; 
                {
                    const float u = baseX / width;
                    const float v = baseY / height;
                    base_ray = camera.get_ray(u, v);
                    baseContritbuion = radiance_gdpt(random, random, ovolume, base_ray, &baseVertices, &basePathPDF, &finalL);
                }
                
                if (basePathPDF > 0)
                    coarse_image[image_index] += baseContritbuion / basePathPDF;
                if (basePathPDF <= 0)
                    continue;
                // 4 sub-path
                const Float3 offset[4] = {
                    Float3( 1, 0, 0),
                    Float3(-1, 0, 0),
                    Float3( 0, 1, 0),
                    Float3( 0,-1, 0)
                };
                std::vector<Vertex> offsetVertices[4];
                for (int offsetDir = 0; offsetDir < 4; ++offsetDir) {
                    const float u = (baseX + offset[offsetDir].x) / width;
                    const float v = (baseY + offset[offsetDir].y) / height;
                    const rt::Ray offset_ray = camera.get_ray(u, v);
                    offsetVertices[offsetDir].reserve(16);

                    enum Result {
                        Naive,
                        Invertible
                    };

                    // Shift
                    rt::Ray current_ray = offset_ray;
                    Result result = Naive;

                    // baseが少なくとも2回は散乱している、有効なパスのときのみ考える
                    Color offsetContribution;
                    float J = 0;
                    if (basePathPDF > 0 && baseVertices.size() >= 3) {
                        const float initial_length = (baseVertices[1].position - base_ray.org).length();

                        const Float3 offset_dir0 = current_ray.dir;
                        const Float3 offset_1 = current_ray.org + initial_length * offset_dir0;
                        
                        if (!ovolume.inside(offset_1)) {
                            // offsetが外に出てしまったので、Naiveにやる
                        } else {
                            result = Invertible;
                            float hitt0 = -1, hitt1 = -1;
                            ovolume.check_intersect(current_ray, &hitt0, &hitt1);

                            const Float3 new_org0 = current_ray.org + (hitt0 + 1e-3f) * offset_dir0;
                            const Float3 new_org1 = current_ray.org + (hitt1 + 1e-3f) * offset_dir0;
                            if (!ovolume.inside(current_ray.org)) {
                                // Volumeの外だったら中からスタート
                                current_ray.org = new_org0;
                            }
                            const Float3 offset_0 = current_ray.org;



                            // Baseの3個目の頂点につなげる
                            const Float3 offset_2 = baseVertices[2].position;

                            // 0 - 1と 1- 2の間のtransmittanceを求める
                            const rt::Ray ray01(offset_0, normalize(offset_1 - offset_0));
                            const rt::Ray ray12(offset_1, normalize(offset_2 - offset_1));
                            const float length01 = (offset_0 - offset_1).length();
                            const float length12 = (offset_1 - offset_2).length();
                            float transmittance01 = 0;
                            float transmittance12 = 0;
                            const int kNumSample = 16;

                            
                            /*
                            {
                                float sum = 0;
                                for (int i = 0; i < kNumSample; ++i) {
                                    const float current = length01 / kNumSample * (i + random.next01());
                                    sum += ovolume.sample_transmittance(ray01.org + current * ray01.dir) / kNumSample;
                                }
                                transmittance01 = exp(-sum);
                            }
                            {
                                float sum = 0;
                                for (int i = 0; i < kNumSample; ++i) {
                                    const float current = length12 / kNumSample * (i + random.next01());
                                    sum += ovolume.sample_transmittance(ray12.org + current * ray12.dir) / kNumSample;
                                }
                                transmittance12 = exp(-sum);
                            }
                            */
                            for (int i = 0; i < kNumSample; ++i) {
                                const float X = ovolume.next_event(ray01, random);
                                if (!(0 <= X && X < length01)) {
                                    transmittance01 += 1.0f / kNumSample;
                                }

                            }
                            /*
                                if (baseVertices[1].transmittance != transmittance01)
                                    printf("%f %f\n", baseVertices[1].transmittance, transmittance01);
                                */
                            for (int i = 0; i < kNumSample; ++i) {
                                const float X = ovolume.next_event(ray12, random);
                                if (!(0 <= X && X < length12)) {
                                    transmittance12 += 1.0f / kNumSample;
                                }
                            }

                            // ヤコビアン計算
                            J = (baseVertices[2].position - baseVertices[1].position).lengthSquared() /
                                (offset_2 - offset_1).lengthSquared();

                            // 残りのTransmittanceはBaseと同じ！
                            float total_transmittance = transmittance01 * transmittance12;
                            for (int i = 3; i < baseVertices.size(); ++i)
                                total_transmittance *= baseVertices[i].transmittance;

                            // その他の要素
                            const float phase = 1.0f / (4.0f * kPI);        
                            const int numScattering = baseVertices.size() - 1;

                            const float throughput = total_transmittance * pow(phase * kAlbedo, numScattering);

                            offsetContribution = throughput * finalL;

                        }
                    }
                    Color accum;
                    if (result == Naive) {
                        /*
                        if (iy > 100)
                        printf("%f %d\n", basePathPDF, baseVertices.size());
                        */

                        float pdf = 1;
                        const Float3 contribution = radiance(random, ovolume, offset_ray, &pdf);
                        if (pdf > 0) {
                            accum = baseContritbuion / basePathPDF - contribution / pdf;
                        } else {
                            accum = baseContritbuion / basePathPDF; 
                        }


                        /*
                        if (ix == 320 && iy == 240) {
                            printf("%f %f %f %f\n", baseContritbuion.x, basePathPDF, contribution.x, pdf);
                            std::cout << accum << std::endl;
                        }
                        */
                    } else {
                        // printf("%f %f %f\n", baseContritbuion.x, offsetContribution.x, J);
                        accum = 0.5f * (baseContritbuion - offsetContribution * J) / basePathPDF;
                        
                    }
                    
                    // 前方差分 or 後方差分
                    if (offsetDir == 1 || offsetDir == 2)
                        diff_image[offsetDir][image_index] += accum;
                    else
                        diff_image[offsetDir][image_index] += -accum;
                }
            }
            coarse_image[image_index] = coarse_image[image_index] / num_sample_per_subpixel;
            for (int i = 0; i < 4; ++i)
                diff_image[i][image_index] = diff_image[i][image_index] / num_sample_per_subpixel;
        }
    }

    std::cout << "Time: " << (timer.end() / 1000.0f) << " sec" << std::endl;
    
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
    }
}

int main(int argc, char *argv[]) {
    using namespace hstd;
    const int width = 640;
    const int height = 480;
    const int sample = 16;


    rt::PerspectiveCamera camera(Float3(-0.2, 1.3, 2.3), Float3(0, 1.2, 0), normalize(Float3(0.05, 1, 0)), 1, (float)height / width, 1.5f, 10000.f);
  
//    rt::PerspectiveCamera camera(Float3(-0.2, 5.3, 5.3), Float3(0, 1.2, 0), normalize(Float3(0.05, 1, 0)), 1, (float)height / width, 1.5f, 10000.f);

    SimpleVolume ovolume("volume.float", 30.0f, kAlbedo,
        Float3(-128 * 0.02 / 2, 0.01, -128 * 0.02 / 2), Float3(0.02, 0.02, 0.02));
    
   //render("vol/reference_8192s.hdr", width, height, 8192, 8, 0, ovolume, camera);

    
    int sample_table[4] = {
        4,
        16,
        64,
        256
    };
    
    
    for (int i = 0; i < 4; ++i) {
        Timer timer;
        timer.begin();
        char buf[256];
        sprintf(buf, "newvol/result2_%03ds", sample_table[i]);
        render_gdpt(buf, width, height, sample_table[i], 8, 0, ovolume, camera);
        
        sprintf(buf, "newvol/result2_%03ds_log.txt", sample_table[i]);
        FILE *fp = fopen(buf, "wt");
        fprintf(fp, "%f sec", timer.end() / 1000.0f);
        fclose(fp);
    }
}