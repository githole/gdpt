#include <vector>

#include "image.h"
#include "vec.h"
#include "timer.h"

typedef gemspt::Vec COLOR;

void load_bin_file(const char* filename, std::vector<COLOR> *image) {
    FILE *fp = fopen(filename, "rb");
    if (fp != NULL) {
        int width, height;

        fread(&width, sizeof(int), 1, fp);
        fread(&height, sizeof(int), 1, fp);

        std::cout << width << " " << height << std::endl;

        image->resize(width * height);

        fread(&(*image)[0], sizeof(COLOR), width * height, fp);

        fclose(fp);
    }
}

void save_hdr_file(const std::string &filename, const COLOR *image, const int width, const int height) {
    hstd::Image nimage(width, height);

    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            nimage.at(ix, iy) = image[iy * width + ix];
        }
    }

    hstd::HDROperator::save(filename.c_str(), &nimage, true);
}

using namespace gemspt;

COLOR sample(std::vector<COLOR>& data, int width, int height, int x, int y) {
    if (x < 0)
        x = 0;
    if (x >= width)
        x = width - 1;
    if (y < 0)
        y = 0;
    if (y >= height)
        y = height - 1;

    return data[y * width + x];
}
double sample(std::vector<double>& data, int width, int height, int x, int y) {
    if (x < 0 || y < 0 || width <= x || height <= y)
        return 0;
    return data[y * width + x];
}

double L2Reconstruction(const char* filename, const char* outputFilename, const char* referenceFileName) {
    std::vector<COLOR> coarse_image;
    std::vector<COLOR> diff[4];
    const int width = 640;
    const int height = 480;

    hstd::Image image(width, height);
    hstd::HDROperator::load(referenceFileName, &image);
    
    char name[4][256] = {
        "__1_0",
        "_-1_0",
        "__0_1",
        "__0-1"
    };
    char buf[256];
    sprintf(buf, "%s_coarse.bin", filename);
    load_bin_file(buf, &coarse_image);
    for (int i = 0; i < 4; ++i) {
        sprintf(buf, "%s_%s.bin", filename, name[i]);
        load_bin_file(buf, &diff[i]);
    }
//    fourierSolve
    
    std::vector<COLOR> out[2];
    out[0].resize(width * height);
    out[1].resize(width * height);
    std::vector<COLOR> dx(width * height);
    std::vector<COLOR> dy(width * height);
    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            int i = iy * width + ix;
            if (ix == 0)
                dx[i] = diff[1][i];
            else
                dx[i] = diff[1][i] + diff[0][i - 1];
            if (iy == 0)
                dy[i] = diff[2][i];
            else
                dy[i] = (diff[2][i] + diff[3][i - width]);
        }
    }
    
    // L2-Norm Reconstruction
    const int kMaxGaussSidel = 500;
    const double kAlpha = 0.2;
    const Vec kLuminance(0.1769, 0.8124, 0.0107);

    for (int i = 0; i < kMaxGaussSidel; ++i) {
        const int current = i % 2;
        const int prev = 1 - current;
        std::cout << i << " ";

        for (int i = 0; i < width * height; ++i)
            out[current][i] = out[prev][i];

        for (int iy = 0; iy < height; ++iy) {
            for (int ix = 0; ix < width; ++ix) {
                const int idx = iy * width + ix;

                out[current][idx] = (1.0 / (4.0 + kAlpha * kAlpha)) * (
                    kAlpha * kAlpha * sample(coarse_image, width, height, ix, iy)
                    + sample(out[current], width, height, ix - 1, iy)
                    + sample(out[current], width, height, ix + 1, iy) 
                    + sample(out[current], width, height, ix, iy - 1)
                    + sample(out[current], width, height, ix, iy + 1)
                    + sample(dx, width, height, ix, iy)
                    - sample(dx, width, height, ix + 1, iy)
                    + sample(dy, width, height, ix, iy)
                    - sample(dy, width, height, ix, iy + 1));
            }
        }

        // Œë·
        double err = 0;
        for (int i = 0; i < width * height; ++i)
            err += sqrt(pow(dot(kLuminance, out[current][i] - out[prev][i]), 2.0)) / (width * height);

        // std::cout << err << " ";

        if (err <= 1e-5 || i == kMaxGaussSidel - 1) {
            char buf[256];
            sprintf(buf, "%s.hdr", outputFilename);
            save_hdr_file(buf, &out[current][0], width, height);

            double refErr = 0;
            for (int iy = 0; iy < height; ++iy) {
                for (int ix = 0; ix < width; ++ix) {
                    refErr += sqrt(pow(dot(kLuminance, (out[current][iy * width + ix] - image.at(ix, iy))), 2.0)) / (width * height);
                }
            }

            return refErr;
        }
    }

    return 0;
}

double L1Reconstruction(const char* filename, const char* outputFilename, const char* referenceFileName) {
    std::vector<COLOR> coarse_image;
    std::vector<COLOR> diff[4];
    const int width = 640;
    const int height = 480;
    
    hstd::Image image(width, height);
    hstd::HDROperator::load(referenceFileName, &image);
    
    char name[4][256] = {
        "__1_0",
        "_-1_0",
        "__0_1",
        "__0-1"
    };
    char buf[256];
    sprintf(buf, "%s_coarse.bin", filename);
    load_bin_file(buf, &coarse_image);
    for (int i = 0; i < 4; ++i) {
        sprintf(buf, "%s_%s.bin", filename, name[i]);
        load_bin_file(buf, &diff[i]);
    }
//    fourierSolve
    
    std::vector<COLOR> out[2];
    out[0].resize(width * height);
    out[1].resize(width * height);
    std::vector<COLOR> dx(width * height);
    std::vector<COLOR> dy(width * height);
    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            int i = iy * width + ix;
            if (ix == 0)
                dx[i] = diff[1][i];
            else
                dx[i] = diff[1][i] + diff[0][i - 1];
            if (iy == 0)
                dy[i] = diff[2][i];
            else
                dy[i] = (diff[2][i] + diff[3][i - width]);
        }
    }

    // L1-Norm Reconstruction
    std::vector<double> coarse_weight(width * height, 1);
    std::vector<double> dx_weight(width * height, 1);
    std::vector<double> dy_weight(width * height, 1);

    const int kMaxL1 = 10;
    const int kMaxGaussSidel = 500;
    const double kAlpha = 0.2;
    const Vec kLuminance(0.1769, 0.8124, 0.0107);
    const double kEPS = 1e-8;

    for (int l1 = 0; l1 < kMaxL1; ++l1) {
        std::cout << l1 << std::endl;
        // Weighted L2-Norm Reconstruction
        int result = 0;
        for (int i = 0; i < kMaxGaussSidel; ++i) {
            const int current = i % 2;
            const int prev = 1 - current;

            for (int i = 0; i < width * height; ++i)
                out[current][i] = out[prev][i];

            for (int iy = 0; iy < height; ++iy) {
                for (int ix = 0; ix < width; ++ix) {
                    const int idx = iy * width + ix;

                    const double Wij_coarse = sample(coarse_weight, width, height, ix, iy);
                    const double Wij_dx = sample(dx_weight, width, height, ix, iy);
                    const double Wi1j_dx = sample(dx_weight, width, height, ix + 1, iy);
                    const double Wij_dy = sample(dy_weight, width, height, ix, iy);
                    const double Wiji_dy = sample(dy_weight, width, height, ix, iy + 1);

                    out[current][idx] =
                        (1.0 / (kAlpha * kAlpha * Wij_coarse + Wij_dx + Wi1j_dx + Wij_dy + Wiji_dy) * (
                        kAlpha * kAlpha * Wij_coarse * sample(coarse_image, width, height, ix, iy)
                        + Wij_dx  * (sample(out[current], width, height, ix - 1, iy) + sample(dx, width, height, ix, iy))
                        + Wi1j_dx * (sample(out[current], width, height, ix + 1, iy) - sample(dx, width, height, ix + 1, iy))
                        + Wij_dy  * (sample(out[current], width, height, ix, iy - 1) + sample(dy, width, height, ix, iy))
                        + Wiji_dy * (sample(out[current], width, height, ix, iy + 1) - sample(dy, width, height, ix, iy + 1))
                        ));
                }
            }

            // Œë·
            double err = 0;
            for (int i = 0; i < width * height; ++i)
                err += sqrt(pow(dot(kLuminance, out[current][i] - out[prev][i]), 2.0)) / (width * height);

            if (err <= 1e-5 || i == kMaxGaussSidel - 1) {
                result = current;
                break;
            }
        }
        // Output
        
        if (l1 == kMaxL1 - 1) {
            char buf[256];
            sprintf(buf, "%s.hdr", outputFilename);
            save_hdr_file(buf, &out[result][0], width, height);

            double refErr = 0;
            for (int iy = 0; iy < height; ++iy) {
                for (int ix = 0; ix < width; ++ix) {
                    refErr += sqrt(pow(dot(kLuminance, (out[result][iy * width + ix] - image.at(ix, iy))), 2.0)) / (width * height);
                }
            }

            return refErr;
        }

        // Reweight
        const double P = 1;
        for (int iy = 0; iy < height; ++iy) {
            for (int ix = 0; ix < width; ++ix) {
                coarse_weight[iy * width + ix] =
                    1.0 / std::max(std::abs(dot(kLuminance, (coarse_image[iy * width + ix] - out[result][iy * width + ix]))), kEPS); 
                
                dx_weight[iy * width + ix] =
                    1.0 / std::max(pow(std::abs(dot(kLuminance, (dx[iy * width + ix] - (sample(out[result], width, height, ix, iy) - sample(out[result], width, height, ix - 1, iy))))), P), kEPS); 
                dy_weight[iy * width + ix] =
                    1.0 / std::max(pow(std::abs(dot(kLuminance, (dy[iy * width + ix] - (sample(out[result], width, height, ix, iy) - sample(out[result], width, height, ix, iy + 1))))), P), kEPS); 
            }
        }

    }

    return 0;
}

int main() {
#if 0
    int sample_table[] = {
        32
    };

    char buf[256];
    char buf2[256];
    for (int i = 0; i < 1; ++i) {
        {
            Timer timer;
            timer.begin();
            sprintf(buf, "../gdpt/output2/result2_%03ds", sample_table[i]); 
            sprintf(buf2, "output2/%03ds_l2", sample_table[i]);
            const double refErr = L2Reconstruction(buf, buf2, "../gdpt/output/reference_8192s.hdr");
            printf("Err: %f\n", refErr);
        }
        
        {
            Timer timer;
            timer.begin();
            sprintf(buf, "../gdpt/output2/result2_%03ds", sample_table[i]); 
            sprintf(buf2, "output2/%03ds_l1", sample_table[i]);
            const double refErr = L1Reconstruction(buf, buf2, "../gdpt/output/reference_8192s.hdr");
            printf("Err: %f\n", refErr);
        }
    }
#endif
    
#if 0
    // coarse vs L1, L2
    int sample_table[] = {
        4,
        16,
        64,
        256,
    };
    const char *dir = "lfreq";
    char buf[256];
    char buf2[256];
    char buf3[256];
    const Vec kLuminance(0.1769, 0.8124, 0.0107);

    
    for (int i = 0; i < 4; ++i) {
        {
            Timer timer;
            timer.begin();
            sprintf(buf, "../gdpt/%s/result2_%03ds_coarse.hdr", dir, sample_table[i]); 
            sprintf(buf3, "../gdpt/%s/reference_8192s.hdr", dir);

            hstd::Image l2Image;
            hstd::Image refImage;

            hstd::HDROperator::load(buf, &l2Image);
            hstd::HDROperator::load(buf3, &refImage);
            
            const int width = l2Image.width();
            const int height = l2Image.height();
            double refErr = 0;
            for (int iy = 0; iy < height; ++iy) {
                for (int ix = 0; ix < width; ++ix) {
                    refErr += sqrt(pow(dot(kLuminance, (l2Image.at(ix, iy)- refImage.at(ix, iy))), 2.0)) / (width * height);
                }
            }

            const double sec = timer.end() / 1000.0f;
            sprintf(buf, "%s/%03ds_corase_log.txt", dir, sample_table[i]);
            FILE *fp = fopen(buf, "wt");
            fprintf(fp, "%f sec\nerr: %f", sec, refErr);
            fclose(fp);
        }
    }
#endif
#if 0
    // volume
    int sample_table[] = {
        4,
        16,
        64,
        256,
    };

    const char *dir = "newvol";

    char buf[256];
    char buf2[256];
    char buf3[256];

    for (int i = 2; i < 4; ++i) {
        {
            Timer timer;
            timer.begin();
            sprintf(buf, "../volume-gdpt/%s/result2_%03ds", dir, sample_table[i]); 
            sprintf(buf2, "%s/%03ds_l2", dir, sample_table[i]);
            sprintf(buf3, "../volume-gdpt/%s/reference_8192s.hdr", dir);
            const double refErr = L2Reconstruction(buf, buf2, buf3);

            const double sec = timer.end() / 1000.0f;
            sprintf(buf, "%s/%03ds_l2_log.txt", dir, sample_table[i]);
            FILE *fp = fopen(buf, "wt");
            fprintf(fp, "%f sec\nerr: %f", sec, refErr);
            fclose(fp);
        }
        
        {
            Timer timer;
            timer.begin();
            sprintf(buf, "../volume-gdpt/%s/result2_%03ds", dir, sample_table[i]); 
            sprintf(buf2, "%s/%03ds_l1", dir, sample_table[i]);
            sprintf(buf3, "../volume-gdpt/%s/reference_8192s.hdr", dir);
            const double refErr = L1Reconstruction(buf, buf2, buf3);

            const double sec = timer.end() / 1000.0f;
            sprintf(buf, "%s/%03ds_l1_log.txt", dir, sample_table[i]);
            FILE *fp = fopen(buf, "wt");
            fprintf(fp, "%f sec\nerr: %f", sec, refErr);
            fclose(fp);
        }
    }
#endif

#if 1
    int sample_table[] = {
        4,
        16,
        64,
        256,
    };

    const char *dir = "output";

    char buf[256];
    char buf2[256];
    char buf3[256];

    for (int i = 0; i < 4; ++i) {
        {
            Timer timer;
            timer.begin();
            sprintf(buf, "../gdpt/%s/result2_%03ds", dir, sample_table[i]); 
            sprintf(buf2, "%s/%03ds_l2", dir, sample_table[i]);
            sprintf(buf3, "../gdpt/%s/reference_8192s.hdr", dir);
            const double refErr = L2Reconstruction(buf, buf2, buf3);

            const double sec = timer.end() / 1000.0f;
            sprintf(buf, "%s/%03ds_l2_log.txt", dir, sample_table[i]);
            FILE *fp = fopen(buf, "wt");
            fprintf(fp, "%f sec\nerr: %f", sec, refErr);
            fclose(fp);
        }
        
        {
            Timer timer;
            timer.begin();
            sprintf(buf, "../gdpt/%s/result2_%03ds", dir, sample_table[i]); 
            sprintf(buf2, "%s/%03ds_l1", dir, sample_table[i]);
            sprintf(buf3, "../gdpt/%s/reference_8192s.hdr", dir);
            const double refErr = L1Reconstruction(buf, buf2, buf3);

            const double sec = timer.end() / 1000.0f;
            sprintf(buf, "%s/%03ds_l1_log.txt", dir, sample_table[i]);
            FILE *fp = fopen(buf, "wt");
            fprintf(fp, "%f sec\nerr: %f", sec, refErr);
            fclose(fp);
        }
    }
#endif

    return 0;
}