#ifndef _HDR_H_
#define _HDR_H_

#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "image.h"

namespace gemspt {


void save_hdr_file(const std::string &filename, const Color *image, const int width, const int height) {
    hstd::Image nimage(width, height);

    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            nimage.at(ix, iy) = image[iy * width + ix];
        }
    }

    hstd::HDROperator::save(filename.c_str(), &nimage, true);
}


void save_hdr_file_positive(const std::string &filename, const Color *image, const int width, const int height) {
    hstd::Image nimage(width, height);

    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            nimage.at(ix, iy) = image[iy * width + ix];

            nimage.at(ix,iy).x = nimage.at(ix, iy).x < 0 ? 0 : nimage.at(ix, iy).x;
            nimage.at(ix,iy).y = nimage.at(ix, iy).y < 0 ? 0 : nimage.at(ix, iy).y;
            nimage.at(ix,iy).z = nimage.at(ix, iy).z < 0 ? 0 : nimage.at(ix, iy).z;
        }
    }

    hstd::HDROperator::save(filename.c_str(), &nimage, true);
}

void save_hdr_file_negative(const std::string &filename, const Color *image, const int width, const int height) {
    hstd::Image nimage(width, height);

    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            nimage.at(ix, iy) = image[iy * width + ix];
            
            nimage.at(ix,iy).x = nimage.at(ix, iy).x > 0 ? 0 : -nimage.at(ix, iy).x;
            nimage.at(ix,iy).y = nimage.at(ix, iy).y > 0 ? 0 : -nimage.at(ix, iy).y;
            nimage.at(ix,iy).z = nimage.at(ix, iy).z > 0 ? 0 : -nimage.at(ix, iy).z;
        }
    }

    hstd::HDROperator::save(filename.c_str(), &nimage, true);
}
};

#endif
