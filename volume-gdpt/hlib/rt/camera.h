#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "../vec3.h"
#include "ray.h"

namespace hstd {

namespace rt {

// ÉJÉÅÉâ
// TODO: ÉÅÉìÉoïœêîÇÃñºëOÇ™Ç§ÇÒÇ±
class Camera {
protected:
    Float3 camera_direction_;
    Float3 camera_lookat_;
    Float3 camera_position_;

    Float3 screen_center_;
    Float3 screen_side_;
    Float3 screen_up_;

    float screen_width_;
    float screen_height_;

    float camera_near_, camera_far_;
    
public:
    Camera(Float3& camera_position, Float3& camera_lookat, Float3& camera_up, float screen_width, float screen_height, float camera_near, float camera_far) : 
    camera_position_(camera_position),
    camera_lookat_(camera_lookat),
    screen_width_(screen_width),
    screen_height_(screen_height),
    camera_near_(camera_near),
    camera_far_(camera_far) {
    }

    float screen_width() const { return screen_width_; }
    float screen_height() const { return screen_height_; }

    Float3 camera_direction() const { return camera_direction_; }
    Float3 screen_up() const { return screen_up_; }
    Float3 screen_side() const { return screen_side_; }

    virtual rt::Ray get_ray(const float u, const float v) const = 0;
};

class PerspectiveCamera : public Camera {
private:
public:
    PerspectiveCamera(Float3& camera_position, Float3& camera_lookat, Float3& camera_up, float screen_width, float screen_height, float camera_near, float camera_far) : 
        Camera(camera_position, camera_lookat, camera_up, screen_width, screen_height, camera_near, camera_far) {
        camera_direction_ = normalize(camera_lookat_ - camera_position_);
        
        const float aspect = (float)screen_width_ / screen_height_;
        screen_side_ = aspect * cross(camera_direction_, camera_up);
        screen_up_ = normalize(cross(camera_direction_, screen_side_));
        screen_center_ = camera_position + camera_near_ * camera_direction_;
    }
    
    rt::Ray get_ray(const float u, const float v) const {
        const float sx = (u * 2.0f) - 1.0f;
        const float sy = (v * 2.0f) - 1.0f; 
        const Float3 pos_on_screen = screen_center_ + sx * screen_side_ + sy * screen_up_;
        rt::Ray ray(camera_position_, normalize(pos_on_screen - camera_position_));

        return ray;
    }
};

} // namespace rt

} // namespace hstd


#endif // _CAMERA_H_