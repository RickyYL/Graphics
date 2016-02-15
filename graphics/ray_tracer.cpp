// c++ -o raytracer -O3 -Wall raytracer.cpp

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>

#if defined __linux__ || defined __APPLE__
#else
#define M_PI     3.141592653589793
#define INFINITY 1e8
#endif

// three dimentional vector
template<typename T>
class Vec3 {
public:
    T x, y, z;
    Vec3(): x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx): x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz): x(xx), y(yy), z(zz) {}
    Vec3<T>  operator+ (const Vec3<T> &v) const { return Vec3<T>(x+v.x, y+v.y, z+v.z); }
    Vec3<T>  operator- (const Vec3<T> &v) const { return Vec3<T>(x-v.x, y-v.y, z-v.z); }
    Vec3<T>  operator- ()                 const { return Vec3<T>(-x,    -y,    -z); }
    Vec3<T>  operator* (const T &f)       const { return Vec3<T>(x*f,   y*f,   z*f); }
    Vec3<T>  operator* (const Vec3<T> &v) const { return Vec3<T>(x*v.x, y*v.y, z*v.z); }
    Vec3<T>  operator/ (const T &f)       const { return Vec3<T>(x/f,   y/f,   z/f); }
    Vec3<T>  operator/ (const Vec3<T> &v) const { return Vec3<T>(x/v.x, y/v.y, z/v.z); }
    Vec3<T>& operator+=(const Vec3<T> &v)       { x += v.x, y += v.y, z += v.z; return *this; }
    Vec3<T>& operator-=(const Vec3<T> &v)       { x -= v.x, y -= v.y, z -= v.z; return *this; }
    Vec3<T>& operator*=(const T &f)             { x *= f,   y *= f,   z *= f;   return *this; }
    Vec3<T>& operator*=(const Vec3<T> &v)       { x *= v.x, y *= v.y, z *= v.z; return *this; }
    Vec3<T>& operator/=(const T &f)             { x /= f,   y /= f,   z /= f;   return *this; }
    Vec3<T>& operator/=(const Vec3<T> &v)       { x /= v.x, y /= v.y, z /= v.z; return *this; }
    T        dot(const Vec3<T> &v)        const { return x * v.x + y * v.y + z * v.z; }
    T        length2()                    const { return x * x + y * y + z * z; }
    T        length()                     const { return sqrt(length2()); }
    Vec3<T>& normalize()                        { if (length2() > 0) *this /= length(); return *this; }
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v) {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

// define vector of reals
typedef Vec3<float> Vec3f;

class Sphere {
public:
    
    Vec3f center;
    float radius, radius2;
    Vec3f surfaceColor, emissionColor;
    float transparency, reflectivity;
    
    Sphere(const Vec3f &center,
           const float &radius,
           const Vec3f &surfaceColor,
           const float &reflectivity = 0,
           const float &transparency = 0,
           const Vec3f &emissionColor = 0)
    : center(center),
    radius(radius),
    radius2(radius * radius),
    surfaceColor(surfaceColor),
    emissionColor(emissionColor),
    transparency(transparency),
    reflectivity(reflectivity)
    {}
    
    // given a ray, determine if the ray hits the sphere
    bool intersect(const Vec3f &rayOrigin, const Vec3f &rayDirection, float &t0, float &t1) const {
        
        Vec3f l = center - rayOrigin;
        float tca = l.dot(rayDirection);
        if (tca < 0) return false;
        float d2 = l.dot(l) - tca * tca;
        if (d2 > radius2) return false;
        
        float thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;
        return true;
    }
};

// This variable controls the maximum recursion depth
#define MAX_RAY_DEPTH 5

float mix(const float &a, const float &b, const float &mix) {
    return b * mix + a * (1 - mix);
}

// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
//
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.

Vec3f trace(const Vec3f &rayOrigin, const Vec3f &rayDirection, const std::vector<Sphere> &spheres, const int &depth) {
    
    // if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
    
    float minDistance = INFINITY;
    const Sphere* sphere = nullptr;
    
    // find intersection of this ray with the sphere in the scene
    for (unsigned i = 0; i < spheres.size(); ++i) {
        float t0 = INFINITY, t1 = INFINITY;
        if (spheres[i].intersect(rayOrigin, rayDirection, t0, t1)) {
            if (t0 < 0)
                t0 = t1;
            if (t0 < minDistance) {
                minDistance = t0;
                sphere = &spheres[i];
            }
        }
    }
    
    // if there's no intersection return black or background color
    if (!sphere)
        return Vec3f(2);
    
    Vec3f surfaceColor = 0; // color of the ray/surface of the object intersected by the ray
    Vec3f pointHit  = rayOrigin + rayDirection * minDistance;
    Vec3f normalHit = pointHit - sphere->center;
    normalHit.normalize();
    
    // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
    float bias = 1e-4; // add some bias to the point from which we will be tracing
    bool inside = false;
    if (rayDirection.dot(normalHit) > 0) {
        normalHit = -normalHit;
        inside = true;
    }
    if ((sphere->transparency > 0 || sphere->reflectivity > 0) && depth < MAX_RAY_DEPTH) {
        float facingratio = -rayDirection.dot(normalHit);
        // change the mix value to tweak the effect
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
        // compute reflection direction (not need to normalize because all vectors
        // are already normalized)
        Vec3f refldir = rayDirection - normalHit * 2 * rayDirection.dot(normalHit);
        refldir.normalize();
        Vec3f reflection = trace(pointHit + normalHit * bias, refldir, spheres, depth + 1);
        Vec3f refraction = 0;
        // if the sphere is also transparent compute refraction ray (transmission)
        if (sphere->transparency) {
            float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
            float cosi = -normalHit.dot(rayDirection);
            float k = 1 - eta * eta * (1 - cosi * cosi);
            Vec3f refrdir = rayDirection * eta + normalHit * (eta *  cosi - sqrt(k));
            refrdir.normalize();
            refraction = trace(pointHit - normalHit * bias, refrdir, spheres, depth + 1);
        }
        // the result is a mix of reflection and refraction (if the sphere is transparent)
        surfaceColor = (reflection * fresneleffect +
                        refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
    } else {
        // it's a diffuse object, no need to raytrace any further
        for (unsigned i = 0; i < spheres.size(); ++i) {
            if (spheres[i].emissionColor.x > 0) {
                // this is a light
                Vec3f transmission = 1;
                Vec3f lightDirection = spheres[i].center - pointHit;
                lightDirection.normalize();
                for (unsigned j = 0; j < spheres.size(); ++j) {
                    if (i != j) {
                        float t0, t1;
                        if (spheres[j].intersect(pointHit + normalHit * bias, lightDirection, t0, t1)) {
                            transmission = 0;
                            break;
                        }
                    }
                }
                surfaceColor += sphere->surfaceColor * transmission *
                std::max(float(0), normalHit.dot(lightDirection)) * spheres[i].emissionColor;
            }
        }
    }
    
    return surfaceColor + sphere->emissionColor;
}

// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of the
// sphere at the intersection point, else we return the background color.

void render(const std::vector<Sphere> &spheres) {
    
    unsigned width = 1920, height = 1080;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.0f);
    
    // Trace rays
    for (unsigned y = 0; y < height; ++y) {
        for (unsigned x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -1);
            raydir.normalize();
            *pixel = trace(Vec3f(0), raydir, spheres, 0);
        }
    }
    
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i)
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
               (unsigned char)(std::min(float(1), image[i].y) * 255) <<
               (unsigned char)(std::min(float(1), image[i].z) * 255);
    ofs.close();
    delete [] image;
}


int main(int argc, char **argv) {
    
    srand48(13);
    std::vector<Sphere> spheres;
    
    // position, radius, surface color, reflectivity, transparency, emission color
    spheres.push_back(Sphere(Vec3f( 0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
    spheres.push_back(Sphere(Vec3f( 0.0,      0, -20),     6, Vec3f(1.00, 0.32, 0.36), 1, 0.5));
    spheres.push_back(Sphere(Vec3f( 5.0,     -1, -15),     2, Vec3f(0.90, 0.76, 0.46), 1, 0.0));
    spheres.push_back(Sphere(Vec3f( 5.0,      0, -25),     3, Vec3f(0.65, 0.77, 0.97), 1, 0.0));
    spheres.push_back(Sphere(Vec3f(-5.5,      0, -15),     4, Vec3f(0.60, 0.30, 0.80), 1, 0.0));
    
    // light
    spheres.push_back(Sphere(Vec3f( 0.0,     20, -30),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
    
    render(spheres);
    
    return 0;
}