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
    Vec3<T>  cross(const Vec3<T> &v)      const { return Vec3<T>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); }
    T        length2()                    const { return x * x + y * y + z * z; }
    T        length()                     const { return sqrt(length2()); }
    Vec3<T>& normalize()                        { if (length2() > 0) *this /= length(); return *this; }
    
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v) {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

// define vector of reals
typedef Vec3<float> vec3f;

class Ray {
public:
    vec3f origin, direction;
    Ray(vec3f _origin = 0, vec3f _direction = 0):
        origin(_origin), direction(_direction) { }
};

class Surface {
public:
    virtual bool isHit(const vec3f &rayorig, const vec3f &raydir, float t0, float t1) const;
};

class Triangle {
public:
    
    vec3f a, b, c;
    vec3f surfaceColor, emissionColor;
    
    Triangle(const vec3f &_a,
             const vec3f &_b,
             const vec3f &_c,
             const vec3f &_surfaceColor,
             const vec3f &_emissionColor = 0) :
    a(_a), b(_b), c(_c),
    surfaceColor(_surfaceColor),
    emissionColor(_emissionColor) { }
    
    bool isHit(const Ray ray, float t0, float &t1, vec3f &hitPoint, vec3f &hitNormal) const {
        
        vec3f rayorig = ray.origin;
        vec3f raydirc = ray.direction;
        
        float d = a.x - b.x, e = a.y - b.y, f = a.z - b.z;
        float g = a.x - c.x, h = a.y - c.y, i = a.z - c.z;
        float j = raydirc.x, k = raydirc.y, l = raydirc.z;
        float m = a.x - rayorig.x, n = a.y - rayorig.y, o = a.z - rayorig.z;
        float M = d*(h*l-k*i) + e*(j*i-g*l) + f*(g*k-h*j);
        
        float t = (- i*(d*n-m*e) - h*(m*f-d*o) - g*(e*o-n*f)) / M;
        if (t < t0 || t > t1)
            return false;
        float gamma = (l*(d*n-m*e) + k*(m*f-d*o) + j*(e*o-n*f)) / M;
        if (gamma < 0 || gamma > 1)
            return false;
        float beta  = (m*(h*l-k*i) + n*(j*i-g*l) + o*(g*k-h*j)) / M;
        if (beta < 0 || beta > 1 - gamma)
            return false;
        
        hitPoint = (ray.origin + ray.direction * t);
        hitNormal = (hitPoint - a).cross(hitPoint - b);
        hitPoint.normalize();
        hitPoint.normalize();
        t1 = t;
        
        return true;
    }
    
};

class Sphere {
public:
    
    vec3f center;
    float radius, radius2;
    vec3f surfaceColor, emissionColor;
    
    Sphere(const vec3f &_center,
           const float &_radius,
           const vec3f &_surfaceColor,
           const vec3f &_emissionColor = 0) :
        center(_center),
        radius(_radius),
        radius2(_radius * _radius),
        surfaceColor(_surfaceColor),
        emissionColor(_emissionColor) { }
    
    bool isHit(const Ray ray, float t0, float &t1, vec3f &hitPoint, vec3f &hitNormal) const {

        vec3f v = center - ray.origin;
        float vd = v.dot(ray.direction);
        float discriminate2 = radius2 - v.dot(v) + vd * vd;
        float intersection = vd - sqrt(discriminate2);
        
        if (discriminate2 < 0)
            return false;
        if (intersection < t0 || intersection > t1)
            return false;
        
        hitPoint = (ray.origin + ray.direction * intersection).normalize();
        hitNormal = (ray.origin + ray.direction * intersection - center).normalize();
        t1 = intersection;
        
        return true;
    }
};

vec3f tracer(const Ray &ray, const std::vector<Sphere> &objects, const vec3f &light, const vec3f &intensity) {
    
    const Sphere *closest = nullptr;
    float minDist = INFINITY;
    vec3f hitPoint, hitNormal;
    
    // search for the first object intersecting with the ray
    for (unsigned i = 0; i < objects.size(); i++) {
        if (objects[i].isHit(ray, 0, minDist, hitPoint, hitNormal))
            closest = &objects[i];
    }
    
    // nothing intersects with the ray
    if (closest == nullptr)
        return 0;
    
    // check if the object is in shadow
    bool inShadow = false;
    
    
    
    
    // shading
    vec3f l = (light - hitPoint).normalize();
    vec3f h = ((ray.origin - ray.direction) + l).normalize();
                          return
/* ambient     shading */ closest->surfaceColor * intensity * 0.1 +
/* Lambertian  shading */ closest->surfaceColor * intensity * std::max(0.0f, hitNormal.dot(l)) +
/* Blinn-Phong shading */ closest->surfaceColor * intensity * std::pow(std::max(0.0f, hitNormal.dot(h)), 3);
}

void render(const std::vector<Sphere> &spheres) {
    
    // initialize
    unsigned width = 1920, height = 1080;
    vec3f *image = new vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.0f);
    vec3f light = vec3f(0,0,20);
    vec3f lightIntensity = vec3f(1,1,1);
    
    // ray tracing and shading
    for (unsigned y = 0; y < height; ++y) {
        for (unsigned x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Ray ray = Ray(vec3f(0), vec3f(xx, yy, -1).normalize());
            *pixel = tracer(ray, spheres, light, lightIntensity);
        }
    }
    
    // write to ppm
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
    
    using namespace std;
    
    vector<Sphere> balls;
    balls.push_back(Sphere(vec3f( 5.0, -1, -15), 2, vec3f(0.90, 0.76, 0.46)));
    balls.push_back(Sphere(vec3f( 0.0,  0, -20), 2, vec3f(1.00, 0.32, 0.36)));
    balls.push_back(Sphere(vec3f( 5.0,  0, -25), 3, vec3f(0.65, 0.77, 0.97)));
    balls.push_back(Sphere(vec3f(-5.5,  0, -15), 2, vec3f(0.60, 0.30, 0.80)));
    
    Ray primRay = Ray(vec3f(1,1,1), vec3f(1,0,0));
    vec3f light = vec3f(0,0,20);
    vec3f lightIntensity = vec3f(1,1,1);
    
    render(balls);
    
    return 0;
}
