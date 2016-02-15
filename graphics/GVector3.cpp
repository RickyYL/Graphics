//
//  GVector3.cpp
//  graphics
//
//  Created by Yuanqi on 2/13/16.
//  Copyright © 2016 Yuanqi. All rights reserved.
//

#include "GVector3.hpp"

Vector3 Vector3::operator+(Vector3 v) {
    return Vector3(x+v.x, y+v.y, z+v.z);
}

Vector3 Vector3::operator-(Vector3 v) {
    return Vector3(x-v.x, y-v.y, z-v.z);
}

Vector3 Vector3::operator*(float n) {
    return Vector3(x*n, y*n, z*n);
}

Vector3 Vector3::operator/(float n) {
    return Vector3(x/n, y/n, z/n);
}

Vector3 abs(Vector3 v) {
    return Vector3(std::fabs(v.x), std::fabs(v.y), std::fabs(v.z));
}

float Vector3::dotProd(Vector3 v) {
    return x*v.x + y*v.y + z*v.z;
}

Vector3 Vector3::crossProd(Vector3 v) {
    return Vector3(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);
}

float Vector3::min() {
    return MIN(x, MIN(y, z));
}

float Vector3::max() {
    return MAX(x, MAX(y, z));
}

float Vector3::norm() {
    return sqrtf(x*x+y*y+z*z);
}

Vector3 Vector3::normalize() {
    return *this / this->norm();
}

float Vector3::distance(Vector3 v) {
    return sqrt((x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z));
}

std::ostream &Vector3::operator<<(std::ostream &os) {
    return os << "x = " << x << " y = " << y << " z = " << z;
}
