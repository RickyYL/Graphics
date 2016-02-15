//
//  GVector3.hpp
//  graphics
//
//  Created by Yuanqi on 2/13/16.
//  Copyright © 2016 Yuanqi. All rights reserved.
//

#ifndef GVector3_hpp
#define GVector3_hpp

#include <iostream>
#include <cmath>
#define MIN(x,y) (x)>(y)?(y):(x)
#define MAX(x,y) (x)>(y)?(x):(y)

class Vector3 {
public:
    
    float x, y, z;
    
    Vector3();
    Vector3(float posX, float posY, float posZ)
        : x(posX), y(posY), z(posZ) {}
    
    Vector3 operator+(Vector3 v);
    Vector3 operator-(Vector3 v);
    Vector3 operator*(float n);
    Vector3 operator/(float n);
    float   dotProd(Vector3 v);
    Vector3 crossProd(Vector3 v);
    float   max();
    float   min();
    float   norm();
    Vector3 normalize();
    float   distance(Vector3 v);
    Vector3 zero() {return Vector3(0,0,0);}
    std::ostream &operator<<(std::ostream &os);
};

Vector3 abs(Vector3 v);





#endif /* GVector3_hpp */
