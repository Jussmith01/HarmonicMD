#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
using namespace std;

//#include <GL/glew.h> // Include glew to get all the required OpenGL headers

//#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtc/type_ptr.hpp>
//#include <glm/gtx/rotate_vector.hpp>


namespace jsm
{
//-----------------------------------------------------------------------------------------//
//**************************************//
//         Quaternion Class             //
//**************************************//
template<class T>
struct quat
{
    //Class Variables
    T w;
    T x;
    T y;
    T z;

    //Default Constructor
    quat()
    {
        this->w = (T)1;
        this->x = (T)0;
        this->y = (T)0;
        this->z = (T)0;
    }

    //Defined Constructor
    quat(T w, T x, T y, T z)
    {
        this->w = (T)w;
        this->x = (T)x;
        this->y = (T)y;
        this->z = (T)z;
    }

    //Class Destructor
    ~quat() {};

    //Class Assignment
    quat& operator=(quat& instance)
    {
        this-> w = instance.w;
        this-> x = instance.x;
        this-> y = instance.y;
        this-> z = instance.z;
        return *this;
    };

    //Class Intrinsic Sum
    quat& operator+=(const quat& instance)
    {
        this->w += instance.w;
        this->x += instance.x;
        this->y += instance.y;
        this->z += instance.z;
        return *this;
    };

    //Class Extrinsic Sum
    quat& operator+(const quat& instance)
    {
        return quat(*this) += instance;
    };

    //Class Intrinsic Product
    quat& operator*=(const quat& instance)
    {
        T w1 = this->w;
        T w2 = instance.w;
        T x1 = this->x;
        T x2 = instance.x;
        T y1 = this->y;
        T y2 = instance.y;
        T z1 = this->z;
        T z2 = instance.z;

        this->w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2;
        this->x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2;
        this->y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2;
        this->z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2;
        return *this;
    };

    //Class Extrinsic Product
    quat& operator*(const quat& instance)
    {
        return quat(*this) *= instance;
    };

    //Class Formatted ostream Printer
    friend std::ostream& operator<<(std::ostream& os, const quat& q)
    {
        stringstream qout;
        qout << "[" << q.w << "," << q.x << "i + " << q.y << "j + " << q.z << "k] ";
        return os << qout.str();
    };
};

//Function calculates the length of the quat class elements
template<typename T>
T quatLength(const quat<T> &q)
{
    T magnitude;
    magnitude = sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    return magnitude;
};

//Function normalizes the given quaternion
template<typename T>
const quat<T> quatNormalize(const quat<T> &q)
{
    quat<T> qr;
    T magnitude = quatLength(q);
    qr.w = q.w / (T)magnitude;
    qr.x = q.x / (T)magnitude;
    qr.y = q.y / (T)magnitude;
    qr.z = q.z / (T)magnitude;
    return qr;
};

//Function Converts the given quaternion to a 4x4 rotation transformation matrix
/*template<typename T>
glm::mat4 quatUnittoMatrix(const quat<T> &q)
{
        glm::mat4 mr;
        mr[0][0] = 1 - 2 * ( q.y * q.y + q.z * q.z );
        mr[0][1] = 2 * ( q.x * q.y - q.z * q.w );
        mr[0][2] = 2 * ( q.x * q.z + q.y * q.w );
        mr[0][3] = 0;
        mr[1][0] = 2 * ( q.x * q.y + q.z * q.w );
        mr[1][1] = 1 - 2 * ( q.x * q.x + q.z * q.z );
        mr[1][2] = 2 * ( q.y * q.z - q.x * q.w );
        mr[1][3] = 0;
        mr[2][0] = 2 * ( q.x * q.z - q.y * q.w );
        mr[2][1] = 2 * ( q.y * q.z + q.x * q.w );
        mr[2][2] = 1 - 2 * ( q.x * q.x + q.y * q.y );
        mr[2][3] = 0; mr[3][0] = 0; mr[3][1] = 0;
        mr[3][2] = 0; mr[3][3] = 1;
        return mr;
};*/

//Function generates the local quaternion of a given glm::vec3 and float (radians)
/*template<typename T>
const quat<T> quatGenLocal( T a,glm::vec3 v )
{
	quat<T> ql;
	ql.w = cos(a/(T)2);
	ql.x = v.x * sin(a/(T)2);
	ql.y = v.y * sin(a/(T)2);
	ql.z = v.z * sin(a/(T)2);
	return ql;
};*/

//Function generates the local quaternion of a given glm::vec3 and float (radians)
/*template<typename T>
glm::vec3 quatRotVec3( quat<T> q,glm::vec3 v3 )
{
        glm::vec4 v4 = glm::vec4(v3,1.0f);
glm::mat4 m4 = quatUnittoMatrix(q);
glm::vec3 rv = glm::vec3(m4 * v4);
        return rv;
};*/
//-----------------------------------------------------------------------------------------//
//**************************************//
//            Vector3 Class             //
//**************************************//
template<class T>
class vec3
{
    //Class Variables
    T xyz[3];

public:
    //Default Constructor
    vec3()
    {
        this->xyz[0] = (T)0;
        this->xyz[1] = (T)0;
        this->xyz[2] = (T)0;
    };

    //Default Constructor
    vec3(T t)
    {
        this->xyz[0] = (T)t;
        this->xyz[1] = (T)t;
        this->xyz[2] = (T)t;
    };

    //Defined Constructor
    vec3(T x, T y, T z)
    {
        this->xyz[0] = (T)x;
        this->xyz[1] = (T)y;
        this->xyz[2] = (T)z;
    };

    //Class Destructor
    ~vec3() {};

    //Class subscript operator
    T &operator[](int i)
    {
        /*if (i<0 || i>2)
        {
            cout << "error: Index out of bounds.\n";
        }*/

        return xyz[i];
    }

    //Class Assignment
    vec3 operator=(const vec3& instance)
    {
        this->xyz[0] = instance.xyz[0];
        this->xyz[1] = instance.xyz[1];
        this->xyz[2] = instance.xyz[2];
        return *this;
    };

    //Class Intrinsic Sum
    vec3 operator+=(const vec3& instance)
    {
        this->xyz[0] += instance.xyz[0];
        this->xyz[1] += instance.xyz[1];
        this->xyz[2] += instance.xyz[2];
        return *this;
    };

    //Class Extrinsic Sum
    vec3 operator+(const vec3& instance)
    {
        return vec3(*this) += instance;
    };

    //Class Intrinsic Sum
    vec3 operator-=( vec3& instance )
    {
        this->xyz[0] -= instance[0];
        this->xyz[1] -= instance[1];
        this->xyz[2] -= instance[2];
        return *this;
    };

    //Class Extrinsic Sum
    vec3 operator-( vec3& instance )
    {
        return vec3(*this) -= instance;
    };

    //Class Intrinsic Product
    vec3 operator*=(T scalar)
    {
        this->xyz[0] *= scalar;
        this->xyz[1] *= scalar;
        this->xyz[2] *= scalar;

        return *this;
    };

    //Class Extrinsic Product
    vec3 operator*(T scalar)
    {
        return vec3(*this) *= scalar;
    };


    //Class Formatted ostream Printer
    friend std::ostream& operator<<(std::ostream& os, const vec3& v)
    {
        stringstream vout;
        vout << " [" << v.xyz[0] << ", " << v.xyz[1] << ", " << v.xyz[2] << "] ";
        return os << vout.str();
    };

    void put(int i,T val)
    {
        xyz[i] = (T)val;
    };

    T get(int i)
    {
        return xyz[i];
    };
};

//Function calculates the dot product of two vec3 class vectors
template<typename T>
inline T dot3( const vec3<T> v1,const vec3<T> v2 )
{
    return v1.xyz[0] * v2.xyz[0] + v1.xyz[1] * v2.xyz[1] + v1.xyz[2] * v2.xyz[2];
};

//Function calculates the cross product of two vec3 class vectors
template<typename T>
inline vec3<T> cross( const vec3<T> va,const vec3<T> vb )
{
    T x = va.xyz[1] * vb.xyz[2] - va.xyz[2] * vb.xyz[1];
    T y = va.xyz[2] * vb.xyz[0] - va.xyz[0] * vb.xyz[2];
    T z = va.xyz[0] * vb.xyz[1] - va.xyz[1] * vb.xyz[0];
    vec3<T> rv(x,y,z);
    return rv;
};

//Function calculates the magnitude of a vec3 class vector
template<typename T>
inline T magnitude( vec3<T> v )
{
    T x2 = v[0] * v[0];
    T y2 = v[1] * v[1];
    T z2 = v[2] * v[2];
    T L = sqrt( x2 + y2 + z2 );
    return L;
};

//Function calculates the magnitude of a vec3 class vector
template<typename T>
inline vec3<T> UniformScalar( vec3<T> v, T s )
{
    vec3<T> vo;
    vo[0] = v[0] * s;
    vo[1] = v[1] * s;
    vo[2] = v[2] * s;
    return vo;
};

//Function calculates the magnitude of a vec3 class vector
template<typename T>
inline vec3<T> zScalar( vec3<T> v, T s )
{
    vec3<T> vo;
    vo[0] = v[0];
    vo[1] = v[1];
    vo[2] = v[2] * s;
    return vo;
};

//Function normalizes the given quaternion
template<typename T>
inline vec3<T> normalize( vec3<T> &v)
{
    vec3<T> vr;
    T mag = (T)magnitude(v);

    T invmag = 1/(T)mag;

    vr[0] = v[0] * invmag;
    vr[1] = v[1] * invmag;
    vr[2] = v[2] * invmag;

    if (vr[0] != vr[0] || vr[1] != vr[1] || vr[2] != vr[2])
    {
        vr[0] = 0;
        vr[1] = 0;
        vr[2] = 0;
    }

    return vr;
};
};
