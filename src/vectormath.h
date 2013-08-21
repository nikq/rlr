/*
 * Vectormath.h
 *
 * implements vectormath
 */


#ifndef __VECTORMATH_H
#define __VECTORMATH_H

#define __NS_VECTORMATH       VECTORMATH
#define __NS_VECTORMATH_BEGIN namespace __NS_VECTORMATH {
#define __NS_VECTORMATH_END   }

#include <math.h>
#include <float.h>
#include <stdio.h>

__NS_VECTORMATH_BEGIN

#define MIN(A,B) ((A<B)?A:B)
#define MAX(A,B) ((A>B)?A:B)

template<typename T>
class Vector3{
public:

  typedef Vector3<T>  self;
  
protected:
private:
public:
  
  T x, y, z;
  
  virtual ~Vector3()  { ; }
  
  Vector3() { ; }
  Vector3( T x_, T y_, T z_ ) : x(x_),y(y_),z(z_){ }
  Vector3( const self& self ) { x = self.x; y = self.y; z = self.z; }

  void dump( void );
  
  inline void set( T x_, T y_, T z_ ) { x = x_; y = y_; z = z_; }
  inline void set( int dim, T value ) {
    switch(dim){
    case 0: x = value; break;
    case 1: y = value; break;
    case 2: z = value; break;
    }
  }
  inline T get ( int dim ) {
    switch(dim){
    case 0:  return x;
    case 1:  return y;
    case 2:  return z;
    default: return T(0);
    }
  }

  inline Vector3 operator+(const Vector3 &b) const {return Vector3(x+b.x, y+b.y, z+b.z);}
  inline Vector3 operator-(const Vector3 &b) const {return Vector3(x-b.x, y-b.y, z-b.z);}
  inline Vector3 operator+(T b)              const {return Vector3(x + b, y + b, z + b);}
  inline Vector3 operator-(T b)              const {return Vector3(x - b, y - b, z - b);}
  inline Vector3 operator*(T b)              const {return Vector3(x * b, y * b, z * b);}
  inline Vector3 mul(const Vector3 &b)       const {return Vector3(x * b.x, y * b.y , z * b.z);}
  inline T       dot(const Vector3 &b)       const {return x * b.x + y * b.y + z * b.z;}
  inline Vector3 operator%(Vector3&b)        const {return Vector3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
  
  inline Vector3 normal(void)                const {return  *this * (1.0 / sqrt(x*x+y*y+z*z)); }
  inline void    normalize(void)                   {*this = *this * (1.0 / sqrt(x*x+y*y+z*z)); }

  inline T       distance2( const Vector3 &b ) const {return (x-b.x)*(x-b.x)+(y-b.y)*(y-b.y)+(z-b.z)*(z-b.z);}
  inline T       distance ( const Vector3 &b ) const {return sqrt( distance2(b)); }
  
  inline T       length2( void ) const {return x*x+y*y+z*z;}
  inline T       length ( void ) const {return sqrt( length2() );}
  
  inline Vector3 min( const Vector3 &b ){ return Vector3( MIN(x, b.x), MIN(y, b.y), MIN(z, b.z)); }
  inline Vector3 max( const Vector3 &b ){ return Vector3( MAX(x, b.x), MAX(y, b.y), MAX(z, b.z)); }
  
  inline Vector3 cross( const Vector3 &b ){
    return Vector3(
      y * b.z - z * b.y,
      z * b.x - x * b.z,
      x * b.y - y * b.x );
  }
};

template class Vector3<double>;
typedef Vector3<double> Vector;


#undef MIN
#undef MAX

template<typename T> void Vector3<T>::dump(void){
  printf("%f,%f,%f\n",x,y,z);
}

__NS_VECTORMATH_END

#endif //__VECTORMATH_H
