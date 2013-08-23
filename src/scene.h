/* -*- coding: shift_jis -*-
 *
 * scene.h
 *
 * Copyright (c) 2013 Hajime UCHIMURA / nikq.
 */


#ifndef __SCENE_H
#define __SCENE_H


//#include <boost/shared_ptr.h>
#include <stdio.h>
#include <stdlib.h>
#include "vectormath.h" // for vectormath
#include "mtseq.h"
#include <vector>

#define __NS_RLR       RLR
#define __NS_RLR_BEGIN namespace __NS_RLR {
#define __NS_RLR_END   }

__NS_RLR_BEGIN


typedef __NS_VECTORMATH::Vector Vector;
typedef __NS_VECTORMATH::Vector Color;      // R,G,B

class Ray{
protected:
private:
public:
  
  Vector pos_;
  Vector dir_;
  
  Ray(){ ; }
  virtual ~Ray(){ ; }
  
  const Vector& origin()    const { return pos_; }
  const Vector& direction() const { return dir_; }
};

class Spectrum{
public:
  
  double waveLength_;
  double power_;

  static double min( void ) { return 380.; }
  static double max( void ) { return 780.; }
  static double rnd(  __NS_MTSEQ::MTSequence &mt ) { return 380 + 400 * mt.genrand_real1() ; }

  Spectrum( ) { ; }
  Spectrum( double w, double p ) : waveLength_(w) , power_( p ) { ; }
  Spectrum absorb( double r ) { return Spectrum( waveLength_, power_ * r ) ; }
  inline Spectrum operator*( double r ) const { return Spectrum( waveLength_, power_ * r ) ; }
  double waveLength(){ return waveLength_ ; }
  
  Color getRGB( void ){
    Color col;
    double w = waveLength_;
    
    if ( w >= 380 && w < 440 ){
      col.x =  -(w - 440.) / (440. - 380.);
      col.y = 0.0;
      col.z = 1.0;
    }else if( w >= 440 && w < 490 ){
      col.x =  0.0;
      col.y = (w - 440.) / (490. - 440.);
      col.z = 1.0;
    }else if( w >= 490 && w < 510){
      col.x =  0.0;
      col.y = 1.0;
      col.z = -(w - 510.) / (510. - 490.);
    } else if( w >= 510 && w < 580 ){
      col.x =  (w - 510.) / (580. - 510.);
      col.y = 1.0;
      col.z = 0.0;
    } else if( w >= 580 && w < 645){
      col.x =  1.0;
      col.y = -(w - 645.) / (645. - 580.);
      col.z = 0.0;
    } else if( w >= 645 && w <= 780 ){
      col.x =  1.0;
      col.y = 0.0;
      col.z = 0.0;
    } else {
      col.x =  0.0;
      col.y = 0.0;
      col.z = 0.0;
    }
    return col * power_;
  }
};


class SceneObject;
class Intersection{
protected:
private:
public:
  
  bool    hit_;
  double  distance_;
  
  Vector point_;
  Vector normal_;
  SceneObject* object_;

  Intersection() : hit_(false), distance_(-1.){;}
  
  void          clear()                 { hit_ = false; }
  
  inline void   setHit     ( bool flag ){ hit_      = flag; }
  inline void   setDistance( double t  ){ distance_ = t;    }
  inline void   setPoint   ( Vector p    ){ point_  = p;    }
  inline void   setNormal  ( Vector norm ){ normal_ = norm; }
  inline void   setObject  ( SceneObject*o){ object_ = o; }
  
  inline bool   hit      ()   const     { return hit_;      }
  inline double distance ()   const     { return distance_; }
  inline Vector getPoint () const { return point_;  }
  inline Vector getNormal() const { return normal_; }
  inline SceneObject* getObject() const{ return object_; }
  
  inline void update( const Intersection& i ){
    if( !i.hit() || i.distance() <= 0. )
      return;
    if( !hit() || distance() > i.distance() ){
      setHit( i.hit() );
      setDistance( i.distance() );
      setPoint( i.getPoint() );
      setNormal( i.getNormal() );
      setObject( i.getObject() );
    }
  }
  
};

class AABB{
protected:
private:
  
  Vector lo_;
  Vector hi_;
  
public:
  
  virtual ~AABB(){;}
  AABB(){ clear(); }

  void dump( void ){
    printf("dumping AABB (%f,%f,%f)>(%f,%f,%f)\n",
           lo_.x, lo_.y, lo_.z,
           hi_.x, hi_.y, hi_.z );
  }
  
  void clear( void ){
    hi_.set( -FLT_MAX, -FLT_MAX, -FLT_MAX );
    lo_.set(  FLT_MAX,  FLT_MAX,  FLT_MAX );
  }
  void expand( Vector v ){ // AABBがvを包含するように拡張する.
    lo_ = lo_.min( v );
    hi_ = hi_.max( v );
  }
  void expand( AABB aabb ){
    lo_ = lo_.min( aabb.lo_ );
    hi_ = hi_.max( aabb.hi_ );
  }
  
  Vector& lo() { return lo_; }
  Vector& hi() { return hi_; }
  
  double distance( Vector& v ) const {
    double dx = 0.f;
    double dy = 0.f;
    double dz = 0.f;
    double d  = 0.f;
    
    if( v.x <= lo_.x )     { dx = lo_.x - v.x; d += dx*dx; }
    else if( hi_.x <= v.x ){ dx = v.x - hi_.x; d += dx*dx; }
    if( v.y <= lo_.y )     { dy = lo_.y - v.y; d += dy*dy; }
    else if( hi_.y <= v.y ){ dy = v.y - hi_.y; d += dy*dy; }
    if( v.z <= lo_.x )     { dz = lo_.z - v.z; d += dz*dz; }
    else if( hi_.z <= v.z ){ dz = v.z - hi_.z; d += dz*dz; }
    
    return sqrt( d );
  }
  
  bool test   ( const Vector& v ){
    return ( lo_.x <= v.x && v.x <= hi_.x &&
             lo_.y <= v.y && v.y <= hi_.y &&
             lo_.z <= v.z && v.z <= hi_.z );
  }

  bool clipRay( const Ray& ray, double& l, double& h ){
    double tmin = 0.f;
    double tmax = FLT_MAX;
    
    for(int axis=0;axis<3;axis++){
      const double d = const_cast<RLR::Vector&>(ray.direction()).get(axis);
      const double p = const_cast<RLR::Vector&>(ray.origin()).get(axis);
      double l = lo_.get(axis);
      double h = hi_.get(axis);

      if( -1e-10 < d && d < 1e-10 ){
        if( p < l || h < p )
          return false;
      }else{
        double ood = 1.0 / d;
        double t1  = (l - p) * ood;
        double t2  = (h - p) * ood;
        double tl,th;
        if( t1 > t2 ){
          tl = t2;
          th = t1;
        }else{
          tl = t1;
          th = t2;
        }
        if( tmin < tl ) tmin = tl;
        if( tmax > th ) tmax = th;
        if( tmin > tmax )
          return false;
      }
    }
    l = tmin;
    h = tmax;
    return true;
  }

  inline bool collide_test( double u0, double u1, double v0, double v1 ){ return !( u1 < v0 || v1 < u0 ); }
  bool collide( const AABB& bb ){
    return ( collide_test( lo_.x, hi_.x, bb.lo_.x, bb.hi_.x ) &&
             collide_test( lo_.y, hi_.y, bb.lo_.y, bb.hi_.y ) &&
             collide_test( lo_.z, hi_.z, bb.lo_.z, bb.hi_.z ) );
  }
};


class SceneMaterial{
protected:
private:
public:
  
  double diffuse_;
  double reflect_;
  double refract_;
  double emitter_;
  double ior_;
  
  Color  color_; // 素材色.

  SceneMaterial():diffuse_(0.),reflect_(0.),refract_(0.),emitter_(0.),ior_(1.){;}
  virtual ~SceneMaterial(){;}
  
  void normalize( void ){
    double sum = diffuse_ + reflect_ + refract_;
    if( sum > 0. ){
      diffuse_ /= sum;
      reflect_ /= sum;
      refract_ /= sum;
    }
  }
  double ior( double w ){ return ior_ - ior_ * 0.24 * (w - 589.3) / 400. ; }
};


class SceneObject{
protected:
private:
public:

  // geometric informations.
  // physical material informations.
  
  typedef SceneObject self;
  SceneMaterial      *material_ptr_;
  
  SceneObject()         : material_ptr_(NULL) { ; }
  virtual ~SceneObject(){ material_ptr_ = NULL;  }
  
  virtual AABB getAABB  ( void ) = 0;               // get bounding of this object.
  virtual bool testAABB ( const AABB& aabb ) = 0;   // test if object is crossing AABB.
  virtual bool testRAY  ( const Ray& ray, Intersection &isect ) = 0;
  virtual Vector getCenter( void ) { return Vector(0,0,0); }
  virtual Vector getRandomPoint( __NS_MTSEQ::MTSequence &mt ) { return Vector(0,0,0); }
  inline SceneMaterial *getMaterial() const { return material_ptr_; }
  void setMaterial( SceneMaterial *mat ) { material_ptr_ = mat; }
};


#define TRIANGLE_EPSILON 1e-10
class Triangle : public SceneObject {
protected:
private:

  Vector vertex_;
  Vector edge0_;
  Vector edge1_;
  Vector normal_;
  
public:
  
  Vector getCenter( void ){
    return vertex_;
  }
  Vector getRandomPoint( __NS_MTSEQ::MTSequence &mt  ) {
    double u,v;
    do{
      u = mt.genrand_real1();
      v = mt.genrand_real1();
    }while( u+v >= 1. );
    return vertex_ + edge0_ * u + edge1_ * v;
  }
  
  void setup( Vector a, Vector b, Vector c ){
    vertex_ = a;
    edge0_  = b - a;
    edge1_  = c - a;
    normal_ = edge1_.cross( edge0_ );
    normal_.normalize();
  }
  
  AABB getAABB  ( void ){
    AABB bb;
    bb.expand( vertex_ );
    bb.expand( vertex_ + edge0_ );
    bb.expand( vertex_ + edge1_ );
    return bb;
  }
  bool testAABB ( const AABB& aabb ){
    return true;
  }
  bool testRAY  ( const Ray& ray, Intersection &isect ) {
    Vector pv, tv, qv;
    double d, u, v, id;

    Vector pos( ray.origin()    );
    Vector dir( ray.direction() );

    //int_count++;
    
    pv = dir.cross( edge1_ );
    d  = edge0_.dot( pv );

    if( d < -TRIANGLE_EPSILON || d > TRIANGLE_EPSILON ) {
      id = 1.0 / d;
      tv = pos - vertex_;
      u  = tv.dot( pv ) * id;
      if( u < 0. || u > 1.0f )
        return false;
      
      qv = tv.cross( edge0_ );
      v  = dir.dot( qv ) * id;
      if( v < 0. || (u+v) > 1.0f )
        return false;
    }else{
      return false;
    }
    double t = edge1_.dot(qv) * id ;
    isect.setHit     ( true );
    isect.setDistance( t );
    
    Vector hit  = pos + dir * t;
    //Vector norm = hit - center_;
    
    isect.setPoint   ( hit );
    isect.setNormal  ( normal_ );
    isect.setObject  ( this );
    
    return true;
  }
};

class Sphere : public SceneObject {
protected:
private:
  
  Vector center_;
  double radius_;
  
public:
  
  void setup( Vector c, double r ){
    center_ = c; radius_ = r;
  }
  
  AABB getAABB  ( void ){
    AABB bb;
    bb.expand( center_ + Vector(radius_, radius_, radius_ )); // 上.
    bb.expand( center_ - Vector(radius_, radius_, radius_ )); // 下.
    return bb;
  }
  
  bool testAABB ( const AABB& aabb ){
    double d = aabb.distance( center_ ); // AABB<>中心の距離
    return d <= radius_;
  }
  
  bool testRAY  ( const Ray& ray, Intersection &isect ){
    const Vector &origin = ray.origin();
    const Vector &dir    = ray.direction();
    
    Vector c = center_ - origin;
    double b    = c.dot( dir );
    double det2 = b * b - c.dot( c ) + radius_*radius_;
    
    if( det2 < 0. )
      return false;
    
    double det = sqrt( det2 );
    double t  = b - det;
    if( t < 0. ){
      t = b + det;
      if( t < 0. )
        return false;
    }
    
    isect.setHit     ( true );
    isect.setDistance( t    );
    
    Vector hit = origin + dir * t;
    Vector norm= hit - center_;
    
    isect.setPoint   ( hit );
    isect.setNormal  ( norm.normal() );
    isect.setObject  ( this );
    return true;
  }
  Vector getCenter( void ){
    return center_;
  }
};


class Scene : public SceneObject {
public:
  
  //typedef boost::shared_ptr< SceneObject > SceneObjectPtr;
  typedef SceneObject *                    SceneObjectPtr;
  typedef std::vector< SceneObjectPtr >    SceneObjectList;
  
protected:
private:
  
  SceneObjectList objects_;
  SceneObjectList emitters_;

public:
  
  void registObject( SceneObjectPtr object ){
    if( object->getMaterial()->emitter_ > 0. )
      emitters_.push_back( object );
    else
      objects_.push_back( object );
  }
  
  AABB getAABB  ( void ) {
    AABB aabb;
    for( SceneObjectList::iterator obj = objects_.begin();
         obj != objects_.end(); ++ obj ){
      aabb.expand( (*obj)->getAABB() );
    }
    return aabb;
  }
  
  bool testAABB ( const AABB& aabb ) {
    return true;
  }
  
  bool testRAY  ( const Ray& ray, Intersection &isect ) {
    Intersection result;
    isect.clear();
    for( SceneObjectList::iterator obj = objects_.begin();
         obj != objects_.end(); ++ obj ){
      (*obj)->testRAY( ray, result );
      isect.update( result );
    }
    return isect.hit();
  }

  SceneObjectList &getObjectList(){  return objects_;  }
  SceneObjectList &getEmitterList(){ return emitters_; }
  Color            getBackground(){ return emitters_.size() == 0 ? Color(1.,1.,1.) : Color(0,0,0); }
};


class Camera{
private:
public:
  Camera(){ zoom_ = 1.0; }
  virtual ~Camera(){;}


  Vector pos_, lookat_, dir_, dirU_, dirV_;
  double zoom_;     // フィルム距離
  double dist_;
  double radius_;
  double focus_offset_;
  void setup(Vector p, Vector l, double z = 1., double r = 0., double f = 1. ){
    Vector v;

    pos_    = p;
    lookat_ = l;
    dir_    = lookat_ - pos_;
    dist_   = dir_.length();
    
    printf("[CAM] position : "); pos_.dump();
    printf("[CAM] lookat   : "); lookat_.dump();
    //printf("[CAM] direction: "); dir_.dump();
    //dist = dir.dot( dir );

    dir_.normalize();
    printf("[CAM] direction: "); dir_.dump();

    v.set( 0.0, 1.0, 0.0 );
    dirU_ = dir_.cross( v );
    dirU_.normalize();
    dirV_ = dirU_.cross( dir_ );
    dirV_.normalize();

    zoom_   = z;
    radius_ = r;
    focus_offset_ = f;
    
    printf("[CAM]  world x : "); dirU_.dump();
    printf("[CAM]  world y : "); dirV_.dump();
    printf("[CAM]  world z : "); dir_.dump();
  }
  void ray( __NS_MTSEQ::MTSequence& mt,double x, double y, Ray& r ){
    double u,v;
    do{
      u = mt.genrand_real1() * 2. - 1.;
      v = mt.genrand_real1() * 2. - 1.;
    }while( (u*u+v*v) > 1. );
    
    r.pos_ = pos_ + (dirU_ * u * radius_) + (dirV_ * v * radius_);
    Vector l = pos_ + dir_ * (dist_ * focus_offset_) + dirU_ * (x * dist_ * focus_offset_ / zoom_) + dirV_ * (y * dist_  * focus_offset_ / zoom_);
    r.dir_ = (l - r.pos_).normal();
  }
};


template< typename SCENE >
class FrameBuffer {
private:
public:
  int width_;
  int height_;
  int sample_;
  int pass_;
  __NS_RLR::Color *film_;
  __NS_RLR::Color *raw_;
  double *weight_;
  double *dx_;
  double *dy_;

  virtual ~FrameBuffer(){ clear(); }
  FrameBuffer() : width_(0), height_(0), sample_(0), pass_(0), film_(NULL), raw_(NULL){;}

  void setup( int w, int h, int s = 1 ) {
    int i;
    float a;
    width_ = w;
    height_= h;
    sample_= s;
    film_ = (__NS_RLR::Color*) malloc( sizeof( __NS_RLR::Color ) * w * h );
    raw_  = (__NS_RLR::Color*) malloc( sizeof( __NS_RLR::Color ) * w * h );
    weight_=(double*)malloc(sizeof(double)*w*h);
    dx_   = new double[ s ];
    dy_   = new double[ s ];

    for( i = 0; i<w*h; i++ ){
      film_[i].set( 0., 0., 0. );
      weight_[i]=0.;
    }

    a = 0.;
    for(i = 0;i<sample_;i++){
      dx_[i] = (double) i / sample_;
      dy_[i] = a - floor( a );
      a      = a + (1. + sqrt(5.))/2.;
    }
  }

  void clear(void)
    {
      width_ = height_ = sample_ = pass_ = 0;
      if( film_ ){ free( film_ ); film_ = NULL; }
      if( raw_  ){ free( raw_ );  raw_ = NULL;  }
      if( dx_   ){ delete [] dx_; dx_ = NULL; }
      if( dy_   ){ delete [] dy_; dy_ = NULL; }
    }
  
  void saturate( void ) {
    int i,j;
    double lo[3] = {FLT_MAX,FLT_MAX,FLT_MAX};
    double hi[3] = {FLT_MIN,FLT_MIN,FLT_MIN};
    for(i=0;i<width_*height_;i++){
      for(j=0;j<3;j++){
        double v = film_[i].get(j);
        if( v > hi[j] ) hi[j] = v;
        if( v < lo[j] ) lo[j] = v;
      }
    }
    printf("saturate %f - %f , %f - %f , %f - %f",hi[0],lo[0],hi[1],lo[1],hi[2],lo[2]);
    for(i=0;i<width_*height_;i++){
      for(j=0;j<3;j++){
        raw_[i].set(j,(film_[i].get(j)-lo[j]) / (hi[j]-lo[j]));
      }
    }
  }
  void normalize( void ) {
    int i,j;
    for(i=0;i<width_*height_;i++){
      double w = 1.;
      if( weight_[i] > 0. )
        w = weight_[i];
      for(j=0;j<3;j++){
        raw_[i].set(j,film_[i].get(j) / w);
      }
    }
  }
  
  int toInt(double x,double g){
    // tone mapping
    return int(pow(1-exp(-x),1/g)*255+.5);
    //return int( x*255. );
  }
  void save_bmp( char *filename, double scale, float gamma ) {
    unsigned char *rgb;
    int pixel;
    int x,y,i,j;
    gamma = 1.f / gamma;
    rgb = (unsigned char*)malloc(width_ * height_ * 3);
    for(y = 0; y < height_; y ++ ){
      for(x = 0;x < width_; x ++ ){
        i = (y*width_)+x;
        j = (((height_-y-1)*width_)+x)*3;
        pixel      = (int) toInt( raw_[i].get(0) * scale, gamma );
        rgb[j ] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);

        pixel      = (int) toInt( raw_[i].get(1) * scale, gamma );
        rgb[j+1] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);

        pixel      = (int) toInt( raw_[i].get(2) * scale, gamma );
        rgb[j+2] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);
      }
    }
    bmp_write( filename, width_, height_, rgb );
    free(rgb);
  }
  void save_hdr( char *filename, double scale, float gamma ) {
    float *rgb;
    int x,y,i,j;
    rgb = (float*)malloc(width_ * height_ * 3 * sizeof(float));
    for(y = 0; y < height_; y ++ ){
      for(x = 0;x < width_; x ++ ){
        i = (y*width_)+x;
        j = i*3;
        rgb[j+0] = raw_[i].get(0);
        rgb[j+1] = raw_[i].get(1);
        rgb[j+2] = raw_[i].get(2);
      }
    }
    hdr_write( filename, width_, height_, rgb );
    free(rgb);
  }


  typedef struct {
    unsigned char  id[2];
    unsigned long  filesize;
    unsigned short reserve1;
    unsigned short reserve2;
    unsigned long  offset; // 26:12,54:40,122:108
    unsigned long  headsize; // 12,40,108
    unsigned long  width;
    unsigned long  height;
    unsigned short plane;
    unsigned short bpp;
    unsigned long  method;
    unsigned long  datasize;
    unsigned long  x_dpm;
    unsigned long  y_dpm;
    unsigned long  palnum;
    unsigned long  imp_pal;
  }BMP_header;

  inline float maxf( float a, float b ){ return (a>b) ? a:b ; }

  // 指数部を返す
  unsigned char exponent( float a ){
    unsigned char *c = (unsigned char*)&a;
    unsigned char  e = ((c[3]<<1)&0xFE) | ((c[2]>>7)&1);
    return e;
  }
  // 同じ指数の最も小さい少数を返す
  float basefloat( float a ){
    float f = a;
    unsigned char *c = (unsigned char*)&f;
    c[2] &= ~1;
    c[1]  =  0;
    c[0]  =  0;
    return f;
  }
  
  int hdr_write(char *name,int width,int height,float *rgb){
    FILE *fp;
    fp = fopen( name, "wb" );
    if(!fp)
      return -1;
    fprintf( fp, "#?RADIANCE\n");
    fprintf( fp, "# Made with RLR\n");
    fprintf( fp, "FORMAT=32-bit_rle_rgbe\n");
    fprintf( fp, "EXPOSURE=          1.0000000000000\n");
    fprintf( fp, "\n");
    fprintf( fp, "-Y %d +X %d\n",height,width);

    //unsigned char magic4[4] = {0x02, 0x02, (width>>8)&0xFF, width&0xFF};
    unsigned char line[ 0x7FFF*4 ];
    
    //for(int y=0;y<height;y++){
    for(int y=height-1;y>=0;y--){
      
      //fwrite( magic4, 1, 4, fp );
      // RGBE化
      for(int x=0;x<width;x++){
        float r = rgb[ (y*width + x)*3 + 0 ];
        float g = rgb[ (y*width + x)*3 + 1 ];
        float b = rgb[ (y*width + x)*3 + 2 ];
        float m = maxf(maxf(r,g),b);
        
        int e;
        int iR, iG, iB, iE;
        
        if(m <= 1e-32) {
          iR = iG = iB = iE = 0;
        } else {
          m = frexp(m, &e) * 255.9999 / m;
          
          iR = (int)(r * m);
          iG = (int)(g * m);
          iB = (int)(b * m);
          iE = (int)(e + 128);
        }
        line[ x * 4 + 0 ] = iR;
        line[ x * 4 + 1 ] = iG;
        line[ x * 4 + 2 ] = iB;
        line[ x * 4 + 3 ] = iE;
      }
      
      fwrite( line, 1, width*4, fp );
      
      // RLE化
      /*
      for(int x=0; x < width ;x++){
        printf("%02x,",line[x]);
      }
      printf("\n");
      
      for(int x=0; x < width ;x++){
        unsigned char run;
        //printf("x:%d\n",x);
        // 連続側
        run = 0;
        while( run < 127 && (x+run) < width && line[x+run] == line[x] ){
          //printf("run %d, %02x,%02x\n",run,line[x],line[x+run]);
          run++;
        }
        //printf("run:%d\n",run);
        if( run > 0 ){
          x   += run;
          run += 0x80;
          fwrite( &run      , 1, 1, fp );
          fwrite( &(line[x]), 1, 1, fp );
        }
        
        
        // 不連続側
        //printf("x:%d\n",x);
        run = 0;
        while( run < 127 && (x+run) < width && line[x+run] != line[x+run+1] ){
          //printf("run %d, %02x,%02x\n",run,line[x+run],line[x+run+1]);
          run++;
        }
        //printf("run:%d\n",run);
        if( run > 0 ){
          fwrite( &run      ,   1, 1, fp );
          fwrite( &(line[x]), run, 1, fp );
          x += run;
        }
      }*/
    }
    
    fclose(fp);
    return 0;
  }
  
  int bmp_write(char *name,int width,int height,unsigned char *rgb) {
    BMP_header bh;
    FILE *fp;
    int x,y,k;
    int r,g,b;
    int mod;

    bh.id[0] = 'B';
    bh.id[1] = 'M';
    bh.filesize = 54 + (width*height*3);
    bh.offset   = 54;
    bh.headsize = 40;
    bh.width    = width;
    bh.height   = height;
    bh.plane    = 1;
    bh.bpp      = 24;
    bh.method   = 0;
    bh.datasize = 0;
    bh.x_dpm    = 3779;
    bh.y_dpm    = 3779;
    bh.palnum   = 0;
    bh.imp_pal  = 0;

    fp = fopen(name,"wb");
    if(!fp) return -1;
    fwrite( bh.id, 1,2 ,fp);
    fwrite( &(bh.filesize), 4,1,fp);
    fwrite( &(bh.reserve1), 2,1,fp);
    fwrite( &(bh.reserve2), 2,1,fp);
    fwrite( &(bh.offset  ), 4,1,fp);
    fwrite( &(bh.headsize), 4,1,fp);
    fwrite( &(bh.width   ), 4,1,fp);
    fwrite( &(bh.height  ), 4,1,fp);
    fwrite( &(bh.plane   ), 2,1,fp);
    fwrite( &(bh.bpp     ), 2,1,fp);
    fwrite( &(bh.method  ), 4,1,fp);
    fwrite( &(bh.datasize), 4,1,fp);
    fwrite( &(bh.x_dpm   ), 4,1,fp);
    fwrite( &(bh.y_dpm   ), 4,1,fp);
    fwrite( &(bh.palnum  ), 4,1,fp);
    fwrite( &(bh.imp_pal ), 4,1,fp);

    mod = (width*3)%4;
    if(mod)
      mod = 4-mod;
    //printf("mod:%d\n",mod);
    for(y=0;y<height;y++){
      for(x=0;x<width;x++){
        k = ((height-y-1)*width + x)*3;

        r = rgb[k  ]; // R
        g = rgb[k+1]; // G
        b = rgb[k+2]; // B
        fputc(b,fp);
        fputc(g,fp);
        fputc(r,fp);
      }
      if(mod>0){
        for(x=0;x<mod;x++){
          fputc(0,fp);
        }
      }
    }
    fclose(fp);
    return 0;
  }

};


__NS_RLR_END


#endif
