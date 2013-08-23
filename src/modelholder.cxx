/* -*- coding: shift_jis -*-
 *
 * MQO file loader.
 *
 * Copyright (c) 2013 Hajime UCHIMURA / nikq.
 */

#include <string.h>
#include <stdio.h>
#include "vectormath.h"
#include "modelholder.h"

__NS_RLRUTIL_BEGIN

// 受け付ける文法："hagehoge"  / hoge( hoge hoge hoge )
//トークンセパレータ：\x20 || \t
char *ModelHolder::tokenize( char *s, char *token )
{
  char *p = (char *)s;
  while( *p && (*p == ' ' || *p == '\t') )
    p++;
  if( *p == 0 )
    return NULL;
  if( *p == '\"' ) {
    p++;
    while(*p != '\"' && *p){
      *token = *p; token++; p++;
    }
    p++;
    *token = 0;
  } else {
    while(*p != ' ' && *p != '\t' && *p != '\r' && *p != '\n' && *p) {
      *token = *p;
      if(*p == '(') {
        token++; p++;
        while( *p != ')' && *p ){
          *token = *p; token++; p++;
        }
        *token = *p;
      }
      token++; p++;
    }
    *token = 0;
  }
  if( *p == '\r' || *p == '\n' || *p == 0 )
    return NULL;
  return p;
}

bool ModelHolder::load( const char *fn )
{
  FILE *fp;
  char line[1024];
  char *p,*p2;
  char token[256],t2[256];

  int nmat;
  int nvert;
  int nface;
  int i,a = 0,b = 0,c = 0,d = 0,m = 0;
  int depth;
  int totalface = 0, currface = 0;

  //VertList vertex;
  //std::vector< __NS_VECTORMATH::Vector > vertices;
  __NS_RLR::Triangle poly;

  double x,y,z;

  fp = fopen(fn,"rb");
  if(fp == NULL)
    return false;
  
  fgets(line,1024,fp);
  if(strcmp(line,"Metasequoia Document\r\n") != 0)
    return false;
  
  fgets(line,1024,fp);
  if(strcmp(line,"Format Text Ver 1.0\r\n") != 0)
    return false;
  
  __NS_RLR::SceneMaterial white;
  
  materials_.clear();
  materials_.push_back( white );
  
  //ロード開始
  while(1) {
    fgets(line,1024,fp);
    //printf("\"%s\"\n",line);
    if(strncmp(line,"Eof",3) == 0 || feof(fp))
      break;

    p = line;
    p = tokenize(p,token);

    if(strcmp(token,"Scene") == 0) {
      __NS_VECTORMATH::Vector pos, look, clook, cpos;
      while(1) {
        fgets(line,1024,fp);
        p = line;
        p = tokenize(p,token);
        if( token[0] == '}' ) {
          // 計算するぞおい
          break;
        }
      }
    }

    if(strncmp(token,"Material",6) == 0) {
      materials_.clear();
      p = tokenize(p,token);
      nmat = atoi(token);
      printf("[MQO] %d materials\n",nmat);
      for(i=0;i<nmat;i++) {
        //マテリアルの読み込み
        fgets(line,1024,fp);
        p = line;
        p = tokenize(p,token);
        printf("[MAT] %2d:%12s :",i,token);
        //終了
        if(token[0] == '}')
          break;

        // 構文サンプル

        __NS_RLR::SceneMaterial m;
        // "mat1"  shader(3) vcol(1) col(1.000 1.000 1.000 1.000) dif(1.000) amb(0.250) emi(0.000) spc(0.000) power(0.00) tex("00tex_master.BMP")
        // "floor" shader(3)         col(1.000 1.000 1.000 1.000) dif(0.800) amb(0.600) emi(0.000) spc(0.000) power(5.00)

        while( 1 ) {
          p = tokenize( p, token );
          if( strncmp( token, "col(" , 4) == 0 ){
            double r,g,b;
            p2 = tokenize(token+4,t2);
            r = atof( t2 );
            p2 = tokenize(p2,t2);
            g = atof( t2 );
            p2 = tokenize(p2,t2);
            b = atof( t2 );
            m.color_ = __NS_RLR::Color( r, g, b );

            // alpha goes to Refract
            p2 = tokenize(p2,t2);
            m.refract_ = 1.0f - atof( t2 );
          }
          if( strncmp( token, "dif(" , 4) == 0 ){
            p2 = tokenize( token+4, t2 );
            m.diffuse_ = atof( t2 );
          }
          if( strncmp( token, "amb(" , 4) == 0 ){
            p2 = tokenize( token+4, t2 );
            //m.ambient = atof( t2 );
          }
          if( strncmp( token, "emi(" , 4) == 0 ){
            //m.light   = atof( token+4 );
            p2 = tokenize( token+4, t2 );
            m.emitter_ = atof( t2 );
          }
          if( strncmp( token, "spc(" , 4) == 0 ){
            //m.reflect = atof( token+4 );
            p2 = tokenize( token+4, t2 );
            m.reflect_ = atof( t2 );
          }
          if( strncmp( token, "power(" , 6) == 0 ){
            p2 = tokenize( token+6, t2 );
            //m.spe_pow = atof( t2 );
            m.ior_ = atof( t2 );
          }
          if( !p )
            break;
        }
        m.normalize();
        printf("C(%5.3f,%5.3f,%5.3f),", m.color_.x, m.color_.y, m.color_.z);
        printf("E %5.3f, D %5.3f,S %5.3f,R %5.3f\n", m.emitter_, m.diffuse_, m.reflect_, m.refract_ );
        materials_.push_back( m );
      }
    }

    if(strncmp(token,"Object",6) == 0) {
      //object
      p = tokenize(p,token);
      printf("[OBJ] %-20s ",token);
      depth = 1;
      currface = totalface;
      
      vertices_.clear();
      indices_.clear();
      
      while(1) {
        fgets(line,1024,fp);
        p = line;
        p = tokenize(p,token);
        if(strcmp(token,"vertex")==0) {
          //vertex
          p = tokenize(p,token);
          nvert = atoi(token);
          printf("%5d,",nvert);
          depth++;
          for(i=0;i<nvert;i++) {
            fgets(line,1024,fp);
            p = line;
            p = tokenize(p,token);
            x = atof(token);
            p = tokenize(p,token);
            y = atof(token);
            p = tokenize(p,token);
            z = atof(token);
            //printf("%d : % 8.3f % 8.3f % 8.3f\n",i,x,y,z);
            
            vertices_.push_back( __NS_VECTORMATH::Vector( x,y,z ) );
          }
        } else if(strcmp(token,"face")==0) {
          //face
          p = tokenize(p,token);
          nface = atoi(token);
          depth++;
          printf("%5d,",nface);
          for(i=0;i<nface;i++) {
            fgets(line,1024,fp);
            p = line;
            p = tokenize( p, token );
            if( token[0] == '3' ) {
              //p = tokenize( p, token );
              m = a = b = c = 0;
              while( 1 ) {
                p = tokenize( p, token );
                if( token[0] == 'V' ) {
                  p2 = tokenize( token+2, t2 );
                  a  = atoi( t2 );
                  p2 = tokenize( p2, t2 );
                  b  = atoi( t2 );
                  p2 = tokenize( p2, t2 );
                  c  = atoi( t2 );
                }
                if( token[0] == 'M' ) {
                  p2 = tokenize( token+2, t2 );
                  m  = atoi( t2 );
                }
                if( !p )
                  break;
              }
              
              indices_.push_back( a );
              indices_.push_back( b );
              indices_.push_back( c );
              
              poly.setup( vertices_[a], vertices_[b], vertices_[c] ); // &_materials[m] );
              poly.setMaterial( &materials_[m] );
              triangles_.push_back( poly );
              
              totalface++;
            } else if(token[0] == '4') {
              while( 1 ) {
                p  = tokenize( p, token );
                if( token[0] == 'V' ){
                  p2 = tokenize( token+2, t2 );
                  a  = atoi( t2 );
                  p2 = tokenize( p2, t2 );
                  b  = atoi( t2 );
                  p2 = tokenize( p2, t2 );
                  c  = atoi( t2 );
                  p2 = tokenize( p2, t2 );
                  d  = atoi( t2 );
                }
                if( token[0] == 'M' ){
                  p2 = tokenize( token+2, t2 );
                  m  = atoi( t2 );
                }
                if( !p )
                  break;
              }
              
              indices_.push_back( a );
              indices_.push_back( b );
              indices_.push_back( c );
              
              indices_.push_back( a );
              indices_.push_back( c );
              indices_.push_back( d );
              
              poly.setup( vertices_[a], vertices_[b], vertices_[c] ); //, &_materials[m] );
              poly.setMaterial( &materials_[m] );
              triangles_.push_back( poly );
              poly.setup( vertices_[a], vertices_[c], vertices_[d] ); //, &_materials[m] );
              poly.setMaterial( &materials_[m] );
              triangles_.push_back( poly );
              totalface++;
              totalface++;
            }
          }
        } else if(strcmp(token,"BVertex")==0) {
          printf("バイナリ形式はサポートされていません.\n");
          break;
        } else if(token[0] == '}') {
          //オブジェクト終端のチェック
          depth--;
          if(depth == 0){
            printf("%6d triangles\n", totalface - currface);
            break;
          }
        }
      }
    }
  }
  printf("[MQO] %d faces,Done.\n",totalface);
  fclose(fp);
  return true;
}

void ModelHolder::registScene( __NS_RLR::Scene &scene ){
  printf("registing faces to scene.\n");
  for( TriangleList::iterator it = triangles_.begin();
       it != triangles_.end();
       ++it ){
    scene.registObject( &(*it) );
  }
}

__NS_RLRUTIL_END

#if 0
#include <deque>
#include <string>
using namespace std;

typedef int Byteoffset;
class ByteStream {
public:

  enum {
    PGL_FLOAT,
    PGL_UNSIGNED_BYTE,
    PGL_SHORT
  };
  
  std::deque<char>  buf;

  int offset( void ){ return static_cast<int>(buf.size()); }
  
  void save( char *fn ){
    int size = buf.size();
    printf("buf.size %d\n",size);
    char *buffer = (char*)malloc( size );
    for(int i=0;i<size;i++)
      buffer[i] = buf[i];
    FILE *fp = fopen( fn, "wb" );
    fwrite(buffer,1,size,fp);
    fclose(fp);
    free(buffer);
  }
  
  // alignキープするためにpadする.
  void align( int a ){
    int s = buf.size();
    int m = (a - (s%a))%a;
    // dprintf(" align overhead %d bytes\n", m);
    for(int i=0;i<m;i++)
      buf.push_back( 0 );
  }
  
  
  void put_string( const char *str ){
    int i;
    int len = strlen( str );
    for(i=0;i<len;i++)
      buf.push_back( str[i] );
    buf.push_back( 0 ); // NULL TERMINATE
    align(16);
  }
  
  void put_uchar( unsigned char value ){
    buf.push_back( value );
  }
  void put_ushort( unsigned short value ){
    char *p = (char*)&value;
    buf.push_back( p[0] );
    buf.push_back( p[1] );
  }
  void put_int( int value ){
    char *p = (char*)&value;
    buf.push_back( p[0] );
    buf.push_back( p[1] );
    buf.push_back( p[2] );
    buf.push_back( p[3] );
  }
  void put_uint( uint32_t value ){
    char *p = (char*)&value;
    buf.push_back( p[0] );
    buf.push_back( p[1] );
    buf.push_back( p[2] );
    buf.push_back( p[3] );
  }
  void put_float( float value, int type = PGL_FLOAT ){
    
    /*
      vertexattribpointerがサポートするタイプ.
    case PGL_UNSIGNED_BYTE:
    case PGL_SHORT:
    case PGL_HALF_FLOAT_RSX:
    case PGL_FLOAT:
    case PGL_X11_Y11_Z10_RSX:
     */
    
    switch( type ){
      
    case PGL_UNSIGNED_BYTE:{
      unsigned char uc = value * 255.f;
      buf.push_back( uc );
    }break;
      
    case PGL_SHORT:{
      signed short ss = (value + 1.0f) * 32767.5f;
      char *p = (char*)&ss;
      buf.push_back( p[0] );
      buf.push_back( p[1] );
    }break;
      
    default:
    case PGL_FLOAT:{
      char *p = (char*)&value;
      buf.push_back( p[0] );
      buf.push_back( p[1] );
      buf.push_back( p[2] );
      buf.push_back( p[3] );
    }break;
      
    }
  }
  Byteoffset put_pointer( void ){
    Byteoffset offs = buf.size();
    put_int( 0 );
    return offs;
  }
  void adjust_pointer( Byteoffset reloc, int offs = -1 ){
    
    if( offs < 0 )
      offs = buf.size();
    else
      offs = 0;
    // dprintf(" relocate %x for offset %x\n", reloc, offs);
    char *p = (char*)&offs;
    buf[ reloc   ] = p[0];
    buf[ reloc+1 ] = p[1];
    buf[ reloc+2 ] = p[2];
    buf[ reloc+3 ] = p[3];
  }
  void adjust_uint( Byteoffset reloc, uint32_t value ){
    unsigned char *p = (unsigned char*)&value;
    buf[ reloc   ] = p[0];
    buf[ reloc+1 ] = p[1];
    buf[ reloc+2 ] = p[2];
    buf[ reloc+3 ] = p[3];
  }
};

/*
  struct ModelBinary {
    
    void     *top_;
    size_t    size_;
    uint32_t  magic_;
    uint32_t  pad0_;
    
    void     *vertex_;
    uint32_t  vertex_count_;
    
    void     *index_;
    uint32_t  index_count_;
    
    void map( void* src ){
      top_    = src;
      vertex_ = PDISTD::ByteOffset( top_, (uint32_t)vertex_ );
      index_  = PDISTD::ByteOffset( top_, (uint32_t)index_  );
    }
  };
 */
#if 0
int main( int argc, char *argv[] ){
  __NS_RLRUTIL::ModelHolder model;
  model.load( argv[1] );
  
  ByteStream buffer;
  
  buffer.put_pointer();
  Byteoffset size_offs = buffer.put_pointer();
  (void)size_offs;
  buffer.put_uint( 0x00000001 ); // magicてきとう.
  buffer.put_uint( 0 ); // pad0

  int vert_count = model.vertices_.size();
  int inde_count = model.indices_.size();
  
  Byteoffset vertex_ptr = buffer.put_pointer();
  buffer.put_uint( vert_count );
  Byteoffset index_ptr = buffer.put_pointer();
  buffer.put_uint( inde_count );

  buffer.adjust_pointer( vertex_ptr );
  for( int v = 0 ; v < vert_count ; v ++ ){
    buffer.put_float( model.vertices_[v].x );
    buffer.put_float( model.vertices_[v].y );
    buffer.put_float( model.vertices_[v].z );
  }
  buffer.align( 16 );
  buffer.adjust_pointer( index_ptr );
  for( int i = 0 ; i < inde_count ; i ++ ){
    buffer.put_uint( model.indices_[i] );
  }
  buffer.save( argv[2] );
}

#endif
#endif
