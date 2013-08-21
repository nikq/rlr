/*
 * modelholder.h
 *
 * MQO file loader.
 * copyrights(c) Hajime UCHIMURA / nikq.
 *
 */
#ifndef __MQOLOAD_H
#define __MQOLOAD_H

#include <stdio.h>
#include <string.h>

#include <vector>

#include "vectormath.h"
#include "scene.h"


#define __NS_RLRUTIL       RLRUTIL
#define __NS_RLRUTIL_BEGIN namespace __NS_RLRUTIL {
#define __NS_RLRUTIL_END   }

__NS_RLRUTIL_BEGIN

class ModelHolder {
//private:
public:

  typedef std::vector< __NS_VECTORMATH::Vector > VertexList;
  typedef std::vector< unsigned int            > IndexList;
  typedef std::vector< __NS_RLR::Triangle      > TriangleList;
  typedef std::vector< __NS_RLR::SceneMaterial > MaterialList;

  VertexList vertices_;
  IndexList  indices_;
  
  TriangleList triangles_; // !! instance holder !!
  MaterialList materials_;
  
public:

  ModelHolder(){ ; }
  virtual ~ModelHolder(){ ; }
  
  char *tokenize( char *s, char *token );
  bool load( const char *fn );          // false for error
  void clear(){ 
    triangles_.clear();
    materials_.clear();
  }
  void registScene( __NS_RLR::Scene &scene );
};

__NS_RLRUTIL_END

#endif
