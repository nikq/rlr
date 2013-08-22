/*
 * bvh.h
 *
 * implements AXIS-ALIGNED Bounding Volume Hyerarchy.
 * copyrights(c) Hajime UCHIMURA / nikq.
 */

#ifndef __BVH_H
#define __BVH_H

#include "vectormath.h" // for vectormath
#include <iterator>
#define __NS_BVH       RLR
#define __NS_BVH_BEGIN namespace __NS_BVH {
#define __NS_BVH_END   }

__NS_BVH_BEGIN

template< typename ObjectList >
class BVHNode{
public:
  BVHNode(){;}
  virtual ~BVHNode(){ leaf_.clear(); }
  
  bool isLeaf_;
  AABB aabb_;
  
  //BVHLeaf< ObjectList > * leaf_;
  ObjectList              leaf_;
  BVHNode< ObjectList > * left_;
  BVHNode< ObjectList > * right_;
};

template< typename ObjectList >
class BVH{
public:
  
  typedef BVH    < ObjectList >  self;
  typedef BVHNode< ObjectList >  Node;
  
private:
  
  Node * root_;
  
public:
  int traverse_count_;
  int intersect_count_;
  int ray_count_;

  BVH( void ) : root_(NULL) { ; }
  virtual ~BVH( void ){
    //dump( root_ );
    release( root_ );
  }
  
  void dump( Node *node ){
    if( !node ){
      printf("NULL\n");
      return;
    }
    if( node->isLeaf_ ){
      printf("%p: leaf, %d objs\n", node, node->leaf_.size());
      return;
    }
    //printf("%p: node, left %p, right %p\n", node, node->left_, node->right_ );
    dump( node->left_  );
    dump( node->right_ );
  }
  
  void release( Node *node ){
    if( !node ){
      return;
    }
    if( !node->isLeaf_ ){
      if( node->left_ ){ release( node->left_ ); delete node->left_;  }
      if( node->right_){ release( node->right_ );delete node->right_; }
    }
  }
  
  Node *buildBvhNode( ObjectList & objects, AABB& bb, int depth ){
    
    if( objects.size() <= 4 || depth >= 16 ){
      // make it leaf.
      Node *node = new Node;
      node->isLeaf_ = true;
      std::copy( objects.begin(), objects.end(), std::back_inserter( node->leaf_ ) );
      node->aabb_.clear();
      for( typename ObjectList::iterator it = objects.begin(); it != objects.end(); ++it ) {
        AABB aabb = (*it)->getAABB();
        node->aabb_.expand( aabb );
      }
      return node;
    }
    
    // check median.
    __NS_VECTORMATH::Vector size = bb.hi() - bb.lo();
    int    axis = 0;
    double width  = size.x;
    double median = bb.lo().x + size.x / 2.;
    if( width <= size.y ){ axis = 1; width = size.y; median = bb.lo().y + size.y / 2.; }
    if( width <= size.z ){ axis = 2; width = size.z; median = bb.lo().z + size.z / 2.; }
    // printf("max axis %d, %f,%f\n",axis, width, median);
    
    ObjectList left, right;
    
    //-------------------------------
    AABB bbleft = bb;//左用
    bbleft.hi().set(axis, median);
    
    AABB bbright = bb;//右用
    bbright.lo().set(axis, median);
    //-------------------------------

    AABB mybb;

    median = 0.;
    double divider = 0.;
    for( typename ObjectList::iterator it = objects.begin(); it != objects.end(); ++it ) {
      AABB aabb = (*it)->getAABB();
      __NS_VECTORMATH::Vector center = (aabb.lo() + aabb.hi()) * 0.5;
      median += center.get( axis );
      divider+= 1.;
    }
    median /= divider;
    for( typename ObjectList::iterator it = objects.begin(); it != objects.end(); ++it ) {
      AABB aabb = (*it)->getAABB();
      mybb.expand( aabb );
      __NS_VECTORMATH::Vector center = (aabb.lo() + aabb.hi()) * 0.5;
      if( center.get(axis) <= median ){
        left.push_back( *it );
      }else{
        right.push_back( *it );
      }
    }
    if( !left.size() ){
      // 右しか無い.
      return buildBvhNode( right, bbright, depth );
    }else if( !right.size() ){
      return buildBvhNode( left, bbleft, depth );
    }else{
      // printf("left %d,right %d\n", left.size(), right.size());
      Node *node = new Node;
      node->isLeaf_ = false;
      node->aabb_   = mybb;
      node->left_   = left.size() ? buildBvhNode( left , bbleft , depth+1 ) : NULL;
      node->right_  = right.size()? buildBvhNode( right, bbright, depth+1 ) : NULL;
      return node;
    }
    /* NOTREACHED */
  }
  
  void build( ObjectList& objects ){
    AABB bb;
    for( typename ObjectList::iterator it = objects.begin(); it != objects.end(); ++it ){
      bb.expand( (*it)->getAABB() );
    }
    
    root_ = buildBvhNode( objects, bb, 0 );
    //dump( root_ );
  }
  
  void traverseNode ( Node* node, const Ray& ray, double tmin, double tmax, Intersection &isect ) {
    
    //printf("traverse %f-%f %d,%f\n",tmin,tmax,isect.hit(),isect.distance());
    
    traverse_count_ ++;
    // NULLなら帰る
    if( !node )
      return;
    
    if( node->isLeaf_ ){
      for( typename ObjectList::iterator it = node->leaf_.begin(); it != node->leaf_.end(); ++it ){
        __NS_RLR::Intersection result;
        (*it)->testRAY( ray, result );
        isect.update( result );
        intersect_count_ ++;
      }
      return;
    }
    
    bool l_hit = false;
    bool r_hit = false;
    bool   hit    = isect.hit();
    double hitmin = isect.distance();
    
    if( hit && hitmin < tmin )
      return;
    
    
    double lo, hi;
    if( !node->aabb_.clipRay( ray, lo, hi ) || lo > tmax || hi < tmin )
      return;
    if( hit && hitmin < lo )
      return;
    
    double l_lo(0.), l_hi(0.);
    double r_lo(0.), r_hi(0.);
    
    if( node->left_  ){
      l_hit = node->left_ ->aabb_.clipRay( ray, l_lo, l_hi );
      if( l_hit ){
        if( hit ){
          if( hitmin < l_lo )
            l_hit = 0;
          else if( l_hi > hitmin)
            l_hi = hitmin;
        }
      }
    }
    if( node->right_ ){
      r_hit = node->right_->aabb_.clipRay( ray, r_lo, r_hi );
      if( r_hit ){
        if( hit ){
          if( hitmin < r_lo )
            r_hit = 0;
          else if( r_hi > hitmin)
            r_hi = hitmin;
        }
      }
    }
    
    if( l_hit && !r_hit ) {
      traverseNode( node->left_ , ray, l_lo, l_hi, isect );
    }else if( r_hit && !l_hit ){
      traverseNode( node->right_, ray, r_lo, r_hi, isect );
    }else if( l_hit && r_hit ){
      if( l_lo < r_lo ){
        traverseNode( node->left_ , ray, l_lo, l_hi, isect );
        if( isect.hit() && r_hi > isect.distance() ) r_hi = isect.distance();
        traverseNode( node->right_, ray, r_lo, r_hi, isect );
      }else{
        traverseNode( node->right_, ray, r_lo, r_hi, isect );
        if( isect.hit() && l_hi > isect.distance() ) l_hi = isect.distance();
        traverseNode( node->left_ , ray, l_lo, l_hi, isect );
      }
    }
  }
  bool testRAY  ( const Ray& ray, Intersection &isect ) {
    isect.clear();
    traverse_count_ = 0;
    intersect_count_= 0;
    ray_count_      ++;
    
    traverseNode( root_, ray, 0, FLT_MAX, isect );
    return isect.hit();
  }

};

__NS_BVH_END

#endif
