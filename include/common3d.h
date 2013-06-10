#ifndef _COMMON3D_H
#define _COMMON3D_H

#include <vector>
#include "g2o/types/slam3d/edge_se3.h"
#include "g2o/types/slam3d/edge_se3_pointxyz.h"
#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/types/slam3d/vertex_pointxyz.h"


#ifdef DEBUGMODE
#define DEBUGMESS(X) std::cout<<X<<std::endl
#else
#define DEBUGMESS(X)
#endif

struct VertexWrapper{
  g2o::OptimizableGraph::Vertex * vertex;
  std::vector<g2o::EdgeSE3PointXYZ *> edges;
};

struct Star3D{
  std::vector<VertexWrapper *> poses;
  std::vector<g2o::EdgeSE3 *> edgesPoses;
  
  int gauge_index;
  
  std::vector<VertexWrapper *> landmarks;
  std::vector<g2o::EdgeSE3PointXYZ *> edgesLandmarks;
  
  std::set<g2o::OptimizableGraph::Edge *> edgesShared;
  std::set<g2o::OptimizableGraph::Edge *> edgesCondensed;
  
  Star3D(){gauge_index = 0;};
  
  bool contains(VertexWrapper * v){
    for(unsigned int i=0; i<landmarks.size(); i++){
      if(landmarks[i]->vertex->id() == v->vertex->id()) return true;
    }
    for(unsigned int i=0; i<poses.size(); i++){
      if(poses[i]->vertex->id() == v->vertex->id()) return true;
    }
    return false;
  }
  
  bool contains(g2o::HyperGraph::Vertex * v){
    for(unsigned int i=0; i<landmarks.size(); i++){
      if(landmarks[i]->vertex->id() == v->id()) return true;
    }
    for(unsigned int i=0; i<poses.size(); i++){
      if(poses[i]->vertex->id() == v->id()) return true;
    }
    return false;
  }
  
  VertexWrapper * getWrapper(g2o::HyperGraph::Vertex * v){
    for(unsigned int i=0; i<landmarks.size(); i++){
      if(landmarks[i]->vertex->id() == v->id()) return landmarks[i];
    }
    for(unsigned int i=0; i<poses.size(); i++){
      if(poses[i]->vertex->id() == v->id()) return poses[i];
    }
    return 0;
  }
  
};



// finds the VertexWrapper pointer in the given list.
// RETURNS 0 if not found
VertexWrapper * findWrapper(std::vector<VertexWrapper *> wrappers, g2o::HyperGraph::Vertex * toFind){
  for(unsigned int i=0; i<wrappers.size(); i++){
    if (wrappers[i]->vertex->id() == toFind->id()) return wrappers[i];
  }
  return 0;

}

#endif // _COMMON3D_H
