#ifndef _COMMON_H
#define _COMMON_H

#include <vector>
#include "g2o/types/slam2d/edge_se2.h"
#include "g2o/types/slam2d/edge_se2_pointxy.h"
#include "g2o/types/slam2d/vertex_se2.h"
#include "g2o/types/slam2d/vertex_point_xy.h"


#ifdef DEBUGMODE
#define DEBUGMESS(X) std::cout<<X<<std::endl
#else
#define DEBUGMESS(X)
#endif

struct VertexWrapper{
  g2o::OptimizableGraph::Vertex * vertex;
  std::vector<g2o::OptimizableGraph::Edge *> edges;
};

struct Star2D{
  std::vector<VertexWrapper *> poses;
  std::vector<g2o::EdgeSE2 *> edgesPoses;
  
  int gauge_index;
  
  std::vector<VertexWrapper *> landmarks;
  std::vector<g2o::OptimizableGraph::Edge *> edgesLandmarks;
  
  std::set<g2o::OptimizableGraph::Edge *> edgesShared;
  std::set<g2o::OptimizableGraph::Edge *> edgesCondensed;
  
  Star2D(){gauge_index = 0;};
  
  void pushState(){
    for(unsigned int i=0; i<this->poses.size(); i++){
      this->poses[i]->vertex->push();
    }
    for(unsigned int i=0; i<this->landmarks.size(); i++){
      this->landmarks[i]->vertex->push();
    }
  }
  
  void popState(){
    // pop all the estimates
    for(unsigned int i=0; i<this->poses.size(); i++){
      this->poses[i]->vertex->pop();
    }
    for(unsigned int i=0; i<this->landmarks.size(); i++){
      this->landmarks[i]->vertex->pop();
    }
  }
  
  void fixGauge(){ 
    // fix the gauge
    this->poses[this->gauge_index]->vertex->setFixed(true);
  }
  
  void unfixGauge(){
    // free the gauge
    this->poses[this->gauge_index]->vertex->setFixed(false);
  }

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

#endif // _COMMON_H

