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



class Pose : public g2o::VertexSE2{
 public:
  std::vector<g2o::EdgeSE2PointXY *> _edges;
    
  Pose(){};
  Pose(double x, double y, double theta);
  
  void addEdge(g2o::EdgeSE2PointXY * toAdd){_edges.push_back(toAdd);};
  
  /* void setEstimateData(double * data){ */
  /*   _estimate[0] = data[0]; */
  /*   _estimate[1] = data[1]; */
  /*   _estimate[2] = data[2]; */
  /* }; */
};

Pose::Pose(double x, double y, double theta){
  _estimate = g2o::SE2(x,y,theta);
}


class Landmark : public g2o::VertexPointXY{
 public:
  bool alreadyInserted;
  
  std::vector<g2o::EdgeSE2PointXY *> _edges;
  
  Landmark(){_estimate[0]=0; _estimate[1]=0; alreadyInserted=false;};
  Landmark(double x, double y){_estimate[0]=x; _estimate[1]=y; alreadyInserted=false;};
  void addEdge(g2o::EdgeSE2PointXY * toAdd){_edges.push_back(toAdd);};
};


class Star{
 public:
  std::vector<Pose *> poses;
  std::vector<g2o::EdgeSE2 *> edgesPoses;
  
  int gauge_index;
  
  std::vector<Landmark *> landmarks;
  std::vector<g2o::EdgeSE2PointXY *> edgesLandmarks;
  
  std::set<g2o::OptimizableGraph::Edge *> edgesShared;
  std::set<g2o::OptimizableGraph::Edge *> edgesCondensed;
  
  Star(){gauge_index = 0;};
  
  bool contains(Landmark * l){
    for(unsigned int i=0; i<landmarks.size(); i++){
      if(landmarks[i]->id() == l->id()) return true;
    }
    return false;
  }
  
  bool contains(Pose *p){
    for(unsigned int i=0; i<poses.size(); i++){
      if(poses[i]->id() == p->id()) return true;
    }
    return false;
  }
  
};


#endif // _COMMON_H
