#include <iostream>
#include "FileReader.cpp"
#include "common.h"
#include "edge_labeler.h"

#include "g2o/types/slam2d/edge_se2.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/slam2d/edge_se2_twopointsxy.h"
#include "g2o/types/slam2d/edge_se2_lotsofxy.h"

#define MAX_IDS 100000
#define STAR_LENGTH 70
#define OPTIMIZATION_STEPS 30
#define FEW_LEVEL 1
#define MANY_LEVEL 2


G2O_REGISTER_TYPE(VERTEX_SE2, Pose);
G2O_REGISTER_TYPE(VERTEX_XY, Landmark);

std::vector<Pose *> poses;
std::vector<Landmark *> landmarks;
std::vector<g2o::EdgeSE2 *> edgesPoses;
std::vector<g2o::EdgeSE2PointXY *> edgesLandmarks;

void * ids[MAX_IDS];

void readPose(std::vector<std::string> &textline){
  int id = atoi(textline[1].c_str());
  double x = atof(textline[2].c_str());
  double y = atof(textline[3].c_str());
  double theta = atof(textline[4].c_str());
  Pose * p = new Pose(x,y,theta);
  p->setId(id);
  poses.push_back(p);
  ids[id] = (void*) p;
}


void readLandmark(std::vector<std::string> &textline){
  int id = atoi(textline[1].c_str());
  double x = atof(textline[2].c_str());
  double y = atof(textline[3].c_str());
  Landmark * l = new Landmark(x,y);
  l->setId(id);
  landmarks.push_back(l);
  ids[id] = (void*) l;
}


void readEdgePoses(std::vector<std::string> &textline){
  g2o::EdgeSE2 * e = new g2o::EdgeSE2();
  
  int id1 = atoi(textline[1].c_str());
  int id2 = atoi(textline[2].c_str());
  
  e->vertices()[0] = (Pose *)ids[id1];
  e->vertices()[1] = (Pose *)ids[id2];
  
  double d_x = atof(textline[3].c_str());
  double d_y = atof(textline[4].c_str());
  double d_t = atof(textline[5].c_str());
  
  g2o::SE2 measurement(d_x, d_y, d_t);
  e->setMeasurement(measurement);
  
  Eigen::Matrix3d info;
  double i00 = atof(textline[6].c_str());
  double i01 = atof(textline[7].c_str());
  double i02 = atof(textline[8].c_str());
  double i11 = atof(textline[9].c_str());
  double i12 = atof(textline[10].c_str());
  double i22 = atof(textline[11].c_str());
  info <<	i00,	i01,	i02,
    		i01,	i11,	i12,
    		i02,	i12,	i22;
  
  e->setInformation(info);
  
  edgesPoses.push_back(e);
}


void readEdgeLandmark(std::vector<std::string> &textline){
  g2o::EdgeSE2PointXY * e = new g2o::EdgeSE2PointXY();
  
  int id1 = atoi(textline[1].c_str());
  int id2 = atoi(textline[2].c_str());
  
  Pose * p = (Pose *)ids[id1];
  Landmark * l = (Landmark *)ids[id2];
  
  p->addEdge(e);
  l->addEdge(e);
  
  e->vertices()[0] = p;
  e->vertices()[1] = l;
  
  double d_x = atof(textline[3].c_str());
  double d_y = atof(textline[4].c_str());
  
  Eigen::Vector2d measurement(d_x, d_y);
  e->setMeasurement(measurement);
  
  Eigen::Matrix2d info;
  double i00 = atof(textline[5].c_str());
  double i01 = atof(textline[6].c_str());
  double i11 = atof(textline[7].c_str());
  info <<	i00,	i01,
		i01,	i11;

  e->setInformation(info);
  
  edgesLandmarks.push_back(e);
}


void readFromFile(std::string filename){
  FileReader fr(filename);
  if(!fr.is_open()){
    std::cout << "cannot read file " << filename << ", quitting." << std::endl;
    exit(1);
  }
  
  std::vector<std::string> textline;
  fr.readLine(&textline);
  while(fr.good()){
    
    // compare returns 0 if the two strings are equal
    if((textline[0].compare(std::string("VERTEX_SE2"))) == 0){
      readPose(textline);
    }
    else if((textline[0].compare(std::string("VERTEX_XY"))) == 0){
      readLandmark(textline);
    }
    else if((textline[0].compare(std::string("EDGE_SE2"))) == 0){
      readEdgePoses(textline);
    }
    else if((textline[0].compare(std::string("EDGE_SE2_XY"))) == 0){
      readEdgeLandmark(textline);
    }
    
    textline.clear();
    fr.readLine(&textline);
  }
}


void init(){
  // for(unsigned int i=0; i<MAX_IDS; i++){
    
  // }
  
}


int main(int argc, char ** argv){
  std::cout << "----CONDENSE----" << std::endl;
  
  DEBUGMESS("<debug mode enabled>");
  
  if(argc < 2){
    std::cout << "Usage: condense <input_file.g2o>" << std::endl;
    return 0;
  }
  
  init();
  
  readFromFile(argv[1]);
  
#ifdef DEBUGMODE
  // list vertices:
  std::cout << "listing poses:" << std::endl;
  for(unsigned int i=0; i<poses.size(); i++){
    std::cout << "Pose " << poses[i]->id() << " = " << poses[i]->estimate()[0] << "\t" << poses[i]->estimate()[1] << "\t" << poses[i]->estimate()[2] << std::endl;
  }
  
  std::cout << "listing landmarks:" << std::endl;
  for(unsigned int i=0; i<landmarks.size(); i++){
    std::cout << "Lmark " << landmarks[i]->id() << " = " << landmarks[i]->estimate()[0] << "\t" << landmarks[i]->estimate()[1] << std::endl;
  }
  
  std::cout << "listing pose-pose edges:" << std::endl;
  for(unsigned int i=0; i<edgesPoses.size(); i++){
    std::cout << "ID1 = " << edgesPoses[i]->vertices()[0]->id() << "\t ID2 = " << edgesPoses[i]->vertices()[1]->id() << std::endl;
  }
  
  std::cout << "listing pose-landmark edges:" << std::endl;
  for(unsigned int i=0; i<edgesLandmarks.size(); i++){
    std::cout << "ID1 = " << edgesLandmarks[i]->vertices()[0]->id() << "\t ID2 = " << edgesLandmarks[i]->vertices()[1]->id() << std::endl;
  }
#endif // DEBUGMODE
  
  std::cout << "Loaded " << poses.size() << " poses" << std::endl;
  std::cout << "Loaded " << landmarks.size() << " landmarks" << std::endl;
  std::cout << "Loaded " << edgesPoses.size() << " edges pose-pose" << std::endl;
  std::cout << "Loaded " << edgesLandmarks.size() << " edges pose-landmark" << std::endl;
  
  // types definition
  typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
  typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
  
  // allocating the optimizer
  g2o::SparseOptimizer * optimizer = new g2o::SparseOptimizer();
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  
  optimizer->setAlgorithm(solver);
  
  // std::cout << "Writing output file " << argv[2] << "...";
  // std::ofstream g2oout(argv[2]);
  
  // populate graph:
  // poses
  for(unsigned int i=0; i<poses.size(); i++){
    // g2oout << "VERTEX_SE2 " << poses[i]->id() << " ";
    // poses[i]->write(g2oout);
    // g2oout << std::endl;
    optimizer->addVertex(poses[i]);
  }
  
  // edges between poses
  for(unsigned int i=0; i<edgesPoses.size(); i++){
    optimizer->addEdge(edgesPoses[i]);
    
    Pose * p1 = (Pose *) (edgesPoses[i]->vertices()[0]);
    Pose * p2 = (Pose *) (edgesPoses[i]->vertices()[1]);
    // g2oout << "EDGE_SE2" << " " << p1->id() << " " << p2->id() << " ";
    // edgesPoses[i]->write(g2oout);
    // g2oout << std::endl;
  }
  
  // landmarks
  for(unsigned int i=0;i<landmarks.size(); i++){
    optimizer->addVertex(landmarks[i]);
  }
  
  // edges pose-landmark
  for(unsigned int i=0; i<edgesLandmarks.size(); i++){
    optimizer->addEdge(edgesLandmarks[i]);
  }
  
  
  
  // start to generate the stars
  std::vector<Star *> stars;
  Star * s;
  unsigned int inserted = 0;
  for(unsigned int i=0; i<poses.size(); i++){
    if(inserted==0){
      s = new Star();
      s->gauge_index = (STAR_LENGTH/2)-1;
      stars.push_back(s);
    }
    
    if(i>0 && inserted==0){ // must add the previous pose/edge
      i--;
    }
    
    Pose * p = poses[i];
    s->poses.push_back(p);
    s->edgesPoses.push_back(edgesPoses[i]);
    inserted = (inserted+1) % STAR_LENGTH;
    
    for(unsigned int e=0; e<p->_edges.size(); e++){
      Landmark * l = (Landmark *) p->_edges[e]->vertices()[1];
      if(!s->contains(l)){
	s->landmarks.push_back(l);
      }
      s->edgesLandmarks.push_back(p->_edges[e]);
    }
    
  }
  std::cout << "generated " << stars.size() << " stars" << std::endl;
  
  g2o::EdgeLabeler labeler(optimizer);
  
  int num_edgesxy = 0;
  int num_triedges = 0;
  
  for(unsigned int star_index=0; star_index<stars.size(); star_index++){
    std::cout << std::endl;
    std::cout << "analyzing star #" << star_index+1 << std::endl;
    
    Star * s = stars[star_index];
    
    std::cout << "poses: " << s->poses.size() << std::endl;
    std::cout << "landmarks: " << s->landmarks.size() << std::endl;
    
    // push all the estimates
    for(unsigned int i=0; i<s->poses.size(); i++){
      s->poses[i]->push();
    }
    for(unsigned int i=0; i<s->landmarks.size(); i++){
      s->landmarks[i]->push();
    }
    
    // fix the gauge
    s->poses[s->gauge_index]->setFixed(true);
    
    // ready to move stuff
    
    // optimize the local map
    g2o::OptimizableGraph::EdgeSet toOptimize;
    for(unsigned int i=0; i<s->edgesPoses.size(); i++){
      if(star_index<stars.size()-1 || i<s->edgesPoses.size()-1)
      toOptimize.insert(s->edgesPoses[i]);
    }
    for(unsigned int i=0;i<s->edgesLandmarks.size(); i++){
      toOptimize.insert(s->edgesLandmarks[i]);
    }
    optimizer->initializeOptimization(toOptimize);
    int optim_result = optimizer->optimize(OPTIMIZATION_STEPS);
    if(optim_result<1){
      std::cout << "!!! optimization failed !!!" <<std::endl;
      // don't label according to the optimized graph.
      
      // go back to the unoptimized state
      // pop all the estimates
      for(unsigned int i=0; i<s->poses.size(); i++){
	s->poses[i]->pop();
      }
      for(unsigned int i=0; i<s->landmarks.size(); i++){
	s->landmarks[i]->pop();
      }
      
      // unfix the gauge
      s->poses[s->gauge_index]->setFixed(false);
      
      continue;
      
    }
    
    // shared variables:
    // the first and the last pose are always shared (except for the first and the last star)
    if(star_index>0){
      g2o::EdgeSE2 * edge = new g2o::EdgeSE2();
      edge->vertices()[0] = s->poses[s->gauge_index];
      edge->vertices()[1] = s->poses[0];
      s->edgesShared.insert(edge);
    }
    if(star_index<stars.size()-1){
      g2o::EdgeSE2 * edge = new g2o::EdgeSE2();
      edge->vertices()[0] = s->poses[s->gauge_index];
      edge->vertices()[1] = s->poses[s->poses.size()-1];
      s->edgesShared.insert(edge);
    }
    
    // look for shared xy points
    std::vector<unsigned int> shared;
    unsigned int count = 0;
    for(unsigned int i=0; i<s->landmarks.size(); i++){
      Landmark * l = s->landmarks[i];
      for(unsigned int e=0; e<l->_edges.size(); e++){
	Pose * p = (Pose *) l->_edges[e]->vertices()[0];
	if(!s->contains(p)){
	  shared.push_back(i);
	  count++;
	  break;
	}
      }
    }
    
    std::cout << "shared XY: " << shared.size() << std::endl;
    
    // create condensed edges for the shared variables
    if(count > 0){
      
      g2o::EdgeSE2LotsOfXY * edge = new g2o::EdgeSE2LotsOfXY();
      edge->setSize(1 + shared.size());
      edge->vertices()[0] = s->poses[s->gauge_index];
      for(unsigned int i=0; i<shared.size(); i++){
	edge->vertices()[1+i] = s->landmarks[shared[i]];
      }
      
      s->edgesShared.insert(edge);
      // std::cout << "edge: ";
      // edge->write(std::cout);
      // std::cout << std::endl;
	  
      
    }
    
    labeler.labelEdges(s->edgesShared);
    
    
    
    // create condensed measurements for the local variables
    for(unsigned int i=0; i<s->poses.size(); i++){
      
      if(i==s->gauge_index) continue;
      
      g2o::EdgeSE2 * edge = new g2o::EdgeSE2();
      if(i<s->gauge_index){
    	edge->vertices()[0] = s->poses[i];
    	edge->vertices()[1] = s->poses[s->gauge_index];
      }
      else{
    	edge->vertices()[0] = s->poses[s->gauge_index];
    	edge->vertices()[1] = s->poses[i];
      }
      
      s->edgesCondensed.insert(edge);
    }
    
        
    // for(unsigned int i=0; i<s->landmarks.size(); i++){
    //   for(unsigned int j=i+1; j<s->landmarks.size(); j++){
    // 	for(unsigned int k=j+1; k<s->landmarks.size(); k++){
    // 	  g2o::EdgeSE2LotsOfXY * edge = new g2o::EdgeSE2LotsOfXY();
    // 	  edge->setSize(4);
    // 	  edge->vertices()[0] = s->poses[s->gauge_index];
    // 	  edge->vertices()[1] = s->landmarks[i];
    // 	  edge->vertices()[2] = s->landmarks[j];
    // 	  edge->vertices()[3] = s->landmarks[k];
    // 	  s->edgesCondensed.insert(edge);
    // 	}
    //   }
    // }
    
    g2o::EdgeSE2LotsOfXY * edge = new g2o::EdgeSE2LotsOfXY();
    edge->setSize(1 + s->landmarks.size());
    edge->vertices()[0] = s->poses[s->gauge_index];
    for(unsigned int i=0; i<s->landmarks.size(); i++){
      edge->vertices()[1 + i] = s->landmarks[i];
    }
    s->edgesCondensed.insert(edge);
    
    labeler.labelEdges(s->edgesCondensed);
    
    // pop all the estimates
    for(unsigned int i=0; i<s->poses.size(); i++){
      s->poses[i]->pop();
    }
    for(unsigned int i=0; i<s->landmarks.size(); i++){
      s->landmarks[i]->pop();
    }
    
    // unfix the gauge
    s->poses[s->gauge_index]->setFixed(false);
    
  }
  
  //std::cout << "created " << num_edgesxy << " binary edges and " << num_triedges << " tri-edges" << std::endl;
  
  std::cout << "generating file few.g2o...";
  
  for(unsigned int i=0; i<stars.size(); i++){
    Star * s = stars[i];
    
    for(std::set<g2o::OptimizableGraph::Edge *>::iterator it=s->edgesShared.begin(); it!=s->edgesShared.end(); it++){
      (*it)->setLevel(FEW_LEVEL);
      optimizer->addEdge((*it));
    }
  }
  
  std::ofstream ofs1("few.g2o");
  optimizer->save(ofs1, FEW_LEVEL);
  
  std::cout << "\tDONE" << std::endl;
  
  std::cout << "generating file many.g2o...";
  
  for(unsigned int i=0; i<stars.size(); i++){
    Star * s = stars[i];
    
    for(std::set<g2o::OptimizableGraph::Edge *>::iterator it=s->edgesCondensed.begin(); it!=s->edgesCondensed.end(); it++){
      (*it)->setLevel(MANY_LEVEL);
      optimizer->addEdge((*it));
    }
  }
  
  std::ofstream ofs2("many.g2o");
  optimizer->save(ofs2, MANY_LEVEL);
  
  std::cout << "\tDONE" << std::endl;
    
  return 0;
}
