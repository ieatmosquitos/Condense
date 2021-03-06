#include <iostream>
#include <string>
//#include <stdlib.h> // atoi
#include <cstdlib> // atoi, itoa

#include "FileReader.cpp"
#include "common.h"
#include "edge_labeler.h"
#include "clustering.cpp"
#include "stability.h"

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/slam2d/types_slam2d.cpp"

#define MAX_IDS 100000
#define FEW_LEVEL 2
#define MANY_LEVEL 3

// types definition
typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

std::vector<VertexWrapper *> poses;
std::vector<VertexWrapper *> landmarks;
std::vector<g2o::EdgeSE2 *> edgesPoses;
std::vector<g2o::OptimizableGraph::Edge *> edgesLandmarks;

void * ids[MAX_IDS];

int _starLength;
int _optimizationSteps;

// studying impact of poses edges on optimizability
bool _createPosesEdges;

// clustering stuff
bool _clusterize; // if false, no clusters are made and ONLY BINARY EDGES ARE CREATED
int _max_clusters;
int _max_landmarks_per_edge; // set lesser than 1 to let the edge be as big as it wants


// robust kernel
g2o::RobustKernel * robust_kernel;


int makeClusters(std::vector<VertexWrapper *> &vertices, int * labels){
  // check special cases
  if(vertices.size() < 1) return 0;
  
  if(vertices.size()<2){
    labels[0] = 0;
    return 1;
  }
  
  if(!_clusterize){
    for(unsigned int i=0; i<vertices.size(); i++){
      labels[i] = 0;
    }
    return 1;
  }
  
  double values[2*vertices.size()];
  double * means = (double *) malloc(sizeof(double));
  
  for(unsigned int i=0; i<vertices.size(); i++){
    g2o::VertexPointXY * l = (g2o::VertexPointXY *) vertices[i]->vertex;
    unsigned int index = i*2;
    Eigen::Vector2d est = l->estimate();
    values[index] = est[0];
    values[index + 1] = est[1];
  }
  
  int clusters = findClusters(values, 2, labels, &means, vertices.size(), _max_clusters);
  return clusters;
}



void getSharedEdges(Star2D * s, std::vector<VertexWrapper *> &shared, std::vector<VertexWrapper *> &local){
  // look for shared xyz points
  for(unsigned int i=0; i<s->landmarks.size(); i++){
    VertexWrapper * l = s->landmarks[i];
    bool isShared = false;
    for(unsigned int e=0; e<l->edges.size(); e++){
      g2o::HyperGraph::Vertex * p =  l->edges[e]->vertices()[0];
      if(!s->contains(p)){
	isShared = true;
	break;
      }
    }
    if(isShared){
      shared.push_back(l);
    }
    else{
      local.push_back(l);
    }
  }
  std::cout << "shared XY: " << shared.size() << std::endl;
  std::cout << "shared+local = " << shared.size() + local.size() << std::endl;
}



// 
void insertSharedEdges(Star2D * s, std::vector<VertexWrapper *> &shared, int clusters, int * labels){
  for(unsigned int c=0; c<clusters; c++){
    std::vector<VertexWrapper *> to;
    for(unsigned int i=0; i<shared.size(); i++){
      if(labels[i] == c){
	to.push_back(shared[i]);
      }
    }
    if(_max_landmarks_per_edge < 1 || _max_landmarks_per_edge > to.size()){
      // create an edge connecting the gauge to all the landmarks in the cluster
      g2o::EdgeSE2LotsOfXY * edge = new g2o::EdgeSE2LotsOfXY();
      edge->setSize(1+to.size());
      edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
      for(unsigned int i=0; i<to.size(); i++){
	edge->vertices()[1+i] = to[i]->vertex;
      }
      s->edgesShared.insert(edge);
      s->edgesCondensed.insert(edge);
    }
    else{
      // create many edges connecting the gauge to subsets of the cluster
      unsigned int count = 0;
      unsigned int last_index = to.size() - 1;
      while(count < to.size()){
	unsigned int so_many = _max_landmarks_per_edge;
	if(count + so_many > to.size()){
	  so_many = to.size() - count;
	}
	
	g2o::EdgeSE2LotsOfXY * edge = new g2o::EdgeSE2LotsOfXY();
	edge->setSize(1 + so_many);
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	for(unsigned int v=0; v<so_many; v++){
	  edge->vertices()[1+v] = to[count+v]->vertex;
	}
	s->edgesShared.insert(edge);
	s->edgesCondensed.insert(edge);
	
	count += so_many;
      }
    }
  }
}



void insertLocalEdges(Star2D * s,  std::vector<VertexWrapper *> &local){
  for(unsigned int i=0; i<local.size(); i++){
    g2o::EdgeSE2LotsOfXY * e = new g2o::EdgeSE2LotsOfXY();
    e->setSize(2);
    e->vertices()[0] = s->poses[s->gauge_index]->vertex;
    e->vertices()[1] = local[i]->vertex;
    s->edgesCondensed.insert(e);
  }
}


void prepareOutputName(std::string original, std::string & fewout, std::string & manyout){
  
  fewout = "few.g2o";
  manyout = "many.g2o";
  
  std::stringstream ss;
  ss << original.substr(0,original.length()-4) << "-sl_" << _starLength << "-mc_" << _max_clusters << "-ml_" << _max_landmarks_per_edge;
  if(!_clusterize){
    ss << "-noclust";
  }
  if(!_createPosesEdges){
    ss << "-nopos";
  }
  
  std::string basename = ss.str();
  //  std::cout << "basename = " << basename << std::endl;

  fewout = basename;
  manyout = basename;
  
  fewout.append("_few.g2o");
  manyout.append("_many.g2o");
  
}


// bool checkIntegrity(Star2D * s){
//   for(unsigned int i=0; i<s->landmarks.size(); i++){
//     VertexWrapper * vw = s->landmarks[i];
    
//     unsigned int count = 0;
    
//     for(unsigned int j=0; j<vw->edges.size(); j++){
//       if(s->contains(vw->edges[j]->vertices()[0])){
// 	count++;
//       }
//     }
    
//     if(count<1) return false;
//   }
//   return true;
// }


void purgeLonelyLandmarks(Star2D * s){
  
  for(unsigned int l=0; l<s->landmarks.size(); l++){
    VertexWrapper * v = s->landmarks[l];
    
    // count how many poses see this landmark within the same star
    unsigned int count=0;
    bool isBearing = false;
    bool areAllBearings = true;
     
    for(unsigned int i = 0; i<v->edges.size(); i++){
      
      g2o::OptimizableGraph::Edge * e = v->edges[i];
      
      if(s->contains(e->vertices()[0])){
 	count++;
 	g2o::EdgeSE2PointXYBearing * b = dynamic_cast<g2o::EdgeSE2PointXYBearing *>(e);
 	if(b!=NULL){
 	  isBearing = true;
 	}
	else{
	  areAllBearings = false;
	}
      }
    }
    
    bool remove_it = false;
     
    // std::cout << "landmark " << v->vertex->id() << " has " << v->edges.size() << " edges" << std::endl;
    // std::cout << "seen by " << count << " poses inside the star" << std::endl;
    if(count < 1){
      std::cerr << "implementation is bugged! check the purgeLonelyLandmarks() function" << std::endl;
      exit(11);
    }
     
    if(count < 2 && isBearing){
      remove_it = true;
    }
    else{ // check if the point is stable enough
      
      if(areAllBearings){
	// build a linear system gathering the informations from the bearings measured from the poses
	
	linSystem linsys;
	for(unsigned int i = 0; i<v->edges.size(); i++){
      	  g2o::OptimizableGraph::Edge * e = v->edges[i];
	  if(s->contains(e->vertices()[0])){
	    g2o::VertexSE2 * pose = dynamic_cast<g2o::VertexSE2 *>(e->vertices()[0]);
	    if(!pose){
	      std::cerr << "watch out! the code is bugged. Look for this sentence in condense.cpp" << std::endl;
	      exit(20);
	    }
	    
	    g2o::EdgeSE2PointXYBearing * b = (g2o::EdgeSE2PointXYBearing *) e;
	    linsys.addConstraint(pose->estimate()[0], pose->estimate()[1], pose->estimate()[2], b->measurement());
	  }
	  
      	}
	
	if(!linsys.checkStability()){
	  remove_it = true;
	}
	 
      }
    }
     
    if(remove_it){
      // should remove this landmark from the star
       
      // search the edge in the star that leads to this landmark
      for(unsigned int i=0; i<s->edgesLandmarks.size(); i++){
    	if(s->edgesLandmarks[i]->vertices()[1] == v->vertex){
    	  s->edgesLandmarks.erase(s->edgesLandmarks.begin()+i);
    	  i--;
    	}
      }
       
      s->landmarks.erase(s->landmarks.begin()+l);
      l--;
       
    }
  }
}


void init(int argc, char** argv){
  
  for(unsigned int i=0; i<MAX_IDS; i++){
    ids[i] = 0;
  }
  
  // initialize options
  _starLength = 70;
  _optimizationSteps = 100;
  _createPosesEdges = true;
  _clusterize = true;
  _max_clusters = 6;
  _max_landmarks_per_edge = 5;
  
  // check options from command line
  for(unsigned int i=2; i<argc; i++){
    std::string option(argv[i]);
    
    bool known = false;
    
    if(option.compare("-os") == 0){
      known = true;
      i++;
      if(i == argc){
	std::cerr << "ERROR: no value specified for option " << option << std::endl;
	exit(1);
      }
      _optimizationSteps = atoi(argv[i]);
    }
    
    else if(option.compare("-sl") == 0){
      known = true;
      i++;
      if(i == argc){
	std::cerr << "ERROR: no value specified for option " << option << std::endl;
	exit(1);
      }
      _starLength = atoi(argv[i]);
    }
    
    else if(option.compare("-nopos") == 0){
      known = true;
      _createPosesEdges = false;
    }
    
    else if(option.compare("-noclust") == 0){
      known = true;
      _clusterize = false;
    }
    
    else if(option.compare("-mc") == 0){
      known = true;
      i++;
      if(i == argc){
	std::cerr << "ERROR: no value specified for option " << option << std::endl;
	exit(1);
      }
      _max_clusters = atoi(argv[i]);
    }
    
    else if(option.compare("-ml") == 0){
      known = true;
      i++;
      if(i == argc){
	std::cerr << "ERROR: no value specified for option " << option << std::endl;
	exit(1);
      }
      _max_landmarks_per_edge = atoi(argv[i]);
    }
    
    if(!known){
      std::cerr << "ERROR: unknown command: " << option << std::endl;
      exit(1);
    }
  }

  
  // robust kernel
  robust_kernel = new g2o::RobustKernelCauchy();
}

int main(int argc, char ** argv){
  
  std::cout << "----CONDENSE 2D----" << std::endl;
  
  // check args
  if(argc<2){
    std::cout << "Usage: condense <input_file.g2o> <options>" << std::endl;
    std::cout << "\toptions:" << std::endl;
    std::cout << "\t-os N\tsets the optimization steps for the local maps to N" << std::endl;
    std::cout << "\t-sl N\tsets the length of the stars to N steps" << std::endl;
    std::cout << "\t-nopos\tinhibits the creation of pose-pose edges" << std::endl;
    std::cout << "\t-noclust\t doesn't clusterize the landmarks" << std::endl;
    std::cout << "\t-mc N\tsets the maximum number of clusters to N" << std::endl;
    std::cout << "\t-ml N\t sets the maximum number of landmarks per edge to N" << std::endl;
    
    return 0;
  }
  
  init(argc,argv);
  
  std::cout << "INFO: _starLength = " << _starLength << std::endl;
  std::cout << "INFO: _optimizationSteps = " << _optimizationSteps << std::endl;
  std::cout << "INFO: _createPosesEdges = " << _createPosesEdges << std::endl;
  std::cout << "INFO: _clusterize = " << _clusterize << std::endl;
  std::cout << "INFO: _max_clusters = " << _max_clusters << std::endl;
  std::cout << "INFO: _max_landmarks_per_edge = " << _max_landmarks_per_edge << std::endl;
  
  // allocate the optimizer
  g2o::SparseOptimizer * optimizer = new g2o::SparseOptimizer();
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  // g2o::OptimizationAlgorithmGaussNewton * solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  g2o::OptimizationAlgorithmLevenberg * solver = new g2o::OptimizationAlgorithmLevenberg(blockSolver);
  
  optimizer->setAlgorithm(solver);
  blockSolver->setOptimizer(optimizer);
  
  // load the graph
  std::ifstream ifs(argv[1]);
  if(!ifs){
    std::cout << "Cannot open file " << argv[1] << std::endl;
    return 1;
  }
  if(!optimizer->load(ifs)){
    std::cout << "Error loading graph" << std::endl;
    return 2;
  }
  std::cout << "Loaded " << optimizer->vertices().size() << " vertices" << std::endl;
  std::cout << "Loaded " << optimizer->edges().size() << " edges" << std::endl;
  
  
  optimizer->initializeOptimization();
  optimizer->computeActiveErrors();
  
  std::cout << "Initial chi2 = " << optimizer->activeChi2() << std::endl;
  
  optimizer->computeInitialGuess();
  
  
  // populate the structures
  for(g2o::HyperGraph::EdgeSet::iterator it=optimizer->edges().begin(); it!=optimizer->edges().end(); it++){
    g2o::HyperGraph::Edge * e = (*it);
    
    bool isEdgePose = false;
    {
      g2o::EdgeSE2 * p = dynamic_cast<g2o::EdgeSE2 *>(e);
      if(p!=NULL){
	isEdgePose = true;
      }
      else{
	g2o::OptimizableGraph::Edge * p = dynamic_cast<g2o::EdgeSE2PointXY *>(e);
	if(p==NULL){
	  p = dynamic_cast<g2o::EdgeSE2PointXYBearing *>(e);
	}
	if(p==NULL){
	  // don't load this edge
	  continue;
	}
      }
    }
    for (std::vector<g2o::HyperGraph::Vertex*>::const_iterator vit = e->vertices().begin(); vit != e->vertices().end(); ++vit) {
      g2o::HyperGraph::Vertex * v = (*vit);
      if(ids[v->id()] != 0){
      	VertexWrapper * prev = (VertexWrapper *) ids[v->id()];
      	if(v->id() != prev->vertex->id()){
      	  std::cout << "ERORE!" << std::endl;
      	  return 3;
      	}
	if(!isEdgePose){
	  prev->edges.push_back((g2o::OptimizableGraph::Edge *)e);
	}
      }
      else{
      	VertexWrapper * vw = new VertexWrapper();
	vw->vertex = (g2o::OptimizableGraph::Vertex *) v;
      	ids[v->id()] = vw;
	if(!isEdgePose){
	  vw->edges.push_back((g2o::OptimizableGraph::Edge *) e);
	}
	g2o::VertexSE2 * p = dynamic_cast<g2o::VertexSE2 *>(v);
	if(p!=NULL){
	  // this is an SE2 vertex
	  poses.push_back(vw);
	}
	else{
	  // this is an XY point
	  landmarks.push_back(vw);
	}
      }
    }
    
    ((g2o::OptimizableGraph::Edge *)e)->setRobustKernel(robust_kernel);
    
    if(isEdgePose) edgesPoses.push_back((g2o::EdgeSE2 *)e);
    else edgesLandmarks.push_back((g2o::OptimizableGraph::Edge *)e);
  }
  
  std::cout << "poses.size() = " << poses.size() << std::endl;
  std::cout << "landmarks.size() = " << landmarks.size() << std::endl;
  std::cout << "edgesPoses.size() = " << edgesPoses.size() << std::endl;
  std::cout << "edgesLandmarks.size() = " << edgesLandmarks.size() << std::endl;
    
  // order poses by id (if necessary);
  for(int i=0; i<poses.size(); i++){
    for(int j=i-1; j>=0; j--){
      if(poses[i]->vertex->id() < poses[j]->vertex->id()){
  	VertexWrapper * tmp = poses[i];
  	poses[i] = poses[j];
  	poses[j] = tmp;
	i--;
      }
    }
  }
  
  int prev = -1;
  for(unsigned int i=0; i<poses.size(); i++){
    int id = poses[i]->vertex->id();
    if(id < prev){
      std::cout << id << " is lower than its previous: " << prev << std::endl;
    }
    prev = id;
    
  }
  
  
  // order edgesPoses such that edge i connects pose i to pose i+1
  for(unsigned int i=0; i<poses.size()-1; i++){
    VertexWrapper * from = poses[i];
    VertexWrapper * to = poses[i+1];
    
    // look for the corresponding edge;
    g2o::EdgeSE2 * edge = 0;
    int e_index;
    for(unsigned int j=0; j<edgesPoses.size(); j++){
      g2o::EdgeSE2 * e = edgesPoses[j];
      if(e->vertices()[0]->id() == from->vertex->id()  &&  e->vertices()[1]->id() == to->vertex->id()){
	edge = e;
	e_index = j;
	break;
      }
    }
    if(edge == 0){
      continue;
    }
    g2o::EdgeSE2 * tmp = edgesPoses[i];
    edgesPoses[i] = edge;
    edgesPoses[e_index] = tmp;
  }
  
  prev = -1;
  for(unsigned int i=0; i<edgesPoses.size(); i++){
    g2o::EdgeSE2 * e = edgesPoses[i];
    int id = e->vertices()[0]->id();
    if(id < prev){
      std::cout << "edge.vertices()[0]->id() = " << id << " is lower than its previous: " << prev << std::endl;
    }
    prev = id;
  }
  
  // optimizer->computeInitialGuess();
  
  // start to generate the stars
  std::vector<Star2D *> stars;
  Star2D * s;
  unsigned int inserted = 0;
  for(unsigned int i=0; i<poses.size(); i++){
    if(inserted==0){
      if(i>0){
	purgeLonelyLandmarks(s);
	s->gauge_index = (s->poses.size()/2);
      }
      s = new Star2D();
      stars.push_back(s);
    }
    
    if(i>0 && inserted==0){ // must add the previous pose/edge
      i--;
    }
    
    VertexWrapper * p = poses[i];
    s->poses.push_back(p);
    s->edgesPoses.push_back(edgesPoses[i]);
    
    inserted = (inserted+1) % _starLength;
    
    for(unsigned int e=0; e<p->edges.size(); e++){
      g2o::HyperGraph::Vertex * l = (g2o::HyperGraph::Vertex *) p->edges[e]->vertices()[1];
      VertexWrapper * vw = s->getWrapper(l);
      if(!vw){
	// look for the right wrapper
	vw = findWrapper(landmarks, l);
	if(!vw){
	  std::cerr << "NNNNNNOOOOOOOOOOOOOO!" << std::endl;
	  exit(12);
	}
	
	s->landmarks.push_back(vw);
	
      }
      s->edgesLandmarks.push_back(p->edges[e]);
    }
    
  }
  purgeLonelyLandmarks(stars[stars.size()-1]);
  
  std::cout << "generated " << stars.size() << " stars" << std::endl;
  
  //  return 0;
  
  g2o::EdgeLabeler labeler(optimizer);
  
  for(unsigned int star_index=0; star_index<stars.size(); star_index++){
    std::cout << std::endl;
    std::cout << "analyzing star #" << star_index+1 << "/" << stars.size() << std::endl;
    
    Star2D * s = stars[star_index];
    
    std::cout << "poses: " << s->poses.size() << std::endl;
    std::cout << "landmarks: " << s->landmarks.size() << std::endl;
    
    // push all the estimates
    s->pushState();
    // fix the gauge
    s->fixGauge();
    
    // ready to move stuff
    
    // move landmarks according to a guess from the edges
    for(unsigned int i=0;i<s->edgesLandmarks.size(); i++){
      g2o::EdgeSE2PointXYBearing * e = dynamic_cast<g2o::EdgeSE2PointXYBearing *>(s->edgesLandmarks[i]);
      if(e!=NULL){
	g2o::OptimizableGraph::VertexSet fixed;
	fixed.insert(e->vertices()[0]);
	e->initialEstimate(fixed, (g2o::OptimizableGraph::Vertex *)e->vertices()[1]);
      }
    }
    
    // optimize the local map
    std::cout << "optimizing the local map..." << std::endl;
    g2o::OptimizableGraph::EdgeSet toOptimize;
    for(unsigned int i=0; i<s->edgesPoses.size(); i++){
      if(star_index<stars.size()-1 || i<s->edgesPoses.size()-1)
      toOptimize.insert(s->edgesPoses[i]);
    }
    for(unsigned int i=0;i<s->edgesLandmarks.size(); i++){
      toOptimize.insert(s->edgesLandmarks[i]);
    }
    optimizer->initializeOptimization(toOptimize);
    optimizer->computeInitialGuess();
    int optim_result = optimizer->optimize(_optimizationSteps);
    
    if(optim_result < 1){
      s->popState();
    }
    
    if(optim_result > 0){
      std::cout << "looking for shared variables" << std::endl;
      std::vector<VertexWrapper *> shared;
      std::vector<VertexWrapper *> local; // non-shared variables
      getSharedEdges(s, shared, local);
    
      if(shared.size() > 0){
	std::cout << "clustering shared landmarks..." <<std::endl;
	int labels[shared.size()];
	int clusters = makeClusters(shared, labels);
	
	std::cout << "found " << clusters << " clusters" << std::endl;
	
	insertSharedEdges(s, shared, clusters, labels);
      }
      insertLocalEdges(s,local);
    }
    
    
    
    // create condensed measurements for the poses
    if(_createPosesEdges){
      for(unsigned int i=0; i<s->poses.size(); i++){
	
	if(i==s->gauge_index) continue;
	
	g2o::EdgeSE2 * edge = new g2o::EdgeSE2();
	
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	edge->vertices()[1] = s->poses[i]->vertex;
	
	
	s->edgesCondensed.insert(edge);
      }
    }
    
    std::cout << "labelling condensed edges" << std::endl;
    labeler.labelEdges(s->edgesCondensed);
    
    // if(optim_result < 1){ // add the original edges instead of the condensed ones
    //   for(unsigned int i=0; i<s->edgesLandmarks.size(); i++){
    // 	s->edgesCondensed.insert(s->edgesLandmarks[i]);
    //   }
    // }
    
    if(optim_result > 0){
      s->popState();
    }
    s->unfixGauge();
    
  }
  
  
  std::string fewname = "questo_no";
  std::string manyname = "questo_no";
  prepareOutputName(std::string(argv[1]), fewname, manyname);
  
  // std::cout << "generating file " << fewname << "...";
  
  // for(unsigned int i=0; i<stars.size(); i++){
  //   Star2D * s = stars[i];
    
  //   for(std::set<g2o::OptimizableGraph::Edge *>::iterator it=s->edgesShared.begin(); it!=s->edgesShared.end(); it++){
  //     (*it)->setLevel(FEW_LEVEL);
  //     optimizer->addEdge((*it));
  //   }
  // }
  
  // std::ofstream ofs1(fewname.c_str());
  // optimizer->save(ofs1, FEW_LEVEL);
  
  // std::cout << "\tDONE" << std::endl;
  
  std::cout << "generating file " << manyname << "...";
  
  for(unsigned int i=0; i<stars.size(); i++){
    Star2D * s = stars[i];
    
    for(std::set<g2o::OptimizableGraph::Edge *>::iterator it=s->edgesCondensed.begin(); it!=s->edgesCondensed.end(); it++){
      (*it)->setLevel(MANY_LEVEL);
      optimizer->addEdge((*it));
    }
  }
  
  std::ofstream ofs2(manyname.c_str());
  optimizer->save(ofs2, MANY_LEVEL);
  
  std::cout << "\tDONE" << std::endl;
  
  
  return 0;
}
