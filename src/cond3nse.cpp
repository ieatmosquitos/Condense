#include <iostream>
#include <string>
#include <stdlib.h> // atoi

#include "FileReader.cpp"
#include "common3d.h"
#include "edge_labeler.h"
#include "clustering.cpp"

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
//#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/slam3d/types_slam3d.cpp"

#define MAX_IDS 300000
#define FEW_LEVEL 2
#define MANY_LEVEL 3

// types definition
typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

std::vector<VertexWrapper *> poses;
std::vector<VertexWrapper *> landmarks;
std::vector<g2o::EdgeSE3 *> edgesPoses;
std::vector<g2o::OptimizableGraph::Edge *> edgesLandmarks;

void * ids[MAX_IDS];

int _starLength;
int _optimizationSteps;

// studying impact of poses edges on optimizability
bool _createPosesEdges; // may be set to false using the -nopos option at launch time

// clustering stuff
bool _clusterize; // if false, no clusters are made and ONLY BINARY EDGES ARE CREATED
int _max_clusters;
int _max_landmarks_per_edge; // set lesser than 1 to let the edge be as big as it wants



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
  
  double values[3*vertices.size()];
  double * means = (double *) malloc(sizeof(double));
  
  for(unsigned int i=0; i<vertices.size(); i++){
    g2o::VertexPointXYZ * l = (g2o::VertexPointXYZ *) vertices[i]->vertex;
    unsigned int index = i*3;
    Eigen::Vector3d est = l->estimate();
    values[index] = est[0];
    values[index + 1] = est[1];
    values[index + 2] = est[2];
  }
  
  int clusters = findClusters(values, 3, labels, &means, vertices.size(), _max_clusters);
  return clusters;
}



void getSharedEdges(Star3D * s, std::vector<VertexWrapper *> &shared){
  // look for shared xyz points
  for(unsigned int i=0; i<s->landmarks.size(); i++){
    VertexWrapper * l = s->landmarks[i];
    for(unsigned int e=0; e<l->edges.size(); e++){
      g2o::HyperGraph::Vertex * p =  l->edges[e]->vertices()[0];
      if(!s->contains(p)){
	shared.push_back(l);
	break;
      }
    }
  }
  std::cout << "shared XYZ: " << shared.size() << std::endl;
}



// 
void insertSharedEdges(Star3D * s, std::vector<VertexWrapper *> &shared, int clusters, int * labels){
  for(unsigned int c=0; c<clusters; c++){
    std::vector<VertexWrapper *> to;
    for(unsigned int i=0; i<shared.size(); i++){
      if(labels[i] == c){
	to.push_back(shared[i]);
      }
    }
    if(_max_landmarks_per_edge < 1 || _max_landmarks_per_edge > to.size()){
      // create an edge connecting the gauge to all the landmarks in the cluster
      g2o::EdgeSE3LotsOfXYZ * edge = new g2o::EdgeSE3LotsOfXYZ();
      edge->setSize(1+to.size());
      edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
      for(unsigned int i=0; i<to.size(); i++){
	edge->vertices()[1+i] = to[i]->vertex;
      }
      s->edgesShared.insert(edge);
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
	
	g2o::EdgeSE3LotsOfXYZ * edge = new g2o::EdgeSE3LotsOfXYZ();
	edge->setSize(1 + so_many);
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	for(unsigned int v=0; v<so_many; v++){
	  edge->vertices()[1+v] = to[count+v]->vertex;
	}
	s->edgesShared.insert(edge);
	
	count += so_many;
      }
    }
  }
}



void insertCondensedEdges(Star3D * s, int clusters, int * labels){
  for(unsigned int c=0; c<clusters; c++){
    std::vector<VertexWrapper *> to;
    for(unsigned int i=0; i<s->landmarks.size(); i++){
      if(labels[i] == c){
	to.push_back(s->landmarks[i]);
      }
    }
    if(_max_landmarks_per_edge < 1 || _max_landmarks_per_edge > to.size()){
      // create an edge connecting the gauge to all the landmarks in the cluster
      g2o::EdgeSE3LotsOfXYZ * edge = new g2o::EdgeSE3LotsOfXYZ();
      edge->setSize(1+to.size());
      edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
      for(unsigned int i=0; i<to.size(); i++){
	edge->vertices()[1+i] = to[i]->vertex;
      }
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
	
	g2o::EdgeSE3LotsOfXYZ * edge = new g2o::EdgeSE3LotsOfXYZ();
	edge->setSize(1 + so_many);
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	for(unsigned int v=0; v<so_many; v++){
	  edge->vertices()[1+v] = to[count+v]->vertex;
	}
	s->edgesCondensed.insert(edge);
	
	count += so_many;
      }
    }
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


bool checkIntegrity(Star3D * s){
  for(unsigned int i=0; i<s->landmarks.size(); i++){
    VertexWrapper * vw = s->landmarks[i];
    
    unsigned int count = 0;
    
    for(unsigned int j=0; j<vw->edges.size(); j++){
      if(s->contains(vw->edges[j]->vertices()[0])){
	count++;
      }
    }
    
    if(count<1) return false;
  }
  return true;
}


void orderLastElement(std::vector<g2o::EdgeSE3 *> &poses){
  int index = poses.size()-1;
  while(index > 0 && poses[index]->vertices()[0]->id() < poses[index-1]->vertices()[0]->id()){
    g2o::EdgeSE3 * tmp = poses[index];
    poses[index] = poses[index-1];
    poses[index-1] = tmp;
    index--;
  }
}


void init(int argc, char** argv){
  for(unsigned int i=0; i<MAX_IDS; i++){
    ids[i] = 0;
  }
  
  // initialize options
  _starLength = 20;
  _optimizationSteps = 50;
  _createPosesEdges = true;
  _clusterize = true;
  _max_clusters = 6;
  _max_landmarks_per_edge = 5;
  
  // check specified options
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
}

int main(int argc, char ** argv){
  
  std::cout << "----CONDENSE 3D----" << std::endl;
  
  // check args
  if(argc<2){
    std::cout << "Usage: cond3nse <input_file.g2o> <options>" << std::endl;
    std::cout << "\toptions:" << std::endl;
    std::cout << "\t-os N\tsets the optimization steps for the local maps to N" << std::endl;
    std::cout << "\t-sl N\tsets the length of the stars to N steps" << std::endl;
    std::cout << "\t-nopos\tinhibits the creation of pose-pose edges" << std::endl;
    std::cout << "\t-noclust\t doesn't clusterize the landmarks" << std::endl;
    std::cout << "\t-mc N\tsets the maximum number of clusters to N" << std::endl;
    std::cout << "\t-ml N\t sets the maximum number of landmarks per edge to N" << std::endl;
    
    return 0;
  }
  
  init(argc, argv);
  
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
  g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  //g2o::OptimizationAlgorithmLevenberg * solver = new g2o::OptimizationAlgorithmLevenberg(blockSolver);
  
  optimizer->setAlgorithm(solver);
  
  
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
      g2o::EdgeSE3 * p = dynamic_cast<g2o::EdgeSE3 *>(e);
      if(p!=NULL){
	isEdgePose = true;
      }
      else{
	g2o::OptimizableGraph::Edge * p = dynamic_cast<g2o::EdgeSE3PointXYZ *>(e);
	if(p==NULL){
	  p = dynamic_cast<g2o::EdgeSE3PointXYZDisparity *>(e);
	  
	  if(p==NULL){
	    p = dynamic_cast<g2o::EdgeSE3PointXYZDepth *>(e);
	    if(p==NULL){
	    
	      // don't load this edge
	      continue;
	    }
	  }
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
	g2o::VertexSE3 * p = dynamic_cast<g2o::VertexSE3 *>(v);
	if(p!=NULL){
	  // this is an SE3 vertex
	  poses.push_back(vw);
	}
	else{
	  // this is an XYZ point
	  landmarks.push_back(vw);
	}
      }  
    }
    
    if(isEdgePose){
      edgesPoses.push_back((g2o::EdgeSE3 *)e);
      orderLastElement(edgesPoses);
    }
    else edgesLandmarks.push_back((g2o::OptimizableGraph::Edge *)e);
  }
  
  std::cout << "poses.size() = " << poses.size() << std::endl;
  std::cout << "landmarks.size() = " << landmarks.size() << std::endl;
  std::cout << "edgesPoses.size() = " << edgesPoses.size() << std::endl;
  std::cout << "edgesLandmarks.size() = " << edgesLandmarks.size() << std::endl;
    
  // order poses by id (if necessary);
  std::cout << "ordering poses..." << std::endl;
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
  
  std::cout << "checking poses order..." << std::endl;
  int prev = -1;
  for(unsigned int i=0; i<poses.size(); i++){
    int id = poses[i]->vertex->id();
    if(id < prev){
      std::cout << id << " is lower than its previous: " << prev << std::endl;
    }
    prev = id;
    
  }
  
  
  // // order edgesPoses such that edge i connects pose i to pose i+1
  // std::cout << "checking edges order..." << std::endl;
  // for(unsigned int i=0; i<poses.size()-1; i++){
  //   VertexWrapper * from = poses[i];
  //   VertexWrapper * to = poses[i+1];
    
  //   // look for the corresponding edge;
  //   g2o::EdgeSE3 * edge = 0;
  //   int e_index;
  //   for(unsigned int j=0; j<edgesPoses.size(); j++){
  //     g2o::EdgeSE3 * e = edgesPoses[j];
  //     if(e->vertices()[0]->id() == from->vertex->id()  &&  e->vertices()[1]->id() == to->vertex->id()){
  // 	edge = e;
  // 	e_index = j;
  // 	break;
  //     }
  //   }
  //   if(edge == 0){
  //     continue;
  //   }
  //   g2o::EdgeSE3 * tmp = edgesPoses[i];
  //   edgesPoses[i] = edge;
  //   edgesPoses[e_index] = tmp;
  // }
  
  
  std::cout << "checking edges order..." << std::endl;
  prev = -1;
  for(unsigned int i=0; i<edgesPoses.size(); i++){
    g2o::EdgeSE3 * e = edgesPoses[i];
    int id = e->vertices()[0]->id();
    if(id < prev){
      std::cout << "edge.vertices()[0]->id() = " << id << " is lower than its previous: " << prev << std::endl;
    }
    prev = id;
  }
  
  
  std::cout << "starting to generate the stars..." << std::endl;

  // start to generate the stars
  std::vector<Star3D *> stars;
  Star3D * s;
  unsigned int inserted = 0;
  for(unsigned int i=0; i<poses.size(); i++){
    if(inserted==0){
      if(i>0){
	// if (!checkIntegrity(s)){
	//   std::cerr << "THESE STARS MAKE ME MAD!" << std::endl;
	//   exit(13);
	// }
	s->gauge_index = (s->poses.size()/2);
      }
      s = new Star3D();
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
	}
	else{
	  s->landmarks.push_back(vw);
	}
      }
      s->edgesLandmarks.push_back(p->edges[e]);
    }
    
  }
  // if (!checkIntegrity(stars[stars.size()-1])){
  //    std::cerr << "THESE STARS MAKE ME MAD!" << std::endl;
  //    exit(13);
  // };
  std::cout << "generated " << stars.size() << " stars" << std::endl;
  
  g2o::EdgeLabeler labeler(optimizer);
  
  for(unsigned int star_index=0; star_index<stars.size(); star_index++){
    std::cout << std::endl;
    std::cout << "analyzing star #" << star_index+1 << "/" << stars.size() << std::endl;
    
    Star3D * s = stars[star_index];
    
    std::cout << "poses: " << s->poses.size() << std::endl;
    std::cout << "landmarks: " << s->landmarks.size() << std::endl;
    
    // push all the estimates
    s->pushState();
    // fix the gauge
    s->fixGauge();
	 
    // ready to move stuff
    
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
    // if(optim_result<1){
    //   std::cout << "!!! optimization failed !!!" <<std::endl;
    //   // don't label according to the optimized graph.
      
    //   // go back to the unoptimized state
    //   // pop all the estimates
    //   for(unsigned int i=0; i<s->poses.size(); i++){
    // 	s->poses[i]->vertex->pop();
    //   }
    //   for(unsigned int i=0; i<s->landmarks.size(); i++){
    // 	s->landmarks[i]->vertex->pop();
    //   }
      
    //   // unfix the gauge
    //   s->poses[s->gauge_index]->vertex->setFixed(false);
      
    //   continue;
      
    // }
    
    
    // shared variables:
    // the first and the last pose are always shared (except for the first and the last star)
    if(_createPosesEdges){
      if(star_index>0){
	g2o::EdgeSE3 * edge = new g2o::EdgeSE3();
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	edge->vertices()[1] = s->poses[0]->vertex;
	s->edgesShared.insert(edge);
      }
      if(star_index<stars.size()-1){
	g2o::EdgeSE3 * edge = new g2o::EdgeSE3();
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	edge->vertices()[1] = s->poses[s->poses.size()-1]->vertex;
	s->edgesShared.insert(edge);
      }
    }
    
    if(optim_result > 0){
      std::cout << "looking for shared variables" << std::endl;
      std::vector<VertexWrapper *> shared;
      getSharedEdges(s, shared);
      
      if(shared.size() > 0){
	std::cout << "clustering shared landmarks..." <<std::endl;
	int labels[shared.size()];
	int clusters = makeClusters(shared, labels);
	
	std::cout << "found " << clusters << " clusters" << std::endl;
	
	insertSharedEdges(s, shared, clusters, labels);
      }
    }
    
    std::cout << "labelling edges to shared variables" << std::endl;
    labeler.labelEdges(s->edgesShared);
    
    
    
    // create condensed measurements for the local variables
    if(_createPosesEdges){
      for(unsigned int i=0; i<s->poses.size(); i++){
      
	if(i==s->gauge_index) continue;
      
	g2o::EdgeSE3 * edge = new g2o::EdgeSE3();
	
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	edge->vertices()[1] = s->poses[i]->vertex;
	
      
	s->edgesCondensed.insert(edge);
      }
    }
    
    if(optim_result>0){
      if(s->landmarks.size() > 0){
	std::cout << "clustering all the landmarks..." <<std::endl;
	int labels[s->landmarks.size()];
	int clusters = makeClusters(s->landmarks, labels);
	
	std::cout << "found " << clusters << " clusters" << std::endl;
	
	insertCondensedEdges(s, clusters, labels);
      }
    }
    
    std::cout << "labelling condensed edges to all the local variables" << std::endl;
    labeler.labelEdges(s->edgesCondensed);
    
    s->popState();
    s->unfixGauge();
    
  }
  
  
   std::string fewname = "questo_no";
  std::string manyname = "questo_no";
  prepareOutputName(std::string(argv[1]), fewname, manyname);
  
  std::cout << "generating file " << fewname << "...";
  
  for(unsigned int i=0; i<stars.size(); i++){
    Star3D * s = stars[i];
    
    for(std::set<g2o::OptimizableGraph::Edge *>::iterator it=s->edgesShared.begin(); it!=s->edgesShared.end(); it++){
      (*it)->setLevel(FEW_LEVEL);
      optimizer->addEdge((*it));
    }
  }
  
  std::ofstream ofs1(fewname.c_str());
  optimizer->save(ofs1, FEW_LEVEL);
  
  std::cout << "\tDONE" << std::endl;
  
  std::cout << "generating file " << manyname << "...";
  
  for(unsigned int i=0; i<stars.size(); i++){
    Star3D * s = stars[i];
    
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
