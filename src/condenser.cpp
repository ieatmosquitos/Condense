#include <iostream>
#include <string>
#include <stdlib.h> // atoi, malloc, free

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
#include "g2o/core/robust_kernel_impl.h"
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

unsigned int _minimum_observations;

// studying impact of poses edges on optimizability
bool _createPosesEdges; // may be set to false using the -nopos option at launch time

// clustering stuff
bool _clusterize; // if false, no clusters are made and ONLY BINARY EDGES ARE CREATED
int _max_clusters;
int _max_landmarks_per_edge; // set lesser than 1 to let the edge be as big as it wants


// robust kernel
g2o::RobustKernel * robust_kernel;


g2o::OptimizationAlgorithmWithHessian* solverWithHessian;


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



void getSharedEdges(Star3D * s, std::vector<VertexWrapper *> &shared, std::vector<VertexWrapper *> &local){
  // look for shared xyz points
  for(unsigned int i=0; i<s->landmarks.size(); i++){
    VertexWrapper * l = s->landmarks[i];
    bool isShared = false;
    for(unsigned int e=0; e<l->edges.size(); e++){
      g2o::HyperGraph::Vertex * p =  l->edges[e]->vertices()[0];
      if(!s->contains(p)){
	isShared=true;
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
	s->edgesShared.insert(edge);
	s->edgesCondensed.insert(edge);
	
	count += so_many;
      }
    }
  }
}

void insertLocalEdges(Star3D * s,  std::vector<VertexWrapper *> &local){
  for(unsigned int i=0; i<local.size(); i++){
    g2o::EdgeSE3LotsOfXYZ * e = new g2o::EdgeSE3LotsOfXYZ();
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


int purgeLandmarks(Star3D * s,g2o::SparseOptimizer  * optimizer){
  unsigned int dropped = 0;
  
  // push the star state
  s->pushState();
  
  g2o::OptimizableGraph::EdgeSet toOptimize;
  for(unsigned int i=0; i<s->edgesLandmarks.size(); i++){
    toOptimize.insert(s->edgesLandmarks[i]);
  }
  for(unsigned int i=0; i<s->edgesPoses.size(); i++){
    toOptimize.insert(s->edgesPoses[i]);
  }
  
  optimizer->initializeOptimization(toOptimize); 
  optimizer->computeInitialGuess();
  
  toOptimize.clear();
  
  for(unsigned int l=0; l<s->landmarks.size(); l++){
    s->pushState();
    
    VertexWrapper * vw = s->landmarks[l];
    bool remove_it = false;
    
    std::vector<g2o::OptimizableGraph::Vertex *> toFix;
    
    // look for poses in the star that see this landmark
    for(unsigned int e=0; e<vw->edges.size(); e++){
      
      if(s->contains(vw->edges[e]->vertices()[0])){
	toOptimize.insert(vw->edges[e]);
	toFix.push_back((g2o::OptimizableGraph::Vertex *) (vw->edges[e]->vertices()[0]));
      }
    }
    
    // fix those poses
    for(unsigned int i=0; i<toFix.size(); i++){
      toFix[i]->setFixed(true);
    }
    
    // solve the landmark
    optimizer->initializeOptimization(toOptimize);
    optimizer->solver()->init();
    if (!solverWithHessian->buildLinearStructure()){
      std::cerr << "FATAL: failure while building linear structure" << std::endl;
      exit(25);
    }
    optimizer->computeActiveErrors();
    solverWithHessian->updateLinearSystem();
    
    vw->vertex->solveDirect();
    
    g2o::OptimizableGraph::Vertex * v = vw->vertex;
    
    Eigen::MatrixXd h(v->dimension(), v->dimension());
    for (int i=0; i<v->dimension(); i++){
      for (int j=0; j<v->dimension(); j++)
	h(i,j)=v->hessian(i,j);
    }
    Eigen::EigenSolver<Eigen::MatrixXd> esolver;
    esolver.compute(h);
    Eigen::VectorXcd ev= esolver.eigenvalues();
    double emin = std::numeric_limits<double>::max();
    double emax = -std::numeric_limits<double>::max();
    for (int i=0; i<ev.size(); i++){
      emin = ev(i).real()>emin ? emin : ev(i).real();
      emax = ev(i).real()<emax ? emax : ev(i).real();
    }
    double d=emin/emax;
    
    if(d<1e-3){
      //      std::cout << "d = " << d << std::endl;
      remove_it = true;
    }
    
    
    if(remove_it){
      dropped++;
      // should remove this landmark from the star
       
      // search the edges in the star that leads to this landmark
      for(unsigned int i=0; i<s->edgesLandmarks.size(); i++){
    	if(s->edgesLandmarks[i]->vertices()[1] == vw->vertex){
    	  s->edgesLandmarks.erase(s->edgesLandmarks.begin()+i);
    	  i--;
    	}
      }
  
      s->landmarks[l]->vertex->pop();
      s->landmarks[l]->vertex->pop();
      s->landmarks.erase(s->landmarks.begin()+l);
      l--;
       
    }
    
    for(unsigned int i=0; i<toFix.size(); i++){
      toFix[i]->setFixed(false);
    }
    
    s->popState();
  }
 
  s->popState();
  return dropped;
}


void init(int argc, char** argv){
  for(unsigned int i=0; i<MAX_IDS; i++){
    ids[i] = 0;
  }
  
  // initialize options
  _starLength = 20;
  _optimizationSteps = 100;
  _createPosesEdges = true;
  _clusterize = true;
  _max_clusters = 6;
  _max_landmarks_per_edge = 5;
  
  _minimum_observations = 1;
  
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
    
    else if(option.compare("-obs") == 0){
      known = true;
      i++;
      if(i == argc){
	std::cerr << "ERROR: no value specified for option " << option << std::endl;
	exit(1);
      }
      _minimum_observations = atoi(argv[i]);
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
  blockSolver->setOptimizer(optimizer);
  
  solverWithHessian = dynamic_cast<g2o::OptimizationAlgorithmWithHessian*>(optimizer->solver());
  
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
    
    ((g2o::OptimizableGraph::Edge *)e)->setRobustKernel(robust_kernel);
    
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
  { // this block requires much memory, that's why we isolate it using brackets
    Star3D * s;
    unsigned int inserted = 0;
    VertexWrapper * candidates [100000];
    int times [100000];
    g2o::OptimizableGraph::Edge ** cand_edges [100000];
    unsigned int cand_index_next = 0;
    unsigned int cand_edges_index_next [100000];
    for(unsigned int i=0; i<poses.size(); i++){
      if(inserted==0){
	unsigned int dropped = cand_index_next;
	unsigned int kept = 0;
	if(i>0){
	  // finalize the star checking which landmarks have to be removed
	  for(unsigned int c=0; c<cand_index_next; c++){
	    if(times[c] >= _minimum_observations){
	      kept ++;
	      // confirm the landmark and its edges
	      s->landmarks.push_back(candidates[c]);
	      for(unsigned int ed=0; ed<cand_edges_index_next[c]; ed++){
		s->edgesLandmarks.push_back(cand_edges[c][ed]);
	      }
	    }
	  }
	  dropped = dropped - kept;
	
	  // clear stuff
	  for(unsigned int ce=0; ce<cand_index_next; ce++){
	    free(cand_edges[ce]);
	  }
	
	  cand_index_next = 0;
	
	
	  dropped = dropped + purgeLandmarks(s,optimizer);
	  std::cout << "dropped: " << dropped << "\t; kept: " << kept << std::endl;
	  s->gauge_index = (s->poses.size()/2);
	}
            
	s = new Star3D();
	stars.push_back(s);
	std::cout << "generating star number " << stars.size() << std::endl;
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

	VertexWrapper * vw;
	int wrapper_index = findWrapperIndex(candidates, cand_index_next, l);
	if(wrapper_index != -1){
	  vw = candidates[wrapper_index];
	  times[wrapper_index] = times[wrapper_index] + 1;
	  cand_edges[wrapper_index][cand_edges_index_next[wrapper_index]]=(p->edges[e]);
	  cand_edges_index_next[wrapper_index] = cand_edges_index_next[wrapper_index] + 1 ;
	}
	else{
	  vw = findWrapper(landmarks, l);
	  if(!vw){
	    std::cerr << "NNNNNNOOOOOOOOOOOOOO!" << std::endl;
	    exit(15);
	  }
	  else{
	    times[cand_index_next] = 1;
	  
	    candidates[cand_index_next] = vw;
	  
	    g2o::OptimizableGraph::Edge ** edges_vect = (g2o::OptimizableGraph::Edge **) malloc(10000*sizeof(g2o::OptimizableGraph::Edge *));
	  
	    //	    std::cout << "adding cand_edges..." << std::endl;
	    cand_edges[cand_index_next] = edges_vect;
	    edges_vect[0] = p->edges[e];
	    cand_edges_index_next[cand_index_next] = 1;
	  
	    cand_index_next ++ ;
	  }
	}
      
      }
    
    }
    
    unsigned int dropped = cand_index_next;
    unsigned int kept = 0;
    // finalize the last star checking which landmarks have to be removed
    for(unsigned int c=0; c<cand_index_next; c++){
      if(times[c] >= _minimum_observations){
	kept ++;
	// confirm the landmark and its edges
	stars[stars.size()-1]->landmarks.push_back(candidates[c]);
	for(unsigned int ed=0; ed<cand_edges_index_next[c]; ed++){
	  stars[stars.size()-1]->edgesLandmarks.push_back(cand_edges[c][ed]);
	}
      }
    }
    dropped = dropped - kept;
	
    // clear stuff
    for(unsigned int ce=0; ce<cand_index_next; ce++){
      free(cand_edges[ce]);
    }
	
    cand_index_next = 0;
	
    dropped = dropped + purgeLandmarks(s,optimizer);
    std::cout << "dropped: " << dropped << "\t; kept: " << kept << std::endl;
    s->gauge_index = (s->poses.size()/2);

  }// end of the highly memory consuming block
  
  purgeLandmarks(stars[stars.size()-1], optimizer);
  
  
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
    
    std::vector<VertexWrapper *> shared;
    std::vector<VertexWrapper *> local;	// non-shared variables
    
    std::cout << "looking for shared variables" << std::endl;
    getSharedEdges(s, shared, local);  
    
    if(optim_result < 1){
      s->popState();
    }
    
    if(optim_result > 0){
      if(shared.size() > 0){
	std::cout << "clustering shared landmarks..." <<std::endl;
	int labels[shared.size()];
	int clusters = makeClusters(shared, labels);
	
	std::cout << "found " << clusters << " clusters" << std::endl;
	
	insertSharedEdges(s, shared, clusters, labels);
      }
    }
    
    insertLocalEdges(s,local);
    
    
    // create condensed measurements for the poses
    if(_createPosesEdges){
      for(unsigned int i=0; i<s->poses.size(); i++){
      
	if(i==s->gauge_index) continue;
      
	g2o::EdgeSE3 * edge = new g2o::EdgeSE3();
	
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
  //   Star3D * s = stars[i];
    
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
