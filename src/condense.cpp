#include <iostream>

#include "FileReader.cpp"
#include "common.h"
#include "edge_labeler.h"
#include "clustering.cpp"

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
//#include "g2o/core/robust_kernel_impl.h"
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

int _starLength = 100;
int _optimizationSteps = 100;

// studying impact of poses edges on optimizability
bool _createPosesEdges = true;

// clustering stuff
bool _clusterize = false; // if false, no clusters are made and ONLY BINARY EDGES ARE CREATED
int _max_clusters = 6;
int _max_landmarks_per_edge = 5; // set lesser than 1 to let the edge be as big as it wants

// clusterizes the edges set in the given star and creates edges according to these clusters
void computeSharedEdges(Star2D * s){
  // look for shared xy points
  std::vector<VertexWrapper *> shared;
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
  std::cout << "shared XY: " << shared.size() << std::endl;
  
  // now, shared contains all the vertexwrappers that wrap landmarks that are shared with other stars
  
  // let's clusterize these landmarks
  
  double values[2*shared.size()];
  double * means = (double *) malloc(sizeof(double));
    
  for(unsigned int i=0; i<shared.size(); i++){
    g2o::VertexPointXY * l = (g2o::VertexPointXY *) shared[i]->vertex;
    unsigned int index = i*2;
    Eigen::Vector2d est = l->estimate();
    values[index] = est[0];
    values[index + 1] = est[1];
  }
  
  if(shared.size() < 1) return;
  
  if(shared.size() < 2){
    g2o::EdgeSE2LotsOfXY * edge = new g2o::EdgeSE2LotsOfXY();
    edge->setSize(2);
    edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
    edge->vertices()[0] = shared[0]->vertex;
    s->edgesShared.insert(edge);
    return;
  }
  
  int labels[shared.size()];
  
  int clusters = findClusters(values, 2, labels, &means, shared.size(), _max_clusters);
  
  std::cout << "found " << clusters << " clusters" << std::endl;
  
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
	
	count += so_many;
      }
    }
  }
  
  
}


void computeBinarySharedEdges(Star2D * s){
  // look for shared xyz points
  std::vector<VertexWrapper *> shared;
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
  std::cout << "shared XY: " << shared.size() << std::endl;
  
  // now, shared contains all the vertexwrappers that wrap landmarks that are shared with other stars
  
  for(unsigned int i=0; i<shared.size(); i++){
    g2o::EdgeSE2LotsOfXY * e = new g2o::EdgeSE2LotsOfXY();
    e->setSize(2);
    e->vertices()[0] = s->poses[s->gauge_index]->vertex;
    e->vertices()[1] = shared[i]->vertex;
    s->edgesShared.insert(e);
  }
}



void computeCondensedEdges(Star2D * s){
  
  if(s->landmarks.size() < 1) return;
  
  if(s->landmarks.size() == 1){
    g2o::EdgeSE2LotsOfXY * edge = new g2o::EdgeSE2LotsOfXY();
    edge->setSize(2);
    edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
    edge->vertices()[1] = s->landmarks[0]->vertex;
    s->edgesCondensed.insert(edge);
    return;
  }
  
  
  double values[s->landmarks.size()*2];
  double * means = (double *) malloc(sizeof(double));
  // populate the values array
  for(unsigned int i=0; i<s->landmarks.size(); i++){
    g2o::VertexPointXY * l = (g2o::VertexPointXY *) s->landmarks[i]->vertex;
    unsigned int index = i*2;
    Eigen::Vector2d est = l->estimate();
    values[index] = est[0];
    values[index+1] = est[1];
  }
  int labels[s->landmarks.size()];
    
  int clusters = findClusters(values, 2, labels, &means, s->landmarks.size(), _max_clusters);
  std::cout << "found " << clusters << " clusters" << std::endl;
  
  for(unsigned int c=0; c<clusters; c++){
    std::vector<VertexWrapper *> to;
    for(unsigned int i=0; i<s->landmarks.size(); i++){
      if(labels[i] == c){
	to.push_back(s->landmarks[i]);
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
	s->edgesCondensed.insert(edge);
	
	count += so_many;
      }
    }
  }
}



void computeBinaryCondensedEdges(Star2D * s){
  // now, shared contains all the vertexwrappers that wrap landmarks that are shared with other stars
  
  for(unsigned int i=0; i<s->landmarks.size(); i++){
    g2o::EdgeSE2LotsOfXY * e = new g2o::EdgeSE2LotsOfXY();
    e->setSize(2);
    e->vertices()[0] = s->poses[s->gauge_index]->vertex;
    e->vertices()[1] = s->landmarks[i]->vertex;
    s->edgesCondensed.insert(e);
  }
}




void init(){
  for(unsigned int i=0; i<MAX_IDS; i++){
    ids[i] = 0;
  }
}

int main(int argc, char ** argv){
  
  std::cout << "----CONDENSE 2D----" << std::endl;
  
  // check args
  if(argc<2){
    std::cout << "Usage: condense <input_file.g2o>" << std::endl;
    return 0;
  }
  
  
  // allocate the optimizer
  g2o::SparseOptimizer * optimizer = new g2o::SparseOptimizer();
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  
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
  


  // start to generate the stars
  std::vector<Star2D *> stars;
  Star2D * s;
  unsigned int inserted = 0;
  for(unsigned int i=0; i<poses.size(); i++){
    if(inserted==0){
      s = new Star2D();
      s->gauge_index = (_starLength/2);
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
  std::cout << "generated " << stars.size() << " stars" << std::endl;
  
  g2o::EdgeLabeler labeler(optimizer);
  
  for(unsigned int star_index=0; star_index<stars.size(); star_index++){
    std::cout << std::endl;
    std::cout << "analyzing star #" << star_index+1 << "/" << stars.size() << std::endl;
    
    Star2D * s = stars[star_index];
    
    std::cout << "poses: " << s->poses.size() << std::endl;
    std::cout << "landmarks: " << s->landmarks.size() << std::endl;
    
    // push all the estimates
    for(unsigned int i=0; i<s->poses.size(); i++){
      s->poses[i]->vertex->push();
    }
    for(unsigned int i=0; i<s->landmarks.size(); i++){
      s->landmarks[i]->vertex->push();
    }
    
    // fix the gauge
    s->poses[s->gauge_index]->vertex->setFixed(true);
	 
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
    int optim_result = optimizer->optimize(_optimizationSteps);
    if(optim_result<1){
      std::cout << "!!! optimization failed !!!" <<std::endl;
      // don't label according to the optimized graph.
      
      // go back to the unoptimized state
      // pop all the estimates
      for(unsigned int i=0; i<s->poses.size(); i++){
  	s->poses[i]->vertex->pop();
      }
      for(unsigned int i=0; i<s->landmarks.size(); i++){
  	s->landmarks[i]->vertex->pop();
      }
      
      // unfix the gauge
      s->poses[s->gauge_index]->vertex->setFixed(false);
      
      continue;
      
    }
    
    // shared variables:
    // the first and the last pose are always shared (except for the first and the last star)
    if(_createPosesEdges){
      if(star_index>0){
	g2o::EdgeSE2 * edge = new g2o::EdgeSE2();
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	edge->vertices()[1] = s->poses[0]->vertex;
	s->edgesShared.insert(edge);
      }
      if(star_index<stars.size()-1){
	g2o::EdgeSE2 * edge = new g2o::EdgeSE2();
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	edge->vertices()[1] = s->poses[s->poses.size()-1]->vertex;
	s->edgesShared.insert(edge);
      }
    }

    if(_clusterize){
      std::cout << "clustering shared landmarks..." <<std::endl;
      computeSharedEdges(s);
    }
    
    else{
      std::cout << "computing binary edges to the shared landmarks" << std::endl;
      computeBinarySharedEdges(s);
    }
    
    std::cout << "labelling edges to shared variables" << std::endl;
    labeler.labelEdges(s->edgesShared);
    
    
    
    // create condensed measurements for the local variables
    if(_createPosesEdges){
      for(unsigned int i=0; i<s->poses.size(); i++){
      
	if(i==s->gauge_index) continue;
      
	g2o::EdgeSE2 * edge = new g2o::EdgeSE2();
	
	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
	edge->vertices()[1] = s->poses[i]->vertex;
	
      
	s->edgesCondensed.insert(edge);
      }
    }
    
    if(_clusterize){
      std::cout << "clustering all the landmarks..." <<std::endl;
      computeCondensedEdges(s);
    }
    else{
      std::cout << "computing binary edges to all the local landmarks" << std::endl;
      computeBinaryCondensedEdges(s);
    }
    
    std::cout << "labelling condensed edges to the local variables" << std::endl;
    labeler.labelEdges(s->edgesCondensed);
    
    // pop all the estimates
    for(unsigned int i=0; i<s->poses.size(); i++){
      s->poses[i]->vertex->pop();
    }
    for(unsigned int i=0; i<s->landmarks.size(); i++){
      s->landmarks[i]->vertex->pop();
    }
    
    // unfix the gauge
    s->poses[s->gauge_index]->vertex->setFixed(false);
    
  }
  
  //std::cout << "created " << num_edgesxy << " binary edges and " << num_triedges << " tri-edges" << std::endl;
  
  std::cout << "generating file few.g2o...";
  
  for(unsigned int i=0; i<stars.size(); i++){
    Star2D * s = stars[i];
    
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
    Star2D * s = stars[i];
    
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
