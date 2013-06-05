#include <iostream>

#include "FileReader.cpp"
#include "common3d.h"
#include "edge_labeler.h"
#include "clustering.cpp"

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
//#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/slam3d/types_slam3d.cpp"

#define MAX_IDS 100000
#define FEW_LEVEL 2
#define MANY_LEVEL 3

// types definition
typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

std::vector<VertexWrapper *> poses;
std::vector<VertexWrapper *> landmarks;
std::vector<g2o::EdgeSE3 *> edgesPoses;
std::vector<g2o::EdgeSE3PointXYZ *> edgesLandmarks;

void * ids[MAX_IDS];

int _starLength = 20;
int _optimizationSteps = 30;
bool _clusterize = true;


// clusterizes the edges set in the given star and creates edges according to these clusters
void computeSharedEdges(Star3D * s){
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
  std::cout << "shared XYZ: " << shared.size() << std::endl;
  
  // now, shared contains all the vertexwrappers that wrap landmarks that are shared with other stars
  
  // let's clusterize these landmarks
  
  double values[3*shared.size()];
  double * means = (double *) malloc(sizeof(double));
    
  for(unsigned int i=0; i<shared.size(); i++){
    g2o::VertexPointXYZ * l = (g2o::VertexPointXYZ *) shared[i]->vertex;
    unsigned int index = i*3;
    Eigen::Vector3d est = l->estimate();
    values[index] = est[0];
    values[index + 1] = est[1];
    values[index + 2] = est[2];
  }
  
  if(shared.size() < 2){
    // just create 1 edge
    // TODO!!!
    return;
  }
  
  int labels[shared.size()];
  
  int max_clusters = 8;
  int clusters = findClusters(values, 3, labels, &means, shared.size(), max_clusters);
  
  std::cout << "found " << clusters << " clusters" << std::endl;
  
  for(unsigned int c=0; c<clusters; c++){
    std::vector<VertexWrapper *> to;
    for(unsigned int i=0; i<shared.size(); i++){
      if(labels[i] == c){
	to.push_back(shared[i]);
      }
    }
    g2o::EdgeSE3LotsOfXYZ * edge = new g2o::EdgeSE3LotsOfXYZ();
    edge->setSize(1+to.size());
    edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
    for(unsigned int i=0; i<to.size(); i++){
      edge->vertices()[1+i] = to[i]->vertex;
    }
    s->edgesShared.insert(edge);
  }
  
  
}


void computeBinarySharedEdges(Star3D * s){
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
  std::cout << "shared XYZ: " << shared.size() << std::endl;
  
  // now, shared contains all the vertexwrappers that wrap landmarks that are shared with other stars
  
  for(unsigned int i=0; i<shared.size(); i++){
    g2o::EdgeSE3LotsOfXYZ * e = new g2o::EdgeSE3LotsOfXYZ();
    e->setSize(2);
    e->vertices()[0] = s->poses[s->gauge_index]->vertex;
    e->vertices()[1] = shared[i]->vertex;
    s->edgesShared.insert(e);
  }
}



void computeCondensedEdges(Star3D * s){
  if(s->landmarks.size() > 1){
    double values[s->landmarks.size()*3];
    double * means = (double *) malloc(sizeof(double));
    // populate the values array
    for(unsigned int i=0; i<s->landmarks.size(); i++){
      g2o::VertexPointXYZ * l = (g2o::VertexPointXYZ *) s->landmarks[i]->vertex;
      unsigned int index = i*3;
      Eigen::Vector3d est = l->estimate();
      values[index] = est[0];
      values[index+1] = est[1];
      values[index+2] = est[2];
    }
    int labels[s->landmarks.size()];
    
    int max_clusters = 8;
    int clusters = findClusters(values, 3, labels, &means, s->landmarks.size(), max_clusters);
    std::cout << "found " << clusters << " clusters" << std::endl;
    
    for(unsigned int c=0; c<clusters; c++){
      std::vector<VertexWrapper *> to;
      for(unsigned int i=0; i<s->landmarks.size(); i++){
	if(labels[i] == c){
	  to.push_back(s->landmarks[i]);
	}
      }
      g2o::EdgeSE3LotsOfXYZ * edge = new g2o::EdgeSE3LotsOfXYZ();
      edge->setSize(1+to.size());
      edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
      for(unsigned int i=0; i<to.size(); i++){
	edge->vertices()[1+i] = to[i]->vertex;
      }
      s->edgesCondensed.insert(edge);
    }
  }
  
}



void computeBinaryCondensedEdges(Star3D * s){
  // now, shared contains all the vertexwrappers that wrap landmarks that are shared with other stars
  
  for(unsigned int i=0; i<s->landmarks.size(); i++){
    g2o::EdgeSE3LotsOfXYZ * e = new g2o::EdgeSE3LotsOfXYZ();
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
  
  std::cout << "----CONDENSE 3D----" << std::endl;
  
  // check args
  if(argc<2){
    std::cout << "Usage: cond3nse <input_file.g2o>" << std::endl;
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
      g2o::EdgeSE3 * p = dynamic_cast<g2o::EdgeSE3 *>(e);
      if(p!=NULL){
	isEdgePose = true;
      }
      else{
	g2o::EdgeSE3PointXYZ * p = dynamic_cast<g2o::EdgeSE3PointXYZ *>(e);
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
	  prev->edges.push_back((g2o::EdgeSE3PointXYZ *)e);
	}
      }
      else{
      	VertexWrapper * vw = new VertexWrapper();
	vw->vertex = (g2o::OptimizableGraph::Vertex *) v;
      	ids[v->id()] = vw;
	if(!isEdgePose){
	  vw->edges.push_back((g2o::EdgeSE3PointXYZ *) e);
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
    
    if(isEdgePose) edgesPoses.push_back((g2o::EdgeSE3 *)e);
    else edgesLandmarks.push_back((g2o::EdgeSE3PointXYZ *)e);
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
    g2o::EdgeSE3 * edge = 0;
    int e_index;
    for(unsigned int j=0; j<edgesPoses.size(); j++){
      g2o::EdgeSE3 * e = edgesPoses[j];
      if(e->vertices()[0]->id() == from->vertex->id()  &&  e->vertices()[1]->id() == to->vertex->id()){
	edge = e;
	e_index = j;
	break;
      }
    }
    if(edge == 0){
      continue;
    }
    g2o::EdgeSE3 * tmp = edgesPoses[i];
    edgesPoses[i] = edge;
    edgesPoses[e_index] = tmp;
  }
  
  prev = -1;
  for(unsigned int i=0; i<edgesPoses.size(); i++){
    g2o::EdgeSE3 * e = edgesPoses[i];
    int id = e->vertices()[0]->id();
    if(id < prev){
      std::cout << "edge.vertices()[0]->id() = " << id << " is lower than its previous: " << prev << std::endl;
    }
    prev = id;
  }
  


  // start to generate the stars
  std::vector<Star3D *> stars;
  Star3D * s;
  unsigned int inserted = 0;
  for(unsigned int i=0; i<poses.size(); i++){
    if(inserted==0){
      s = new Star3D();
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
    
    Star3D * s = stars[star_index];
    
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
    if(star_index>0){
      g2o::EdgeSE3 * edge = new g2o::EdgeSE3();
      edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
      edge->vertices()[1] = s->poses[0]->vertex;
      //s->edgesShared.insert(edge);
    }
    if(star_index<stars.size()-1){
      g2o::EdgeSE3 * edge = new g2o::EdgeSE3();
      edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
      edge->vertices()[1] = s->poses[s->poses.size()-1]->vertex;
      //s->edgesShared.insert(edge);
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
    for(unsigned int i=0; i<s->poses.size(); i++){
      
      if(i==s->gauge_index) continue;
      
      g2o::EdgeSE3 * edge = new g2o::EdgeSE3();
      if(i<s->gauge_index){
    	edge->vertices()[0] = s->poses[i]->vertex;
    	edge->vertices()[1] = s->poses[s->gauge_index]->vertex;
      }
      else{
    	edge->vertices()[0] = s->poses[s->gauge_index]->vertex;
    	edge->vertices()[1] = s->poses[i]->vertex;
      }
      
      //s->edgesCondensed.insert(edge);
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
    Star3D * s = stars[i];
    
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
    Star3D * s = stars[i];
    
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
