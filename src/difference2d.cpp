/* 
   DIFFERENCE2D
   This program takes two g2o graph files and, assuming that they contain the same vertices, looks for the alignment that maximizes the overlap between the two graphs, and computes the total error, that is the sum of all the distances between the two estimates of each vertex.
 */
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
//#include "g2o/core/robust_kernel_impl.h"
//#include "g2o/types/slam3d/types_slam3d.cpp"
#include "g2o/types/slam2d/types_slam2d.cpp"
#include <vector>

// types definition
typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

// global variables
double _weights[2] = {1, 0};
pthread_mutex_t _mutex_error;
pthread_mutex_t _mutex_vertices;
pthread_mutex_t _mutex_cout;
double _best_error = -1;
int _num_threads = 4;


struct thread_struct{
  g2o::HyperGraph::VertexIDMap *map1;
  g2o::HyperGraph::VertexIDMap *map2;
  std::vector<g2o::OptimizableGraph::Vertex *> *vertices;
  bool * visited;
};


struct visualizer_struct{
  bool * visited;
  unsigned int size;
};



// returns the absolute value of the given element
template<class T> T abs(T a){
  if(a<0) return -a;
  return a;
}


// computes the difference between two angles, managing the problems around PI
double computeAnglesDifference(double ang1, double ang2){
  if(abs(ang1-ang2) > M_PI){	// this passes through the discontinuous point π
    if(ang1 > 0){	// ang2 is negative
      return -(2*M_PI - (ang1-ang2));
    }
    else{	// ang1 negative, ang2 positive
      return (2*M_PI - (ang2-ang1));
    }
  }
  return (ang1-ang2);
}


void computeVector3(double * p1, double * p2, double * v){
  double dx = p2[0] - p1[0];
  double dy = p2[1] - p1[1];
  v[0] = p1[0];
  v[1] = p1[1];
  v[2] = atan2(dy, dx);
}


// computes the transformation matrix that, applied to the second vector, brings it to overlap with the first vector
// the vectors are considered to be a minimal representation of the postion-orientation of an object in the 2d plane
// E.G. [1 2 PI/4] stands for: posizion=[1, 2] -- orientation=PI/4
Eigen::Matrix3d getTransformation(double * p1, double * p2){
  double t[2] = {p1[0] - p2[0], p1[1]-p2[1]};
  double r = computeAnglesDifference(p1[2], p2[2]);
  
  double sinr = sin(r);
  double cosr = cos(r);
  
  Eigen::Matrix3d ret;
  ret <<
    cosr, -sinr, t[0],
    sinr, cosr , t[1],
    0   , 0    , 1   ;
  return ret;
}



void getError(g2o::HyperGraph::VertexIDMap &map1, g2o::HyperGraph::VertexIDMap &map2, g2o::SE2 transf, double * error){
  for(g2o::HyperGraph::VertexIDMap::iterator it=map1.begin(); it!=map1.end(); it++){
    g2o::OptimizableGraph::Vertex * v1 = (g2o::OptimizableGraph::Vertex *) it->second;
    
    g2o::VertexSE2 * p1 = dynamic_cast<g2o::VertexSE2 *>(v1);
    
    if(p1!=NULL){ // it's a pose
      g2o::VertexSE2 * p2 = (g2o::VertexSE2 *) map2[p1->id()];
      
      g2o::SE2 diff = transf*p2->estimate();
      error[0] += sqrt(pow(diff[0]-p1->estimate()[0], 2) + pow(diff[1]-p1->estimate()[1], 2));
      error[1] += abs(computeAnglesDifference(diff[2],p1->estimate()[2]));
    }
    else{ // it is a landmark
      g2o::VertexPointXY * p1 = dynamic_cast<g2o::VertexPointXY *>(v1);
      
      if(p1==NULL){
	std::cerr << "ERROR: unknown tipe of vertex" << std::endl;
	exit(7);
      }
      
      g2o::VertexPointXY * p2 = (g2o::VertexPointXY *) map2[p1->id()];
      
      Eigen::Vector2d diff = transf*p2->estimate();
      error[0] += sqrt(pow(diff[0]-p1->estimate()[0], 2) + pow(diff[1]-p1->estimate()[1], 2));
    }
    
  }
}


void checkBestError(double error, int id){
  pthread_mutex_lock(&_mutex_error);
  if(error < _best_error || _best_error == -1){
    _best_error = error;
    pthread_mutex_lock(&_mutex_cout);
    std::cout << std::endl << id << ": new record: " << _best_error << std::endl;
    pthread_mutex_unlock(&_mutex_cout);
  }
  pthread_mutex_unlock(&_mutex_error);
}


void analyze(g2o::HyperGraph::VertexIDMap &map1, g2o::HyperGraph::VertexIDMap &map2, std::vector<g2o::OptimizableGraph::Vertex *> &vertices, bool*visited){
  unsigned int index = 0;
  
  while(index < vertices.size()){
    pthread_mutex_lock(&_mutex_vertices);
    while(index<vertices.size() && visited[index]){
      index++;
    }
    if(index >= vertices.size()){
      pthread_mutex_unlock(&_mutex_vertices);
      continue;
    }
    
    visited[index] = true;
    pthread_mutex_unlock(&_mutex_vertices);
    
    g2o::OptimizableGraph::Vertex * v1 =  vertices[index];
    
    int dim = v1->minimalEstimateDimension();
    // if dim is 3 it is a pose, if dim is 2 it is a landmark
    
    // get the corresponding vertex in the other graph
    g2o::OptimizableGraph::Vertex * v2 = (g2o::OptimizableGraph::Vertex *) map2[v1->id()];
    
    
    if(dim == 3){
      // if it is a pose, we don't need anything else
      
      g2o::VertexSE2 * p1 = dynamic_cast<g2o::VertexSE2 *>(v1);
      if(p1==NULL){
	std::cerr << "Vertex " << v1->id() << " has dimension 3 but it is not a VertexSE2. Cannot handle it. I'm out!" << std::endl;
	exit(8);
      }
      
      g2o::VertexSE2 * p2 = (g2o::VertexSE2 *) v2;
      
      g2o::SE2 transf = (p2->estimate().inverse() * p1->estimate());
      
      double error[2] = {0, 0};
      getError(map1, map2, transf, error);
      
      double condensed_error = error[0] * _weights[0] + error[1] * _weights[1];
      checkBestError(condensed_error, p1->id());
      
    }
    else{
      // it is a landmark
      g2o::VertexPointXY * l1 = (g2o::VertexPointXY *) v1;
      if(l1==NULL){
	std::cerr << "Vertex " << v1->id() << " has dimension 2 but it is not a VertexPointXY. Cannot handle it. I'm out!" << std::endl;
	exit(8);
      }
      
      g2o::VertexPointXY * l2 = (g2o::VertexPointXY *) v2;
      
      // we need something else to allineate the two graphs:
      // we will take another vertex in the first graph and its correspondant in the other graph,
      // and use it to give an "orientation" to the landmark.
      // Since we don't know which should we choose, we have to try all of them.
      
      for(g2o::HyperGraph::VertexIDMap::iterator it=map1.begin(); it!=map1.end(); it++){
	g2o::OptimizableGraph::Vertex * w1 = (g2o::OptimizableGraph::Vertex *) it->second;
	if(w1->id() == l1->id()) continue;	// can't align with self
	
	// pthread_mutex_lock(&_mutex_cout);
	// std::cout << std::endl << "\taligning " << l1->id() << " with " << w1->id();
	// pthread_mutex_unlock(&_mutex_cout);
	
	g2o::OptimizableGraph::Vertex * w2 = (g2o::OptimizableGraph::Vertex *) map2[w1->id()];
	
	// we don't actually care whether w is 3 or 2 dims, but this information will be needed to instantiate the right pointer (there is no superclass with the estimate() function that does not require the dimension information)
	double est_w1[3];
	double est_w2[3];
	
	if(w1->minimalEstimateDimension() == 3){
	  g2o::VertexSE2 * e1 = dynamic_cast<g2o::VertexSE2 *>(w1);
	  
	  if(e1 == NULL){
	    std::cerr << "Vertex " << w1->id() << " has dimension 3 but it is not a VertexSE2. Cannot handle it. I'm out!" << std::endl;
	    exit(8);
	  }
	  
	  g2o::VertexSE2 * e2 = (g2o::VertexSE2 *)(w2);
	  
	  est_w1[0] = e1->estimate()[0];
	  est_w1[1] = e1->estimate()[1];
	  est_w2[0] = e2->estimate()[0];
	  est_w2[1] = e2->estimate()[1];
	}
	else{
	  g2o::VertexPointXY * e1 = dynamic_cast<g2o::VertexPointXY *>(w1);
	  
	  if(e1 == NULL){
	    std::cerr << "Vertex " << w1->id() << " has dimension 2 but it is not a VertexPointXY. Cannot handle it. I'm out!" << std::endl;
	    exit(8);
	  }
	  
	  g2o::VertexPointXY * e2 = (g2o::VertexPointXY *)(w2);
	  
	  est_w1[0] = e1->estimate()[0];
	  est_w1[1] = e1->estimate()[1];
	  est_w2[0] = e2->estimate()[0];
	  est_w2[1] = e2->estimate()[1];
	}
	
	
	double d1[2] = {est_w1[0] - l1->estimate()[0], est_w1[1] - l1->estimate()[1]};
	double d2[2] = {est_w2[0] - l2->estimate()[0], est_w2[1] - l2->estimate()[1]};
	
	if( (d1[0] == 0 && d1[1] == 0) || (d2[0] == 0 && d2[1] == 0)){ // cannot determine an alignment
	  continue;
	}
	
	double angle1 = atan2(d1[1], d1[0]);
	double angle2 = atan2(d2[1], d2[0]);
	
	g2o::SE2 p1(l1->estimate()[0], l1->estimate()[1], angle1);
	g2o::SE2 p2(l2->estimate()[0], l2->estimate()[1], angle2);
	
	g2o::SE2 transf = (p2.inverse() * p1);
	
	double error[2] = {0, 0};
	getError(map1, map2, transf, error);
	
	double condensed_error = error[0] * _weights[0] + error[1] * _weights[1];     
	checkBestError(condensed_error, l1->id());
      }
    }
  }
}


void * launch_analyze(void * t_struct){
  thread_struct * arg = (thread_struct *) t_struct;
  analyze(*(arg->map1), *(arg->map2), *(arg->vertices), arg->visited);
}


void visualize(bool * visited, unsigned int size){
  double advanced = 0;
  
  unsigned int index = 0;
  bool evolved = true;
  while(index < size){
    
    while(!visited[index]){
      evolved = true;
      sleep(2);
    }
    index++;
    if(evolved){
      std::cout << "progress: " << floor(((double)index/(double)size)*10000+0.5)/100 << " %" << std::endl;
    }
    evolved = false;
  }
}


void * launch_visualize(void * v_struct){
  visualizer_struct * arg = (visualizer_struct *) v_struct;
  visualize(arg->visited, arg->size);
}


void init(){
  if(pthread_mutex_init(&_mutex_cout, NULL) == -1){
    std::cout << "errors initializing mutexes" << std::endl;
    exit(9);
  }
  if(pthread_mutex_init(&_mutex_error, NULL) == -1){
    std::cout << "errors initializing mutexes" << std::endl;
    exit(9);
  }
  if(pthread_mutex_init(&_mutex_vertices, NULL) == -1){
    std::cout << "errors initializing mutexes" << std::endl;
    exit(9);
  }
}


int main(int argc, char ** argv){
  
  // check input
  if(argc < 3){
    std::cout << "Usage: difference <file1.g2o> <file2.g2o>" << std::endl;
    return 0;
  }
  
  init();
  
   // allocate the optimizers
  g2o::SparseOptimizer * optimizer1 = new g2o::SparseOptimizer();
  g2o::SparseOptimizer * optimizer2 = new g2o::SparseOptimizer();
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  
  optimizer1->setAlgorithm(solver);
  optimizer2->setAlgorithm(solver);
  
  // load the graphs
  std::ifstream ifs1(argv[1]);
  std::ifstream ifs2(argv[2]);
  if(!ifs1){
    std::cout << "Cannot open file " << argv[1] << std::endl;
    return 1;
  }
  if(!ifs2){
    std::cout << "Cannot open file " << argv[2] << std::endl;
    return 1;
  }
  if(!optimizer1->load(ifs1)){
    std::cout << "Error loading graph #1" << std::endl;
    return 2;
  }
  if(!optimizer2->load(ifs2)){
    std::cout << "Error loading graph #2" << std::endl;
    return 2;
  }
  std::cout << "graph #1:";
  std::cout << "\tLoaded " << optimizer1->vertices().size() << " vertices" << std::endl;
  std::cout << "\tLoaded " << optimizer1->edges().size() << " edges" << std::endl;
  
  std::cout << "graph #2:";
  std::cout << "\tLoaded " << optimizer1->vertices().size() << " vertices" << std::endl;
  std::cout << "\tLoaded " << optimizer1->edges().size() << " edges" << std::endl;
  
  
  
  g2o::HyperGraph::VertexIDMap map1 = optimizer1->vertices();
  g2o::HyperGraph::VertexIDMap map2 = optimizer2->vertices();
  
  // safety checks
  
  std::cout << "safety checks...";
  if(optimizer1->vertices().size() != optimizer2->vertices().size()){
    std::cout << "!!! the two graphs don't have the same number of vertices !!!" << std::endl;
    return 3;
  }
  
  for(g2o::HyperGraph::VertexIDMap::iterator it=map1.begin(); it!=map1.end(); it++){
    
    g2o::OptimizableGraph::Vertex * v1 = (g2o::OptimizableGraph::Vertex *) it->second;
    
    int dim = v1->minimalEstimateDimension();
    
    if(dim!=3 & dim!=2){
      std::cerr << "Vertex #" << v1->id() << " has strange dimension: " << dim << std::endl;
      return 4;
    }
    
    // get the corresponding vertex in the other graph
    g2o::OptimizableGraph::Vertex * v2 = (g2o::OptimizableGraph::Vertex *) map2[v1->id()];
    assert(v2->id() == v1->id());
    if(v2->minimalEstimateDimension() != dim){
      std::cerr << "ERROR: vertex " << v2->id() << " has different dimension in the two graphs" << std::endl;
      return 5;
    }
  }
  
  std::cout << "done" << std::endl;
  
  std::vector<g2o::OptimizableGraph::Vertex*> vertices;
  
  for(g2o::HyperGraph::VertexIDMap::iterator it=map1.begin(); it!=map1.end(); it++){
    vertices.push_back((g2o::OptimizableGraph::Vertex *)it->second);
  }
  
  std::cout << "vertices vector filled" << std::endl;
  
  bool visited[vertices.size()];
  for(unsigned int i=0; i<vertices.size(); i++){
    visited[i] = false;
  }
  
  thread_struct t_struct;
  t_struct.map1 = &map1;
  t_struct.map2 = &map2;
  t_struct.vertices = &vertices;
  t_struct.visited = visited;
  
  pthread_t threads[_num_threads + 1];
  
  for(unsigned int i=0; i<_num_threads; i++){
    pthread_create(threads+i, NULL, launch_analyze, &t_struct);
  }
  
  visualizer_struct v_struct;
  v_struct.visited = visited;
  v_struct.size = vertices.size();
  
  pthread_create(threads+_num_threads, NULL, launch_visualize, &v_struct);

  std::cout << "threads launched" << std::endl;
  
  for(unsigned int i=0; i<_num_threads+1; i++){
    pthread_join(threads[i], NULL);
  }
  
  std::cout << "minimum error = " << _best_error;
}
