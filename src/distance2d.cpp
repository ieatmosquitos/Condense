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
pthread_mutex_t _mutex_error;
pthread_mutex_t _mutex_vertices;
pthread_mutex_t _mutex_cout;
double _total_error = 0;
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



double getError(g2o::HyperGraph::VertexIDMap &map1, g2o::HyperGraph::VertexIDMap &map2, g2o::OptimizableGraph::Vertex * v){
  
  // v is the fixed point wrt which all the distances are computed.
  // w is the point that is measured
  // the number (E.G. v1, v2) indicates the graph, so v1 and v2 are the same point v but in the first or second graph
  
  double v1_est[2];
  double v2_est[2];
  double w1_est[2];
  double w2_est[2];
  
  // understand if it is a pose/landmark (needed to call estimate())
  g2o::VertexSE2 * v1 = dynamic_cast<g2o::VertexSE2 *>(v);
  
  if(v1 == NULL){
    // it's a landmark
    g2o::VertexPointXY * v1 = dynamic_cast<g2o::VertexPointXY *>(v);
    if(v1==NULL){
      std::cerr << "Vertex " << v->id() << " isn't a VertexSE2 nor a VertexPointXY. Cannot handle it. I'm out!" << std::endl;
      exit(8);
    }
    g2o::VertexPointXY * v2 = (g2o::VertexPointXY *)map2[v1->id()];
    v1_est[0] = v1->estimate()[0];
    v1_est[1] = v1->estimate()[1];
    v2_est[0] = v2->estimate()[0];
    v2_est[1] = v2->estimate()[1];
  }
  else{
    g2o::VertexSE2 * v2 = (g2o::VertexSE2 *)map2[v1->id()];
    v1_est[0] = v1->estimate()[0];
    v1_est[1] = v1->estimate()[1];
    v2_est[0] = v2->estimate()[0];
    v2_est[1] = v2->estimate()[1];
  }
  
  double error = 0;
  
  for(g2o::HyperGraph::VertexIDMap::iterator it=map1.begin(); it!=map1.end(); it++){
    g2o::OptimizableGraph::Vertex * w1 = (g2o::OptimizableGraph::Vertex *) it->second;
        
    g2o::VertexSE2 * p1 = dynamic_cast<g2o::VertexSE2 *>(w1);
    
    if(p1!=NULL){ // it's a pose
      g2o::VertexSE2 * p2 = (g2o::VertexSE2 *) map2[p1->id()];
      
      w1_est[0] = p1->estimate()[0];
      w1_est[1] = p1->estimate()[1];
      w2_est[0] = p2->estimate()[0];
      w2_est[1] = p2->estimate()[1];
    }
    else{ // it is a landmark
      g2o::VertexPointXY * p1 = dynamic_cast<g2o::VertexPointXY *>(w1);
      
      if(p1==NULL){
	std::cerr << "ERROR: unknown tipe of vertex" << std::endl;
	exit(7);
      }
      
      g2o::VertexPointXY * p2 = (g2o::VertexPointXY *) map2[p1->id()];
      
      w1_est[0] = p1->estimate()[0];
      w1_est[1] = p1->estimate()[1];
      w2_est[0] = p2->estimate()[0];
      w2_est[1] = p2->estimate()[1];
    }
    
    double dist1 = sqrt(pow(w1_est[0] - v1_est[0], 2) + pow(w1_est[1] - v1_est[1], 2));
    double dist2 = sqrt(pow(w2_est[0] - v2_est[0], 2) + pow(w2_est[1] - v2_est[1], 2));
    error += abs(dist2 - dist1);
  }
  
  return error;
}


// void checkBestError(double error, int id){
//   pthread_mutex_lock(&_mutex_error);
//   if(error < _best_error || _best_error == -1){
//     _best_error = error;
//     pthread_mutex_lock(&_mutex_cout);
//     std::cout << std::endl << id << ": new record: " << _best_error << std::endl;
//     pthread_mutex_unlock(&_mutex_cout);
//   }
//   pthread_mutex_unlock(&_mutex_error);
// }


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
    
    g2o::OptimizableGraph::Vertex * v =  vertices[index];
    
    double error = getError(map1, map2, v);
    
    pthread_mutex_lock(&_mutex_error);
    _total_error += error;
    pthread_mutex_unlock(&_mutex_error);
    
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
  
  std::cout << "average error = " << _total_error/vertices.size() << std::endl;
}
