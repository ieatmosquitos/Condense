#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
//#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/slam3d/types_slam3d.cpp"
#include "g2o/types/slam2d/types_slam2d.cpp"
#include <vector>

// types definition
typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;


int main(int argc, char ** argv){
  // check input
  if(argc < 3){
    std::cout << "Usage: compare <file1.g2o> <file2.g2o>" << std::endl;
    std::cout << "vertices are taken from file1, edges from file2" << std::endl;
    return 0;
  }
  
  // allocate the optimizers
  g2o::SparseOptimizer * optimizer1 = new g2o::SparseOptimizer();	// vertices taken from here
  g2o::SparseOptimizer * optimizer2 = new g2o::SparseOptimizer();	// edges taken from here
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  
  optimizer1->setAlgorithm(solver);
  optimizer2->setAlgorithm(solver);
  //optimizer3->setAlgorithm(solver);
  
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
  std::cout << "graph #1:" << std::endl;
  std::cout << "\tLoaded " << optimizer1->vertices().size() << " vertices" << std::endl;
  std::cout << "\tLoaded " << optimizer1->edges().size() << " edges" << std::endl;
  
  std::cout << "graph #2:" << std::endl;
  std::cout << "\tLoaded " << optimizer2->vertices().size() << " vertices" << std::endl;
  std::cout << "\tLoaded " << optimizer2->edges().size() << " edges" << std::endl;
  
  // safety checks
  
  std::cout << "safety checks...";
  if(optimizer1->vertices().size() != optimizer2->vertices().size()){
    std::cout << "!!! the two graphs don't have the same number of vertices !!!" << std::endl;
    return 3;
  }
  std::cout << std::endl;
  
  
  g2o::HyperGraph::VertexIDMap map1 = optimizer1->vertices();
  g2o::HyperGraph::VertexIDMap map2 = optimizer2->vertices();
  
  double estimate[100];
  for(g2o::HyperGraph::VertexIDMap::iterator it=map1.begin(); it!=map1.end(); it++){
    g2o::OptimizableGraph::Vertex * v1 = (g2o::OptimizableGraph::Vertex *) it->second;
    // look for the same vertex in the other map
    g2o::OptimizableGraph::Vertex * v2 = (g2o::OptimizableGraph::Vertex *) map2[v1->id()];
    
    v1->getEstimateData(estimate);
    v2->setEstimateData(estimate);
    
  }
  
  
  optimizer2->initializeOptimization();
  optimizer2->computeActiveErrors();
  std::cout << "chi^2 = " << optimizer2->activeChi2() << std::endl;
  std::cout << std::endl;
  
}
