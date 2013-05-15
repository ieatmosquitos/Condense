#include <iostream>
#include "FileReader.cpp"

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/slam2d/types_slam2d.cpp"

int main(int argc, char**argv){
  if(argc < 2){
    std::cout << "Usage: listedges <filename.g2o>" << std::endl;
    exit(0);
  }
  
  
  // types definition
  typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
  typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

  
  g2o::SparseOptimizer * optimizer = new g2o::SparseOptimizer();
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  
  optimizer->setAlgorithm(solver);
  
  std::ifstream ifs(argv[1]);
  if(!ifs){
    std::cout << "cannot read file " << argv[1] << std::endl;
    exit(1);
  }
  
  
  
  if(!optimizer->load(ifs)){
    std::cout << "error loading graph from file " << argv[1] << std::endl;
    exit(1);
  };
  
  for (g2o::HyperGraph::EdgeSet::const_iterator it = optimizer->edges().begin(); it != optimizer->edges().end(); ++it) {
    g2o::OptimizableGraph::Edge* e = static_cast<g2o::OptimizableGraph::Edge*>(*it);
    
    std::cout << "vertices: " << e->vertices().size() << std::endl;
    std::cout << "ids: ";
    for(unsigned int i=0; i<e->vertices().size(); i++){
      std::cout << "\t" << e->vertices()[i]->id();
    }
    std::cout << std::endl;
    
    e->write(std::cout);
    std::cout << std::endl;
    std::cout << std::endl;
  }
  
  return 0;
}
