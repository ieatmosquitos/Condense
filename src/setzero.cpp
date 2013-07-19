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


std::string computeOutputName(std::string input_name){
  std::string name = input_name;
  if(name.substr(name.length()-4, name.length()).compare(".g2o") == 0 ){
    name = name.substr(0, name.length()-4); 
  }
  
  name.append("_tozero.g2o");
  
  return name;
}


int main(int argc, char ** argv){
  if(argc < 2){
    std::cout << "Usage: compare <file1.g2o>" << std::endl;
    std::cout << "all the vertices from the input graph are set to 0, and the output is written to a file named 'tozero.g2o'" << std::endl;
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

  g2o::HyperGraph::VertexIDMap vertices = optimizer->vertices();
  
  double estimate[7] = {0, 0, 0, 0, 0, 0, 1}; // this should be big enough for every kind of vertex
  
  for(g2o::HyperGraph::VertexIDMap::iterator it=vertices.begin(); it!=vertices.end(); it++){
    g2o::OptimizableGraph::Vertex * v = (g2o::OptimizableGraph::Vertex *) it->second;
    
    v->setEstimateData(estimate);
    
  }
  
  std::string outname = computeOutputName(std::string(argv[1]));
  
  std::ofstream ofs(outname.c_str());
  optimizer->save(ofs);
  
  std::cout << "output saved to '" << outname << "'" << std::endl;
}
