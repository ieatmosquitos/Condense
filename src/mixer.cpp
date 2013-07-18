#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
//#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/slam3d/types_slam3d.cpp"
#include "g2o/types/slam2d/types_slam2d.cpp"
#include <vector>
#include <stdlib.h>

// types definition
typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;


double chi2;
//double ideal_chi2;
double initial_chi2;
bool optim_fail;

int max_clusters = -1;
int landmarks_per_edge = -1;
bool make_clusters = true;
int star_length = -1;

int preoptim = -1;

bool diff_number = false;

void setParameters(char * filename){
  std::string s(filename);
  
  unsigned long sl = s.find("-sl");
  if(sl==std::string::npos) return;
  unsigned long mc = s.find("-mc");
  if(mc==std::string::npos) return;
  unsigned long ml = s.find("-ml");
  if(ml==std::string::npos) return;
  
  
  unsigned long nocl = s.find("-noclust");
  if(nocl!=std::string::npos){
    make_clusters = false;
  }
  else{
    nocl = s.find("_", ml+4);
    
    if(nocl==std::string::npos){
      make_clusters = false;
      nocl = s.find("-", ml+4);
    }
  }
  
  std::stringstream ss;
  unsigned long _i = s.find("-i_");
  if(_i!=std::string::npos){
    unsigned long i_ = s.find("_", _i+3);
    ss << s.substr(_i+3, i_ - (_i+3)) << ' ';
  }
  
  ss << s.substr(sl+4, mc- (sl+4)) << ' ' << s.substr(mc+4, ml-(mc+4)) << ' ' << s.substr(ml+4, nocl - (ml+4));
  
  if(_i!=std::string::npos){
    ss >> preoptim;
  }
  
  ss >> star_length >> max_clusters >> landmarks_per_edge;
  
  // std::cout << "star_length = " << star_length << std::endl;
  
  // std::cout << "max_clusters = " << max_clusters << std::endl;
  
  // std::cout << "landmarks_per_edge = " << landmarks_per_edge << std::endl;
  
  // if(!make_clusters){
  //   std::cout << "don't make clusters!" << std::endl;
  // }

  

}


int main(int argc, char ** argv){
  // check input
  if(argc < 4){
    std::cout << "Usage: analyze <file1.g2o> <file2.g2o> <outfile.g2o>" << std::endl;
    std::cout << "vertices are taken from file1, edges from file2" << std::endl;
    std::cout << "the mixed graph is stored in outfile" << std::endl;
    return 0;
  }
  
  setParameters(argv[1]);
  
  // allocate the optimizers
  g2o::SparseOptimizer * optimizer1 = new g2o::SparseOptimizer();	// vertices taken from here
  g2o::SparseOptimizer * optimizer2 = new g2o::SparseOptimizer();	// edges taken from here
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(blockSolver);
  
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
  
  std::ofstream ofs(argv[3]);
  
  // std::cout << "graph #1:" << std::endl;
  // std::cout << "\tLoaded " << optimizer1->vertices().size() << " vertices" << std::endl;
  // std::cout << "\tLoaded " << optimizer1->edges().size() << " edges" << std::endl;
  
  // std::cout << "graph #2:" << std::endl;
  // std::cout << "\tLoaded " << optimizer2->vertices().size() << " vertices" << std::endl;
  // std::cout << "\tLoaded " << optimizer2->edges().size() << " edges" << std::endl;
  
  // safety checks
  
  // std::cout << "safety checks...";
  if(optimizer1->vertices().size() != optimizer2->vertices().size()){
    std::cout << "!!! the two graphs don't have the same number of vertices !!!" << std::endl;
    
    diff_number = true;
  }
  // std::cout << std::endl;
  
  
  g2o::HyperGraph::VertexIDMap map1 = optimizer1->vertices();
  g2o::HyperGraph::VertexIDMap map2 = optimizer2->vertices();
  
  //  optimizer2->initializeOptimization();
  //  optimizer2->optimize(10);
  //  optimizer2->computeActiveErrors();
  //  ideal_chi2 = optimizer2->activeChi2();
  
  
  double estimate[100];
  for(g2o::HyperGraph::VertexIDMap::iterator it=map1.begin(); it!=map1.end(); it++){
    g2o::OptimizableGraph::Vertex * v1 = (g2o::OptimizableGraph::Vertex *) it->second;
    // look for the same vertex in the other map
    g2o::OptimizableGraph::Vertex * v2 = (g2o::OptimizableGraph::Vertex *) map2[v1->id()];
    
    v1->getEstimateData(estimate);
    v2->setEstimateData(estimate);
    
  }
  
  optimizer2->save(ofs);
  
  std::cout << "files mixed, output saved to " << argv[3] << std::endl;
  
  // if(optim_result){
  //   optimizer2->computeActiveErrors();
  //   ofs << "chi^2 after 10 LM optimization steps: " << optimizer2->activeChi2() << std::endl;
  // }
  // else{
  //   ofs << "optimization ERROR, cannot compute final chi^2" << std::endl;
  // }
  return 0;
}
