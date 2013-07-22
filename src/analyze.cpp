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
  if(argc < 5){
    std::cout << "Usage: analyze <file1.g2o> <file2.g2o> <temp> <log.txt>" << std::endl;
    std::cout << "vertices are taken from file1, edges from file2" << std::endl;
    std::cout << "<temp> is a file used for intermediate processing" << std::endl;
    std::cout << "informations on the optimization effects are added to the log file" << std::endl;
    return 0;
  }
  
  setParameters(argv[1]);
  
  // allocate the optimizers
  g2o::SparseOptimizer * optimizer1 = new g2o::SparseOptimizer();	// vertices taken from here
  g2o::SparseOptimizer * optimizer2 = new g2o::SparseOptimizer();	// edges taken from here
  g2o::SparseOptimizer * optimizer3 = new g2o::SparseOptimizer();	// used to compute the error
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(blockSolver);
  //g2o::OptimizationAlgorithmGaussNewton * solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  
  optimizer1->setAlgorithm(solver);
  optimizer2->setAlgorithm(solver);
  optimizer3->setAlgorithm(solver);
  
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
  
  std::ofstream ofs;
  ofs.open(argv[4], std::ios_base::app);
  
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
  g2o::HyperGraph::VertexIDMap map3 = optimizer3->vertices();
  
  double estimate[100];
  int removedVert = 0;
  int removedEdg = 0;
  for(g2o::HyperGraph::EdgeSet::iterator it=optimizer2->edges().begin(); it!=optimizer2->edges().end(); it++){
    g2o::HyperGraph::Edge * e = (*it);
    
    g2o::OptimizableGraph::Edge * pose = dynamic_cast<g2o::EdgeSE2 *>(e);
    if(pose == NULL){
      pose = dynamic_cast<g2o::EdgeSE3 *>(e);
    }
    
    // if(pose != NULL){//it's an odometry edge, both vertices must be added (if not present) in the 3rd graph
    //   g2o::HyperGraph::Vertex * v1 = e->vertices()[0];
    //   g2o::HyperGraph::Vertex * v2 = e->vertices()[1];
      
    //   if(map3[v1->id()]==0){
    // 	optimizer3->addVertex(v1);
    //   }
    //   else{
    // 	std::cout << "it was already there!" << std::endl;
    //   }
    //   if(map3[v2->id()]==0){
    // 	optimizer3->addVertex(v2);
    //   }
    //   optimizer3->addEdge(e);
    // }
    
    if(pose != NULL){
      g2o::OptimizableGraph::Vertex * v1 = (g2o::OptimizableGraph::Vertex *) e->vertices()[0];
      g2o::OptimizableGraph::Vertex * v2 = (g2o::OptimizableGraph::Vertex *) e->vertices()[1];
      
      g2o::OptimizableGraph::Vertex * c1 = (g2o::OptimizableGraph::Vertex *) map1[v1->id()];	// corresponding vertex in the other graph
      if(c1 != 0){
	c1->getEstimateData(estimate);
	v1->setEstimateData(estimate);
      }
      else{
	optimizer2->removeVertex(v1);
	removedVert++;
      }
      
      g2o::OptimizableGraph::Vertex * c2 = (g2o::OptimizableGraph::Vertex *) map1[v2->id()];	// corresponding vertex in the other graph
      if(c2 != 0){
	c2->getEstimateData(estimate);
	v2->setEstimateData(estimate);
      }
      else{
	optimizer2->removeVertex(v2);
	removedVert++;
      }
      
      if(c1 == 0 || c2 == 0){
	optimizer2->removeEdge(e);
	removedEdg++;
      }
      
      continue;
    }
    
    g2o::HyperGraph::Vertex * v = e->vertices()[1];
    if(map1[v->id()] == 0){
      optimizer2->removeVertex(v);
      removedVert++;
      optimizer2->removeEdge(e);
      removedEdg++;
    }
    else{
      g2o::OptimizableGraph::Vertex * v = (g2o::OptimizableGraph::Vertex *) e->vertices()[1];
      g2o::OptimizableGraph::Vertex * c = (g2o::OptimizableGraph::Vertex *) map1[v->id()];
      
      c->getEstimateData(estimate);
      v->setEstimateData(estimate);
    }
  }
  
  std::cout << "removed " << removedVert << " vertices and " << removedEdg << " edges" << std::endl;
  
  //  optimizer2->initializeOptimization();
  //  optimizer2->optimize(10);
  //  optimizer2->computeActiveErrors();
  //  ideal_chi2 = optimizer2->activeChi2();
  
  std::ofstream outfs(argv[3]);
  optimizer2->save(outfs);
  outfs.close();
  
  optimizer2->clear();
  
  std::ifstream reloader(argv[3]);
  optimizer3->load(reloader);
  
  optimizer3->initializeOptimization();
  optimizer3->computeActiveErrors();
  
  // ofs << "estimates from: " << argv[1] << std::endl;
  // ofs << "edges from: " << argv[2] << std::endl;
  // ofs << "initial chi^2: " << optimizer2->activeChi2() << std::endl;
  
  initial_chi2 = optimizer3->activeChi2();
  
  int optim_result = optimizer3->optimize(10);
  
  if(!optim_result){
    optim_fail = true;
  }
  else{
    optim_fail = false;
  }
  
  optimizer3->computeActiveErrors();
  chi2 = optimizer3->activeChi2();
  
  ofs << argv[1] << "," << argv[2] << "," << optimizer1->vertices().size() << "," << optimizer2->vertices().size() << "," << star_length << "," << max_clusters << "," << landmarks_per_edge << "," << (make_clusters? " " : "NO CLUST") << "," << initial_chi2 << "," << preoptim << "," << chi2 << "," << (optim_fail? "fail" : "success");
  if(diff_number){
    ofs << "!!! DIFF_VERTICES_NUMBER !!!";
  }
  ofs << std::endl;

  
  // if(optim_result){
  //   optimizer2->computeActiveErrors();
  //   ofs << "chi^2 after 10 LM optimization steps: " << optimizer2->activeChi2() << std::endl;
  // }
  // else{
  //   ofs << "optimization ERROR, cannot compute final chi^2" << std::endl;
  // }
  return 0;
}
