#include <iostream>
#include "FileReader.cpp"
#include "common.h"

#include "g2o/types/slam2d/edge_se2.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/slam2d/edge_se2_twopointsxy.h"

#define MAX_IDS 100000

G2O_REGISTER_TYPE(VERTEX_SE2, Pose);
G2O_REGISTER_TYPE(VERTEX_XY, Landmark);

std::vector<Pose *> poses;
std::vector<Landmark *> landmarks;
std::vector<g2o::EdgeSE2 *> edgesPoses;
std::vector<g2o::EdgeSE2PointXY *> edgesLandmarks;

std::vector<g2o::EdgeSE2TwoPointsXY *> triEdges;

void * ids[MAX_IDS];

void readPose(std::vector<std::string> &textline){
  int id = atoi(textline[1].c_str());
  double x = atof(textline[2].c_str());
  double y = atof(textline[3].c_str());
  double theta = atof(textline[4].c_str());
  Pose * p = new Pose(x,y,theta);
  p->setId(id);
  poses.push_back(p);
  ids[id] = (void*) p;
}


void readLandmark(std::vector<std::string> &textline){
  int id = atoi(textline[1].c_str());
  double x = atof(textline[2].c_str());
  double y = atof(textline[3].c_str());
  Landmark * l = new Landmark(x,y);
  l->setId(id);
  landmarks.push_back(l);
  ids[id] = (void*) l;
}


void readEdgePoses(std::vector<std::string> &textline){
  g2o::EdgeSE2 * e = new g2o::EdgeSE2();
  
  int id1 = atoi(textline[1].c_str());
  int id2 = atoi(textline[2].c_str());
  
  e->vertices()[0] = (Pose *)ids[id1];
  e->vertices()[1] = (Pose *)ids[id2];
  
  double d_x = atof(textline[3].c_str());
  double d_y = atof(textline[4].c_str());
  double d_t = atof(textline[5].c_str());
  
  g2o::SE2 measurement(d_x, d_y, d_t);
  e->setMeasurement(measurement);
  
  Eigen::Matrix3d info;
  double i00 = atof(textline[6].c_str());
  double i01 = atof(textline[7].c_str());
  double i02 = atof(textline[8].c_str());
  double i11 = atof(textline[9].c_str());
  double i12 = atof(textline[10].c_str());
  double i22 = atof(textline[11].c_str());
  info <<	i00,	i01,	i02,
    		i01,	i11,	i12,
    		i02,	i12,	i22;
  
  e->setInformation(info);
  
  edgesPoses.push_back(e);
}


void readEdgeLandmark(std::vector<std::string> &textline){
  g2o::EdgeSE2PointXY * e = new g2o::EdgeSE2PointXY();
  
  int id1 = atoi(textline[1].c_str());
  int id2 = atoi(textline[2].c_str());
  
  Pose * p = (Pose *)ids[id1];
  Landmark * l = (Landmark *)ids[id2];
  
  p->addEdge(e);
  l->addEdge(e);
  
  e->vertices()[0] = p;
  e->vertices()[1] = l;
  
  double d_x = atof(textline[3].c_str());
  double d_y = atof(textline[4].c_str());
  
  Eigen::Vector2d measurement(d_x, d_y);
  e->setMeasurement(measurement);
  
  Eigen::Matrix2d info;
  double i00 = atof(textline[5].c_str());
  double i01 = atof(textline[6].c_str());
  double i11 = atof(textline[7].c_str());
  info <<	i00,	i01,
		i01,	i11;

  e->setInformation(info);
  
  edgesLandmarks.push_back(e);
}


void readFromFile(std::string filename){
  FileReader fr(filename);
  if(!fr.is_open()){
    std::cout << "cannot read file " << filename << ", quitting." << std::endl;
    exit(1);
  }
  
  std::vector<std::string> textline;
  fr.readLine(&textline);
  while(fr.good()){
    
    // compare returns 0 if the two strings are equal
    if((textline[0].compare(std::string("VERTEX_SE2"))) == 0){
      readPose(textline);
    }
    else if((textline[0].compare(std::string("VERTEX_XY"))) == 0){
      readLandmark(textline);
    }
    else if((textline[0].compare(std::string("EDGE_SE2"))) == 0){
      readEdgePoses(textline);
    }
    else if((textline[0].compare(std::string("EDGE_SE2_XY"))) == 0){
      readEdgeLandmark(textline);
    }
    
    textline.clear();
    fr.readLine(&textline);
  }
}


void init(){
  // for(unsigned int i=0; i<MAX_IDS; i++){
    
  // }
  
}


int main(int argc, char ** argv){
  std::cout << "----CONDENSE----" << std::endl;
  
  DEBUGMESS("<debug mode enabled>");
  
  if(argc < 3){
    std::cout << "Usage: condense <input_file.g2o> <output_file.g2o>" << std::endl;
    return 0;
  }
  
  init();
  
  readFromFile(argv[1]);
  
#ifdef DEBUGMODE
  // list vertices:
  std::cout << "listing poses:" << std::endl;
  for(unsigned int i=0; i<poses.size(); i++){
    std::cout << "Pose " << poses[i]->id() << " = " << poses[i]->estimate()[0] << "\t" << poses[i]->estimate()[1] << "\t" << poses[i]->estimate()[2] << std::endl;
  }
  
  std::cout << "listing landmarks:" << std::endl;
  for(unsigned int i=0; i<landmarks.size(); i++){
    std::cout << "Lmark " << landmarks[i]->id() << " = " << landmarks[i]->estimate()[0] << "\t" << landmarks[i]->estimate()[1] << std::endl;
  }
  
  std::cout << "listing pose-pose edges:" << std::endl;
  for(unsigned int i=0; i<edgesPoses.size(); i++){
    std::cout << "ID1 = " << edgesPoses[i]->vertices()[0]->id() << "\t ID2 = " << edgesPoses[i]->vertices()[1]->id() << std::endl;
  }
  
  std::cout << "listing pose-landmark edges:" << std::endl;
  for(unsigned int i=0; i<edgesLandmarks.size(); i++){
    std::cout << "ID1 = " << edgesLandmarks[i]->vertices()[0]->id() << "\t ID2 = " << edgesLandmarks[i]->vertices()[1]->id() << std::endl;
  }
#endif // DEBUGMODE
  
  std::cout << "Loaded " << poses.size() << " poses" << std::endl;
  std::cout << "Loaded " << landmarks.size() << " landmarks" << std::endl;
  std::cout << "Loaded " << edgesPoses.size() << " edges pose-pose" << std::endl;
  std::cout << "Loaded " << edgesLandmarks.size() << " edges pose-landmark" << std::endl;
  
  // find mergeable couples
  std::vector<Landmark *> first_vertex;
  std::vector<Landmark *> second_vertex;
  
  for (unsigned int i=0; i<landmarks.size(); i++){
    Landmark * l1 = landmarks[i];
    // store l1 associated poses
    int l1_poses[l1->_edges.size()];
    for (unsigned int obs=0; obs<l1->_edges.size(); obs++){
      Pose * p = (Pose *) (l1->_edges[obs]->vertices()[0]);
      l1_poses[obs] = p->id();
    }
    
    // check all the other landmarks
    for(unsigned int j = i+1; j<landmarks.size(); j++){
      Landmark * l2 = landmarks[j];
      // store l2 associated poses
      int l2_poses[l2->_edges.size()];
      for (unsigned int obs=0; obs<l2->_edges.size(); obs++){
	Pose * p = (Pose *) (l2->_edges[obs]->vertices()[0]);
	l2_poses[obs] = p->id();
      }
      
      // compare int sets
#ifdef DEBUGMODE
      std::cout << "comparing sets: " << std::endl;
      std::cout << "[";
      for(unsigned int obs=0; obs<l1->_edges.size(); obs++){
	std::cout << " " << l1_poses[obs];
      }
      std:: cout << " ]" << std::endl;
      std::cout << "[";
      for(unsigned int obs=0; obs<l2->_edges.size(); obs++){
	std::cout << " " << l2_poses[obs];
      }
      std:: cout << " ]" << std::endl;
#endif	// DEBUGMODE
      
      int * set1 = l1_poses;
      unsigned int set1_size = l1->_edges.size();
      int * set2 = l2_poses;
      unsigned int set2_size = l2->_edges.size();
      if(l1->_edges.size() > l2->_edges.size()){
	set1 = l2_poses;
	set2 = l1_poses;
	set1_size = l2->_edges.size();
	set2_size = l1->_edges.size();
      }
      
      bool found_all = true;
      for(unsigned int id1=0; id1<set1_size && found_all; id1++){
	bool found = false;
	for(unsigned int id2=0; id2<set2_size; id2++){
	  if (set1[id1] == set2[id2]){
	    found = true;
	    break;
	  }
	}
	if(!found){
	  found_all = false;
	  //break;
	}
      }
      if (found_all){
	first_vertex.push_back(l1);
	second_vertex.push_back(l2);
#ifdef DEBUGMODE
	std::cout << "OK!" << std::endl;
#endif	// DEBUGMODE
      }
    }
  }
  
#ifdef DEBUGMODE
  std::cout << "first_vertex.size() = " << first_vertex.size() << std::endl;
  std::cout << "second_vertex.size() = " << second_vertex.size() << std::endl;
#endif	// DEBUGMODE
  
#ifdef DEBUGMODE
  std::cout << std::endl << "mergeable vertex couples:" << std::endl;
  for(unsigned int fv = 0; fv<first_vertex.size(); fv++){
    std::cout << first_vertex[fv]->id() << "\t-->\t" << second_vertex[fv]->id() << std::endl;
    std::cout << first_vertex[fv]->id() << " seen by [";
    for(unsigned int i=0; i<first_vertex[fv]->_edges.size(); i++){
      std::cout << " " << first_vertex[fv]->_edges[i]->vertices()[0]->id();
    }
    std::cout << " ]" << std::endl;
    std::cout << second_vertex[fv]->id() << " seen by [";
    for(unsigned int i=0; i<second_vertex[fv]->_edges.size(); i++){
      std::cout << " " << second_vertex[fv]->_edges[i]->vertices()[0]->id();
    }
    std::cout << " ]" << std::endl;
    std::cout << std::endl;
    
  }
#endif	// DEBUGMODE
  
  // types definition
  typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
  typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
  
  
  std::cout << "look for the biggest shared set of landmarks" << std::endl;
  Pose * has_bigger_shared_set_1;
  Pose * has_bigger_shared_set_2;
  std::vector<Landmark *> bigger_shared_set;
  for(unsigned int i = 0; i<poses.size(); i++){
    Pose * p = poses[i];
    
#ifdef DEBUGMODE
    std::cout << "analyzing pose #" << i << " (id " << p->id() << ")" << std::endl;
#endif	// DEBUGMODE

    std::vector<Landmark *> p_sees;	// landmarks seen by p
    std::vector<g2o::EdgeSE2PointXY *> p_edges;	  // measurements for seen landmarks
    
    std::vector<Pose *> adjacent_poses;	// an adjacent pose is a pose that "sees" at least one of the landmarks seen by this one
    
    // populate the p_sees and adjacent_poses vectors
    for(unsigned int m = 0; m<p->_edges.size(); m++){
      Landmark * l = (Landmark *) (p->_edges[m]->vertices()[1]);
      p_sees.push_back(l);
      p_edges.push_back(p->_edges[m]);
      
      for(unsigned int from=0; from<l->_edges.size(); from++){
	Pose * adj = (Pose *) l->_edges[from]->vertices()[0];
	if(adj->id() == p->id()) continue;	// p is not adjacent to itself (else it gets the higher overlapping value)
	bool already_added = false;
	for(unsigned int j=0; j<adjacent_poses.size(); j++){
	  if(adjacent_poses[j] == adj){
	    already_added = true;
	    break;
	  }
	}
	if(!already_added) adjacent_poses.push_back(adj);
      }
    }
    
    DEBUGMESS("p_sees populated");
    
    // now for every adjacent pose must compute the "overlap" value. (number of seen vertices in common)
    for(unsigned int a=0; a<adjacent_poses.size(); a++){
#ifdef DEBUGMODE
      std::cout << "adj #" << a << std::endl;
#endif	// DEBUGMODE
      Pose * adj = adjacent_poses[a];
      
      std::vector<Landmark *> adj_sees;	// landmarks seen by adj
      
      // populate adj_sees
      for(unsigned int u=0; u<adj->_edges.size(); u++){
	adj_sees.push_back((Landmark*)(adj->_edges[u]->vertices()[1]));
      }
      
      std::vector<Landmark *> shared_vertices;
      bool consider_this[p_sees.size()];
      for(unsigned int c=0; c<p_sees.size(); c++){
	consider_this[c] = false;
      }
      
      // compare the two sets
      for(unsigned int s1=0; s1<p_sees.size(); s1++){
	for(unsigned int s2=0; s2<adj_sees.size(); s2++){
	  if(p_sees[s1] == adj_sees[s2]) {
	    shared_vertices.push_back((Landmark *) (p_sees[s1]));
	    consider_this[s1] = true;
	  }
	}
      }
      
#ifdef DEBUGMODE
      std::cout << "shared landmarks: " << shared_vertices.size() << std::endl;
#endif	// DEBUGMODE
      
      for(unsigned int i=0; i<p_sees.size(); i++){
#ifdef DEBUGMODE
	std::cout << "\ti: " << i << " ";
#endif	// DEBUGMODE
	if(!consider_this[i]) {
#ifdef DEBUGMODE
	  std::cout << "NOT considered" << std::endl;
#endif	// DEBUGMODE
	  continue;
	}
#ifdef DEBUGMODE
	std::cout << std::endl;
#endif	// DEBUGMODE
	
	Eigen::Vector2d meas1 = p_edges[i]->measurement();
	
	for(unsigned int j=i+1; j<p_sees.size(); j++){
	  if(!consider_this[j]) {
#ifdef DEBUGMODE
	    std::cout << "\t\tj: " << j << std::endl;
#endif	// DEBUGMODE
	    continue;
	  }
	  
#ifdef DEBUGMODE
	  std::cout << "\t\tj: " << "\tgenerate an edge " << std::endl;
#endif	// DEBUGMODE
	  // create a tri-edge
	  Eigen::Vector2d meas2 = p_edges[j]->measurement();
	  g2o::EdgeSE2TwoPointsXY * edge = new g2o::EdgeSE2TwoPointsXY();
	  edge->vertices()[0] = p;
	  edge->vertices()[1] = p_sees[i];
	  edge->vertices()[2] = p_sees[j];
	  Eigen::Vector4d meas(meas1[0], meas1[1], meas2[0], meas2[1]);
	  edge->setMeasurement(meas);
	  
	  Eigen::Matrix2d info1 = p_edges[i]->information();
	  Eigen::Matrix2d info2 = p_edges[j]->information();
	  Eigen::Matrix4d tri_info;
	  tri_info <<
	    info1(0,0),	info1(0,1),	0,		0,
	    info1(1,0),	info1(1,1),	0,		0,
	    0,		0,		info2(0,0),	info2(0,1),
	    0,		0,		info2(1,0),	info2(1,1);
	  
	  // tri_info <<
	  //   10, 2, 2, 2,
	  //   2, 10, 2, 2,
	  //   2, 2, 10, 2,
	  //   2, 2, 2, 10;
	  
	  edge->setInformation(tri_info);
	  triEdges.push_back(edge);
	}
      }      
      
      if(shared_vertices.size() > bigger_shared_set.size()){
	bigger_shared_set = shared_vertices;
	has_bigger_shared_set_1 = p;
	has_bigger_shared_set_2 = adj;
      }
    }
  }
  
  std::cout << "bigger shared landmarks set: " << bigger_shared_set.size() << " elements" << std::endl << "<";
  for(unsigned int i=0; i<bigger_shared_set.size(); i++){
    std::cout << " " << bigger_shared_set[i]->id();
  }
  std::cout << " >" << std::endl;
  std::cout << "seen by poses " << has_bigger_shared_set_1->id() << " and " << has_bigger_shared_set_2->id() << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "generated " << triEdges.size() << " tri-edges" << std::endl;
  
  std::cout << std::endl;
  
  // allocating the optimizer
  g2o::SparseOptimizer * optimizer = new g2o::SparseOptimizer();
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  
  optimizer->setAlgorithm(solver);
  
  std::cout << "Writing output file " << argv[2] << "...";
  std::ofstream g2oout(argv[2]);
  
  // populate graph:
  // poses
  for(unsigned int i=0; i<poses.size(); i++){
    // g2oout << "VERTEX_SE2 " << poses[i]->id() << " ";
    // poses[i]->write(g2oout);
    // g2oout << std::endl;
    optimizer->addVertex(poses[i]);
  }
  
  // edges between poses
  for(unsigned int i=0; i<edgesPoses.size(); i++){
    optimizer->addEdge(edgesPoses[i]);
    
    Pose * p1 = (Pose *) (edgesPoses[i]->vertices()[0]);
    Pose * p2 = (Pose *) (edgesPoses[i]->vertices()[1]);
    // g2oout << "EDGE_SE2" << " " << p1->id() << " " << p2->id() << " ";
    // edgesPoses[i]->write(g2oout);
    // g2oout << std::endl;
  }
  
  std::vector<g2o::EdgeSE2TwoPointsXY *> already_added;
  
  // insert tri-edges in the graph
  for(unsigned int i=0; i<triEdges.size(); i++){
    g2o::EdgeSE2TwoPointsXY * e = triEdges[i];
    
    bool skip_this = false;
    // check if it is already in the graph
    for(unsigned int j=0; j<already_added.size(); j++){
      g2o::EdgeSE2TwoPointsXY * compare = already_added[j];
      if(e->vertices()[0]->id() == compare->vertices()[0]->id() && e->vertices()[1]->id() == compare->vertices()[1]->id() && e->vertices()[2]->id() == compare->vertices()[2]->id()) skip_this = true;
      if(e->vertices()[0]->id() == compare->vertices()[0]->id() && e->vertices()[1]->id() == compare->vertices()[2]->id() && e->vertices()[2]->id() == compare->vertices()[1]->id()) skip_this = true;
    }
    if(skip_this) continue;
    
    Landmark * l1 = (Landmark *) (e->vertices()[1]);
    Landmark * l2 = (Landmark *) (e->vertices()[2]);
    
    // add the landmarks to the graph
    if(!l1->alreadyInserted){
      l1->alreadyInserted = true;
      optimizer->addVertex(l1);
      // g2oout << "VERTEX_XY " << l1->id() << " ";
      // l1->write(g2oout);
      // g2oout << std::endl;
    }
    if(!l2->alreadyInserted){
      l2->alreadyInserted = true;
      optimizer->addVertex(l2);
      // g2oout << "VERTEX_XY " << l2->id() << " ";
      // l2->write(g2oout);
      // g2oout << std::endl;
    }
    
    already_added.push_back(e);
    
    optimizer->addEdge(e);
    
    // Eigen::Vector4d meas = e->measurement();
    // Eigen::Matrix4d info = e->information();
    
    // Eigen::Vector2d meas1(meas[0], meas[1]);
    // Eigen::Vector2d meas2(meas[2], meas[3]);
    
    // Eigen::Matrix2d info1;
    // Eigen::Matrix2d info2;
    // info1 << info(0,0), info(0,1), info(1,0), info(1,1);
    // info2 << info(2,2), info(2,3), info(3,2), info(3,3);
  
    // g2o::EdgeSE2PointXY * e1 = new g2o::EdgeSE2PointXY();
    // g2o::EdgeSE2PointXY * e2 = new g2o::EdgeSE2PointXY();
    
    // e1->vertices()[0] = e->vertices()[0];
    // e2->vertices()[0] = e->vertices()[0];
    // e1->vertices()[1] = l1;
    // e2->vertices()[1] = l2;
    
    // e1->setMeasurement(meas1);
    // e2->setMeasurement(meas2);
    
    // e1->setInformation(info1);
    // e2->setInformation(info2);
    
    // optimizer->addEdge(e1);
    // optimizer->addEdge(e2);
    
    // g2oout << "EDGE_SE2_TWOPOINTSXY " << e->vertices()[0]->id() << " " << l1->id() << " " << l2->id() << " ";
    // e->write(g2oout);
    // g2oout << std::endl;
  }
  
  optimizer->save(g2oout);
  
  std::cout << "done" << std::endl;
  
  std::cout << already_added.size() << " different edges were actually created" << std::endl << std::endl;
  
  // try the optimizer
  optimizer->setVerbose(true);
  optimizer->vertex(0)->setFixed(true);
  optimizer->initializeOptimization();
  std::cout << "trying to optimize" << std::endl;
  optimizer->optimize(1);
  
  return 0;
}
