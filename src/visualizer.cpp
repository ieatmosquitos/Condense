#include <iostream>
#include <vector>
// #include "tools.cpp"
// #include "Retta2d.cpp"
// #include <list>
#include <iostream>
#include <cstdlib>
// #include <opencv2/opencv.hpp>
// #include <opencv2/highgui/highgui.hpp>
#include <drawer.h>
#include "common.h"
#include "FileReader.cpp"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/types/slam2d/vertex_se2.h"
#include "g2o/types/slam2d/edge_se2.h"
#include "g2o/types/slam2d/vertex_point_xy.h"
#include "g2o/types/slam2d/edge_se2_pointxy.h"
// #include "g2o/types/slam2d/types_slam2d.h"

#define DRAWER_MULTIPLIER 10
#define MAX_STEPS 7000

rdrawer::RobotDrawer * _drawer;	// used for drawing on the screen
g2o::SparseOptimizer * optimizer;	// optimizer declaration

bool free_running;
bool next_step;
bool draw;

// types definition
typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;


#define MAX_IDS 100000

void * ids[MAX_IDS];


  
G2O_REGISTER_TYPE(VERTEX_SE2, Pose);
G2O_REGISTER_TYPE(VERTEX_XY, Landmark);

std::vector<Pose *> poses;
std::vector<Landmark *> landmarks;
std::vector<g2o::EdgeSE2 *> edgesPoses;
std::vector<g2o::EdgeSE2PointXY *> edgesLandmarks;

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



void handleEvents(sf::RenderWindow * window){
  sf::Event event;
  if(window->GetEvent(event)){
    switch(event.Type){
    case sf::Event::KeyPressed:
      switch(event.Key.Code){
      case sf::Key::N:
      case sf::Key::Return:
	next_step = true;
	break;
      case sf::Key::F:
      case sf::Key::Space:
	free_running = !free_running;
	break;
      case sf::Key::P:
	_drawer->zoom(1.1);
	break;
      case sf::Key::L:
	_drawer->zoom(0.9);
	break;
      // case sf::Key::R:
      // 	next_run = true;
      // 	break;
      // case sf::Key::O:
      // 	optimization_active = !optimization_active;
      // 	break;
      case sf::Key::D:
	draw = !draw;
	break;
      case sf::Key::Escape:
	exit(0);
	break;
      default:
	break;
      }	// end of switch event.key.code
      break;
    case sf::Event::MouseWheelMoved:
      if(event.MouseWheel.Delta<0) _drawer->zoom(1.1f);
      else _drawer->zoom(0.9f);
      break;
    default:
      break;
    }	// end of switch event.type
  }
}


void init(){

  next_step = false;
  free_running = false;
  draw = true;
  
  // creating the drawer
  _drawer = new rdrawer::RobotDrawer();

  // allocating the optimizer
  optimizer = new g2o::SparseOptimizer();
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
  
  optimizer->setAlgorithm(solver);
  
}

int main(int argc, char** argv){
  std::cout << "Visualizer" << std::endl;
  std::cout << "==================" << std::endl;
  
  // check the input
  if(argc < 2){
    std::cout << "Usage: visualize <input_file>" << std::endl;
    exit(0);
  }
  
  init();
  
  readFromFile(argv[1]);
  
  std::cout << "Loaded " << poses.size() << " poses" << std::endl;
  std::cout << "Loaded " << landmarks.size() << " landmarks" << std::endl;
  std::cout << "Loaded " << edgesPoses.size() << " edges pose-pose" << std::endl;
  std::cout << "Loaded " << edgesLandmarks.size() << " edges pose-landmark" << std::endl;
  

  std::vector<Landmark *> seen_landmarks;
  
  unsigned int last_step = MAX_STEPS;
  if(poses.size() < MAX_STEPS){
    last_step = (unsigned int) poses.size();
  }
  
  for(unsigned int i=0; i<last_step; i++){
    
    _drawer->clearAll();
    
    std::cout << "step " << i << std::endl;
    
    next_step = false;
    
    
    for(unsigned int p=0; p<i; p++){
      if(draw){
	_drawer->addTrajectoryStep(poses[p]->estimate()[0] * DRAWER_MULTIPLIER, poses[p]->estimate()[1] * DRAWER_MULTIPLIER, poses[p]->estimate()[2]);
      }
    }
    
    for(unsigned int e=0; e<poses[i]->_edges.size(); e++){
      Landmark * l = (Landmark *) (poses[i]->_edges[e]->vertices()[1]);
      bool found = false;
      for(unsigned int seen=0; seen<seen_landmarks.size(); seen++){
	if(seen_landmarks[seen]->id() == l->id()){
	  found = true;
	}
      }
      if(!found){
	seen_landmarks.push_back(l);
      }
    }
    
    for(unsigned int s=0; s<seen_landmarks.size(); s++){
      if(draw){
	_drawer->addAndCreateLandmark(seen_landmarks[s]->estimate()[0] * DRAWER_MULTIPLIER, seen_landmarks[s]->estimate()[1] * DRAWER_MULTIPLIER);
      }
    
    
      rdrawer::RobotPose drawer_robot_pose;
      drawer_robot_pose.x = poses[i]->estimate()[0] * DRAWER_MULTIPLIER;
      drawer_robot_pose.y = poses[i]->estimate()[1] * DRAWER_MULTIPLIER;
      drawer_robot_pose.theta = poses[i]->estimate()[2];
      _drawer->setRobotPose(drawer_robot_pose);
    }
    
    

    if(free_running){
      handleEvents(_drawer->getWindow());
      if(draw) _drawer->draw();
      usleep(20);
    }
  
    while(!next_step && !free_running){
      handleEvents(_drawer->getWindow());
      _drawer->draw();
      usleep(20);
    }
    
  }
  
  _drawer->clearAll();
  
  // add stuff to the optimizer
  for(unsigned int i=0; i<last_step; i++){
    optimizer->addVertex(poses[i]);
    
    if(i>0){
      optimizer->addEdge(edgesPoses[i-1]);
    }
    
    for(unsigned int e=0; e<poses[i]->_edges.size(); e++){
      Landmark * l = (Landmark *) (poses[i]->_edges[e]->vertices()[1]);
      bool found = false;
      for(unsigned int seen=0; seen<seen_landmarks.size(); seen++){
	if(seen_landmarks[seen]->id() == l->id()){
	  found = true;
	}
      }
      if(found){
	if(!l->alreadyInserted){
	  l->alreadyInserted = true;
	  optimizer->addVertex(l);
	}
	optimizer->addEdge(poses[i]->_edges[e]);
      }
    }
  }
  
  std::ofstream preopt("preoptim.g2o");
  optimizer->save(preopt);
  
  g2o::OptimizableGraph::EdgeSet edgesToOptimize;
  
  //unsigned int first_index = poses.size()-500;
  unsigned int first_index = 0;
  
  for(unsigned int i=first_index; i<poses.size(); i++){
    if (i>0) edgesToOptimize.insert(edgesPoses[i-1]);
    Pose * p = poses[i];
    for(unsigned int e=0; e<p->_edges.size(); e++){
      edgesToOptimize.insert(p->_edges[e]);
    }
  }
  
  poses[first_index]->setFixed(true);
  optimizer->initializeOptimization(edgesToOptimize);
  optimizer->optimize(10);
  
  for(unsigned int p=0; p<last_step; p++){
    _drawer->addTrajectoryStep(poses[p]->estimate()[0] * DRAWER_MULTIPLIER, poses[p]->estimate()[1] * DRAWER_MULTIPLIER, poses[p]->estimate()[2]);
  }
  
  for(unsigned int s=0; s<seen_landmarks.size(); s++){
    _drawer->addAndCreateLandmark(seen_landmarks[s]->estimate()[0] * DRAWER_MULTIPLIER, seen_landmarks[s]->estimate()[1] * DRAWER_MULTIPLIER);
  }
  
  rdrawer::RobotPose drawer_robot_pose;
  drawer_robot_pose.x = poses[last_step-1]->estimate()[0] * DRAWER_MULTIPLIER;
  drawer_robot_pose.y = poses[last_step-1]->estimate()[1] * DRAWER_MULTIPLIER;
  drawer_robot_pose.theta = poses[last_step-1]->estimate()[2];
  _drawer->setRobotPose(drawer_robot_pose);
    
  
  while(true){
    handleEvents(_drawer->getWindow());
    _drawer->draw();
    usleep(200);
  }
  _drawer->clearAll();
  
  return 0;
}
