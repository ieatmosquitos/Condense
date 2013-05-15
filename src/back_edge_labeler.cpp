#include <Eigen/Dense>
#include "edge_labeler.h"
#include "g2o/stuff/unscented.h"

namespace g2o {

  using namespace std;

  typedef SigmaPoint<VectorXd> MySigmaPoint;

  EdgeLabeler::EdgeLabeler(SparseOptimizer* optimizer) {
    _optimizer = optimizer;
  }
  static int nummero = 0;
  int EdgeLabeler::labelEdges(std::set<OptimizableGraph::Edge*>& edges){
    // assume the system is "solved"
 
    // compute the sparse pattern of the inverse
    std::set<std::pair<int, int> > pattern;
    for (std::set<OptimizableGraph::Edge*>::iterator it=edges.begin(); it!=edges.end(); it++){
      augmentSparsePattern(pattern, *it);
    }
    

    SparseBlockMatrix<MatrixXd> spInv;

    bool result = computePartialInverse(spInv, pattern);
    //cerr << "partial inverse computed = " << result << endl;
    //cerr << "non zero blocks" << spInv.nonZeroBlocks() << endl;

    if (! result ){
      return -1;
    }
    int count=0;
    for (std::set<OptimizableGraph::Edge*>::iterator it=edges.begin(); it!=edges.end(); it++){
      count += labelEdge(spInv, *it) ? 1 : 0;
    }
    return count;
  }

  void EdgeLabeler::augmentSparsePattern(std::set<std::pair<int, int> >& pattern, OptimizableGraph::Edge* e){
    for (size_t i=0; i<e->vertices().size(); i++){
      const OptimizableGraph::Vertex* v=(const OptimizableGraph::Vertex*) e->vertices()[i];
      int ti=v->hessianIndex();
      if (ti==-1)
	continue;
      for (size_t j=i; j<e->vertices().size(); j++){
	const OptimizableGraph::Vertex* v=(const OptimizableGraph::Vertex*) e->vertices()[j];
	int tj = v->hessianIndex();
	if (tj==-1)
	  continue;
	if (ti<tj){
	  pattern.insert(std::make_pair(ti, tj));
	} else {
	  pattern.insert(std::make_pair(tj, ti));
	}
      }
    }
  }
  
  bool EdgeLabeler::computePartialInverse(SparseBlockMatrix<MatrixXd>& spinv, const std::set<std::pair<int,int> >& pattern){
    std::vector<std::pair<int, int> > blockIndices(pattern.size());
    // Why this does not work???
    //std::copy(pattern.begin(),pattern.end(),blockIndices.begin());

    int k=0;
    for(std::set<std::pair<int, int> >::const_iterator it= pattern.begin(); it!=pattern.end(); it++){
      blockIndices[k++]=*it;
    }

    //cerr << "sparse pattern contains " << blockIndices.size() << " blocks" << endl;
    return _optimizer->computeMarginals(spinv, blockIndices);
  }

  bool EdgeLabeler::labelEdge( const SparseBlockMatrix<MatrixXd>& spinv, OptimizableGraph::Edge* e){

    Eigen::Map<MatrixXd> info(e->informationData(), e->dimension(), e->dimension());
    // cerr << "original information matrix" << endl;
    // cerr << info << endl;

    int maxDim=0;
    for (size_t i=0; i<e->vertices().size(); i++){
      const OptimizableGraph::Vertex* v=(const OptimizableGraph::Vertex*) e->vertices()[i];
      int ti=v->hessianIndex();
      if (ti==-1)
	continue;
      maxDim+=v->minimalEstimateDimension();
    }


    //cerr << "maxDim= " << maxDim << endl;
    MatrixXd cov(maxDim, maxDim);
    int cumRow=0;
    bool somethingWrong = false;
    for (size_t i=0; i<e->vertices().size(); i++){
      const OptimizableGraph::Vertex* vr=(const OptimizableGraph::Vertex*) e->vertices()[i];
      int ti=vr->hessianIndex();
      if (ti>-1) {
	int cumCol=0;
  
	for (size_t j=0; j<e->vertices().size(); j++){
	  const OptimizableGraph::Vertex* vc=(const OptimizableGraph::Vertex*) e->vertices()[j];
	  int tj = vc->hessianIndex();
	  // cerr << "ti = " << ti << "\ttj = " << tj << endl;
	  if (tj>-1){
	    // cerr << "tj>-1" << endl;
	    // cerr << "ti=" << ti << " tj=" << tj 
	    //    << " cumRow=" << cumRow << " cumCol=" << cumCol << endl;
	    if (ti<=tj){
              // cerr << "ti<tj" << endl;
	      // cerr << "spinv.block(ti,tj)" << endl;
	      // cerr << spinv.block(ti, tj) << endl;
	      // cerr << "spinv.block(tj,ti)" << endl;
	      // cerr << spinv.block(tj, ti) << endl;
		
	      // cerr << "cblock_ptr" << spinv.block(ti, tj) << endl;
	      // cerr << "cblock.size=" << spinv.block(ti, tj)->rows() << "," << spinv.block(ti, tj)->cols() << endl;
	      // cerr << "cblock" << endl;
	      // cerr << *spinv.block(ti, tj) << endl;
        
	      int r;
	      int c;
	      bool mustTranspose = false;
	      if(!spinv.block(ti,tj)){
		r = tj;
		c = ti;
		mustTranspose = true;
	      }
	      else{
		r=ti;
		c=tj;
	      }
	      assert(spinv.block(r, c));
	
	      if(mustTranspose){
		cerr << "AAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
		somethingWrong = true;
		cov.block(cumRow, cumCol, vr->minimalEstimateDimension(), vc->minimalEstimateDimension()) = spinv.block(r,c)->transpose();
	      }
	      else{
		cov.block(cumRow, cumCol, vr->minimalEstimateDimension(), vc->minimalEstimateDimension()) = *(spinv.block(r,c));
	      }
	    } else { // ti > tj

	      // cerr << "tj<ti" << endl;
	      // cerr << "spinv.block(ti,tj)" << endl;
	      // cerr << spinv.block(ti, tj) << endl;
	      // cerr << "spinv.block(tj,ti)" << endl;
	      // cerr << spinv.block(tj, ti) << endl;
	      int r;
	      int c;
	      bool mustTranspose = false;
	      if(!spinv.block(ti,tj)){
		r = tj;
		c = ti;
		mustTranspose = true;
	      }
	      else{
		r=ti;
		c=tj;
	      }
	      assert(spinv.block(r, c));
	
	      // cerr << "cblock.size=" << spinv.block(tj, ti)->cols() << "," << spinv.block(tj, ti)->rows() << endl;
	      // cerr << "cblock" << endl;
	      // cerr << spinv.block(tj, ti)->transpose() << endl;
	
	      if(mustTranspose){
		cov.block(cumRow, cumCol, vr->minimalEstimateDimension(), vc->minimalEstimateDimension()) = spinv.block(r,c)->transpose();
	      }
	      else{
		cerr << "BUAAAAAAAA" << endl;
		somethingWrong = true;
		cov.block(cumRow, cumCol, vr->minimalEstimateDimension(), vc->minimalEstimateDimension()) = *(spinv.block(r,c));
	      }
        
	    }
	    cumCol += vc->minimalEstimateDimension();
	  }
	}
	cumRow += vr->minimalEstimateDimension();
      }
    }
    if (somethingWrong){
      char buf[1024]; sprintf(buf, "irina_%05d.dat",nummero);
      nummero++;
      spinv.writeOctave(buf, true);
    }
    // cerr << "covariance assembled" << endl;
    // cerr << cov << endl;
    // now cov contains the aggregate marginals of the state variables in the edge
    VectorXd incMean(maxDim);
    incMean.fill(0);
    std::vector<MySigmaPoint, Eigen::aligned_allocator<MySigmaPoint> > incrementPoints;
    if (! sampleUnscented(incrementPoints, incMean, cov)){
      cerr << "sampleUnscented fail" << endl;
      return false;
    }
    // now determine the zero-error measure by applying the error function of the edge
    // with a zero measurement
    // TODO!!!
    bool smss = e->setMeasurementFromState();
    if (! smss) {
      cerr << "FATAL: Edge::setMeasurementFromState() not implemented" << endl;
    }
    assert(smss && "Edge::setMeasurementFromState() not implemented");



    //std::vector<MySigmaPoint> globalPoints(incrementPoints.size());
    std::vector<MySigmaPoint, Eigen::aligned_allocator<MySigmaPoint> > errorPoints(incrementPoints.size());

    // for each sigma point, project it to the global space, by considering those variables
    // that are involved
    //cerr << "sigma points are extracted, remapping to measurement space" << endl;
    for (size_t i=0; i<incrementPoints.size(); i++) {
      int cumPos=0;
      //VectorXd globalPoint(maxDim);

      // push all the "active" state variables
      for (size_t j=0; j<e->vertices().size(); j++){
        OptimizableGraph::Vertex* vr=(OptimizableGraph::Vertex*) e->vertices()[j];
        int tj=vr->hessianIndex();
        if (tj==-1)
          continue;
        vr->push();
      }

      
      for (size_t j=0; j<e->vertices().size(); j++){
        OptimizableGraph::Vertex* vr=(OptimizableGraph::Vertex*) e->vertices()[j];
        int tj=vr->hessianIndex();
        if (tj==-1)
          continue;
        vr->oplus(&incrementPoints[i]._sample[cumPos]);
        //assert(vr->getMinimalEstimateData(&globalPoint[cumPos]) && "Vertex::getMinimalEstimateData(...) not implemented");
        cumPos+=vr->minimalEstimateDimension();
      }

      // construct the sigma point in the global space
      // globalPoints[i]._sample=globalPoint;
      // globalPoints[i]._wi=incrementPoints[i]._wi;
      // globalPoints[i]._wp=incrementPoints[i]._wp;

      // construct the sigma point in the error space
      e->computeError();
      Map<VectorXd> errorPoint(e->errorData(),e->dimension());

      errorPoints[i]._sample=errorPoint;
      errorPoints[i]._wi=incrementPoints[i]._wi;
      errorPoints[i]._wp=incrementPoints[i]._wp;
      
      // pop all the "active" state variables
      for (size_t j=0; j<e->vertices().size(); j++){
        OptimizableGraph::Vertex* vr=(OptimizableGraph::Vertex*) e->vertices()[j];
        int tj=vr->hessianIndex();
        if (tj==-1)
          continue;
        vr->pop();
      }

    }

    // reconstruct the covariance of the error by the sigma points 
    MatrixXd errorCov(e->dimension(), e->dimension());
    VectorXd errorMean(e->dimension()); 
    reconstructGaussian(errorMean, errorCov, errorPoints);
    info=errorCov.inverse();
    
    // cerr << "remapped information matrix" << endl;
    // cerr << info << endl;
    return true;
  }

}
