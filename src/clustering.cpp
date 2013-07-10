#include <stdlib.h> // rand
#include <math.h>

template<typename T> T absolute(T a){
  if(a<0) return -a;
  return a;
}



// classifies the input R-D points in K clusters using the K-means algorithm
// returns FALSE if it was not possible to clusterize the data, else returns TRUE
// inputs:	• v - array containing the values
// 		• R - dimension
//		• label - array where labels will be put
//		• means - means will be put here. This is supposed to be an array of dimension R*K
//		• N - the number of samples
//		• K - the number of clusters to generate
// NOTE ON THE INPUT:	the values are supposed to be on a single array, the R parameter will specify the size of each sample.
// 			E.G: 3D samples - < 1x, 1y, 1z, 2x, 2y, 2z, ....> (R here is 3)
bool Kmeans(double * v, unsigned int R, int * label, double * means, unsigned int N, unsigned int K){
  
  if(K>N) return false;
  
  for(unsigned int i=0; i<N; i++){
    label[i] = -1;
  }
  
  // std::cout << "initializing the means" << std::endl;
  
  // initialize the means using points in the set
  {
    int chosen[K];
    
    
    // choose wisely: choose points that maximize the mutual distance
    // first choose one with the biggest/smallest value on one of the dimensions
    double mins[R];
    int mins_ind[R];
    double maxs[R];
    int maxs_ind[R];
    for(unsigned int r=0; r<R; r++){
      mins[r] = v[r];
      mins_ind[r] = 0;
      maxs[r] = v[r];
      maxs_ind[r] = 0;
    }
    for(unsigned int n=1; n<N; n++){
      unsigned int index = n*R;
      for(unsigned int r=0; r<R; r++){
	if(v[index+r] < mins[r]){
	  mins[r] = v[index+r];
	  mins_ind[r] = n;
	}
	else if(v[index+r] > maxs[r]){
	  maxs[r] = v[index+r];
	  maxs_ind[r] = n;
	}
      }
    }
    if(rand()%2){
      int rnd = rand();
      chosen[0] = mins_ind[rnd%R];
    }
    else{
      int rnd = rand();
      chosen[0] = maxs_ind[rnd%R];
    }
    int m0_index = chosen[0]*R;
    for(unsigned int r=0; r<R; r++){
      means[r] = v[m0_index+r];
    }
    
    for(unsigned int i=1; i<K; i++){
      chosen[i] = -1;
      
      double max_dist=0;
      unsigned int m = 0;	// index (out of N) of the point that maximizes the distance from the already setted means.
      
      for(unsigned int n=0; n<N; n++){
	unsigned int index = n*R;
	double distance = 0;
	for(unsigned int j=0; j<i; j++){ // for each of the already computed means
	  unsigned int mindex = j*R;
	  for(unsigned int r=0; r<R; r++){
	    distance += pow(means[mindex+r] - v[index+r],2);
	  }
	}
	// distance = sqrt(distance)/i; // not needed,   x > g ---> f(x) > f(g)
	if(distance > max_dist){
	  max_dist = distance;
	  m = n;
	}
	
      }
      
      chosen[i] = m;
      
      int index = chosen[i]*R;
      for(unsigned int r=0; r<R; r++){
	means[i*R+r] = v[index+r];
      }
    }
  }
  
  int iterations = 0;
  
  bool something_changed = true;
  while(something_changed && iterations<1000){
    iterations++;
    something_changed = false;
    
    // std::cout << "itaration #" << iterations << std::endl;
    // std::cout << "current means:" << std::endl;
    // for(unsigned int k=0; k<K; k++){
    //   std::cout << "\t[";
    //   for(unsigned int i=0; i<R; i++){
    // 	std::cout << " " << means[k][i];
    //   }
    //   std::cout << " ]" << std::endl;
    // }
    
    double new_m[K][R];	// accumulators
    double cardinality[K];

    for(unsigned int k=0; k<K; k++){
      for(unsigned int i=0; i<R; i++){
	new_m[k][i] = 0;
      }
      cardinality[k] = 0;
    }
    
    for(unsigned int n=0; n<N; n++){	// for every element
      int vindex = n*R;
      
      int prev_label = label[n];
      double best_distance = 0;
      for(unsigned int r=0; r<R; r++){
	best_distance += pow((means[r] - v[vindex+r]), 2);
      }
      best_distance = sqrt(best_distance);
      
      label[n] = 0;
      
      // compute the distance from all the means and choose the smaller one
      for(unsigned int k=1; k<K; k++){
	double distance = 0;
	unsigned int kindex = k*R;
	for(unsigned int r=0; r<R; r++){
	  distance += pow((means[kindex+r] - v[vindex+r]), 2);
	}
	distance = sqrt(distance);
	if(distance < best_distance){
	  best_distance = distance;
	  label[n] = k;
	}
      }
      
      if(label[n] != prev_label){
	something_changed = true;
      }
      
      // update accumulators
      for(unsigned int r=0; r<R; r++){
	new_m[label[n]][r] += v[vindex+r];
      }
      
      cardinality[label[n]] = cardinality[label[n]] + 1;
    }
    
    // compute new centroids
    for(unsigned int k=0; k<K; k++){
      double prev[R];
      unsigned int kindex = k*R;
      if(cardinality[k] == 0){ // just leave this poor mean where it is
	continue;
      }
      for(unsigned int r=0; r<R; r++){
	prev[r] = means[kindex+r];
	means[kindex+r] = new_m[k][r]/cardinality[k];
	if(means[kindex+r] != prev[r]){
	  something_changed = true;
	}
      }
    }
  }
  
  // std::cout << "performed " << iterations << " iterations" << std::endl;
  // std::cout << "labels: [";
  // for(unsigned int i=0; i<N; i++){
  //   std::cout << " " << label[i];
  // }
  // std::cout << " ]" << std::endl;
    
  return true;
}





// runs Kmeans more than once, and stops when two consecutive iterations bring to the same clusterization
bool stableKmeans(double * v, unsigned int R, int * label, double * means, unsigned int N, unsigned int K){
  
  int run_count = 0;	// counts the number of times Kmeans is called
  bool is_stable = false;	// set to true when the solution is considered "stable"
  
  int previous_labels[N];
  
  while(!is_stable){
    run_count++;
    is_stable = true;
    
    // store previous solution
    if(run_count>1){
      for(unsigned int i=0; i<N; i++){
	previous_labels[i] = label[i];
      }
    }
    else{
      for(unsigned int i=0; i<N; i++){
	previous_labels[i] = -1;
      }
    }
    
    if(!Kmeans(v,R,label,means,N,K)){
      return false;
    }
    if(run_count < 2){
      is_stable = false;
      continue;
    }
    // else
    // check if something changed
    int first_el[K];
    int corresponding[K];
    for(unsigned int k=0; k<K; k++){	// for every cluster
      first_el[k] = -1;
      // look for the first element in that cluster
      for(unsigned int i=0; i<N; i++){
	if(label[i] == k){
	  first_el[k] = i;
	  break;
	}
      }
      if(first_el[k] != -1){
	corresponding[k] = previous_labels[first_el[k]];
      }
      // else no sample is now close to this mean
      
    }
    
    // now check all the correspondancies
    for(unsigned int n=0; n<N; n++){
      int k = label[n];
      if(previous_labels[n] != corresponding[k]){
	is_stable = false;
	break;
      }
    }
  }
  
  //std::cout << "runs count: " << run_count << std::endl;
}




// same as above, but also stores the means at the given pointer. NOTE: this creates the means array!
// moreover, the user can specify the maximum number of clusters to be created
// classifies the input R-D points in K clusters using the K-means algorithm
// increments K while <a certain condition> is verified
// returns -1 if it was not possible to clusterize the data, else returns the number of clusters
// inputs:	• v - array containing the values
// 		• R - dimension
//		• label - array where labels will be put
//		• N - the number of samples
// NOTE ON THE INPUT:	the values are supposed to be on a single array, the R parameter will specify the size of each sample.
// 			E.G: 3D samples - < 1x, 1y, 1z, 2x, 2y, 2z, ....> (R here is 3)
int findClusters(double * input, unsigned int R, int * label, double ** ret_means, unsigned int N, unsigned int max_clusters){
  
  // normalize the input
  double normalizer = absolute(input[0]);
  unsigned int inputsize = R*N;
  for(unsigned int i=0; i<inputsize; i++){
    double value = absolute(input[i]);
    if( value > normalizer) normalizer = value;
  }
  double v[inputsize];
  for(unsigned int i=0; i<inputsize; i++){
    v[i] = input[i] / normalizer;
  }
  
  unsigned int K = 2;
  
  double best_score;
  int best_labels[N];
  int best_K;  
  bool first_iteration = true;
  
  while(K<=N && K<=max_clusters){
    
    //std::cout << std::endl << "K = " << K << std::endl;
    
    double * means = (double*) malloc(K*R*sizeof(double));
    
    stableKmeans(v,R,label,means,N,K);
    
    // std::cout << "means:";
    // for(unsigned int k=0; k<K; k++){
    //   std::cout << " <";
    //   unsigned int kindex = k*R;
    //   for(unsigned int r=0; r<R; r++){
    // 	std::cout << " " << means[kindex+r];
    //   }
    //   std::cout << " >";
    // }
    // std::cout << std::endl;
    
        
    // compute the score for the current asset
    double errors[K];
    double cardin[K];
    
    for(unsigned int k=0; k<K; k++){
      unsigned int kindex = k*R;
      errors[k] = 0;
      cardin[k] = 0;
      
      // compute the average error
      for(unsigned int n=0; n<N; n++){
	if(label[n] == k){
	  double distance = 0;
	  unsigned int vindex = n*R;
	  for(unsigned int r=0; r<R; r++){
	    distance += pow(means[kindex+r] - v[vindex+r],2);
	  }
	  distance = sqrt(distance);
	  errors[k] += distance;
	  cardin[k]++;
	}
      }
      if(cardin[k]>0){
	//errors[k] = errors[k]/(double)cardin[k];
      }
    }
    
    double total_score = 0;
    double kscores[K];
    
    for(unsigned int k=0; k<K; k++){
      kscores[k] = -1;
      unsigned int kindex = k*R;
      
      for(unsigned int j=0; j<K; j++){
	if(j==k) continue;
	unsigned int jindex = j*R;
	
	double distance = 0;
	for(unsigned int r=0; r<R; r++){
	  distance += pow(means[kindex+r] - means[jindex+r], 2);
	}
	distance = sqrt(distance);
	
	
	double score;
	if(cardin[k] != 0  &&  cardin[k] != 0){
	  score = (errors[k] + errors[j]) / distance;
	}
	else{
	  score = 1000;
	}
	
	if(kscores[k] == -1 || score > kscores[k]){
	  kscores[k] = score;
	}
	
      }
      
      total_score += kscores[k];
    }
    
    double asset_score = total_score/K;
    
    //std::cout << "score for K = " << K << ": " << asset_score << std::endl;
    
    if(first_iteration || asset_score < best_score){
      //std::cout << "\tit's a new record!!!" << std::endl;
      first_iteration = false;
      best_score = asset_score;
      best_K = K;
      for(unsigned int n=0; n<N; n++){
	best_labels[n] = label[n];
      }
      free(*ret_means);
      *ret_means = means;
    }
    else{
      free(means);
    }
    
    K++;
    
  }
  
  for(unsigned int n=0; n<N; n++){
    label[n] = best_labels[n];
  }
  
  // unnormalize means
  for(unsigned int i=0; i<best_K*R; i++){
    (*ret_means)[i] = (*ret_means)[i] * normalizer;
  }
  
  return best_K;
  
}





int findClusters(double * input, unsigned int R, int * label, unsigned int N){
  double * means = (double *) malloc(sizeof(double));
  
  return findClusters(input, R, label, &means, N, 8);
}



int stableFindClusters(double * input, unsigned int R, int * label, double ** ret_means, unsigned int N, unsigned int max_clusters){
  
  int run_count = 0;
  int prev_labels[N];
  int prev_K = -1;
  
  bool stable = false;
  
  while(!stable && run_count < 10){
    run_count++;
    stable = true;
    
    if(run_count>1){
      for(unsigned int i=0; i<N; i++){
	prev_labels[i] = label[i];
      }
    }
    else{
      for(unsigned int i=0; i<N; i++){
	prev_labels[i] = -1;
      }      
    }
    
    int fc = findClusters(input, R, label, ret_means, N, max_clusters);
    
    if(run_count < 2  ||  fc != prev_K){
      stable = false;
      prev_K = fc;
      continue;
    }
    // else
    // check if something changed
    
    prev_K = fc;
    
    int first_el[fc];
    int corresponding[fc];
    for(unsigned int k=0; k<fc; k++){	// for every cluster
      first_el[k] = -1;
      // look for the first element in that cluster
      for(unsigned int i=0; i<N; i++){
	if(label[i] == k){
	  first_el[k] = i;
	  break;
	}
      }
      if(first_el[k] != -1){
	corresponding[k] = prev_labels[first_el[k]];
      }
      // else no sample is now close to this mean
      
    }
    
    // now check all the correspondancies
    for(unsigned int n=0; n<N; n++){
      int k = label[n];
      if(prev_labels[n] != corresponding[k]){
	stable = false;
	break;
      }
    }
  }
  
  return prev_K;
    
}





int stableFindClusters(double * input, unsigned int R, int * label, unsigned int N){
  double * means = (double *) malloc(sizeof(double));
  
  return stableFindClusters(input, R, label, &means, N, 8);
}
