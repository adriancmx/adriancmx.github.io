
######################################
#            Multiple Try            #
######################################

# This is an example of a simple multiple try structure.

x=rep(0,5000)
x0=-1
i=1
while (i<5001) {
  #we will use normal distribution as the proposal distribution
  #as k=2, we generate two independent samples from the normal distributions.
  y1<-rnorm(1,x0,1)
  y2<-rnorm(1,x0,1) 
  
  ##the weight function is calculated，
  ##lambda should be some symmetric non-negative function, 
  ##here we could choose lambda function to be absolution function.
  w1<-density_f(y1)*dnorm(x0,y1,1)*abs(x0-y1) 
  w2<-density_f(y2)*dnorm(x0,y2,1)*abs(x0-y2)
  
  ##choose y between y1 and y2.
  u<-runif(1,0,1)
  s<-w1/(w1+w2) 
  if(u<s){
    y=y1 
  }else{ 
    y=y2
  }
  
  ##choose reference x
  x1<-rnorm(1,y,1)
  x2<-x0 
  w3<-density_f(x1)*dnorm(y,x1,1)*abs(x1-y) 
  w4<-density_f(x2)*dnorm(y,x1,1)*abs(x2-y) 
  
  ##acceptance probability is decided for y.
  prob<-(w1+w2)/(w3+w4) 
  u<-runif(1,0,1)
  if(u<min(1,prob)){
    x[i]=y
    i=i+1
    x0=y
  }
}


######################################
#            Likelihood              #
######################################

# This is matlab version of the Gaussian Mixture model likelihood

#############################
# Part 1 likelihood for merge
#############################

# % calculate Likelihood
# Model1 = priorModel.copy();
# term1 = 0;
# for j = 1:size(c_i_split_idx, 2)
# if j == 1
# term1 = term1 + Model1.GetPrior(X, c_i_split_idx(j));
# Model1 = Model1.AddOne(X, c_i_split_idx(j));
# else
#   term1 = term1 + Model1.GetPosteriorPred(X, c_i_split_idx(j));
# Model1 = Model1.AddOne(X, c_i_split_idx(j));
# end
# end
# Model2 = priorModel.copy();
# term2 = 0;
# for j = 1:size(c_j_split_idx, 2)
# if j == 1
# term2 = term2 + Model2.GetPrior(X, c_j_split_idx(j));
# Model2 = Model2.AddOne(X, c_j_split_idx(j));
# else
#   term2 = term2 + Model2.GetPosteriorPred(X, c_j_split_idx(j));
# Model2 = Model2.AddOne(X, c_j_split_idx(j));
# end
# end
# Model3 = priorModel.copy();
# term3 = 0;
# for j = 1:size(S_with_ij_idx, 2)
# if j == 1
# term3 = term3 + Model3.GetPrior(X, S_with_ij_idx(j));
# Model3 = Model3.AddOne(X, S_with_ij_idx(j));
# else
#   term3 = term3 + Model3.GetPosteriorPred(X, S_with_ij_idx(j));
# Model3 = Model3.AddOne(X, S_with_ij_idx(j));
# end
# end
# Likelihood = term1+term2-term3;


#############################
# Part 2 likelihood for split
#############################


# % calculate Likelihood
# Model1 = priorModel.copy();
# term1 = 0;
# for j = 1:size(S_with_ij_idx, 2)
# if j == 1
# term1 = term1 + Model1.GetPrior(X, S_with_ij_idx(j));
# Model1 = Model1.AddOne(X, S_with_ij_idx(j));
# else
#   term1 = term1 + Model1.GetPosteriorPred(X, S_with_ij_idx(j));
# Model1 = Model1.AddOne(X, S_with_ij_idx(j));
# end
# end
# Model2 = priorModel.copy();
# term2 = 0;
# for j = 1:size(c_i_idx, 2)
# if j == 1
# term2 = term2 + Model2.GetPrior(X, c_i_idx(j));
# Model2 = Model2.AddOne(X, c_i_idx(j));
# else
#   term2 = term2 + Model2.GetPosteriorPred(X, c_i_idx(j));
# Model2 = Model2.AddOne(X, c_i_idx(j));
# end
# end
# Model3 = priorModel.copy();
# term3 = 0;
# for j = 1:size(c_j_idx, 2)
# if j == 1
# term3 = term3 + Model3.GetPrior(X, c_j_idx(j));
# Model3 = Model3.AddOne(X, c_j_idx(j));
# else
#   term3 = term3 + Model3.GetPosteriorPred(X, c_j_idx(j));
# Model3 = Model3.AddOne(X, c_j_idx(j));
# end
# end
# Likelihood = term1-term2-term3;