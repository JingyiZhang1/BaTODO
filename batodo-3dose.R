#######################################################################################################
#             BaTODO: A Novel Multi-Dose Randomized Basket Trial Design for Dose Testing              #
#                          and Optimization in Oncology Drug Development                              # 
#                                                                                                     #
#                                     rlin@mdanderson.org                                             #
#######################################################################################################
# ----------------------------------------------------------------------------------------------------#
#        Main input variables                                                                         #
#        ntrial            -->  number of replications                                                #
#        nsample           -->  Maximum sample size for each dose arm                                 #
#        n1                -->  Stage 1 sample size  for each dose arm                                #
#        peff              -->  True toxicity rate                                                    #
#        ptox              -->  True efficacy rate                                                    #
#        theta0T           -->  Null toxicity rate                                                    #   
#        theta0            -->  Null efficacy rate                                                    #
#        weightu           -->  penalty of toxicity in the utility function                           #
#        delta1            -->  maximum acceptable difference in utility                              #
#        cutoff_utility    -->  cutoff value for identifying optimal doses                            #
# ----------------------------------------------------------------------------------------------------#


batodo_3dose <- function(ntrial, peff, ptox, n1, nsample, theta0, theta0T, weightu, delta1, cutoff_utility){
  
  library(R2jags)
  library(truncnorm)
  
  set.seed(6)
  
  nmodel <- 15
  ## possible model
  model <- list()
  model[[1]] <- c(1,1,1,1); model[[2]] <- c(1,2,2,2)
  model[[3]] <- c(2,1,2,2); model[[4]] <- c(2,2,1,2)
  model[[5]] <- c(2,2,2,1); model[[6]] <- c(1,1,2,2)
  model[[7]] <- c(1,2,1,2); model[[8]] <- c(1,2,2,1)
  model[[9]] <- c(1,1,2,3); model[[10]] <- c(1,2,1,3)
  model[[11]] <- c(1,2,3,1); model[[12]] <- c(2,1,1,3)
  model[[13]] <- c(2,1,3,1); model[[14]] <- c(2,3,1,1)
  model[[15]] <- c(1,2,3,4)
  
  ## all clusters
  model_cluster <- list()
  model_cluster[[1]] <- c(1,2,3,4); model_cluster[[2]] <- c(1,2,3)
  model_cluster[[3]] <- c(1,2,4);   model_cluster[[4]] <- c(1,3,4)
  model_cluster[[5]] <- c(2,3,4);   model_cluster[[6]] <- c(1,2)
  model_cluster[[7]] <- c(1,3);     model_cluster[[8]] <- c(1,4)
  model_cluster[[9]] <- c(2,3);     model_cluster[[10]] <- c(2,4)
  model_cluster[[11]] <- c(3,4);    model_cluster[[12]] <- 1
  model_cluster[[13]] <- 2;         model_cluster[[14]] <- 3
  model_cluster[[15]] <- 4
  
  
  zeta.mu <- 2
  zeta.tau <- 2
  tau.theta1 <- 10
  nclusters <- 15
  nsampling <- 10000
  nburnin <- 1000
  cutoff_interim <- cutoff_interimT <- 0.2505
  cutoff_final <- cutoff_finalT <- 0.844
  prior_model_weight <- rep(1,nmodel)/nmodel
  
  nks <- dim(peff)[1] # number of tumor types
  nds <- dim(peff)[2] # number of doses
  
  ## identify the true optimal and promising doses
  opt_u <- function(peff,ptox,theta0,theta0T,delta1,weightu){
    uu <- round(peff-ptox*weightu,4)
    uu[which(peff<=theta0 | ptox>=theta0T)] <- -100 ## exclude overly toxic or futile doses
    
    opt <- promising <- NULL
    if(sum(uu>-100)>0){
      opt <- which(uu==max(uu))
      promising <- which(uu>=(max(uu)-delta1))
    }
    
    if(length(opt)>0){return(list(opt=opt,promising=promising))}else{return(0)}
  }
  
  opt_dose_uu <- opt_dose_uugo <- matrix(0, nrow = nks, ncol = nds)
  
  for(ik in 1:nks){
    temp <- opt_u(peff[ik,],ptox[ik,],theta0,theta0T,delta1,weightu)
    if(is.list(temp)){
      opt_dose_uu[ik,temp$opt] <- 1
      opt_dose_uugo[ik,temp$promising] <- 1
    }
  }
  
  which.is.max <- function(a){ return(which(a==max(a))) }
  
  ## parameter initialization
  data.IA <- dataE.IA <- dataT.IA <- matrix(n1,nrow=nks,ncol=nds)
  ia.postE  <- ia_pE_bma <- matrix(0,nrow=nks,ncol=nds)
  ia.postT  <- ia_pT_bma <- matrix(0,nrow=nks,ncol=nds)
  
  data.FA <- dataE.FA <- dataT.FA <- matrix(n1,nrow=nks,ncol=nds)
  fa.postE <- fa_pE_bma <- re_sampleE <- matrix(0,nrow=nks,ncol=nds)
  fa.postT <- fa_pT_bma <- re_sampleT <- matrix(0,nrow=nks,ncol=nds)
  
  sel <- sel_opt <- matrix(0,nrow=nks,ncol=nds)
  
  
  data.IA_ntrial <- dataE.IA_ntrial <- dataT.IA_ntrial <- array(n1,c(nks,nds,ntrial))
  ia.postE_ntrial  <- ia_pE_bma_ntrial <- array(0,c(nks,nds,ntrial))
  ia.postT_ntrial  <- ia_pT_bma_ntrial <- array(0,c(nks,nds,ntrial))
  
  data.FA_ntrial <- dataE.FA_ntrial <- dataT.FA_ntrial <-array(n1,c(nks,nds,ntrial))
  fa.postE_ntrial <- fa_pE_bma_ntrial <- array(0,c(nks,nds,ntrial))
  fa.postT_ntrial <- fa_pT_bma_ntrial <- array(0,c(nks,nds,ntrial))
  
  sel_ntrial <- sel_opt_ntrial <- array(0,c(nks,nds,ntrial))
  
  ## function: calculate marginal likelihood
  bma_ndlm_marginl <- function(datar, datan, nks, nds, nsampling,model,nmodel, theta0,zeta.mu,zeta.tau,tau.theta1){
    
    marginlik <- rep(0,nmodel) ## marginal likelihood for each model
    
    rmu11 <- rmu21 <- rmu31 <- matrix(0,nsampling,nmodel) ## randomly generated mu for tumor type 1
    rmu12 <- rmu22 <- rmu32 <- matrix(0,nsampling,nmodel) ## randomly generated mu for tumor type 2
    rmu13 <- rmu23 <- rmu33 <- matrix(0,nsampling,nmodel) ## randomly generated mu for tumor type 3
    rmu14 <- rmu24 <- rmu34 <- matrix(0,nsampling,nmodel) ## randomly generated mu for tumor type 4
    
    ##### model
    
    ub <- 10
    
    for(imodel in c(1,6,7,8)){
      rsigma.mu <- cbind(runif(nsampling, 0, zeta.mu),runif(nsampling, 0, zeta.mu))
      
      rtheta.j1 <- cbind(rnorm(nsampling, qnorm(theta0), tau.theta1),rnorm(nsampling, qnorm(theta0), tau.theta1))
      rtau.theta2 <- cbind(runif(nsampling, 0, zeta.tau),runif(nsampling, 0, zeta.tau))
      rtheta.j2 <- cbind(rnorm(nsampling, rtheta.j1[,1], rtau.theta2[,1]),rnorm(nsampling, rtheta.j1[,2], rtau.theta2[,2]))
      rtheta.j3 <- cbind(rnorm(nsampling, rtheta.j2[,1], rtau.theta2[,1]),rnorm(nsampling, rtheta.j2[,2], rtau.theta2[,2]))
      
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1[,model[[imodel]][1]], sd = rsigma.mu[,model[[imodel]][1]])
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2[,model[[imodel]][1]], sd = rsigma.mu[,model[[imodel]][1]])
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3[,model[[imodel]][1]], sd = rsigma.mu[,model[[imodel]][1]])
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1[,model[[imodel]][2]], sd = rsigma.mu[,model[[imodel]][2]])
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2[,model[[imodel]][2]], sd = rsigma.mu[,model[[imodel]][2]])
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3[,model[[imodel]][2]], sd = rsigma.mu[,model[[imodel]][2]])
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1[,model[[imodel]][3]], sd = rsigma.mu[,model[[imodel]][3]])
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2[,model[[imodel]][3]], sd = rsigma.mu[,model[[imodel]][3]])
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3[,model[[imodel]][3]], sd = rsigma.mu[,model[[imodel]][3]])
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1[,model[[imodel]][4]], sd = rsigma.mu[,model[[imodel]][4]])
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2[,model[[imodel]][4]], sd = rsigma.mu[,model[[imodel]][4]])
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3[,model[[imodel]][4]], sd = rsigma.mu[,model[[imodel]][4]])
    }
    
    for(imodel in 2){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu11[,imodel], sd = rtau.theta2)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu21[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 3){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu12[,imodel], sd = rtau.theta2)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu22[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 4){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu13[,imodel], sd = rtau.theta2)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu23[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 5){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu14[,imodel], sd = rtau.theta2)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu24[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 9){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu13[,imodel], sd = rtau.theta2)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu23[,imodel], sd = rtau.theta2)
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu14[,imodel], sd = rtau.theta2)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu24[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 10){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu12[,imodel], sd = rtau.theta2)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu22[,imodel], sd = rtau.theta2)
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu14[,imodel], sd = rtau.theta2)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu24[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 11){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu12[,imodel], sd = rtau.theta2)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu22[,imodel], sd = rtau.theta2)
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu13[,imodel], sd = rtau.theta2)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu23[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 12){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu11[,imodel], sd = rtau.theta2)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu21[,imodel], sd = rtau.theta2)
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu14[,imodel], sd = rtau.theta2)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu24[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 13){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu11[,imodel], sd = rtau.theta2)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu21[,imodel], sd = rtau.theta2)
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu13[,imodel], sd = rtau.theta2)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu23[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 14){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu11[,imodel], sd = rtau.theta2)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu21[,imodel], sd = rtau.theta2)
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu12[,imodel], sd = rtau.theta2)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu22[,imodel], sd = rtau.theta2)
      
      rsigma.mu <- runif(nsampling, 0, zeta.mu)
      rtheta.j1 <- rnorm(nsampling, qnorm(theta0))
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rtheta.j2 <- rnorm(nsampling, rtheta.j1, rtau.theta2)
      rtheta.j3 <- rnorm(nsampling, rtheta.j2, rtau.theta2)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j1, sd = rsigma.mu)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j2, sd = rsigma.mu)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rtheta.j3, sd = rsigma.mu)
    }
    
    for(imodel in 15){
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu11[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu21[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu11[,imodel], sd = rtau.theta2)
      rmu31[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu21[,imodel], sd = rtau.theta2)
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu12[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu22[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu12[,imodel], sd = rtau.theta2)
      rmu32[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu22[,imodel], sd = rtau.theta2)
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu13[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu23[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu13[,imodel], sd = rtau.theta2)
      rmu33[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu23[,imodel], sd = rtau.theta2)
      
      rtau.theta2 <- runif(nsampling, 0, zeta.tau)
      rmu14[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = qnorm(theta0), sd = tau.theta1)
      rmu24[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu14[,imodel], sd = rtau.theta2)
      rmu34[,imodel] <- rtruncnorm(nsampling, a=-ub, b=ub, mean = rmu24[,imodel], sd = rtau.theta2)
    }
    
    
    
    loglik  <-datar[1,1]*log(pnorm(rmu11)) + (datan[1,1]-datar[1,1])*log(1-pnorm(rmu11)) +
      datar[1,2]*log(pnorm(rmu21)) + (datan[1,2]-datar[1,2])*log(1-pnorm(rmu21)) +
      datar[1,3]*log(pnorm(rmu31)) + (datan[1,3]-datar[1,3])*log(1-pnorm(rmu31)) + 
      datar[2,1]*log(pnorm(rmu12)) + (datan[2,1]-datar[2,1])*log(1-pnorm(rmu12)) +
      datar[2,2]*log(pnorm(rmu22)) + (datan[2,2]-datar[2,2])*log(1-pnorm(rmu22)) +
      datar[2,3]*log(pnorm(rmu32)) + (datan[2,3]-datar[2,3])*log(1-pnorm(rmu32)) + 
      datar[3,1]*log(pnorm(rmu13)) + (datan[3,1]-datar[3,1])*log(1-pnorm(rmu13)) +
      datar[3,2]*log(pnorm(rmu23)) + (datan[3,2]-datar[3,2])*log(1-pnorm(rmu23)) +
      datar[3,3]*log(pnorm(rmu33)) + (datan[3,3]-datar[3,3])*log(1-pnorm(rmu33)) + 
      datar[4,1]*log(pnorm(rmu14)) + (datan[4,1]-datar[4,1])*log(1-pnorm(rmu14)) +
      datar[4,2]*log(pnorm(rmu24)) + (datan[4,2]-datar[4,2])*log(1-pnorm(rmu24)) +
      datar[4,3]*log(pnorm(rmu34)) + (datan[4,3]-datar[4,3])*log(1-pnorm(rmu34))
    
    for(imodel in 1:nmodel){ marginlik[imodel] <- mean(exp(loglik[which(is.finite(loglik[,imodel])),imodel])) }
    
    return(marginlik/sum(marginlik))
  }
  
  ## function: Bayesian analysis for IA or FA
  ## time_analysis = 'interim' or 'final'
  bma_analysis <- function(datar,datan,theta0,time_analysis){
    
    post_model <- rep(0,nmodel)
    
    post_pr <- post_p <- array(0,c(nmodel,nks,nds))
    post_p_bma <- post_pr_bma <- matrix(0,nrow=nks,ncol=nds)
    
    post_samples_c <- post_samples <- list()
    
    for(imodel in 1:nmodel){
      post_samples[[imodel]] <- list()
      for(id in 1:nds){ post_samples[[imodel]][[id]] <- matrix(0,nrow=nsampling-nburnin,ncol=nks) }
    }
    
    ## model for 1 tumor type
    model_str_single_subgroup <- " model {
          for(id in 1:nds){
            y[id] ~ dbin(p[id],n[id])
            probit(p[id]) <- mu[id]
          }

          mu[1] ~ dnorm(theta,1/tau.theta1/tau.theta1)
          mu[2] ~ dnorm(mu[1],1/tau.theta2/tau.theta2)
          mu[3] ~ dnorm(mu[2],1/tau.theta2/tau.theta2)
          
          tau.theta2 ~ dunif(0,zeta.tau)
          
        }
        "
    ## model for 2 tumor types
    model_str2 <-  "model{
              for(id in 1:nds){
                for(ik in 1:nks){
                  y[ik,id] ~ dbin(p[ik,id],n[ik,id])
                  probit(p[ik,id]) <- mu[ik,id]
                  mu[ik,id] ~ dnorm(theta.j[id],1/sigma.mu/sigma.mu)
                }
              }
              sigma.mu ~ dunif(0,zeta.mu)
    
              theta.j[1] ~ dnorm(theta,1/tau.theta1/tau.theta1)
              theta.j[2] ~ dnorm(theta.j[1],1/tau.theta2/tau.theta2)
              theta.j[3] ~ dnorm(theta.j[2],1/tau.theta2/tau.theta2)
              
              tau.theta2 ~ dunif(0,zeta.tau)
              
        }"
    
    ## model for 3 tumor types
    model_str3 <-  "model{
              for(id in 1:nds){
                for(ik in 1:nks){
                  y[ik,id] ~ dbin(p[ik,id],n[ik,id])
                  probit(p[ik,id]) <- mu[ik,id]
                  mu[ik,id] ~ dnorm(theta.j[id],1/sigma.mu/sigma.mu)
                }
              }
              sigma.mu ~ dunif(0,zeta.mu)
    
              theta.j[1] ~ dnorm(theta,1/tau.theta1/tau.theta1)
              theta.j[2] ~ dnorm(theta.j[1],1/tau.theta2/tau.theta2)
              theta.j[3] ~ dnorm(theta.j[2],1/tau.theta2/tau.theta2)
              
              tau.theta2 ~ dunif(0,zeta.tau)
              
        }"
    
    ## model for 4 tumor types
    model_str4 <-  "model{
              for(id in 1:nds){
                for(ik in 1:nks){
                  y[ik,id] ~ dbin(p[ik,id],n[ik,id])
                  probit(p[ik,id]) <- mu[ik,id]
                  mu[ik,id] ~ dnorm(theta.j[id],1/sigma.mu/sigma.mu)
                }
              }
              sigma.mu ~ dunif(0,zeta.mu)
    
              theta.j[1] ~ dnorm(theta,1/tau.theta1/tau.theta1)
              theta.j[2] ~ dnorm(theta.j[1],1/tau.theta2/tau.theta2)
              theta.j[3] ~ dnorm(theta.j[2],1/tau.theta2/tau.theta2)
              
              tau.theta2 ~ dunif(0,zeta.tau)
              
        }"
    
    
    ##### 
    for(icluster in 1:11){
      iks <- model_cluster[[icluster]]
      n_iks <- length(iks)
      data_jags <- list(y=datar[iks,],n=datan[iks,],nks=n_iks,nds=nds,
                        theta=qnorm(theta0), zeta.mu=zeta.mu, zeta.tau=zeta.tau,tau.theta1=tau.theta1)
      model_str <- model_str2
      if(length(iks)==3){model_str <- model_str3}
      if(length(iks)==4){model_str <- model_str4}
      fit <- jags(data=data_jags, parameters.to.save = c('p'), model.file = textConnection(model_str),
                  n.chains = 1, n.iter = nsampling,n.burnin = nburnin, n.thin = 1)
      
      post_samples_c[[icluster]] <- list()
      post_samples_c[[icluster]][[1]] <- fit$BUGSoutput$sims.list$p[,,1]
      post_samples_c[[icluster]][[2]] <- fit$BUGSoutput$sims.list$p[,,2]
      post_samples_c[[icluster]][[3]] <- fit$BUGSoutput$sims.list$p[,,3]
    };
    
    for(icluster in 12:15){
      iks <- model_cluster[[icluster]]
      data_jags <- list(y=datar[iks,],n=datan[iks,],nds=nds, theta=qnorm(theta0),zeta.tau=zeta.tau,tau.theta1=tau.theta1)
      fit <- jags(data=data_jags, parameters.to.save = c('p'), model.file = textConnection(model_str_single_subgroup),
                  n.chains = 1, n.iter = nsampling,n.burnin = nburnin, n.thin = 1)
      
      post_samples_c[[icluster]] <- list()
      post_samples_c[[icluster]][[1]] <- as.matrix(fit$BUGSoutput$sims.list$p[,1])
      post_samples_c[[icluster]][[2]] <- as.matrix(fit$BUGSoutput$sims.list$p[,2])
      post_samples_c[[icluster]][[3]] <- as.matrix(fit$BUGSoutput$sims.list$p[,3])
    };
    
    ## calculate posterior weight
    post_model <- bma_ndlm_marginl(datar, datan, nks, nds, nsampling, model, nmodel, theta0, zeta.mu, zeta.tau, tau.theta1)
    
    ## merge the posterior samples
    for(id in 1:nds){
      
      post_samples[[1]][[id]] <- post_samples_c[[1]][[id]]
      post_samples[[2]][[id]][,which(model[[2]]==1)] <- post_samples_c[[12]][[id]];post_samples[[2]][[id]][,which(model[[2]]==2)] <- post_samples_c[[5]][[id]]
      post_samples[[3]][[id]][,which(model[[3]]==1)] <- post_samples_c[[13]][[id]];post_samples[[3]][[id]][,which(model[[3]]==2)] <- post_samples_c[[4]][[id]]
      post_samples[[4]][[id]][,which(model[[4]]==1)] <- post_samples_c[[14]][[id]];post_samples[[4]][[id]][,which(model[[4]]==2)] <- post_samples_c[[3]][[id]]
      post_samples[[5]][[id]][,which(model[[5]]==1)] <- post_samples_c[[15]][[id]];post_samples[[5]][[id]][,which(model[[5]]==2)] <- post_samples_c[[2]][[id]]
      post_samples[[6]][[id]][,which(model[[6]]==1)] <- post_samples_c[[6]][[id]];post_samples[[6]][[id]][,which(model[[6]]==2)] <- post_samples_c[[11]][[id]]
      post_samples[[7]][[id]][,which(model[[7]]==1)] <- post_samples_c[[7]][[id]];post_samples[[7]][[id]][,which(model[[7]]==2)] <- post_samples_c[[10]][[id]]
      post_samples[[8]][[id]][,which(model[[8]]==1)] <- post_samples_c[[8]][[id]];post_samples[[8]][[id]][,which(model[[8]]==2)] <- post_samples_c[[9]][[id]]
      post_samples[[9]][[id]][,which(model[[9]]==1)] <- post_samples_c[[6]][[id]];post_samples[[9]][[id]][,which(model[[9]]==2)] <- post_samples_c[[14]][[id]];
      post_samples[[9]][[id]][,which(model[[9]]==3)] <- post_samples_c[[15]][[id]]
      post_samples[[10]][[id]][,which(model[[10]]==1)] <- post_samples_c[[7]][[id]];post_samples[[10]][[id]][,which(model[[10]]==2)] <- post_samples_c[[13]][[id]];
      post_samples[[10]][[id]][,which(model[[10]]==3)] <- post_samples_c[[15]][[id]]
      post_samples[[11]][[id]][,which(model[[11]]==1)] <- post_samples_c[[8]][[id]];post_samples[[11]][[id]][,which(model[[11]]==2)] <- post_samples_c[[13]][[id]];
      post_samples[[11]][[id]][,which(model[[11]]==3)] <- post_samples_c[[14]][[id]]
      post_samples[[12]][[id]][,which(model[[12]]==1)] <- post_samples_c[[9]][[id]];post_samples[[12]][[id]][,which(model[[12]]==2)] <- post_samples_c[[12]][[id]];
      post_samples[[12]][[id]][,which(model[[12]]==3)] <- post_samples_c[[15]][[id]]
      post_samples[[13]][[id]][,which(model[[13]]==1)] <- post_samples_c[[10]][[id]];post_samples[[13]][[id]][,which(model[[13]]==2)] <- post_samples_c[[12]][[id]];
      post_samples[[13]][[id]][,which(model[[13]]==3)] <- post_samples_c[[14]][[id]]
      post_samples[[14]][[id]][,which(model[[14]]==1)] <- post_samples_c[[11]][[id]];post_samples[[14]][[id]][,which(model[[14]]==2)] <- post_samples_c[[12]][[id]];
      post_samples[[14]][[id]][,which(model[[14]]==3)] <- post_samples_c[[13]][[id]]
      post_samples[[15]][[id]][,which(model[[15]]==1)] <- post_samples_c[[12]][[id]];post_samples[[15]][[id]][,which(model[[15]]==2)] <- post_samples_c[[13]][[id]]
      post_samples[[15]][[id]][,which(model[[15]]==3)] <- post_samples_c[[14]][[id]];post_samples[[15]][[id]][,which(model[[15]]==4)] <- post_samples_c[[15]][[id]]
      for(imodel in 1:nmodel){
        post_pr[imodel,,id] <- colMeans(post_samples[[imodel]][[id]]>theta0)  ## posterior probability using model i
        post_p[imodel,,id] <- colMeans(post_samples[[imodel]][[id]]) ## posterior mean estimate using model i
      }
      post_pr_bma[,id] <- colSums(matrix(rep(post_model,nks),ncol=nks)*post_pr[,,id]) ## averaged posterior probability
      post_p_bma[,id] <- colSums(matrix(rep(post_model,nks),ncol=nks)*post_p[,,id])  ## averaged posterior mean estimate
    }
    
    
    
    if(is.na(post_pr_bma[1])){ return(0) }else{
      if(time_analysis == 'interim'){
        return(list(post_pr_bma=post_pr_bma,post_p_bma=post_p_bma,
                    post_model=post_model, post_pr=post_pr,post_p=post_p))
      }
      
      if(time_analysis == 'final'){
        
        weighted_sample <- array(0, c(nsampling-nburnin,nds,nks))
        re_sample <- array(0, c(nsampling-nburnin,nds,nks))
        for(ik in 1:nks){
          for(id in  1:nds){
            for(imodel in 1:nmodel){ weighted_sample[,id,ik] <- weighted_sample[,id,ik] + post_samples[[imodel]][[id]][,ik] * post_model[imodel] }
            
            n_re_sample <- as.vector(rmultinom(1, size = nsampling-nburnin, prob = post_model))
            re_sample_jk <- NULL
            for(imodel in 1:nmodel){
              re_sample_jk <- c(re_sample_jk, post_samples[[imodel]][[id]][sample(1:(nsampling-nburnin), n_re_sample[imodel]),ik])
            }
            re_sample[,id,ik] <- re_sample_jk
            
          }
        }
        
        return(list(post_pr_bma=post_pr_bma,post_p_bma=post_p_bma,
                    post_model=post_model, post_pr=post_pr,post_p=post_p,
                    weighted_sample=weighted_sample,re_sample=re_sample))
      }
    }
  }
  
  
  ## function: decision making
  ## selecting both optimal and promising doses
  decision_sel <- function(cutoff_utility, ik, ia.postE, ia.postT, fa.postE, fa.postT, uufa){
    
    sel_ik <- rep(0,nds)
    
    setII <- which(ia.postE[ik,]>cutoff_interim & ia.postT[ik,]>cutoff_interimT & fa.postE[ik,]>cutoff_final & fa.postT[ik,]>cutoff_finalT)
    
    if(length(setII)==1){
      
      sel[ik,setII] <- 1
      
    }else if (length(setII)==2){
      
      ## identify the optimal dose
      ## pr(utility[j]=max(utility)|data)
      max_count <- apply(uufa[,setII,ik],1,which.is.max)
      prob_max <- c(sum(max_count==1),sum(max_count==2))
      max_eff_dose <- which.is.max(prob_max)
      
      sel_ik[setII[max_eff_dose]] <- 1
      
      ## identify the promising doses
      ## pr(max(utility[j]) - utility_optimal < delta1|data)
      candi_dose <- which.is.max(-prob_max)
      if(mean(uufa[,setII[max_eff_dose],ik] - uufa[,setII[candi_dose],ik]<delta1)>cutoff_utility){
        sel_ik[setII[candi_dose]] <- 1 }
      
    }else if(length(setII)==3){
      
      max_count <- apply(uufa[,,ik],1,which.is.max)
      prob_max <- c(sum(max_count==1),sum(max_count==2),sum(max_count==3))
      max_eff_dose <- which.is.max(prob_max)
      
      sel_ik[max_eff_dose] <- 1
      
      candi_dose <- setII[-max_eff_dose]
      for(ican in candi_dose){
        if(mean(uufa[,max_eff_dose,ik] - uufa[,ican,ik]<delta1)>cutoff_utility){
          sel_ik[ican] <- sel_ik[ican] + 1 }
      }
      
    }
    
    return(sel_ik)
    
  }
  
  
  ## function: decision making
  ## selecting optimal
  decision_sel_opt <- function(ik, ia.postE, ia.postT, fa.postE, fa.postT, uufa){
    
    sel_opt_ik <- rep(0,nds)
    
    setII <- which(ia.postE[ik,]>cutoff_interim & ia.postT[ik,]>cutoff_interimT & fa.postE[ik,]>cutoff_final & fa.postT[ik,]>cutoff_finalT)
    
    if(length(setII)==1){
      
      sel_opt_ik[setII] <- 1
      
    }else if (length(setII)==2){
      
      max_count <- apply(uufa[,setII,ik],1,which.is.max)
      prob_max <- c(sum(max_count==1),sum(max_count==2))
      max_eff_dose <- which.is.max(prob_max)
      
      sel_opt_ik[setII[max_eff_dose]] <- 1
      
    }else if(length(setII)==3){
      
      max_count <- apply(uufa[,,ik],1,which.is.max)
      prob_max <- c(sum(max_count==1),sum(max_count==2),sum(max_count==3))
      max_eff_dose <- which.is.max(prob_max)
      
      sel_opt_ik[max_eff_dose] <- 1
      
    }
    
    return(sel_opt_ik)
    
  }
  
  for(itrial in 1:ntrial){
    
    ##### data generation and analysis #####
    
    ## interim data generation
    for(ik in 1:nks){
      for(id in 1:nds){
        dataE.IA[ik,id] <- sum(rbinom(1,n1,peff[ik,id]))
        dataT.IA[ik,id] <- sum(rbinom(1,n1,ptox[ik,id]))} }
    
    ## interim analysis for efficacy
    result_iaE <- bma_analysis(datar=dataE.IA,datan=data.IA,theta0=theta0, time_analysis = 'interim')
    if(!is.list(result_iaE)){ return(0) }
    ia.postE <- result_iaE$post_pr_bma # posterior probability
    ia_pE_bma <- result_iaE$post_p_bma # posterior mean estimate of efficacy
    post_model_iaE <- result_iaE$post_model # posterior model weight
    
    ## interim analysis for toxicity
    result_iaT <- bma_analysis(datar=dataT.IA,datan=data.IA,theta0=theta0T, time_analysis = 'interim')
    if(!is.list(result_iaT)){ return(0) }
    ia.postT <- 1-result_iaT$post_pr_bma  # posterior probability
    ia_pT_bma <- result_iaT$post_p_bma  # posterior mean estimate of toxicity
    post_model_iaT <- result_iaT$post_model # posterior model weight
    
    # stage 2 data generation
    for(ik in 1:nks){
      for(id in 1:nds){
        if(ia.postE[ik,id]<cutoff_interim | ia.postT[ik,id]<cutoff_interimT){ # ineffective
          dataE.FA[ik,id] <- dataE.IA[ik,id]
          dataT.FA[ik,id] <- dataT.IA[ik,id]
        }else{  ## continue the trial for acceptable dose arms
          dataE.FA[ik,id] <- dataE.IA[ik,id] + sum(rbinom(nsample - n1,1,peff[ik,id]))
          dataT.FA[ik,id] <- dataT.IA[ik,id] + sum(rbinom(nsample - n1,1,ptox[ik,id]))
          data.FA[ik,id] <- nsample }
      }
    }
    
    
    ## final analysis and decision making
    if(sum(data.FA==nsample)>0){
      result_faE <- bma_analysis(datar=dataE.FA,datan=data.FA,theta0=theta0, time_analysis = 'final')
      if(!is.list(result_faE)){ return(0) }
      fa.postE <- result_faE$post_pr_bma # posterior probability
      fa_pE_bma <- result_faE$post_p_bma # posterior mean estimate of efficacy
      post_model_faE <- result_faE$post_model  # posterior model weight
      re_sampleE <- result_faE$re_sample
      
      result_faT <- bma_analysis(datar=dataT.FA,datan=data.FA,theta0=theta0T, time_analysis = 'final')
      if(!is.list(result_faT)){ return(0) }
      fa.postT <- 1-result_faT$post_pr_bma # posterior probability
      fa_pT_bma <- result_faT$post_p_bma # posterior mean estimate of toxicity
      post_model_faT <- result_faT$post_model  # posterior model weight
      re_sampleT <- result_faT$re_sample
      
      ## posterior utility
      uufa <- array(0,c(dim(re_sampleT)[1],nds,nks))
      
      for(ik in 1:nks){
        
        uufa[,,ik] <- re_sampleE[,,ik] - weightu*re_sampleT[,,ik]
        
        ik_decision_opt <- decision_sel_opt(ik, ia.postE, ia.postT, fa.postE, fa.postT, uufa)
        sel_opt[ik,] <- ik_decision_opt
        ik_decision <- decision_sel(cutoff_utility=cutoff_utility, ik, ia.postE, ia.postT, fa.postE, fa.postT, uufa)
        sel[ik,] <- ik_decision
      }
    }
    
    data.IA_ntrial[,,itrial] <- data.IA
    dataE.IA_ntrial[,,itrial] <- dataE.IA
    dataT.IA_ntrial[,,itrial] <- dataT.IA
    
    ia.postE_ntrial[,,itrial]  <- ia.postE
    ia_pE_bma_ntrial[,,itrial] <- ia_pE_bma
    ia.postT_ntrial[,,itrial]  <- ia.postT
    ia_pT_bma_ntrial[,,itrial] <- ia_pT_bma
    
    if(sum(data.FA==nsample)==0){    
      data.FA_ntrial[,,itrial] <- data.IA
      dataE.FA_ntrial[,,itrial] <- dataE.IA
      dataT.FA_ntrial[,,itrial] <- dataT.IA
      
      fa.postE_ntrial[,,itrial] <- ia.postE
      fa_pE_bma_ntrial[,,itrial] <- ia_pE_bma
      fa.postT_ntrial[,,itrial] <- ia.postT
      fa_pT_bma_ntrial[,,itrial] <- ia_pT_bma 
    }else{
      data.FA_ntrial[,,itrial] <- data.FA
      dataE.FA_ntrial[,,itrial] <- dataE.FA
      dataT.FA_ntrial[,,itrial] <- dataT.FA
      
      fa.postE_ntrial[,,itrial] <- fa.postE
      fa_pE_bma_ntrial[,,itrial] <- fa_pE_bma
      fa.postT_ntrial[,,itrial] <- fa.postT
      fa_pT_bma_ntrial[,,itrial] <- fa_pT_bma
    }
    
    
    sel_ntrial[,,itrial] <- sel
    sel_opt_ntrial[,,itrial] <- sel_opt
    
  }

  
  
  interim_data <- round(cbind(apply(dataE.IA_ntrial,c(1,2),mean), apply(dataT.IA_ntrial,c(1,2),mean), apply(data.IA_ntrial,c(1,2),mean)),2)
  rownames(interim_data) <- c('Tumor type 1', 'Tumor type 2', 'Tumor type 3', 'Tumor type 4')
  colnames(interim_data) <- c('Eff-dose1', 'Eff-dose2', 'Tox-dose1', 'Tox-dose2', 'N-dose1', 'N-dose2')
  
  interim_estmate <- round(cbind(apply(ia_pE_bma_ntrial,c(1,2),mean), apply(ia_pT_bma_ntrial,c(1,2),mean)),2)
  rownames(interim_estmate) <- c('Tumor type 1', 'Tumor type 2', 'Tumor type 3', 'Tumor type 4')
  colnames(interim_estmate) <- c('Eff-dose1', 'Eff-dose2', 'Tox-dose1', 'Tox-dose2')
  
  final_data <- cbind(apply(dataE.FA_ntrial,c(1,2),mean),apply(dataT.FA_ntrial,c(1,2),mean),apply(data.FA_ntrial,c(1,2),mean))
  rownames(final_data) <- c('Tumor type 1', 'Tumor type 2', 'Tumor type 3', 'Tumor type 4')
  colnames(final_data) <- c('Eff-dose1', 'Eff-dose2', 'Tox-dose1', 'Tox-dose2', 'N-dose1', 'N-dose2')
  
  final_estmate <- round(cbind(apply(fa_pE_bma_ntrial,c(1,2),mean), apply(fa_pT_bma_ntrial,c(1,2),mean)),2)
  rownames(final_estmate) <- c('Tumor type 1', 'Tumor type 2', 'Tumor type 3', 'Tumor type 4')
  colnames(final_estmate) <- c('Eff-dose1', 'Eff-dose2', 'Tox-dose1', 'Tox-dose2')
  
  
  sel_opt_sum <- round(apply(sel_opt_ntrial,c(1,2),mean),2)
  rownames(sel_opt_sum) <- c('Tumor type 1', 'Tumor type 2', 'Tumor type 3', 'Tumor type 4')
  colnames(sel_opt_sum) <- c('Dose1', 'Dose2')
  
  sel_sum <- round(apply(sel_ntrial,c(1,2),mean),2)
  rownames(sel_sum) <- c('Tumor type 1', 'Tumor type 2', 'Tumor type 3', 'Tumor type 4')
  colnames(sel_sum) <- c('Dose1', 'Dose2')
  
  
  re.list <- list(interim_data = interim_data, interim_estmate = interim_estmate, 
                  final_data = final_data, final_estmate = final_estmate, 
                  optimal_selection = sel_opt_sum, 
                  promising_selection = sel_sum)
  
  return(re.list)
}

## example
# ntrial <- 2
# peff <- rbind(c(0.20, 0.20, 0.20), c(0.20, 0.20, 0.20), c(0.20, 0.20, 0.20), c(0.20, 0.20, 0.20))
# ptox <- rbind(c(0.40, 0.40, 0.40), c(0.40, 0.40, 0.40), c(0.40, 0.40, 0.40), c(0.40, 0.40, 0.40))
# n1 <- 10
# nsample <- 20
# theta0 <- 0.2
# theta0T <- 0.4
# weightu <- 0.25
# delta1 <- 0.05
# cutoff_utility <- 0.4
# 
# batodo_3dose(ntrial, peff, ptox, n1, nsample, theta0, theta0T, weightu, delta1, cutoff_utility)

