
### ::: All for simulation ::: ###
# #- DON'T RUN EXCEPT CHTC #################### CHTC - Starts ##################

args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
}

arguments = commandArgs(trailingOnly=TRUE)
n_iter = as.numeric(arguments[[1]])
n_rep = as.numeric(arguments[[2]])
ni = as.numeric(arguments[[3]])
nj = as.numeric(arguments[[4]])

# ##################### CHTC - Ends ########################### - DON'T RUN ENDS




###### ====== All for 2 slopes model ====== ######

# Library
{
  library(lme4)
  library(brms)
  library(tidyverse)
  options(mc.cores = 8)
}

# FOR arguments
# args = commandArgs(trailingOnly = T)




### :::::: 2 slopes :::::: ###

### :::::: Data Generation :::::: ###

dgf <- function(ni, nj, gamma00, gamma01, gamma02, gamma03, gamma04, gamma05, gamma06, gamma07, u_0, u_1, w_0, w_1){
  # RVs set up
  x = matrix(99, nrow = ni, ncol = 4) 
  z = matrix(99, nrow = nj, ncol = 3) 
  y = matrix(99, nrow = ni*nj, ncol = 1) 
  r = matrix(99, nrow = ni*nj, ncol = 1) 
  sim = matrix(99, nrow = ni*nj, ncol = 15)
  U = matrix(99, nrow = ni, ncol = 3) # 1 for random intercept, more for random slope
  W = matrix(99, nrow = nj, ncol = 3) # 1 for random intercept, more for random slope
  
  # Random effects set up
  # Correlation matrix for subject
  corU <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow =3, ncol =3,byrow = T)
  sdU <- matrix(c(u_0, 0, 0, 0, u_1, 0, 0, 0, u_2), nrow =3, ncol =3,byrow = T)
  covU <- sdU%*%corU%*%sdU
  I <- covU
  
  # Correlation matrix for item
  corW <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow =3, ncol =3,byrow = T)
  sdW <- matrix(c(w_0, 0, 0, 0, w_1, 0, 0, 0, w_2), nrow =3, ncol =3,byrow = T)
  covW <- sdW%*%corW%*%sdW
  J <- covW
  
  # For subject
  for (i in 1:ni){
    # Random effects
    U[i,1] <- I[1,1]*rnorm(1, 0, 1)
    U[i,2] <- I[1,2]*U[i,1] + I[2,2]*rnorm(1) 
    U[i,3] <- I[1,3]*U[i,1] + I[1,2]*U[i,2] + I[3,3]*rnorm(1)
    
    
    # RVs
    x[i,1] <- rnorm(1,0,2)
    x[i,2] <- rnorm(1,0,3)
    x[i,3] <- rnorm(1,0,1.25)
    x[i,4] <- rbinom(1, 1, .1)
    
    # For stimulus
    for (j in 1:nj) {
      # Random effects
      W[j,1]<- J[1,1]*rnorm(1)
      W[j,2]<- J[1,2]*W[j,1] + J[2,2]*rnorm(1)
      W[j,3] <- J[1,3]*W[j,1] + J[1,2]*W[j,2] + J[3,3]*rnorm(1)
      
      # RVs
      z[j,1] <- rnorm(1,0,1.75)
      z[j,2] <- rnorm(1,0,2)
      z[j,3] <- rnorm(1,0,2.25)
    }
  }
  ind <- 1 # for index
  for (i in 1:ni){
    for(j in 1:nj){
      r[ind, 1] <- rnorm(1,0,1)
      y[ind, 1] = gamma00 + gamma01*x[i,1] + gamma02*z[j,1] + gamma03*x[i,2] + gamma04*z[j,2] + gamma05*x[i,3] + gamma06*z[j,3] + gamma07*x[i,4] + U[i,1] + U[i,2]*x[i,1] + W[j,1] + W[j,2]*z[j,1] + r[ind, 1]
      tmp <- c(i,j,y[ind, 1],x[i,1],z[j,1],x[i,2],z[j,2],x[i,3],z[j,3],x[i,4],U[i,1],U[i,2],W[j,1],W[j,2],r[ind, 1])
      sim[ind, ] <- tmp
      colnames(sim) <- c("i", "j", "y", "x1", "z1", "x2", "z2", "x3", "z3", "x4", "u0", "u1", "w0", "w1", "r")
      ind <- ind + 1
    }
  }
  
  # write.table(sim, file = paste("dt_", "4slopes_", ni, "_", nj, ".csv", sep = ""))
  
  return(sim)
  
}


### :::::: Estimation :::::: ###
alz <- function(df, ni, nj){
  
  ### :::::: Set up :::::: ###
  {
    lm_ml <- list()
    lm_reml <- list()
    bms <- list()
    converge <- data.frame(ml1 = NA, ml2 = NA, reml1 = NA, reml2 = NA)
  }
  
  ### Specify the priors ###
  pri = c(
    set_prior('normal(1,10)', class = 'Intercept'),
    set_prior('normal(-2,10)', coef = 'x1',class='b'),
    set_prior('normal(2,10)', coef = 'z1',class='b'),
    set_prior('normal(.3,10)', coef = 'x2',class='b'),
    set_prior('normal(.2,10)', coef = 'z2',class='b'),
    set_prior('normal(.5,10)', coef = 'x3',class='b'),
    set_prior('normal(1.5,10)', coef = 'z3',class='b'),
    set_prior('normal(1,10)', coef = 'x4',class='b'),
    set_prior('cauchy(0,5)', coef = 'Intercept', class = 'sd', group = 'i'),
    set_prior('cauchy(0,5)', coef = 'Intercept', class = 'sd', group = 'j'))
  
  
  
  ### For ML ###
  
  # For NM
  lm_ml[[1]] <- lmer(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (1|j), data = df, REML=F, control = lmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 100000)))
  lm_ml[[2]] <- lmer(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (0+x1|i) + (1|j) + (0+z1|j), data = df, REML=F, control=lmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=100000)))
  
  # For BOBYQA
  
  lm_ml[[3]] <- lmer(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (1|j), data = df, REML=F, control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  lm_ml[[4]] <- lmer(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (0+x1|i) + (1|j) + (0+z1|j), data = df, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
  
  ### For REML
  # For NM
  lm_reml[[1]]<- lmer(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (1|j), data = df, REML=T, control = lmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 100000)))
  lm_reml[[2]] <- lmer(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (0+x1|i) + (1|j) + (0+z1|j), data = df, REML=T, control=lmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=100000)))
  
  # For BOBYQA
  lm_reml[[3]] <- lmer(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (1|j), data = df, REML=T, control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  lm_reml[[4]] <- lmer(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (0+x1|i) + (1|j) + (0+z1|j), data = df, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
  
  ### For Bayes model ###
  bms[[1]] <- brm(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (1|j), data = df, prior = pri, iter = 10000, thin = 10, chains = 4)
  bms[[2]] <- brm(y ~  x1 + z1 + x2 + z2 + x3 + z3 + x4 + (1|i) + (0+x1|i) + (1|j) + (0+z1|j), data = df, prior = pri, iter = 10000, thin = 10, chains = 4)
  
  
  ### Summarize the results ###
  # Convergence
  converge[, "ml1"] <- lm_ml[[1]]@optinfo$conv$opt #converged = 0 (if 1, not converged)
  converge[, "ml2"] <- lm_ml[[2]]@optinfo$conv$opt
  converge[, "reml1"] <- lm_reml[[1]]@optinfo$conv$opt #converged = 0 (if 1, not converged)
  converge[, "reml2"] <- lm_reml[[2]]@optinfo$conv$opt
  
  # ML
  ml_f1 <- summary(lm_ml[[1]])$coefficients
  ml_r1 <- summary(lm_ml[[1]])$varcor
  ml_f2 <- summary(lm_ml[[2]])$coefficients
  ml_r2 <- summary(lm_ml[[2]])$varcor
  
  # REML
  reml_f1 <- summary(lm_reml[[1]])$coefficients
  reml_r1 <- summary(lm_reml[[1]])$varcor
  reml_f2 <- summary(lm_reml[[2]])$coefficients
  reml_r2 <- summary(lm_reml[[2]])$varcor
  
  # For brms model
  b1 <- rbind(summary(bms[[1]])$fixed, summary(bms[[1]])$random$i, summary(bms[[1]])$random$j)
  b2 <- rbind(summary(bms[[2]])$fixed, summary(bms[[2]])$random$i, summary(bms[[2]])$random$j)
  out_b <- as.data.frame(rbind(b1, b2))
  ml_f <- as.data.frame(rbind(ml_f1, ml_f2))
  reml_f <- as.data.frame(rbind(reml_f1, reml_f2))
  
  
  out_all <- list("out_b" = out_b, "ml_f" = ml_f, "reml_f" = reml_f, "ml_r1" = ml_r1, "ml_r2" = ml_r2, "converge" = converge, 
                  "lm_ml" = lm_ml, "lm_reml" = lm_reml, "bms" = bms)
  
  # saveRDS(out_all, file = "out_all.rds")
  
  return(out_all)
  
}


### :::::: Performance and replication:::::: ###


perf <- function(out_all, gamma00, gamma01, gamma02, gamma03, gamma04, gamma05, gamma06, gamma07,
                 u_0, u_1, w_0, w_1){
  
  ## Set up
  true_val <- data.frame(true = c(gamma00, gamma01, gamma02, gamma03, gamma04, gamma05, gamma06, gamma07, u_0, w_0,
                                  gamma00, gamma01, gamma02, gamma03, gamma04, gamma05, gamma06, gamma07,
                                  u_0, u_1, w_0, w_1))
  dim2 <- nrow(true_val)
  
  # Parameter recovery
  
  bias <- out_all$out_b$Estimate - true_val
  rhat_each <- out_all$out_b$Rhat
  
  rhat_p1 <- sum(round(out_all$out_b$Rhat[1:10], 2) > 1)/10
  rhat_p2 <- sum(round(out_all$out_b$Rhat[10:dim2], 2) > 1)/(dim2-10)
  
  ci_width <- as.data.frame(out_all$out_b$`u-95% CI` - out_all$out_b$`l-95% CI`)
  
  
  perf_list <- list("bias" = bias, "rhat" = rhat_each, "rhat_p1"  = rhat_p1, "rhat_p2"  = rhat_p2, "ci_width" = ci_width)
  
  # saveRDS(perf_list, file = "perf_list.rds")
  
  return(perf_list)
}


### :::::: Replication :::::: ###


run_sim <- function(n_rep, ni, nj, gamma00, gamma01, gamma02, gamma03, gamma04, gamma05, gamma06, gamma07,
                    u_0, u_1, w_0, w_1){
  
  # Set up
  seed <- sample(1:1000, size = n_rep, replace = F)
  
  df <- list()
  out_all <- list()
  perf<- list()
  
  
  for (h in 1:n_rep){
    
    set.seed(seed[h])
    # setwd("/Users/mingya/Research/CREM_Carolyn_2019/Thesis/PHD_R") # For local laptop
    # setwd("/home/mhuang233/ccrem_carolyn_19") # for CHTC
    
    #newdir <- paste(h, "_4slopes_", "seed", seed[h], "_", ni, "_", nj, sep = "")
    
    #dir.create(newdir)
    #setwd(newdir)
    
    df[[h]] <- as.data.frame(dgf(ni = ni, nj = nj, gamma00 = gamma00, gamma01 = gamma01, gamma02 = gamma02, gamma03 = gamma03,
                                 gamma04 = gamma04, gamma05 = gamma05, gamma06 = gamma06, gamma07 = gamma07,
                                 u_0 = u_0 , u_1 = u_1, w_0 = w_0, w_1 = w_1))
    out_all[[h]] <- alz(df = df[[h]], ni = ni, nj = nj)
    perf[[h]] <- perf(out_all = out_all[[h]], gamma00, gamma01, gamma02, gamma03, gamma04, gamma05, gamma06, gamma07,
                      u_0, u_1, w_0, w_1)
    
  }
  
  # setwd("/Users/mingya/Research/CREM_Carolyn_2019/Thesis/PHD_R") # for local laptop
  # setwd("/home/mhuang233/ccrem_carolyn_19") # for CHTC
  
  #saveRDS(perf, file = "perf_all2_100.rds")
  return(perf)
  
}



### Specify the parameter ###
{
  n_iter = n_iter
  n_rep = n_rep
  ni <- ni
  nj <- nj
  gamma00 <- 1
  gamma01 <- -2
  gamma02 <- 2
  gamma03 <- .3
  gamma04 <- .2
  gamma05 <- .5
  gamma06 <- 1.5
  gamma07 <- 1
  u_0 <- .4
  u_1 <- .2
  u_2 <- .6
  w_0 <- .7
  w_1 <- .4
  w_2 <- .9
}


### RUN ###
out2 <- run_sim(n_rep = n_rep, ni = ni, nj = nj, gamma00 = gamma00, gamma01 = gamma01, gamma02 = gamma02, gamma03 = gamma03,
                    gamma04 = gamma04, gamma05 = gamma05, gamma06 = gamma06, gamma07 = gamma07,
                    u_0 = u_0 , u_1 = u_1, w_0 = w_0, w_1 = w_1)


# Pull out the Results
#rnames <- c("1gamma00", "1gamma01", "1gamma02", "1gamma03", "1gamma04", "1gamma05", "1gamma06", "1gamma07", "1u_0", "1w_0",
#            "2gamma00", "2gamma01", "2gamma02", "2gamma03", "2gamma04", "2gamma05", "2gamma06", "2gamma07",
#            "2u_0", "2u_1", "2w_0", "2w_1")


#bias_est <- as.data.frame(sapply(out_list, "[[", "bias"), row.names = rnames)
#ciw_est <- as.data.frame(sapply(out_list, "[[", "ci_width"), row.names = rnames)
#rhat_each <- as.data.frame(sapply(out_list, "[[", "rhat"), row.names = rnames)
#names(ciw_est) <- NULL


#rhat1_est <- sapply(out_list, "[[", "rhat_p1")
#rhat2_est <- sapply(out_list, "[[", "rhat_p2")
#rhat_est <- as.data.frame(rbind(rhat1_est, rhat2_est))

#saveRDS(bias_est, file = "bias_est2_100.rds")
#saveRDS(ciw_est, file = "ciw_est2_100.rds")
#saveRDS(rhat_est, file = "rhat_est2_100.rds")
#saveRDS(rhat_each, file = "rhat_each2_100.rds")

#saveRDS(out_list, file = paste0(n_iter, "_", n_rep, "2out_", ni, "_", nj, ".rds"), compress = T)


# my.file.rename <- function(from, to) {
#  todir <- dirname(to)
#  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
#  file.rename(from = from,  to = to)
#}

# my.file.rename(from = paste0("2out_", ni, "_", nj, ".rds"),
#               to = paste0("/staging/mhuang233/2out_", ni, "_", nj, ".rds"))


saveRDS(out2, file = paste0("cc", n_iter, "_", n_rep, "_2_", ni, "_", nj, ".rds"), compress = T)

# file.copy(from = paste0("cc", n_iter, "_", n_rep, "_2_", ni, "_", nj, ".rds"),
#          to = paste0("/staging/mhuang233/cc", n_iter, "_", n_rep, "_2_", ni, "_", nj, ".rds"))

# file.remove(paste0("cc", n_iter, "_", n_rep, "_2_", ni, "_", nj, ".rds"))



