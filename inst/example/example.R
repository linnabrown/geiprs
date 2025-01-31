library(devtools)

install_github("junyangq/glmnetPlus")
install_github("chrchang/plink-ng", subdir="/2.0/cindex")
install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
install_github("dajmcdon/sparsegl")


unlink(".RData")
rm(list=ls())
#########################################

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ROCR))
suppressMessages(library(utils))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))
suppressMessages(library(stats))
suppressMessages(library(cindex))
suppressMessages(library(pgenlibr))
suppressMessages(library(glmnet))
suppressMessages(library(magrittr))
suppressMessages(library(gglasso))
suppressMessages(library(sparsegl))
############################################

testing=F
if(testing){
  genotype_file="/SFS/archive/data1/bardstms/shenjud/2023SummerIntern/simulation/simulated_genotype/N_40000/pgxsimQCed.trainingANDvalidation"
  # genotype_file="/SFS/project/comp/BARDS/PGx/PGx_Projects/2023SummerIntern/snpnet_ge/data/genotype/pgxsimQCed.trainingANDvalidation"
  phenotype_file="/SFS/archive/data1/bardstms/shenjud/2023SummerIntern/simulation/simulated_phenotype/dis1001/sim_phenotype_adjusted_forSNPNET_seed_100001.txt"
  env="envir"
  tau=0.9
  phenotype="resid"
  split="split"
  memory=40000
  result_dir="/SFS/archive/data1/bardstms/shenjud/2023SummerIntern/simulation/simulation_results/PRS/snpnet-ge/models/dis1001/seed_100060_tau_0.9"
}else{
  args = commandArgs(trailingOnly = TRUE)
  genotype_file=args[1]
  phenotype_file=args[2]
  env=args[3]
  tau=as.numeric(args[4])
  phenotype=args[5]
  split=args[6]
  memory=as.numeric(args[7])
  result_dir=args[8]
}

source('R/geiprs.R')
source('R/functions.R')


configs <- list(
  verbose=TRUE,
  num.groups.batch = 500,
  nCores = 8,
  verbose=FALSE,
  prevIter = 0,
  memory = memory,
  save = TRUE,
  results.dir = result_dir,
  save.dir="result",
  plink2.path="/SFS/project/comp/BARDS/PGx/PGx_Projects/2023SummerIntern/snpnet_ge/snpnet_new/plink2"
)

if(testing==T){
  genotype.pfile = genotype_file
  phenotype.file = phenotype_file
  phenotype = phenotype 
  env = env 
  tau = tau 
  family = "gaussian"
  covariates = NULL
  nlambda = 100 
  lambda.min.ratio = 1e-2
  p.factor = NULL 
  split.col = split
  configs = configs
  lambda=NULL
  status.col = NULL
  mem = NULL
  alpha=1

}

output<- geiprs(
  genotype.pfile = genotype_file, 
  phenotype.file = phenotype_file, 
  phenotype = phenotype, 
  env = env, 
  tau = tau, 
  family = "gaussian", 
  covariates = NULL,
  nlambda = 100, 
  lambda.min.ratio = 1e-2,
  p.factor = NULL, 
  split.col = split,
  configs = configs)

saveRDS(output, file = file.path(result_dir, "final_output.rds" ))