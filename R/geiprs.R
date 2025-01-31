#' Fit the Lasso/Elastic-Net for Large Phenotype-Genotype and Genotype-Environment interaction Datasets
#'
#' Fit the entire lasso, elastic-net, group lasso, and sparse group lasso solution path using the Batch Screening Iterative Lasso (BASIL) algorithm
#' on large phenotype-genotype and Genotype-Environment interaction datasets.
#'
#' REFERENCE:
#' Le Huang#, Wujuan Zhong#, Song Zhai, and Judong Shen. "GEiPRS: A Fast and Powerful Machine Learning Method 
#' for Polygenic Risk Score Prediction by Leveraging Genotype-Environment Interactions". 
#' bioRxiv (2025): xxx
#'
#' @usage geiprs(genotype.pfile, phenotype.file, phenotype, family = NULL, covariates = NULL, alpha
#'   = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04), lambda = NULL,
#'   split.col = NULL, p.factor = NULL, status.col = NULL, mem = NULL, configs = NULL)
#'
#' @param genotype.pfile the PLINK 2.0 pgen file that contains genotype.
#'                       We assume the existence of genotype.pfile.{pgen,pvar.zst,psam}.
#' @param phenotype.file the path of the file that contains the phenotype values and can be read as
#'                       as a table. There should be FID (family ID) and IID (individual ID) columns
#'                       containing the identifier for each individual, and the phenotype column(s).
#'                       (optional) some covariate columns and a column specifying the
#'                       training/validation split can be included in this file.
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file.
#' @param env the name of the environment. Must be the same as the corresponding column name in
#'                  the phenotype file. It represents group information
#' @param tau the sparse group lasso mixing parameter, where the penalty is defined as 
#'.           \eqn{tau * ||beta||_1 + (1-tau)\sum w_g||beta_g||_2}  If tau = 0, it is group lasso penalty; if
#'            tau =1, it is lasso penalty. Only environment != NULL can consider this parameter
#' @param family the type of the phenotype: "gaussian", "binomial", or "cox". If not provided or NULL,
#'               it will be detected based on the number of levels in the response.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#' @param alpha the elastic-net mixing parameter, where the penalty is defined as
#'              alpha * ||beta||_1 + (1-alpha)/2 * ||beta||_2^2. alpha = 1 corresponds to the lasso penalty,
#'              while alpha = 0 corresponds to the ridge penalty.
#' @param nlambda the number of lambda values - default is 100.
#' @param lambda.min.ratio smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value,
#'                         i.e. the smallest value for which all coefficients are zero. The default
#'                         depends on the sample size nobs relative to the number of actual variables
#'                         nvars (after QC filtering). If nobs > nvars, the default is 0.0001, close to zero.
#'                         If nobs < nvars, the default is 0.01. A very small value of lambda.min.ratio
#'                         will lead to a saturated fit in the nobs < nvars case.
#' @param lambda one can specify the full lambda list on which the lasso/elastic-net will be solved.
#'               Once provided, `lambda` and `lambda.min.ratio` will be ignored. It can be used for refitting
#'               after the optimal parameter is selected by validation.
#' @param split.col the column name in the phenotype file that specifies the membership of individuals to
#'                  the training or the validation set. The individuals marked as "train" and "val" will
#'                  be treated as the training and validation set, respectively. When specified, the
#'                  model performance is evaluated on both the training and the validation sets.
#' @param p.factor a named vector of separate penalty factors applied to each coefficient. This is
#' a number that multiplies lambda to allow different shrinkage. If not provided, default is 1
#' for all variables. Otherwise should be complete and positive for all variables.
#' @param status.col the column name for the status column for Cox proportional hazards model.
#'                   When running the Cox model, the specified column must exist in the phenotype file.
#' @param mem Memory (MB) available for the program. It tells PLINK 2.0 the amount of memory it can
#' harness for the computation. IMPORTANT if using a job scheduler.
#' @param configs a list of other config parameters.
#'                \describe{
#'                 \item{missing.rate}{variants are excluded if the missing rate exceeds this level. Default is 0.1.}
#'                 \item{MAF.thresh}{variants are excluded if the minor allele frequency (MAF) is lower
#'                                than this level. Default is 0.001.}
#'                 \item{nCores}{the number of cores used for computation. You may use the maximum number
#'                            of cores available on the computer. Default is 1, single core.}
#'                 \item{num.snps.batch}{the number of variants added to the strong set in each iteration. Default is 1000.}
#'                 \item{num.groups.batch}{the number of groups added to the strong set in each iteration. Default is 1000. 
#'                  This parameter will be applied to sparse group lasso/ group lasso.}
#'                 \item{niter}{The number of maximum iteration in the algorithm. Note that each iteration
#'                              may be able to find solutions for more than one lambda value. The default is 50}
#'                 \item{prevIter}{if non-zero, it indicates the last successful iteration in the procedure so that
#'                              we can restart from there. niter should be no less than prevIter.}
#'                 \item{save}{a logical value whether to save the intermediate results (e.g. in case of job failure and restart).}
#'                 \item{results.dir}{the path to the directory where meta and intermediate results are saved.}
#'                 \item{meta.dir}{the relative path to the subdirectory used to store the computed
#'                              summary statistics, e.g. mean, missing rate, standard deviation (when `standardization = TRUE`).
#'                              Needed when `save = TRUE`. Default is `"meta.dir/`.}
#'                 \item{save.dir}{the relative path to the subdirectory used to store the intermediate
#'                              results so that we may look into or recover from later.
#'                              Needed when `save = TRUE`. Default is `"results/`.}
#'                 \item{excludeSNP}{character vector containing genotype names to exclude from
#'                                  the analysis}
#'                 \item{nlams.init}{the number of lambdas considered in the first iteration.
#'                              Default 10 is a reasonable number to start with.}
#'                 \item{nlams.delta}{the length of extended lambdas down the sequence when there are few
#'                              left in the current sequence (remember we don't fit all lambdas
#'                              every iteration, only extend when most of the current ones have been completed and validated). Default is 5.}
#'                 \item{glmnet.thresh}{the convergence threshold used in glmnet/glmnetPlus.}
#'                 \item{keep}{one may specify keep file in plink format to focus on a subset of individuals.}
#'                 \item{use.glmnetPlus}{a logical value whether to use glmnet with warm start, if
#'                              the glmnetPlus package is available. Currently only "gaussian" family is supported.}
#'                 \item{early.stopping}{a logical value indicating whether early stopping based on validation metric is desired.}
#'                 \item{stopping.lag}{a parameter for the stopping criterion such that the procedure stops after
#'                              this number of consecutive decreases in the validation metric.}
#'                 \item{verbose}{a logical value indicating if more detailed messages should be printed.}
#'                 \item{KKT.verbose}{a logical value indicating if details on KKT check should be printed.}
#'                 \item{increase.size}{the increase in batch size if the KKT condition fails often in recent iterations.
#'                              Default is half of the batch size.}
#'                 \item{plink2.path}{the user-specified path to plink2 (default: plink2)}
#'                 \item{zstdcat.path}{the user-specified path to zstdcat (default: zstdcat)}
#'                 \item{zcat.path}{the user-specified path to zcat (to read a zcat compressed phenotype file) (default: zcat)}
#'                 \item{rank}{if TRUE, then the smallest lambda indices when each variable enters the model are recorded}
#'                }
#' @return A list containing the solution path, the metric evaluated on training/validation set and others.
#'            \describe{ metric.train = metric.train, metric.val = metric.val, fit.results = fit.results,
#'              full.lams = full.lams, a0 = a0, beta = beta, configs = configs, var.rank=var.rank,
#'              lambda.min.ratio = lambda.min.ratio, stats = stats
#'  
#'             }        
#'
#' @importFrom data.table ':='
#'
#' @export
#' @examples
#' \dontrun{
#' configs <- list(
#'   verbose = TRUE,
#'   num.groups.batch = 500,
#'   nCores = 8,
#'   memory = 40000,
#'   save = TRUE,
#'   results.dir = "results",
#'   save.dir = "intermediate",
#'   plink2.path = "/path/to/plink2"
#' )
#'
#' output <- geiprs(
#'   genotype.pfile = "path/to/genotype",
#'   phenotype.file = "path/to/phenotype",
#'   phenotype = "trait",
#'   env = "environment",
#'   tau = 0.9,
#'   family = "gaussian",
#'   configs = configs
#' )
#' }

#env variable is a string, which indicate the column name of the phenotype file
#we also need to adjust enviroment
geiprs <- function(genotype.pfile, phenotype.file, phenotype, env = NULL, tau = NULL, family = NULL, covariates = NULL,
                   alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04),
                   lambda = NULL, split.col = NULL, p.factor = NULL, status.col = NULL, mem = NULL,
                   configs = NULL) {
  ##############################
  #Do not annotate it; When project is finished, please annotate it and delete the 
  #executable file plink2
  # configs[['plink2.path']] = './plink2'
  ###############################
  
  # Whether it is Sparse Group Lasso (SGL). Note This does not include group lasso 
  isSGL = !is.null(env) && !is.null(tau) && tau>0 && tau < 1 
  #snpnetlogger(isSGL)
  # Whether it is Group Lasso (GL)
  isGL = !is.null(env) && tau==0
  ID <- ALT <- NULL
  validation <- (!is.null(split.col))
  time.start <- Sys.time()
  snpnetLogger('Start snpnet-ge', log.time = time.start)
  snpnetLogger('Preprocessing start..')
  
  ### --- Read genotype IDs --- ###
  ids <- list(); phe <- list()
  # read sample id, format: FID_IID
  
  #PSAM format
  #FID	IID	SEX	PHENO1
  # per0	per0	2	-0.832928
  # per1	per1	2	-1.01316
  # per2	per2	2	1.05199
  ids[['psam']] <- readIDsFromPsam(paste0(genotype.pfile, '.psam'))
  ### --- combine the specified configs with the default values --- ###
  if (!is.null(lambda)) nlambda <- length(lambda)
  configs <- setupConfigs(configs, genotype.pfile, phenotype.file, phenotype, env, covariates, alpha, nlambda, split.col, p.factor, status.col, mem)
  if (configs[['prevIter']] >= configs[['niter']]) stop("prevIter is greater or equal to the total number of iterations.")
  
  ### --- Read phenotype file --- ###
  #phe[['master']] includes phenotype, covariates, and enviroment(if existing), which has the same order as ids[['psam']]
  
  phe[['master']] <- readPheMaster(phenotype.file, ids[['psam']], family, env, covariates, phenotype, status.col, split.col, configs)
  if (configs[['verbose']]) snpnetLogger(sprintf("nrow of phe is %d", nrow(phe[['master']]))) 
  ### --- infer family and update the configs --- ###
  # family is data family, bionomial, gausian, etc.  status.col used for cox model
  if (is.null(family)) family <- inferFamily(phe[['master']], phenotype, status.col)
  configs <- updateConfigsWithFamily(configs, family)
  # print config
  if (configs[['verbose']]) print(configs)
  
  ### --- Check whether to use glmnet or glmnetPlus --- ###
  # fdev: minimum fractional change in deviance for stopping path
  # devmax: minimum fractional change in deviance for stopping path
  if (configs[['use.glmnetPlus']]) {
    glmnet.settings <- glmnetPlus::glmnet.control()
    on.exit(do.call(glmnetPlus::glmnet.control, glmnet.settings))
    glmnetPlus::glmnet.control(fdev = 0, devmax = 1)
  } else {
    glmnet.settings <- glmnet::glmnet.control()
    on.exit(do.call(glmnet::glmnet.control, glmnet.settings))
    glmnet::glmnet.control(fdev = 0, devmax = 1)
  }
  
  ### --- Process phenotypes --- ###
  if (family == "binomial"){
    # The input binary phenotype is coded as 2/1 (case/control)
    # For glmnet, we map this to 1/0 (case/control)
    # The following expression will replace -9 (missing) with -10, but
    # the set of individuals with no-missing values are already computed.
    if (min(phe[['master']][[phenotype]], na.rm = T) >= 1 && max(phe[['master']][[phenotype]], na.rm = T) <= 2) {
      phe[['master']][[phenotype]] <- phe[['master']][[phenotype]] - 1
    }
  }
  
  ### --- Define the set of individual IDs for training (and validation) set(s) --- ###
  #If not indicate, train all of the dataset
  if(configs[['verbose']]) snpnetLogger(sprintf("rowcount of phe[['master']] is %d", nrow(phe[['master']])))
  if(is.null(split.col)){
    splits <- c('train')
    ids[['train']] <- phe[['master']]$ID
  }else{
    splits <- c('train', 'val')
    for(s in splits){
      ids[[s]] <- phe[['master']]$ID[ phe[['master']][[split.col]] == s ]
      if(configs[['verbose']]) snpnetLogger(sprintf("sample size or rowcount of ids[['%s']] is %d", s, length(ids[[s]])))
    }
  }
  
  
  
  ### --- Prepare the feature matrix --- ###
  # s: train, validation
  # features['train'] and features['valid'] insert environment and covariates first
  features <- list()
  for(s in splits){
    phe[[s]] <- phe[['master']][match(ids[[s]], phe[['master']]$ID), ]
    rownames(phe[[s]]) <- phe[[s]]$ID
    if (length(c(env, covariates)) > 0) {
      if(!is.null(env)){
        phe[[s]][, c(env)] <- scale(phe[[s]][, c(env), with=F])
      }
      features[[s]] <- as.matrix(phe[[s]][, c(env, covariates), with = F])
    } else {
      features[[s]] <- NULL
    }
    if(configs[['verbose']]) snpnetLogger(sprintf("The number of individuals in %s set: %d", s, dim(phe[[s]])[1]))
  }
  
  ### --- Prepare the response --- ###
  response <- list() ; status <- list() ; surv <- list() ; pred <- list()
  for(s in splits){
    response[[s]] <- phe[[s]][[phenotype]]
    if (family == "cox") {
      status[[s]] <- phe[[s]][[status.col]]
      surv[[s]] <- survival::Surv(response[[s]], status[[s]])
    }
  }
  
  ### --- Read genotypes --- ###
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ', paste0(genotype.pfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  if(configs[['verbose']]) snpnetLogger(sprintf("Total snps count in pgen file: %d", length(vars)))
  configs[["excludeSNP"]] <- base::intersect(configs[["excludeSNP"]], vars)
  pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, '.pvar.zst')) #PVAR format, an extension of .bim format.
  pgen <- list()
  
  for(s in splits) {
    if(configs[['verbose']]) snpnetLogger(sprintf("ids[['psam']] count: %d, ids[['%s']] count: %d", length(ids[['psam']]), s, length(ids[[s]])))
    pgen[[s]] <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), pvar=pvar, sample_subset=match(ids[[s]], ids[['psam']]))
  }
  pgenlibr::ClosePvar(pvar)
  
  #calculate genotype counts
  
  stats <- computeStats(genotype.pfile, phe[['train']]$ID, configs = configs)
  
  ### --- Keep track of the lambda index at which each variant is first added to the model, if required --- ###
  if (configs[['rank']]){
    var.rank <- rep(configs[['nlambda']]+1, length(vars))
    names(var.rank) <- vars
  } else{
    var.rank = NULL
  }
  
  ### --- End --- ###
  snpnetLoggerTimeDiff("Preprocessing end.", time.start, indent=1)
  
  # get the residual after adjusting covariates. Let covariates combine environment
  # This step only needs glmnet
  # the features[['train']] already included environment and covariates in cox regression
  if (configs[['prevIter']] == 0) {
    #This is the first iteration, which will calculate the residual adjusted for
    # covariates and environment, then calculate the product of residual and X.
    snpnetLogger("Iteration 0")
    if (family == "cox"){
      glmmod <- glmnet::glmnet(as.matrix(features[['train']]), surv[['train']], family="cox", standardize=F, lambda=c(0))
      residual <- computeCoxgrad(stats::predict(glmmod, newx=as.matrix(features[['train']])), response[['train']], status[['train']])
    } else {
      glmmod <- stats::glm(
        stats::as.formula(paste(phenotype, " ~ ", paste(c(1, env, covariates), collapse = " + "))),
        data = phe[['train']], family = family
      )
      residual <- matrix(stats::residuals(glmmod, type = "response"), ncol = 1)
    }
    rownames(residual) <- rownames(phe[['train']]) #sample size id
    colnames(residual) <- c('0') # 0 represents 0th residual, if it is not the initial status, it will be lambda index
    
    if (configs[['verbose']]) snpnetLogger("  Start computing inner product for initialization ...")
    
    time.prod.init.start <- Sys.time()
    #compute the inner product. row is snp, column is lambda and lambda with E.
    prod.full <- computeProduct(residual, genotype.pfile, vars, phe[['train']], env, stats, configs, iter=0) / nrow(phe[['train']])

    
    
    if (isSGL | isGL){
      prod.full<-rearrange_prod(prod.full)
      score <- getScore(prod.full)[, 1] # it will become 1xp vector with names; each column is group name (or snp id)
      # if (!is.null(p.factor)){score <- score/p.factor[names(score)]} # Divide the score by the penalty factor
      prod.full.names <- rownames(prod.full)
      # sanity check
      stopifnot(all(prod.full.names[seq(1,length(prod.full.names), 2)] == names(score)))
      lambda.max <- max(score, na.rm = T) / (1 - tau)
    } else{
      score <- abs(prod.full[, 1])
      if (!is.null(p.factor)){score <- score/p.factor[names(score)]} # Divide the score by the penalty factor
      score <- score / max(alpha, 1e-3)
    }

    
    if (configs[['verbose']]) snpnetLoggerTimeDiff("  End computing inner product for initialization.", time.prod.init.start)
    
    nobs <- nrow(phe[['train']]) # sample size
    nvars <- length(vars)-length(stats[["excludeSNP"]]) # number of variants
    if(isSGL | isGL){
      nvars <- 2 * nvars
    }
    ##########compute Lambdas#################
    if (is.null(lambda)) {
      configs[['lambda.min.ratio']] <- lambda.min.ratio
      full.lams <- computeLambdas(lambda.max, configs[['nlambda']], configs[['lambda.min.ratio']])
    } else {
      configs['lambda.min.ratio'] <- list(NULL)
      full.lams <- lambda
    }
    
    lambda.idx <- 1
    num.lams <- configs[["nlams.init"]] #number of lambdas consider first iteration. Default count is 10.
    features.to.keep <- names(glmmod$coefficients[-1]) # Those features currently are environment and covars.
    
    prev.beta <- NULL # Record previous beta
    num.new.valid <- NULL  # track number of new valid solutions every iteration, to adjust length of current lambda seq or size of additional variables
    
    metric.train <- rep(NA, length(full.lams))
    metric.val <- rep(NA, length(full.lams))
    
    increase.group.size <- FALSE
    fit.results <- list()
    beta <- list()
    a0 <- list()
    prev.max.valid.idx <- 0
    earlyStopNow <- FALSE
  } else {
    time.load.start <- Sys.time()
    snpnetLogger(paste0("Recover iteration ", configs[['prevIter']]))
    current.configs <- configs
    load(file.path(configs[['results.dir']], configs[["save.dir"]], paste0("output_iter_", configs[['prevIter']], ".RData")))
    configs <- current.configs
    chr.to.keep <- setdiff(features.to.keep, c(covariates, env)) 
    for(s in splits){
      if (!is.null(features[[s]])) {
        #features[[s]][, (chr.to.keep) := prepareFeatures(pgen[[s]], vars, chr.to.keep, stats)]
        tmp.mat = prepareFeatures(pgen[[s]], vars, chr.to.keep, stats)
        common_cols <- intersect(colnames(features[[s]]), colnames(tmp.mat))
        tmp.mat <- tmp.mat[, !common_cols, with = FALSE, drop = FALSE]
        features[[s]] = cbind( as.matrix(features[[s]]), as.matrix(tmp.mat))
        rm(tmp.mat)
      } else {
        features[[s]] <- as.matrix(prepareFeatures(pgen[[s]], vars, chr.to.keep, stats))
      }
    }
    prev.max.valid.idx <- max.valid.idx
    snpnetLoggerTimeDiff("Time elapsed on loading back features", time.load.start)
    earlyStopNow <- (validation && checkEarlyStopping(metric.val, max.valid.idx, configs[['prevIter']], configs))
  }
  cat("\n")
  # end of pre-processing
  
  # start screening
  if(! earlyStopNow){
    for (iter in (configs[['prevIter']]+1):configs[['niter']]) {
      time.iter.start <- Sys.time()
      snpnetLogger(paste0("Iteration ", iter), log.time=time.iter.start)
      
      # num.lams means that the number of lambdas in this iteration.
      # each iteration, num.lams will increase by nlams.delta (default is 5)
      num.lams <- min(num.lams + ifelse(lambda.idx >= num.lams-configs[["nlams.delta"]]/2, configs[["nlams.delta"]], 0),
                      configs[['nlambda']])   ## extend lambda list if necessary. Extended by nlams.delta(default=5)
      #max(c(utils::tail(num.new.valid, 3), 1)). num.new.valid: list: track the number of new valid solutions in each iteration.
      #utils::tail(vector, 3) return lat 3 elements of object vector
      #max() used fo the the maximum number from the vector[:-3] and 1.
      num.lams <- min(num.lams, lambda.idx + ifelse(is.null(num.new.valid), Inf, max(c(utils::tail(num.new.valid, 3), 1))))
      if(configs[['verbose']]) snpnetLogger(sprintf("number of lam for this iteration: %d", num.lams))
      ### --- Update the feature matrix --- ###
      ### [SGL/GL] select top M groups###
      if (configs[['verbose']]) snpnetLogger("Start updating feature matrix ...", indent=1)
      time.update.start <- Sys.time()
      ### Union with active set ###
      if (iter > 1) {
        # features.to.keep: active set and covariates
        # features.to.discard: non-active set
        # set features.to.discard as NULL in  order to only keep active set.
        #features.to.keep should have the whole group (G and GxE)

        #select non active set. features.to.discard should have the whole group (G and GxE)
        features.to.discard <- setdiff(colnames(features[['train']]), features.to.keep)
        if (length(features.to.discard) > 0) {
          #for(s in splits) features[[s]][, (features.to.discard) := NULL]
          for(s in splits) {
            features[[s]] = (features[[s]]) [, !colnames( features[[s]] ) %in% features.to.discard, drop = FALSE]
          }
        }

        features.to.keep.nocovar = features.to.keep[!(features.to.keep %in% c(env, covariates))] 
        group.uniq.feat <- unique(gsub('_E', '', features.to.keep.nocovar))
        #judge whether 
        # stopifnot(length(group.uniq.feat) == 1/2 * length(features.to.keep.nocovar))
        if(length(group.uniq.feat) == 1/2 * length(features.to.keep.nocovar)){
          current_variable=list(group.uniq.feat=group.uniq.feat, features.to.keep.nocovar=features.to.keep.nocovar)
          #saveRDS(current_variable, file="/SFS/project/comp/BARDS/PGx/PGx_Projects/2023SummerIntern/snpnet_ge/snpnet_new/analysis_script/checkpoint.rds")
        }
        which.in.model <- which(names(score) %in% group.uniq.feat)
        if(length(which.in.model) > 0){
          score[which.in.model] <- NA
        }
      }
      if (!is.null(p.factor)) {score <- score/p.factor[names(score)]}
      #ranking the score (of snps or groups) with decreasing order
      sorted.score <- sort(score, decreasing = T, na.last = NA)
      if(configs[['verbose']]){
        snpnetLogger("print top 10 group scores:")
        print(sorted.score[1:10])
      }
      
      
      if (length(sorted.score) > 0) {
        # select top M features/ groups and union with active set to construct strong set. 
        if (!(isSGL) & !(isGL)){
          features.to.add <- names(sorted.score)[1:min(configs[['num.snps.batch']], length(sorted.score))]
          for(s in splits){
            tmp.features.add <- prepareFeatures(pgen[[s]], vars, features.to.add, stats)
            if (!is.null(features[[s]])) {
              features[[s]][, colnames(tmp.features.add) := tmp.features.add]
            } else {
              features[[s]] <- tmp.features.add
            }
            rm(tmp.features.add)
          }
        } else{
          # When the penalty of the model is group lasso or sparse group lasso
          groups.to.add <- names(sorted.score)[1:min(configs[['num.groups.batch']], length(sorted.score))]
          features.to.add <- c(groups.to.add, paste0(groups.to.add, "_E"))
          if(configs[['verbose']]){
            snpnetLogger(sprintf("%d of groups or %d of features to be added", length(groups.to.add), length(features.to.add)))
          }
          for(s in splits){
            # splits: train or valid
            tmp.features.add <- prepareFeaturesEnv(pgen[[s]], vars, groups.to.add, stats, env, phe[[s]])
            if (!is.null(features[[s]])){
              #features[[s]][, colnames(tmp.features.add) := tmp.features.add]
              common_cols <- intersect(colnames(features[[s]]), colnames(tmp.features.add))
              tmp.features.add <- tmp.features.add[, !common_cols, with = FALSE, drop = FALSE]
              features[[s]] = cbind( as.matrix(features[[s]]), as.matrix(tmp.features.add) )
            }else{
              features[[s]] <- as.matrix(tmp.features.add)
            }
            rm(tmp.features.add)
            if(configs[['verbose']]) snpnetLogger(sprintf("The dim of features[[%s]] is: %d %d", s, dim(features[[s]])[1], dim(features[[s]])[2]))
          }
        }
        
      } else {
        break
      }
      if (increase.group.size)  # increase batch size when no new valid solution is found in the previous iteration, but after another round of adding new variables
        configs[['num.groups.batch']] <- configs[['num.groups.batch']] + configs[['increase.size']]
      if (configs[['verbose']]) snpnetLoggerTimeDiff("End updating feature matrix.", time.update.start, indent=2)
      if (configs[['verbose']]) {
        snpnetLogger(paste0("- # ever-active variables: ", length(features.to.keep), "."), indent=2)
        snpnetLogger(paste0("- # newly added variables: ", length(features.to.add), "."), indent=2)
        snpnetLogger(paste0("- Total # variables in the strong set: ", ncol(features[['train']]), "."), indent=2)
      }
      
      ### --- Fit models --- ###
      if (configs[['verbose']]){
        if(configs[['use.glmnetPlus']]){
          snpnetLogger("Start fitting Glmnet with glmnetPlus ...", indent=1)
        }else{
          snpnetLogger("Start fitting Glmnet ...", indent=1)
        }
      }
      
      # generate membership group
      if(isGL | isSGL){
        # covariates: environment, covariates, 
        # get the number of environment and covariates
        cov_len <- length(c(env, covariates))
        # set the group.info
        g_and_ge_len <- ncol(features[['train']]) - cov_len
        all.feature.names <- colnames(features[['train']])
        feature.names <- all.feature.names[-(1:cov_len)]
        prefixes <- str_replace(feature.names, "_E$", "")
        #membership vector
        membership_vector <- match(prefixes, unique(prefixes)) + 1
        group <- c(rep(1, each = cov_len), membership_vector)
        #number of unique groups
        num_of_groups<-length(unique(prefixes)) # except the covariates 
        if(configs[['verbose']]){
          snpnetLogger(sprintf("# of groups: %d, # of G and GxE: %d", num_of_groups, g_and_ge_len))
        }
      }
      # penalty.factor
      if (is.null(p.factor)){
        if(isGL | isSGL){
          penalty.factor <- c(0, rep(sqrt(2), num_of_groups))
        } else{
          penalty.factor <- rep(1, ncol(features[['train']]))
          # do not penalize covariate! feature[['train']] includes covariates
          penalty.factor[seq_len(length(c(covariates, env)))] <- 0
        }
        
      } else {
        # in sgl and gl, the penalty factor will be null. They won't go to this branch
        penalty.factor <- c(rep(0, length(c(covariates, env))), p.factor[colnames(features[['train']])[-(1:length(c(env, covariates)))]])
      }
      current.lams <- full.lams[1:num.lams]
      current.lams.adjusted <- full.lams[1:num.lams] * sum(penalty.factor) / length(penalty.factor)  # adjustment to counteract penalty factor normalization in glmnet
      time.fit.start <- Sys.time()
      
      if (configs[['use.glmnetPlus']]) {
        snpnetLogger("run glmnetPlus")
        start.lams <- lambda.idx   # start index in the whole lambda sequence
        if (!is.null(prev.beta)) {
          beta0 <- rep(1e-20, ncol(features[['train']]))
          beta0[match(names(prev.beta), colnames(features[['train']]))] <- prev.beta
        } else {
          beta0 <- prev.beta
        }
        if(family == "cox"){
          glmfit <- glmnetPlus::glmnet(
            features[['train']], surv[['train']], family = family, alpha = alpha,
            lambda = current.lams.adjusted[start.lams:num.lams], penalty.factor = penalty.factor,
            standardize = configs[['standardize.variant']], thresh = configs[['glmnet.thresh']], beta0 = beta0
          )
          pred.train <- stats::predict(glmfit, newx = features[['train']])
          residual <- computeCoxgrad(pred.train, response[['train']], status[['train']])
        } else {
          
          # we need to add sparse group lasso and group lasso
          #beta0 is the beta from the last lambda
          glmfit <- glmnetPlus::glmnet(
            features[['train']], response[['train']], family = family, alpha = alpha,
            lambda = current.lams.adjusted[start.lams:num.lams], penalty.factor = penalty.factor,
            standardize = configs[['standardize.variant']], thresh = configs[['glmnet.thresh']],
            type.gaussian = "naive", beta0 = beta0
          )
          if(family=="gaussian"){
            residual <- glmfit$residuals
            pred.train <- response[['train']] - residual
          }else{
            pred.train <- stats::predict(glmfit, newx = as.matrix(features[['train']]), type = "response")
            residual <- response[['train']] - pred.train
          }
          
        }
        
      } else { # configs[['use.glmnetPlus']] == FALSE
        snpnetLogger("Not use glmnetplus")
        # start.lams <- 1
        start.lams <- lambda.idx
        #tmp.features.matrix <- as.matrix(features[['train']])
        # tmp.features.matrix = Matrix(as.matrix(features[['train']]), sparse=TRUE)
        if(configs[['verbose']]){
          snpnetLogger(sprintf("dimension of feature matrix %d %d", dim(features[['train']])[1], dim(features[['train']])[2]))
        }
        
        if(family=="cox"){
          # glmfit <- glmnet::glmnet(
          #     features[['train']], surv[['train']], family = family, alpha = alpha,
          #     lambda = current.lams.adjusted[start.lams:num.lams], penalty.factor = penalty.factor,
          #     standardize = configs[['standardize.variant']], thresh = configs[['glmnet.thresh']]
          # )
          # pred.train <- stats::predict(glmfit, newx = features[['train']])
          # residual <- computeCoxgrad(pred.train, response[['train']], status[['train']])
        }else{
          if (isGL){
            snpnetLogger("run GL")
            # Group Lasso REFERENCE: https://cran.r-project.org/web/packages/gglasso/gglasso.pdf
            if(family == "gaussian"){
              loss = "ls"
            }else{ #binomial
              loss = "logit"
            }
            glfit = gglasso::gglasso(
              x=features[['train']], y=response[['train']], group=group, 
              loss=loss, lambda = current.lams[start.lams:num.lams], 
              eps = 1e-8, maxit = 1e+06, pf = penalty.factor
            )
            pred.train <- predict(glfit, newx = features[['train']], type = "link")
            residual <- response[['train']] - pred.train
            snpnetLogger("The current lambdas:")
            print(head(current.lams[start.lams:num.lams]))
            snpnetLogger("The residual of group lasso:")
            print(head(residual))
            fitted.model<-glfit
            
          }else if(isSGL){
            snpnetLogger("run SGL")
            #SGL: Reference: https://cran.r-project.org/web/packages/SGL/SGL.pdf
            
            if(configs[['verbose']]){
              snpnetLogger("memory usage")
              print(gc())
            }
            
            if(configs[['verbose']]){
              snpnetLogger("objects' memory usage (Top 10)")
              objs <- ls()
              obj_sizes <- sapply(objs, function(x) object.size(get(x)))
              obj_sizes_mb <- round(obj_sizes / (1024^2), 1)  # Convert to megabytes
              obj_memory_usage <- data.frame(Object=objs, Memory_MB=obj_sizes_mb)
              obj_memory_usage <- obj_memory_usage[order(-obj_memory_usage$Memory_MB), ]
              rownames(obj_memory_usage) = NULL 
              print(obj_memory_usage[1:10,])
            }
            
            sglfit <- sparsegl_modified(features[['train']], y=response[['train']], 
                                       group = group, maxit=1e+6, eps=1e-8, 
                                       lambda= current.lams[start.lams:num.lams], 
                                       standardize=FALSE, pf_group = penalty.factor, family= family, asparse=tau)
            #pred.train <- predict(sglfit, newx = features[['train']])
            
            pred.train <- predict_sparsegl(sglfit, newx = features[['train']])
            residual <- response[['train']] - pred.train

            fitted.model<-sglfit
            
          }else{
            glmfit <- glmnet::glmnet(
              features[['train']], response[['train']], family = family, alpha = alpha,
              lambda = current.lams.adjusted[start.lams:num.lams], penalty.factor = penalty.factor,
              standardize = configs[['standardize.variant']], thresh = configs[['glmnet.thresh']],
              type.gaussian = "naive"
            )
            pred.train <- stats::predict(glmfit, newx = features[['train']], type = "response")
            residual <- response[['train']] - pred.train
            fitted.model<-glmfit
            
          }
          
        }
        #rm(tmp.features.matrix) # save memory
      }
      
      fit.results[[iter]] <- fitted.model
      rownames(residual) <- rownames(phe[['train']]) # assign sample ids to residual as row
      colnames(residual) <- start.lams:num.lams # assign lambda indices to residual as col. 
      # Note that this is not lambda
      if (configs[['verbose']]) snpnetLoggerTimeDiff("End fitting Models.", time.fit.start, indent=2)
      
      ### --- KKT Check --- ###
      if (configs[['verbose']]) snpnetLogger("Start checking KKT condition ...", indent=1)
      time.KKT.start <- Sys.time()
      
      check.obj <- KKT.check(
        residual, genotype.pfile, vars, phe[['train']], env, nrow(phe[['train']]),
        current.lams[start.lams:num.lams], 1,
        stats, fitted.model, configs, iter, p.factor, alpha, tau, isSGL, isGL
      )
      snpnetLogger("KKT check obj done ...", indent=1)
      
      # update the max valid index in the whole lambda sequence
      max.valid.idx <- check.obj[["max.valid.idx"]] + (start.lams - 1)
      lambda.idx <- max.valid.idx + 1
      
      # Update the lambda index of variants added
      if (configs[['rank']] && check.obj[["max.valid.idx"]] > 0){
        tmp <- 1
        for (lam.idx in start.lams:max.valid.idx){
          
          #looking for the current active features to the each valid lambda
          current_active <- setdiff(names(which(fitted.model$beta[, tmp] != 0)), c(env, covariates))
          tmp <- tmp + 1
          var.rank[current_active] = pmin(var.rank[current_active], lam.idx)
        }
      }
      
      if (configs[['use.glmnetPlus']] && check.obj[["max.valid.idx"]] > 0) {
        prev.beta <- glmfit$beta[, check.obj[["max.valid.idx"]]]
        prev.beta <- prev.beta[prev.beta != 0]
      }
      num.new.valid[iter] <- check.obj[["max.valid.idx"]]
      # if (configs[['use.glmnetPlus']]) {
      #   num.new.valid[iter] <- check.obj[["max.valid.idx"]]
      # } else {
      #   # num.new.valid[iter] <- check.obj[["max.valid.idx"]] - ifelse(iter > 1, num.new.valid[iter-1], 0)
      #   num.new.valid[iter] <- check.obj[["max.valid.idx"]]
      # }
      
      if ( prev.max.valid.idx == max.valid.idx ) {
        # there is no valid solution in this iteration
        features.to.keep <- union(features.to.keep, features.to.add)
        increase.group.size <- TRUE
      } else {
        for (j in 1:check.obj[["max.valid.idx"]]) {
          #a0 is the intercept value
          if(isGL | isSGL){
            a0[[j + (start.lams - 1)]] <- as.numeric(fitted.model$b0[j])
          }else{
            a0[[j + (start.lams - 1)]] <- as.numeric(fitted.model$a0[j])
          }
          
          beta[[j + (start.lams - 1)]] <- fitted.model$beta[, j]
        }
        if (validation) {
          time.val.pred.start <- Sys.time()
          if (family == "cox") {
            pred.val <- stats::predict(glmfit, newx = as.matrix(features[['val']]), lambda = current.lams.adjusted[start.lams:max.valid.idx])
          } else if (configs[['use.glmnetPlus']]) {
            pred.val <- glmnetPlus::predict.glmnet(glmfit, newx = as.matrix(features[['val']]), lambda = current.lams.adjusted[start.lams:max.valid.idx], type = "response")
          } else {
            if(isSGL){
              #pred.val <- predict(fitted.model, newx = Matrix(as.matrix(features[['val']], sparse=TRUE)), s = current.lams[start.lams:max.valid.idx], type="response")
              pred.val <- predict_sparsegl(fitted.model, newx = Matrix(as.matrix(features[['val']], sparse=TRUE)), s = current.lams[start.lams:max.valid.idx], type="response")
            }else if(isGL){
              
              # pred.val <- predict(fitted.model, newx = Matrix(as.matrix(features[['val']], sparse=TRUE)), s = current.lams[start.lams:max.valid.idx], type="link")
              pred.val <- predict(fitted.model, newx = features[['val']], s = current.lams[start.lams:max.valid.idx], type="link")

            }else{
              pred.val <- glmnet::predict.glmnet(fitted.model, newx = as.matrix(features[['val']]), lambda = current.lams.adjusted[start.lams:max.valid.idx], type = "response")
            }
          }
          snpnetLoggerTimeDiff("Time of prediction on validation matrix", time.val.pred.start, indent=2)
        }
        
        # compute metric
        if (family == "cox") {
          metric.train[start.lams:max.valid.idx] <- computeMetric(pred.train[, 1:check.obj[["max.valid.idx"]], drop = F], surv[['train']], configs[['metric']])
          if (validation) metric.val[start.lams:max.valid.idx] <- computeMetric(pred.val, surv[['val']], configs[['metric']])
        } else {
          snpnetLogger('metric train')
          metric.train[start.lams:max.valid.idx] <- computeMetric(pred.train[, 1:check.obj[["max.valid.idx"]], drop = F], response[['train']], configs[['metric']])
          if (validation){
            snpnetLogger('metric val.')
            metric.val[start.lams:max.valid.idx] <- computeMetric(pred.val, response[['val']], configs[['metric']])
          }
        }
        
        score <- check.obj[["score"]]
        # create Active group set
        # active variables or not
        is.ever.active <- apply(fitted.model$beta[, 1:check.obj[["max.valid.idx"]], drop = F], 1, function(x) any(x != 0))
        # remove environment or covariates if they exist
        if((!is.null(env)) | (!is.null(covariates))){
          is.ever.active <- is.ever.active[! (names(is.ever.active) %in% c(env, covariates))]
        }
        # get active groups: active group: beta_g or beta_ge is nonzero
        active.features<-rownames(fitted.model$beta)[is.ever.active]
        if(length(active.features)>0){
          active.group <- unique(gsub("_E", "", active.features))
          # get active features from those groups
          active.feat <- rep(active.group, each = 2)
          active.feat[seq(2, length(active.feat), 2)] = paste0(active.feat[seq(2, length(active.feat), 2)], '_E')
          features.to.keep <- union(active.feat, features.to.keep)
        }
        increase.group.size <- FALSE
        
      }
      
      if (configs[['verbose']]) snpnetLoggerTimeDiff("End checking KKT condition.", time.KKT.start, indent=2)
      
      if (configs[['save']]) {
        save(metric.train, metric.val, fit.results, full.lams, a0, beta, prev.beta, max.valid.idx,
             features.to.keep, num.lams, lambda.idx, score, num.new.valid, increase.group.size, configs,
             file = file.path(configs[['results.dir']], configs[["save.dir"]], paste0("output_iter_", iter, ".RData")))
      }
      
      if ( prev.max.valid.idx < max.valid.idx ) {
        # there are valid solution(s) for at least one lambda
        if (validation) {
          snpnetLogger('Training and validation metric:', indent=1)
        }else{
          snpnetLogger('Training metric:', indent=1)
        }
        for (klam in (prev.max.valid.idx+1):max.valid.idx) {
          if (validation) {
            snpnetLogger(paste0("- Lambda idx ", klam, ". Training: ", metric.train[klam], ". Validation: ", metric.val[klam]), indent=1)
          } else {
            snpnetLogger(paste0("- Lambda idx ", klam, ". Training: ", metric.train[klam], ". "), indent=1)
          }
        }
        prev.max.valid.idx <- max.valid.idx
      }
      time.iter.end <- Sys.time()
      snpnetLoggerTimeDiff(paste0("End iteration ", iter, '.'), time.iter.start, time.iter.end, indent=1)
      snpnetLoggerTimeDiff("The total time since start.", time.start, time.iter.end, indent=2)
      
      ### --- Check stopping criteria --- ####
      if (max.valid.idx == configs[['nlambda']]) break
      if (validation && checkEarlyStopping(metric.val, max.valid.idx, iter, configs)) break
    }
  }
  snpnetLoggerTimeDiff("End snpnet.", time.start)
  if(! configs[['save']]) cleanUpIntermediateFiles(configs)
  if(configs[['verbose']]) print(gc())
  
  out <- list(metric.train = metric.train, metric.val = metric.val, fit.results = fit.results,
              full.lams = full.lams, a0 = a0, beta = beta, configs = configs, var.rank=var.rank,
              lambda.min.ratio = lambda.min.ratio, stats = stats)
  out
}

