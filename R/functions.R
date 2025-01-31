#' Predict from the Fitted Object or File
#'
#' @usage predict_geiprs(fit = NULL, saved_path = NULL, new_genotype_file, new_phenotype_file,
#'   phenotype, gcount_path = NULL, meta_dir = NULL, meta_suffix = ".rda", covariate_names = NULL,
#'   split_col = NULL, split_name = NULL, idx = NULL, family = NULL, snpnet_prefix = "output_iter_",
#'   snpnet_suffix = ".RData", snpnet_subdir = "results", configs = list(zstdcat.path = "zstdcat",
#'   zcat.path='zcat'))
#'
#' @param fit Fitted object returned from the snpnet function. If not specified, `saved_path` has to
#'   be provided.
#' @param saved_path Path to the file that saves the fit object. The full path is constructed as
#'   ${saved_path}/${snpnet_subdir}/${snpnet_prefix}ITER${snpnet_suffix}, where ITER will be the
#'   maximum index found in the snpnet subdirectory. If not specified, `fit` has to be provided.
#' @param new_genotype_file Path to the new suite of genotype files. new_genotype_file.{pgen, psam,
#'   pvar.zst}.
#'   must exist.
#' @param new_phenotype_file Path to the phenotype. The header must include FID, IID. Used for extracting covariates and computing metrics.
#' @param phenotype Name of the phenotype for which the fit was computed.
#' @param env Name of the env for the adjustment variable and interaction term.
#' @param gcount_path Path to the saved gcount file on which the meta statistics can be computed. Only if `saved_path` is specified.
#' @param meta_dir (Depreciated) Path to the saved meta statistics object. The full path is constructed as ${meta_dir}/${STAT}${meta_suffix}, where such files should exist for STAT = pnas, means and optionally sds. Only if `saved_path` is specified.
#' @param meta_suffix (Depreciated) Extension suffix of the meta statistics files. Only if `saved_path` is specified.
#' @param covariate_names Character vector of the names of the adjustment covariates.
#' @param split_col Name of the split column. If NULL, all samples will be used.
#' @param split_name Vector of split labels where prediction is to be made. Should be a combination of "train", "val", "test".
#' @param idx Vector of lambda indices on which the prediction is to be made. If not provided, will predict on all lambdas found.
#' @param family Type of the phenotype: "gaussian" for continuous phenotype and "binomial" for binary phenotype.
#' @param snpnet_prefix Prefix of the snpnet result files used to construct the full path. Only if `saved_path` is specified.
#' @param snpnet_suffix Extension suffix of the snpnet result files used to construct the full path. Only if `saved_path` is specified.
#' @param snpnet_subdir Name of the snpnet result subdirectory holding multiple result files for one phenotype. Only if `saved_path` is specified.
#' @param configs Additional list of configs including path to either zstdcat or zcat.
#'
#' @return A list containing the prediction and the resopnse for which the prediction is made.
#'
#' @export
predict_geiprs <- function(fit = NULL, saved_path = NULL, new_genotype_file, new_phenotype_file, phenotype, env,
                           gcount_path = NULL, meta_dir = NULL, meta_suffix = ".rda",
                           covariate_names = NULL, split_col = NULL, split_name = NULL, idx = NULL,
                           family = NULL,
                           snpnet_prefix = "output_iter_", snpnet_suffix = ".RData", snpnet_subdir = "results",
                           configs = list(zstdcat.path = "zstdcat", zcat.path='zcat')) {

  
  ID <- original_ID <- NULL  # to deal with "no visible binding for global variable"
  ALT <- MISSING_CT <- OBS_CT <- HAP_ALT_CTS <- HET_REF_ALT_CTS <- TWO_ALT_GENO_CTS <- NULL
  stats_msts <- stats_means <- stats_pNAs <- stats_SDs <- NULL

  if (is.null(fit) && is.null(saved_path)) {
    stop("Either fit object or file path to the saved object should be provided.\n")
  }
  if (is.null(fit)) {
    phe_dir <- file.path(saved_path, snpnet_subdir)
    files_in_dir <- list.files(phe_dir)
    result_files <- files_in_dir[startsWith(files_in_dir, snpnet_prefix) & endsWith(files_in_dir, snpnet_suffix)]
    max_iter <- max(as.numeric(gsub(snpnet_suffix, "", gsub(pattern = snpnet_prefix, "", result_files))))

    e <- new.env()
    load(file.path(saved_path, snpnet_subdir, paste0(snpnet_prefix, max_iter, snpnet_suffix)), envir = e)
    a0 <- e$a0
    beta <- e$beta

    stats <- list()
    if (!is.null(gcount_path)) {
      gcount_df <-
        data.table::fread(gcount_path) %>%
        dplyr::rename(original_ID = ID) %>%
        dplyr::mutate(
          ID = paste0(original_ID, '_', ALT),
          stats_pNAs  = MISSING_CT / (MISSING_CT + OBS_CT),
          stats_means = (HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_msts  = (HAP_ALT_CTS + HET_REF_ALT_CTS + 4 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_SDs   = stats_msts - stats_means * stats_means
        )
      stats[["pnas"]]  <- gcount_df %>% dplyr::pull(stats_pNAs)
      stats[["means"]] <- gcount_df %>% dplyr::pull(stats_means)
      stats[["sds"]]   <- gcount_df %>% dplyr::pull(stats_SDs)
      for(key in names(stats)){
        names(stats[[key]]) <- gcount_df %>% dplyr::pull(ID)
      }
    } else {
      stats[["pnas"]] <- readRDS(file.path(meta_dir, paste0("pnas", meta_suffix)))
      stats[["means"]] <- readRDS(file.path(meta_dir, paste0("means", meta_suffix)))
      if (file.exists(file.path(meta_dir, paste0("sds", meta_suffix)))) {
        stats[["sds"]] <- readRDS(file.path(meta_dir, paste0("sds", meta_suffix)))
      }
    }
  } else {
    a0 <- fit$a0
    beta <- fit$beta
    stats <- fit$stats
  }

  #5462
  feature_names <- unique(unlist(sapply(beta, function(x) names(x[x != 0]))))
  feature_names <- setdiff(feature_names, c(covariate_names, env)) # features names now include E, G and GxE
  if (is.null(idx)) idx <- seq_along(a0)

  ids <- list()
  #67443
  ids[["psam"]] <- readIDsFromPsam(paste0(new_genotype_file, '.psam'))

  #get phenotype data
  #61278

  phe_master <- readPheMaster(new_phenotype_file, ids[['psam']], family, env, covariate_names, phenotype, NULL, split_col, configs)

  #61278
  if (length(c(env,covariate_names)) > 0) {
    cov_master <- as.matrix(phe_master[, c(env, covariate_names), with = F])
    cov_no_missing <- apply(cov_master, 1, function(x) all(!is.na(x)))
    phe_master <- phe_master[cov_no_missing, ]
  }

  if (is.null(family)) family <- inferFamily(phe_master, phenotype, NULL)
  if (is.null(configs[["metric"]])) configs[["metric"]] <- setDefaultMetric(family)
  #61278
  if (is.null(split_col)) {
    split_name <- "test"
    ids[["test"]] <- phe_master$ID
  } else {
    for (split in split_name) {

      ids[[split]] <- phe_master$ID[phe_master[[split_col]] == split]
      if (length(ids[[split]]) == 0) {
        warning(paste("Split", split, "doesn't exist in the phenotype file. Excluded from prediction.\n"))
        split_name <- setdiff(split_name, split)
      }
    }
  }


  phe <- list()
  for (split in split_name) {
    ids_loc <- match(ids[[split]], phe_master[["ID"]])
    phe[[split]] <- phe_master[ids_loc]
  }

  covariates <- list() # the covariates also include env

  for (split in split_name) {
    if (length(c(env, covariate_names)) > 0) {
      covariates[[split]] <- phe[[split]][, c(env, covariate_names), with = FALSE]
    } else {
      covariates[[split]] <- NULL
    }
  }
  #612767
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(configs[["zstdcat.path"]], ' ', paste0(new_genotype_file, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  pvar <- pgenlibr::NewPvar(paste0(new_genotype_file, '.pvar.zst'))
  chr <- list()
  for (split in split_name) {
    chr[[split]] <- pgenlibr::NewPgen(paste0(new_genotype_file, '.pgen'), pvar = pvar, sample_subset = match(ids[[split]], ids[["psam"]]))
  }
  pgenlibr::ClosePvar(pvar)

  features <- list()

  for (split in split_name) {

    if (!is.null(covariates[[split]])) {
      # 61278 1
      features[[split]] <- data.table::data.table(scale(covariates[[split]]))
      g_and_gxe<-prepareFeaturesEnv2(chr[[split]], vars, feature_names, stats, env, phe[[split]])
      feature_names_new <- colnames(g_and_gxe)
      features[[split]][, (feature_names_new) := g_and_gxe]
      rm(g_and_gxe) # gc
    } else {
      features[[split]] <- prepareFeatures(chr[[split]], vars, feature_names, stats)
    }
  }

  pred <- list()
  metric <- list()
  metric.prs <- list()
  metric.prs.logpval <- list()
  for (split in split_name) {
    pred[[split]] <- array(dim = c(nrow(features[[split]]), length(idx)),
                           dimnames = list(ids[[split]], paste0("s", idx-1)))
    metric[[split]] <- rep(NA, length(idx))
    metric.prs[[split]] <- rep(NA, length(idx))
    metric.prs.logpval[[split]] <- rep(NA, length(idx))
    names(metric[[split]]) <- paste0("s", idx-1)
    names(metric.prs[[split]]) <- paste0("s", idx-1)
    names(metric.prs.logpval[[split]]) <- paste0("s", idx-1)
  }

  pred_prs_g <- pred
  pred_prs_ge <- pred
  

  response <- list()
  for (split in split_name) {
    response[[split]] <- phe[[split]][[phenotype]]
  }
  
  for (split in split_name) {
    for (i in idx) {
      active_names <- names(beta[[i]])[beta[[i]] != 0]
      pred[[split]] <- getPred(active_names=active_names, feat=features[[split]], bias=a0[[i]], beta_i=beta[[i]], i=i, pred_data=pred[[split]])
      
      active_names_g <- active_names[(!grepl("_E$", active_names)) & (!active_names %in% covariates[[split]]) ]
      pred_prs_g[[split]] <- getPRS(active_names_g, features[[split]], a0[[i]], beta[[i]], i, pred_prs_g[[split]])
      
      active_names_ge <- active_names[(grepl("_E$", active_names)) & (!active_names %in% covariates[[split]]) ]
      pred_prs_ge[[split]] <- getPRSGE(active_names_ge, features[[split]], a0[[i]], beta[[i]], i, pred_prs_ge[[split]])
    }
    metric[[split]] <- computeMetric(pred[[split]], response[[split]], configs[["metric"]])
    result_prs <- computeMetricPRS(features[[split]][[env]], pred_prs_g[[split]], pred_prs_ge[[split]], response[[split]], configs[["metric"]])
    metric.prs[[split]] <- result_prs[1]
    metric.prs.logpval[[split]] <- result_prs[2]
  }

  list(prediction = pred, prs_g=pred_prs_g, prs_ge=pred_prs_ge, response = response, env_val=features[[split]][[env]], metric = metric, metric.prs=metric.prs, metric.prs.logpval=metric.prs.logpval)
}
#get prediction or prs
getPred <- function(active_names, feat, bias, beta_i, i, pred_data){
  # if it is environment, make sure let environment_beta, and snps exist at the same time
  if (length(active_names) > 0) {
        active_names <- Reduce(intersect, list(v1= active_names, v2= colnames(feat)))
        features_single <- as.matrix(feat[, active_names, with = F])
  } else {
        features_single <- matrix(0, nrow(feat), 0)
  }
  pred_single <- bias + features_single %*% beta_i[active_names]
  pred_data[, c(paste0("s", i-1))] <- as.matrix(pred_single)
  pred_data
}

#get prediction or prs
getPRS <- function(active_names, feat, bias, beta_i, i, pred_data){
  # if it is environment, make sure let environment_beta, and snps exist at the same time
  if (length(active_names) > 0) {
        active_names <- Reduce(intersect, list(v1= active_names, v2= colnames(feat)))
        features_single <- as.matrix(feat[, active_names, with = F])
  } else {
        features_single <- matrix(0, nrow(feat), 0)
  }
  pred_single <- features_single %*% beta_i[active_names]
  pred_data[, c(paste0("s", i-1))] <- as.matrix(pred_single)
  pred_data
}

#get prs_ge
getPRSGE <- function(active_GE_names, feat, bias, beta_i, i, pred_data){
  #this feat should be feature matrix
  # if it is environment, make sure let environment_beta, and snps exist at the same time
  if (length(active_GE_names) > 0) {
        
        # get active G names
        active_G_names <- sub("_E$", "", active_GE_names) 
        
        #get intersection with all features. In order to get dosage
        active_G_names.inter <- Reduce(intersect, list(v1= active_G_names, v2= colnames(feat)))
        
        #extract genotype dosage
        features_single <- as.matrix(feat[, active_G_names.inter, with = F])

        #rename active_G_names.inter to active_GE_names.inter by add "_E"
        active_GE_names.inter <- paste0(active_G_names.inter, "_E")
  } else {
        features_single <- matrix(0, nrow(feat), 0)
  }
  pred_single <- features_single %*% beta_i[active_GE_names.inter]
  pred_data[, c(paste0("s", i-1))] <- as.matrix(pred_single)
  pred_data
}
#' @importFrom data.table set as.data.table
#' @importFrom magrittr %>%
#' @importFrom dplyr n
prepareFeatures <- function(pgen, vars, names, stat) {
  buf <- pgenlibr::ReadList(pgen, match(names, vars), meanimpute=F)
  features.add <- as.data.table(buf)
  colnames(features.add) <- names
  for (j in 1:length(names)) {
    set(features.add, i=which(is.na(features.add[[j]])), j=j, value=stat[["means"]][names[j]])
  }
  features.add
}

#' @importFrom data.table set as.data.table
#' @importFrom magrittr %>%
#' @importFrom dplyr n
# try a small data to check whether it is correct
prepareFeaturesEnv<- function(pgen.one, vars, names, stats, env, phenome) {
  
  buf <- pgenlibr::ReadList(pgen.one, match(names, vars), meanimpute=F)
  features.add <- as.data.table(buf) # this table only include the buf data
  colnames(features.add) <- names
  # SNP1, SNP2, ..., SNPn
  # 1, 1, 2
  # 1, 2, 1
  for (j in 1:length(names)) {
    # this is used for filling the na values with the mean values.
    # in table featurs.add, look for the multiple i rows and j th columns, to fill
    # the values stat[['means']] of the mean dosage value of a variant.
    # set() function extract each column to update the NA value into the mean values.
    set(features.add, i=which(is.na(features.add[[j]])), j=j, value=stats[["means"]][names[j]])
  }
  #let's cbind the datatable with the env * datatable
  phenome <- as.data.frame(phenome)
  env.vec <- phenome[, c(env), drop=TRUE] # the rownames of phenotype is sample ID
  #multiply the env and features.add
  features.add.gxe <- sweep(features.add, 1, env.vec, "*")
  colnames(features.add.gxe) <- paste0(names, "_E")
  features.add.g_and_gxe <- cbind(features.add, features.add.gxe)
  #output
  features.add.g_and_gxe
}

#' @importFrom data.table set as.data.table
#' @importFrom magrittr %>%
#' @importFrom dplyr n
# prepare Features Env (consider multi alleic problem)
prepareFeaturesEnv2<- function(pgen.one, vars, names, stats, env, phenome) {
  #names: 5461; vars:612767
  #simulation: 3129, vars:40000
  names.G = sub("_E$", "", names)
  # remove the duplicate genotypes 4701; 2894
  names.G.uniq = unique(names.G)
  #match those genotype data 4701; 2894
  match_results <- match(names.G.uniq, vars)
  # It has some NA due to allele flip. We can directly remove those G that does
  # 4691; 2894
  names.G.uniq.intersect = Reduce(intersect, list(v1 = names.G.uniq, v2 = vars))
  # get flipped alleles
  # fliped.allles<-setdiff(names.G.uniq, names.G.uniq.intersect)

  buf <- pgenlibr::ReadList(pgen.one, match(names.G.uniq.intersect, vars), meanimpute=F)
  #61278  4691; 10000 2894
  features.add <- as.data.table(buf) # this table only include the buf data
  colnames(features.add) <- names.G.uniq.intersect
  for (j in 1:length(names.G.uniq.intersect)) {
    # this is used for filling the na values with the mean values.
    # in table featurs.add, look for the multiple i rows and j th columns, to fill
    # the values stat[['means']] of the mean dosage value of a variant.
    # set() function extract each column to update the NA value into the mean values.
    set(features.add, i=which(is.na(features.add[[j]])), j=j, value=stats[["means"]][names.G.uniq.intersect[j]])
  }
  #let's cbind the datatable with the env * datatable
  phenome <- as.data.frame(phenome)
  env.vec <- phenome[, c(env), drop=TRUE] # the rownames of phenotype is sample ID
  #scale the env
  env.vec <- scale(env.vec)
  #multiply the env and features.add
  # features.add.gxe <- sweep(features.add, 1, env.vec, "*")
  features.add.gxe <- features.add * env.vec
  colnames(features.add.gxe) <- paste0(names.G.uniq.intersect, "_E")
  features.add.g_and_gxe <- cbind(features.add, features.add.gxe)
  #output
  #Lastly, filter the feature.add.gxe with names. names:5449 features.add.g_and_gxe: 61278  9382
  # names.inter <- Reduce(intersect, list(v1=names, v2=colnames(features.add.g_and_gxe)))
  # features.add.g_and_gxe2 61278  5449
  # features.add.g_and_gxe2 <- features.add.g_and_gxe[, names.inter, with=FALSE]
  features.add.g_and_gxe
}

#names might contain G and GxE
  #first, we only maintain the G
 

computeLambdas <- function(lambda.max, nlambda, lambda.min.ratio) {
  # lambda.max <- max(score, na.rm = T)
  
  lambda.min <- lambda.max * lambda.min.ratio
  # exponential decay for lambda values
  full.lams <- exp(seq(from = log(lambda.max), to = log(lambda.min), length.out = nlambda))
  full.lams
}

inferFamily <- function(phe, phenotype, status){
    if (all(unique(phe[[phenotype]] %in% c(0, 1, 2, -9)))) {
        family <- "binomial"
    } else if(!is.null(status) && (status %in% colnames(phe))) {
        family <- "cox"
    } else {
        family <- "gaussian"
    }
    family
}

readIDsFromPsam <- function(psam){
    FID <- IID <- NULL  # to deal with "no visible binding for global variable"
    df <- data.table::fread(psam)
    # print(head(df))
    if (!('#FID' %in% colnames(df))) {
      if ('IID' %in% colnames(df)) {
        warning("#FID column not found in the psam file. Assume and use FID = IID.")
        df['#FID'] <- df['IID']
      } else {
        stop('IID column not found in the psam file.')
      }
    }
    df <- df %>%
    dplyr::rename('FID' = '#FID') %>%
    dplyr::mutate(ID = paste(FID, IID, sep='-'))
    # print(head(df))
    df$ID
}

cat_or_zcat <- function(filename, configs=list(zstdcat.path='zstdcat', zcat.path='zcat')){
    if(stringr::str_ends(basename(filename), '.zst')){
        return(configs[['zstdcat.path']])
    }else if(stringr::str_ends(basename(filename), '.gz')){
        return(configs[['zcat.path']])
    }else{
        return('cat')
    }
}

readPlinkKeepFile <- function(keep_file){
    ID <- NULL  # to deal with "no visible binding for global variable"
    keep_df <- data.table::fread(keep_file, colClasses='character', stringsAsFactors=F)
    keep_df$ID <- paste(keep_df$V1, keep_df$V2, sep='-')
    keep_df %>% dplyr::pull(ID)
}

#' Read covariates, phenotype(s), and env(if existing) from the provided file path
#'
#' Read covariates, phenotype(s), and env(if existing) from the provided file path. Exclude individuals that contain
#' any missing value in the covariates, miss all phenotype values or do not have corresponding
#' genotypes.
#'
#' @param phenotype.file the path of the file that contains the phenotype values and can be read as
#'                       as a table. There should be FID (family ID) and IID (individual ID) columns
#'                       containing the identifier for each individual, and the phenotype column(s).
#'                       (optional) some covariate columns and a column specifying the
#'                       training/validation split can be included in this file.
#' @param psam.ids a vector of ids read from the psam file.
#' @param family the type of the phenotype: "gaussian", "binomial", or "cox".
#' @param env the string of env name used for adjustment and interaction term.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file.
#' @param status the column name for the status column for Cox proportional hazards model.
#'               When running the Cox model, the specified column must exist in the phenotype file.
#' @param split.col the column name in the phenotype file that specifies the membership of individuals to
#'                  the training or the validation set. The individuals marked as "train" and "val" will
#'                  be treated as the training and validation set, respectively. When specified, the
#'                  model performance is evaluated on both the training and the validation sets.
#' @param configs a list of other config parameters. See more description in the `snpnet` function.
#'
#' @return a data.table including the requested columns.
#'
#' @export
readPheMaster <- function(phenotype.file, psam.ids, family, env, covariates, phenotype, status, split.col, configs){
  # > phenotype.file = new_phenotype_file
  # > psam.ids = ids[['psam']]
  # > covariates
  # Error: object 'covariates' not found
  # > covariates<-covariate_names
  # > status=NULL
  # > split.col
  # Error: object 'split.col' not found
  # > split.col<-split_col
  sort_order <- . <- ID <- NULL  # to deal with "no visible binding for global variable"

  if(!is.null(family) && family == 'cox'){
        selectCols <- c("FID", "IID", env, covariates, phenotype, status, split.col)
    } else{
        selectCols <- c("FID", "IID", env, covariates, phenotype, split.col)
    }
    # print("selectCols")
    # print(selectCols)

    # sed -e "s/^#//g" used for remove the pound symbol ("#") for each line
    phe.master.unsorted <- data.table::fread(
      cmd=paste(cat_or_zcat(phenotype.file, configs), phenotype.file, ' | sed -e "s/^#//g"'),
      colClasses = c("FID" = "character", "IID" = "character"), select = selectCols
    )
    #define ID column with the name of "FID_IID"
    phe.master.unsorted$ID <- paste(phe.master.unsorted$FID, phe.master.unsorted$IID, sep='-')

    # make sure the phe.master has the same individual ordering as in the genotype data
    # so that we don't have error when opening pgen file with sample subset option.
    phe.master <- phe.master.unsorted %>%
      dplyr::left_join(
        data.frame(ID = psam.ids, stringsAsFactors=F) %>%
          dplyr::mutate(sort_order = 1:n()),
        by='ID'
      ) %>%
      dplyr::arrange(sort_order) %>% dplyr::select(-sort_order) %>%
      data.table::as.data.table()
    rownames(phe.master) <- phe.master$ID

    # set the missing value to from -9 to NA
    for (name in c(covariates, env, phenotype)) {
      set(phe.master, i = which(phe.master[[name]] == -9), j = name, value = NA) # missing phenotypes are encoded with -9
    }

    # focus on individuals with complete covariates values
    if (is.null(covariates) && is.null(env) ) {
      phe.no.missing <- phe.master
    } else if (is.null(covariates) && !is.null(env)) {
      phe.no.missing <- phe.master %>%
        dplyr::filter_at(dplyr::vars(env), dplyr::all_vars(!is.na(.)))
    } else if (!is.null(covariates) && is.null(env)) {
      phe.no.missing <- phe.master %>%
        dplyr::filter_at(dplyr::vars(covariates), dplyr::all_vars(!is.na(.)))
    } else {
      phe.no.missing <- phe.master %>%
        dplyr::filter_at(dplyr::vars(c(env, covariates)), dplyr::all_vars(!is.na(.)))
    }
    # print("no missing 3")
    # print(nrow(phe.no.missing))
    # print(head(phe.no.missing %>% dplyr::filter(split %in% c('val'))))
    # focus on individuals with at least one observed phenotype values
    phe.no.missing <- phe.no.missing %>%
      dplyr::filter_at(dplyr::vars(phenotype), dplyr::any_vars(!is.na(.))) %>%
      dplyr::filter(ID %in% psam.ids) # check if we have genotype
    # print(tail(phe.no.missing))
    phe.no.missing.IDs <- phe.no.missing$ID
    
    # print("no missing 4")
    # print(head(phe.no.missing %>% dplyr::filter(split %in% c('val'))))
    # print(head(phe.no.missing))
    if(!is.null(split.col)){
        # focus on individuals in training and validation set
        phe.no.missing.IDs <- intersect(
            phe.no.missing.IDs,
            phe.master$ID[ (phe.master[[split.col]] %in% c('train', 'val', 'test')) ]
        )
        # print("no missing id 5")
        # print( length(phe.master$ID[ (phe.master[[split.col]] %in% c('train', 'val', 'test')) ]))
        print(length(phe.master$ID[ (phe.master[[split.col]] %in% c('train', 'val', 'test')) ]))
    }
    if(!is.null(configs[['keep']])){
        # focus on individuals in the specified keep file
        phe.no.missing.IDs <- intersect(phe.no.missing.IDs, readPlinkKeepFile(configs[['keep']]))
    }
    # print("length of no missing")
    # print(length(phe.no.missing.IDs))
    checkMissingPhenoWarning(phe.master, phe.no.missing.IDs)
    # print("small function")
    # print(head(phe.master))
    # print("no missing id")
    # print(head(phe.no.missing.IDs))
    # 
    phe.master[ phe.master$ID %in% phe.no.missing.IDs, ]
}

checkMissingPhenoWarning <- function(phe.master, phe.no.missing.IDs){
  # Show warning message if there are individuals (in phe file)
  # that have (genotype or phenotype) missing values.
    phe.missing.IDs <- phe.master$ID[ ! phe.master$ID %in% phe.no.missing.IDs ]
    if(length(phe.missing.IDs) > 0){
        warning(sprintf(
          'We detected missing values for %d individuals (%s ...).\n',
          length(phe.missing.IDs),
          paste(utils::head(phe.missing.IDs, 5), collapse=", ")
        ))
    }
}
# HAP_ALT_CTS: The number of haploid alternate alleles observed for a given SNP.
# HET_REF_ALT_CTS: The number of heterozygous genotypes observed for a given SNP, where one allele is the reference allele and the other allele is the alternate allele.
# TWO_ALT_GENO_CTS: The number of homozygous alternate genotypes observed for a given SNP.
# OBS_CT: The total number of non-missing genotypes observed for a given SNP.
computeStats <- function(pfile, ids, configs) {
  ID <- original_ID <- NULL  # to deal with "no visible binding for global variable"
  ALT <- MISSING_CT <- OBS_CT <- HAP_ALT_CTS <- HET_REF_ALT_CTS <- TWO_ALT_GENO_CTS <- NULL
  stats_msts <- stats_means <- stats_pNAs <- stats_SDs <- NULL

  keep_f       <- paste0(configs[['gcount.full.prefix']], '.keep')
  gcount_tsv_f <- paste0(configs[['gcount.full.prefix']], '.gcount.tsv')

  dir.create(dirname(configs[['gcount.full.prefix']]), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(gcount_tsv_f)) {
      gcount_df <- data.table::fread(gcount_tsv_f)
  } else {
      # To run plink2 --geno-counts, we write the list of IDs to a file
      # print("ids")
      # print(head(ids))
      dd <- data.frame(ID = ids) %>%
      tidyr::separate(ID, into=c('FID', 'IID'), sep='-')
      # print("head of dd")
      # print(head(dd))
      #write table to keep_f file
      dd%>%data.table::fwrite(keep_f, sep='\t', col.names=F)

      # Run plink2 --geno-counts
      cmd_plink2 <- paste(
          configs[['plink2.path']],
          '--threads', configs[['nCores']],
          '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
          '--keep', keep_f,
          '--out', configs[['gcount.full.prefix']],
          '--memory', configs[['memory']],
          '--geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs'
      )
     
      if (!is.null(configs[['mem']])) cmd_plink2 <- paste(cmd_plink2, '--memory', configs[['mem']])
      #if intern=F, the return is the error code. If intern=T, the return is an R object
      system(cmd_plink2, intern=F, wait=T)

      # read the gcount file
      # stats_means are dosage data (count of AA + count of AB + )
      gcount_df <-
        data.table::fread(paste0(configs[['gcount.full.prefix']], '.gcount')) %>%
        dplyr::rename(original_ID = ID) %>%
        dplyr::mutate(
          ID = paste0(original_ID, '_', ALT),
          stats_pNAs  = MISSING_CT / (MISSING_CT + OBS_CT),
          stats_means = (HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_msts  = (HAP_ALT_CTS + HET_REF_ALT_CTS + 4 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_SDs   = stats_msts - stats_means * stats_means
        )
  }

  out <- list()
  out[["pnas"]]  <- gcount_df %>% dplyr::pull(stats_pNAs)
  out[["means"]] <- gcount_df %>% dplyr::pull(stats_means)
  out[["sds"]]   <- gcount_df %>% dplyr::pull(stats_SDs)

  for(key in names(out)){
      names(out[[key]]) <- gcount_df %>% dplyr::pull(ID)
  }
  out[["excludeSNP"]] <- names(out[["means"]])[(out[["pnas"]] > configs[["missing.rate"]]) | (out[["means"]] < 2 * configs[["MAF.thresh"]])]
  out[["excludeSNP"]] <- out[["excludeSNP"]][ ! is.na(out[["excludeSNP"]]) ]
  out[["excludeSNP"]] <- base::unique(c(configs[["excludeSNP"]], out[["excludeSNP"]]))

  if (configs[['save']]){
      gcount_df %>% data.table::fwrite(gcount_tsv_f, sep='\t')
      saveRDS(out[["excludeSNP"]], file = file.path(dirname(configs[['gcount.full.prefix']]), "excludeSNP.rda"))
  }

  out
}

readBinMat <- function(fhead, configs){
    # This is a helper function to read binary matrix file (from plink2 --variant-score zs bin)
    rows <- data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ', fhead, '.vars.zst'), head=F)$V1
    cols <- data.table::fread(paste0(fhead, '.cols'), head=F)$V1
    bin.reader <- file(paste0(fhead, '.bin'), 'rb')
    M = matrix(
        readBin(bin.reader, 'double', n=length(rows)*length(cols), endian = configs[['endian']]),
        nrow=length(rows), ncol=length(cols), byrow = T
    )
    close(bin.reader)
    colnames(M) <- cols
    rownames(M) <- rows
    if (! configs[['save.computeProduct']]) system(paste(
        'rm', paste0(fhead, '.cols'), paste0(fhead, '.vars.zst'),
        paste0(fhead, '.bin'), sep=' '
    ), intern=F, wait=T)
    M
}

# phe: phenotype data includes env
# env: the string of env
# phe might contain the GxE columns and env column
# vars: a list of snps
# residual: col: sample_id, col: lambda
computeProduct <- function(residual, pfile, vars, phenome, env, stats, configs, iter) {
  ID <- NULL  # to deal with "no visible binding for global variable"
  time.computeProduct.start <- Sys.time()
  snpnetLogger('Start computeProduct()', indent=2, log.time=time.computeProduct.start)

  gc_res <- gc() # garbage collection, used for report memory usage of variables
  if(configs[['KKT.verbose']]) print(gc_res)

  snpnetLogger('Start plink2 --variant-score', indent=3, log.time=time.computeProduct.start)
  dir.create(file.path(configs[['results.dir']], configs[["save.dir"]]), showWarnings = FALSE, recursive = T)
  print(file.path(file.path(configs[['results.dir']], configs[["save.dir"]])))
  
  # file name of residuals
  residual_f <- file.path(configs[['results.dir']], configs[["save.dir"]], paste0("residuals_iter_", iter, ".tsv"))

  # write residuals to a file
  residual_df <- data.frame(residual)
  #the column name of residual_df will lambda_idx_1, lambda_idx_2, ..., ; lambda_idx_n
  colnames(residual_df) <- paste0('lambda_idx_', colnames(residual))

  # if env exists, rbind residual * env vector to the right of the residual data.frame
  # Just stack the new residual dataframe to the old one.
  if(!is.null(env)){
    #convert data.table into data.frame
   
    #append res_dot_e with res 
    phenome.df <- as.data.frame(phenome)
    env.vec <- phenome.df[, c(env), drop=F] 
    rownames(env.vec) <- phenome.df$ID # row_name is sample ID (format: xxx_xxx) which is the same as residual_df.
    # residual column is lambda_idx_1, lambda_idx_2
    # merge env vector with residudal_df (sample_size * lambda) by rowname
    residual_df.with_ids <- residual_df
    residual_df.with_ids$ID <- rownames(residual_df)
    env.vec.with_ids <- env.vec
    env.vec.with_ids$ID <- rownames(env.vec)
    res_and_e <- left_join(residual_df.with_ids, env.vec.with_ids, by="ID")
    rownames(res_and_e) <- res_and_e$ID
    res_and_e <- subset(res_and_e, select = -c(ID))
    # multiply residual and env
    # extract env, sample_size x 1
    env_vec  <- res_and_e[, c(env)]
    # extract residual, sample_size x num.lams
    res_df = res_and_e[, !(colnames(res_and_e) %in% c(env))]
    
    # this is a vector (numeric) or matrix(data.frame)
    res_dot_e = env_vec * res_df
    #just directly append, so the membership is 11112222
    # [r(l)^T (r(l) * E)^T]
    #res_dot_e: n x lambda.num
    # res_and_e.colnames<-colnames(res_and_e)
    # lam.ids <- res_and_e.colnames[1:(length(res_and_e.colnames)-1)]
    if(is.data.frame(res_dot_e)){
      rownames(res_dot_e) <- rownames(res_and_e) #sample ids
      colnames(res_dot_e) <- paste0(colnames(residual_df), '_E')
      res_and_res_dot_e<-cbind(res_df, res_dot_e)
    }else{
      # if res_dot_e is numeric
      # names(res_dot_e) <- paste0(lam.ids, '_E')
      # assign row names to res_dot_e
      names(res_dot_e) <- rownames(res_and_e)
      res_and_res_dot_e<-cbind(res_df, res_dot_e)
      colnames(res_and_res_dot_e) <- c(colnames(residual_df), paste0(colnames(residual_df), '_E'))
    }
    res_and_res_dot_e <- as.data.frame(res_and_res_dot_e)
    #output res and res x E data to residual_f
    res_and_res_dot_e %>%
      tibble::rownames_to_column("ID") %>%
      tidyr::separate(ID, into=c('#FID', 'IID'), sep='-') %>%
      data.table::fwrite(residual_f, sep='\t', col.names=T)
  }else{
    residual_df %>%
    tibble::rownames_to_column("ID") %>%
    tidyr::separate(ID, into=c('#FID', 'IID'), sep='-') %>%
    data.table::fwrite(residual_f, sep='\t', col.names=T)
  }

  # Run plink2 --geno-counts
  #This place should consider gxe term, which will be replaced
  cmd_plink2 <- paste(
    configs[['plink2.path']],
    '--threads', configs[['nCores']],
    '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
    '--read-freq', paste0(configs[['gcount.full.prefix']], '.gcount'),
    '--keep', residual_f,
    '--memory', configs[['memory']],
    '--out', stringr::str_replace_all(residual_f, '.tsv$', ''),
    '--variant-score', residual_f, 'zs', 'bin'
  )

  if (!is.null(configs[['mem']])) {
    cmd_plink2 <- paste(cmd_plink2, '--memory', as.integer(configs[['mem']]) - ceiling(sum(as.matrix(gc_res)[,2])))
  }

  system(cmd_plink2, intern=F, wait=T)
  
  # each row is a snp and each col is lambda or lambda * E
  prod.full <- readBinMat(stringr::str_replace_all(residual_f, '.tsv$', '.vscore'), configs)
  if (! configs[['save.computeProduct']] ) system(paste(
      'rm', residual_f, stringr::str_replace_all(residual_f, '.tsv$', '.log'), sep=' '
  ), intern=F, wait=T)

  snpnetLoggerTimeDiff('End plink2 --variant-score.', time.computeProduct.start, indent=4)

  #col represents residuals(lambda) (and residuals(lambda)*E if it is SGL or GL), row represents variants.
  rownames(prod.full) <- vars
  if (configs[["standardize.variant"]]) {
      for(residual.col in 1:ncol(residual)){
        # prod.full[, residual.col] <- apply(prod.full[, residual.col], 1, "/", stats[["sds"]])
        prod.full[, residual.col] <- prod.full[, residual.col] / stats[['sds']]
      }
  }
  prod.full[stats[["excludeSNP"]], ] <- NA
  snpnetLoggerTimeDiff('End computeProduct().', time.computeProduct.start, indent=3)
  prod.full
}

rearrange_prod <- function(prod.full){
  # rearrange product matrix from p by 2*m into 2*n by m
  #row.names are group names.
  row.names<-rownames(prod.full)
  
  # Repeat each element twice
  row.names.new <- rep(row.names, each = 2)
  
  # Append '_E' to every second element
  row.names.new[seq(2, length(row.names.new), 2)] <- paste0(row.names.new[seq(2, length(row.names.new), 2)], '_E')
  
  # transpose start
  p <- nrow(prod.full)
  m2 <- ncol(prod.full)
  m <- m2 / 2
  left <- prod.full[, 1:m]
  right <- prod.full[, (m+1):m2]
  prod.full.trans <- matrix(NA, nrow = 2*p, ncol = m)
  prod.full.trans[seq(1, 2*p, 2), ] = left # odd
  prod.full.trans[seq(2, 2*p, 2), ] = right
  # prod.full.trans <- matrix(prod.full, nrow = 2*p, ncol = m, byrow = TRUE)
  rownames(prod.full.trans) <- row.names.new
  colnames(prod.full.trans) <- colnames(prod.full)[seq(1, m, 1)]
  # prod.full.trans = as.data.frame(prod.full.trans)
  prod.full.trans
}

soft_threshold <- function(x, tau, lambda){
  # x: row: G and GxE, col: lambdas
  # lambda: The lambda value(s) corresponds to x.
  # ncol(x) should be the same as number of lambda
  sgn <- sign(x) #sign of x
  abs_x <- abs(x) #absolute value of x
  #pmax() is used to take the element-wise maximum between abs_x - tau and 0.
  #(x)+ = max{x, 0}
  lambda.prop <- lambda * tau
  diff <- sweep(x = abs_x, MARGIN=2, STATS=lambda.prop, FUN='-')
  return (sgn* pmax(diff, 0))
}

# Get the score calculated by the left part of the strong rule inequation
#lambda.last: the lambda value(s) from the last iteration
getScore <- function(prod.full){
  # We hope the row of prod.full is: rsid1, rsid1_E, rsid2, rsid2_E
  row.names <- rownames(prod.full)
  m2 <- nrow(prod.full)
  if(!identical(row.names[seq(2,m2,2)], paste0(row.names[seq(1, m2, 2)], "_E"))){
    v_sort <- sub("_E", "", row.names)
    v_sort.uniq <- unique(v_sort)
    # Repeat each element twice
    v_new <- rep(v_sort.uniq, each = 2)
    # Append '_E' to every second element
    
    v_new[seq(2, length(v_new), 2)] <- paste0(v_new[seq(2, length(v_new), 2)], '_E')
    prod.full <- prod.full[v_new, , drop=FALSE]
  } 
  #calculate L2 norms
  # The order might not the same
  prod.full.lambdas <- prod.full[seq(1, m2, 2), , drop=FALSE]
  prod.full.lambdaEs <- prod.full[seq(2, m2, 2), , drop=FALSE]
  
  row.name = rownames(prod.full.lambdas)
  col.name = colnames(prod.full.lambdas)
  # if(class(prod.full) == "data.frame"){
  #   prod.full.lambdas = as.matrix(prod.full.lambdas)
  #   prod.full.lambdaEs = as.matrix(prod.full.lambdaEs)
  # }
  # Ensure prod.full.lambdas and prod.full.lambdaEs have the same dimensions:
  score = NULL
  if(identical(dim(prod.full.lambdas), dim(prod.full.lambdaEs))){
    score <- sqrt(prod.full.lambdaEs^2 + prod.full.lambdas^2)
    rownames(score) = row.name
    colnames(score) = col.name
  }else{
    if (configs[['verbose']]) snpnetLogger("  [Error]: Matrices do not have the same dimensions ...")
  }
  score 
}

getScore2 <- function(prod.full, tau, lambda){
  # We hope the row of prod.full is: rsid1, rsid1_E, rsid2, rsid2_E
  row.names <- rownames(prod.full)
  m2 <- nrow(prod.full)
  if(!identical(row.names[seq(2,m2,2)], paste0(row.names[seq(1, m2, 2)], "_E"))){
    v_sort <- sub("_E", "", row.names)
    v_sort.uniq <- unique(v_sort)
    # Repeat each element twice
    v_new <- rep(v_sort.uniq, each = 2)
    # Append '_E' to every second element
    
    v_new[seq(2, length(v_new), 2)] <- paste0(v_new[seq(2, length(v_new), 2)], '_E')
    prod.full <- prod.full[v_new, , drop=FALSE]
  } 
  #soft thresholding the products
  prod.full <- soft_threshold(prod.full, tau, lambda)
  #calculate L2 norms
  # The order might not the same
  prod.full.lambdas <- prod.full[seq(1, m2, 2), , drop=FALSE]
  prod.full.lambdaEs <- prod.full[seq(2, m2, 2), , drop=FALSE]
  
  row.name = rownames(prod.full.lambdas)
  # print(row.name[1:10])
  col.name = colnames(prod.full.lambdas)
  # if(class(prod.full) == "data.frame"){
  #   prod.full.lambdas = as.matrix(prod.full.lambdas)
  #   prod.full.lambdaEs = as.matrix(prod.full.lambdaEs)
  # }
  # Ensure prod.full.lambdas and prod.full.lambdaEs have the same dimensions:
  score = NULL
  if(identical(dim(prod.full.lambdas), dim(prod.full.lambdaEs))){
    score <- sqrt(prod.full.lambdaEs^2 + prod.full.lambdas^2)
    # print("kkkkkkkkkkkkkkkkkkkkkkkk_print_score_shape1.2.2")
    # print(dim(score))
    # print(length(score))
    # print(class(score))
    rownames(score) = row.name
    colnames(score) = col.name
    # print(head(score))
  }else{
    if (configs[['verbose']]) snpnetLogger("  [Error]: Matrices do not have the same dimensions ...")
  }
  score
}
# train is the training data, it is always phe[['train']]

KKT.check <- function(residual, pfile, vars, train, env, n.train, 
                      current.lams, prev.lambda.idx, 
                      stats, fitted.model, configs, iter, p.factor=NULL, alpha = NULL, tau = NULL, isSGL=FALSE, isGL=FALSE) {
  time.KKT.check.start <- Sys.time()
  if (is.null(alpha)) alpha <- 1
  if (is.null(tau)) tau <- 1
  if (configs[['KKT.verbose']]) snpnetLogger('Start KKT.check()', indent=1, log.time=time.KKT.check.start)
  
  prod.full <- computeProduct(residual, pfile, vars, train, env, stats, configs, iter) / n.train
  prod.full <- rearrange_prod(prod.full)
  scores <- getScore2(prod.full, tau, current.lams) 
  
  if(!is.null(p.factor)){
    prod.full <- sweep(prod.full, 1, p.factor, FUN="/")
  }
  
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- computeProduct.', indent=2, start.time=time.KKT.check.start)
  num.lams <- length(current.lams)
  
  # get strong features (features: G and GxE)
  # strong features: strong group strong features
  if (length(configs[["covariates"]]) + length(configs[['env']]) > 0) {
    strong.vars <- match(rownames(fitted.model$beta[-(1:(length(configs[["env"]]) + length(configs[['covariates']]))), , drop = FALSE]), rownames(prod.full))
    strong.vars.names <- rownames(fitted.model$beta[-(1:(length(configs[["env"]]) + length(configs[['covariates']]))), , drop = FALSE])
    
  } else {
    strong.vars <- match(rownames(fitted.model$beta), rownames(prod.full))
    strong.vars.names <- rownames(fitted.model$beta)
  }
  strong.groups.names <- unique(gsub('_E', '', strong.vars.names))
  
  # get strong group weak features. Strong group weak features
  # are the features whose prefix only appear once in fitted models
  
  strong.groups.weak.feat.idx <- NULL
  strong.groups.weak.feat.names <- NULL
  freq <- table(gsub("_E", "", strong.vars.names))
  if(length(names(freq)[freq == 1]) > 0){
    strong.groups.weak.feat.names <- grep(pattern = paste(names(freq)[freq == 1], collapse = "|"), 
                                          strong.vars.names, 
                                          value = TRUE)
    sgwf.withE <- endsWith(strong.groups.weak.feat.names, "_E") # BOOLEAN
    sgwf.withoutE <- as.logical(FALSE ^ sgwf.withE)
    strong.groups.weak.feat.names[sgwf.withE] <- gsub('_E', '', strong.groups.weak.feat.names[sgwf.withE])
    strong.groups.weak.feat.names[sgwf.withoutE] <- paste0(strong.groups.weak.feat.names[sgwf.withoutE], '_E')
    # idx of strong group weak features in prod.full
    strong.groups.weak.feat.idx <- match(strong.groups.weak.feat.names, rownames(prod.full))
  }
  
  # get names of all features
  all.vars.names <- rownames(prod.full)
  # get names of all groups
  all.groups.names <- unique(gsub('_E', '', all.vars.names))
  
  # get names of weak groups
  weak.groups.names <- setdiff(all.groups.names, strong.groups.names)
  
  # get names of weak group weak features
  weak.groups.weak.feat.names <- rep(weak.groups.names, each = 2)
  weak.groups.weak.feat.names[seq(2, length(weak.groups.weak.feat.names), 2)] <- paste0(weak.groups.weak.feat.names[seq(2, length(weak.groups.weak.feat.names), 2)], '_E')
  
  # idx of weak group weak features in prod.full
  weak.groups.weak.feat.idx <- match(weak.groups.weak.feat.names, all.vars.names) # the length this is 2 * length of weak group
  
  # idx of weak features. Weak features are from strong groups or weak groups.
  weak.vars <- setdiff(1:nrow(prod.full), strong.vars)
  
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- strong.vars.', indent=2, start.time=time.KKT.check.start)
  
  # get the coefficients from strong variables.
  if (length(configs[["covariates"]] + length(configs[['env']])) > 0) {
      strong.coefs <- fitted.model$beta[-(1:(length(configs[['env']]) + length(configs[["covariates"]]))), , drop = FALSE]
  } else {
      strong.coefs <- fitted.model$beta
  }
  # if (!(isSGL) & !(isGL)){
  #   prod.full[strong.vars, ] <- prod.full[strong.vars, , drop = FALSE] - (1-alpha) * as.matrix(strong.coefs) *
  #     matrix(current.lams, nrow = length(strong.vars), ncol = length(current.lams), byrow = T)
  # }
  
  # construct comparison matrix, which is the upperbound of
  # strong rules inequation for KKT checking
  
  if (configs[['KKT.check.aggressive.experimental']]) {
      # An approach to address numerial precision issue.
      # We do NOT recommended this procedure
    # prod.strong <- prod.full[strong.vars, , drop = FALSE]
    # max.abs.prod.strong <- apply(abs(prod.strong), 2, max, na.rm = T)
    # mat.cmp <- matrix(max.abs.prod.strong, nrow = length(weak.vars), ncol = length(current.lams), byrow = T)
  } else {
    
      # upperbound for weak group level checking; lambda * (1-tau) * sqrt(w); w: group length
      mat.cmp.weakGroup <- matrix(current.lams * (1 - tau) * sqrt(2), 
                              nrow = length(weak.groups.names), 
                              ncol = length(current.lams), 
                              byrow = T)  
      # upperbound for weak group's weak features checking, lambda * tau
      mat.cmp.weakGroupWeakFeat  <-  matrix(current.lams * tau, 
                                           nrow = length(weak.groups.weak.feat.idx), 
                                           ncol = length(current.lams), 
                                           byrow = T) 
      # upperbound fro strong group's weak features checking, lambda * tau
      if(length(strong.groups.weak.feat.idx) > 0){
          mat.cmp.strongGroupWeakFeat <- matrix(current.lams * tau, 
                                                nrow = length(strong.groups.weak.feat.idx), 
                                                ncol=length(current.lams), 
                                                byrow=T)
      }
    
    }
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- mat.cmp.', indent=2, start.time=time.KKT.check.start)
  
  # <1>. Construct group level weak features violation checking matrix
  # Left part of weak group level checking
  weakGroup.score <- getScore2(prod.full[weak.groups.weak.feat.idx, , drop=FALSE], tau, current.lams)
  snpnetLogger(paste0("dimension of weak group score", str(dim(weakGroup.score))))
  # weak group level checking; Count violations
  violates.weakGroups <- weakGroup.score - mat.cmp.weakGroup > 0 # violate=TRUE, not violate=FALSE
  rownames(violates.weakGroups) <- rownames(weakGroup.score)
  # extend violation matrix from group level to G and GxE level
  violates.GxE <- violates.weakGroups
  rownames(violates.GxE) <- paste0(rownames(violates.weakGroups), "_E")
  # construct group level weak features violation checking matrix
  violates.features.weakGroupLevel <- rbind(violates.weakGroups, violates.GxE)
  
  # <2>. Construct weak group weak features (wgwf) violations.
  violates.features.wgwf <- (abs(prod.full[weak.groups.weak.feat.idx, , drop=FALSE]) - mat.cmp.weakGroupWeakFeat > 0 )
  rownames(violates.features.wgwf) <- weak.groups.weak.feat.names
  
  # <3>. Construct strong group weak features (sgwf) violations.
  if(length(strong.groups.weak.feat.idx) > 0){
    violates.features.sgwf <- (abs(prod.full[strong.groups.weak.feat.idx, , drop=FALSE]) - mat.cmp.strongGroupWeakFeat > 0 )
    rownames(violates.features.sgwf) <- strong.groups.weak.feat.names
    
    # <4>. Construct strong group violations.
    # If sgwf violates, its group also violates 
    # This violation is for strong group level
    strong.group.withWeakFeat <- gsub('_E', '', strong.groups.weak.feat.names)
    strong.group.withWeakFeatByE <- paste0(strong.group.withWeakFeat, '_E')
    violates.strongGroup.1 <- violates.features.sgwf
    rownames(violates.strongGroup.1) <- strong.group.withWeakFeat
    violates.strongGroup.2 <- violates.features.sgwf
    rownames(violates.strongGroup.2) <- strong.group.withWeakFeatByE
    violates.features.strongGroupLevel <- rbind(violates.strongGroup.1, violates.strongGroup.2)
  }
  
  
  # Put those matrix together. Construct prod.full level BOOOLEAN matrix
  violate.group <- matrix(FALSE, nrow=nrow(prod.full), ncol=length(current.lams))
  violate.feat <- matrix(FALSE, nrow=nrow(prod.full), ncol=length(current.lams))
  
  violate.group[match(rownames(violates.features.weakGroupLevel), rownames(prod.full)),] = violates.features.weakGroupLevel
  violate.feat[match(rownames(violates.features.wgwf), rownames(prod.full)),] = violates.features.wgwf
  
  if(length(strong.groups.weak.feat.idx) > 0){
    violate.group[match(rownames(violates.features.strongGroupLevel), rownames(prod.full)),] = violates.features.strongGroupLevel
    violate.feat[match(rownames(violates.features.sgwf), rownames(prod.full)),] = violates.features.sgwf
  }
  
  # violates.features.temp1 <-  violates.features.temp1[match(rownames(violates.features.temp2),
  #                                                           rownames(violates.features.temp1)),]
  # The matrix follows the following logical operation
  # TRUE, TRUE -> TRUE
  # FALSE,TRUE->FALSE
  # TRUE,FALSE->FALSE
  # FALSE,FALSE->FALSE
  
  #num.violates represent number of features violated, not group.
  
  num.violates <- apply(violate.group & violate.feat, 2, function(x) sum(x, na.rm = T))
  print("print count of violations")
  print(num.violates)
  idx.violation <- which((num.violates != 0) & ((1:num.lams) >= prev.lambda.idx))
  print("idx.violation")
  print(idx.violation)
  # look for max valid lambda idx
  max.valid.idx <- ifelse(length(idx.violation) == 0, num.lams, min(idx.violation) - 1)
  print("max.valid.idx")
  print(max.valid.idx)
 
  if (max.valid.idx > 0) {
    if(isSGL | isGL){
      # print("print max.valid.idx")
      # print(max.valid.idx)
      # print("print prod.full dimension")
      # print(dim(prod.full))
      # print("values of the matrix")
      # print(prod.full[, max.valid.idx, drop=FALSE][1:10,])
      # currently prod.full is just a vector
      score.tmp <- getScore2(prod.full[, max.valid.idx, drop=FALSE], tau, current.lams[max.valid.idx])
      score <- score.tmp[, 1]
      names(score) <- rownames(score.tmp)
      # print(head(score))
      
    }
    else{
      score <- abs(prod.full[, max.valid.idx])
    }
   
  } else {
    score <- NULL
  }
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- score.', indent=2, start.time=time.KKT.check.start)

  out <- list(max.valid.idx = max.valid.idx, score = score)

  if (configs[['KKT.verbose']]) {
    # this is just used for output something
    gene.names <- rownames(prod.full) # if model is GL/SGL, rownames are gene names and GxE names
    strong.names <- rownames(strong.coefs) # Strong features (include gene or GxE), which are used for fitting the model
    active <- matrix(FALSE, nrow(prod.full), num.lams)
    active[match(strong.names, gene.names), ] <- as.matrix(strong.coefs != 0) # strong.coefs row is G+GxE, col is lambda
    inactive <- matrix(FALSE, nrow(prod.full), num.lams)
    inactive[match(strong.names, gene.names), ] <- as.matrix(strong.coefs == 0)
    # active matrix and inactive matrix are inverse.

    prod.strong <- prod.full[strong.vars, , drop = FALSE] #strong.vars: position of strong vars (var: G or GxE)
    prod.weak <- prod.full[weak.vars, , drop = FALSE]

    min.abs.prod.active <- apply(abs(prod.full*active), 2, function(x) min(x[x > 0], na.rm = T))
    max.abs.prod.active <- apply(abs(prod.full*active), 2, max, na.rm = T)
    max.abs.prod.inactive <- apply(abs(prod.full*inactive), 2, max, na.rm = T)
    max.abs.prod.strong <- apply(abs(prod.strong), 2, max, na.rm = T)
    max.abs.prod.weak <- apply(abs(prod.weak), 2, max, na.rm = T)

    print(data.frame(
      lambda = current.lams,
      num.active = apply(active, 2, sum, na.rm = T),
      min.abs.prod.active = min.abs.prod.active,
      max.abs.prod.active = max.abs.prod.active,
      num.inactive = apply(inactive, 2, sum, na.rm = T),
      max.abs.prod.inactive = max.abs.prod.inactive,
      max.abs.prod.strong = max.abs.prod.strong,
      max.abs.prod.weak = max.abs.prod.weak,
      num.violates = num.violates
    ))
  }
  out
}

setDefaultMetric <- function(family){
    if (family == "gaussian") {
        metric <- 'r2'
    } else if (family == "binomial") {
        metric <- 'auc'
    } else if (family == "cox") {
        metric <- 'C'
    } else {
        stop(paste0('The specified family (', family, ') is not supported!'))
    }
    metric
}

#start here...
computeMetric <- function(pred, response, metric.type) {
    if (metric.type == 'r2') {
        metric <- 1 - apply((response - pred)^2, 2, sum) / sum((response - mean(response))^2)
    } else if (metric.type == 'auc') {
        metric <- apply(pred, 2, function(x) {
            pred.obj <- ROCR::prediction(x, factor(response))
            auc.obj <- ROCR::performance(pred.obj, measure = 'auc')
            auc.obj@y.values[[1]]
        })
    } else if (metric.type == 'd2') {
        d0 <- glmnet::coxnet.deviance(NULL, response)
        metric <- apply(pred, 2, function(p) {
            d <- glmnet::coxnet.deviance(p, response)
            1 - d/d0
        })
    } else if (metric.type == 'C'){
      metric <- apply(pred, 2, function(p) {
        cindex::CIndex(p, response[,1], response[,2])
      })
    }
    metric
}
computeMetricPRS <- function(env_col, prs_g, prs_ge, response, metric.type) {
  model<-lm(response ~ env_col + prs_g + prs_ge:env_col)
  fit <- summary(model)
  stderr<-fit$coefficients["env_col:prs_ge", "Std. Error"]
  tval <- fit$coefficients["env_col:prs_ge","t value"]
  #This is loge!!! I need to double check whether this value is correct or not next time!!
  # negLog10P <- -2 * pt(-abs(tval), df = model$df.residual, log.p = TRUE)
  negLog10P = -(log(2) + pt(abs(coef(fit)["env_col:prs_ge","t value"]), model$df.residual, lower.tail = FALSE, log.p = TRUE))/log(10)
  # print(log_pval)
  metric <- fit$r.squared
  # print(fit$coefficients)
  # pval <- -log10(fit$coefficients[4,4])
  return(c(metric, negLog10P))
}
checkEarlyStopping <- function(metric.val, max.valid.idx, iter, configs){
    if (max.valid.idx <= configs[['stopping.lag']]) {
      earlyStop <- FALSE
    } else {
      max.valid.idx.lag <- max.valid.idx-configs[['stopping.lag']]
      max.val.1 <- max(metric.val[1:(max.valid.idx.lag)])
      max.val.2 <- max(metric.val[(max.valid.idx.lag+1):max.valid.idx])
      snpnetLogger(sprintf('stopping lag=%g, max.val.1=%g max.val.2=%g', max.valid.idx.lag, max.val.1, max.val.2))
      if (
          (configs[['early.stopping']]) &&
          (max.valid.idx > configs[['stopping.lag']]) &&
          (max.val.1 > max.val.2)
      ) {
          snpnetLogger(sprintf(
              "Early stopped at iteration %d (Lambda idx=%d ) with validation metric: %.14f.",
              iter, which.max(metric.val), max(metric.val, na.rm = T)
          ))
          snpnetLogger(paste0(
              "Previous ones: ",
              paste(metric.val[(max.valid.idx-configs[['stopping.lag']]+1):max.valid.idx], collapse = ", "),
              "."
          ), indent=1)
          earlyStop <- TRUE
      } else {
          earlyStop <- FALSE
      }
    }
    earlyStop
}

cleanUpIntermediateFiles <- function(configs){
    for(subdir in c(configs[["save.dir"]], configs[["meta.dir"]])){
        system(paste(
            'rm', '-rf', file.path(configs[['results.dir']], subdir), sep=' '
        ), intern=F, wait=T)
    }
}

computeCoxgrad <- function(glmfits, time, d){
    apply(glmfits, 2, function(f){coxgrad(f,time,d,w=rep(1,length(f)))})
}

simplifyList_Col <- function(x) {
  sample <- x[[1L]]
  if (is.matrix(sample)) {
    x <- do.call(rbind, x)
  } else {
    x <- unlist(x)
  }
  return(x)
}

checkGlmnetPlus <- function(use.glmnetPlus, family) {
    if (!requireNamespace("glmnet") && !requireNamespace("glmnetPlus"))
        stop("Please install at least glmnet or glmnetPlus.")
    if(is.null(use.glmnetPlus))
        use.glmnetPlus <- (family == "gaussian")
    if(use.glmnetPlus){
        if (!requireNamespace("glmnetPlus")) {
            warning("use.glmnetPlus was set to TRUE but glmnetPlus not found... Revert back to glmnet.")
            use.glmnetPlus <- FALSE
        }
    }
    use.glmnetPlus
}

#revised by Le Huang, add env param in the code
setupConfigs <- function(configs, genotype.pfile, phenotype.file, phenotype, env, covariates, alpha, nlambda, split.col, p.factor, status.col, mem){
    out.args <- as.list(environment())
    defaults <- list(
        missing.rate = 0.1,
        MAF.thresh = 0.001,
        nCores = 1,
        glmnet.thresh = 1e-07,
        nlams.init = 10,
        nlams.delta = 5,
        num.snps.batch = 1000,
        num.groups.batch=500,
        vzs=TRUE, # geno.pfile vzs
        increase.size = NULL,
        standardize.variant = FALSE,
        early.stopping = TRUE,
        stopping.lag = 2,
        niter = 50,
        keep = NULL,
        lambda.min.ratio = NULL,
        KKT.verbose = FALSE,
        use.glmnetPlus = FALSE,
        save = FALSE,
        save.computeProduct = FALSE,
        prevIter = 0,
        results.dir = NULL,
        meta.dir = 'meta',
        save.dir = 'results',
        verbose = FALSE,
        KKT.check.aggressive.experimental = FALSE,
        gcount.basename.prefix = 'snpnet.train',
        gcount.full.prefix=NULL,
        endian="little",
        metric=NULL,
        plink2.path='plink2',
        zstdcat.path='zstdcat',
        zcat.path='zcat',
        rank = TRUE,
        excludeSNP = NULL,
        memory = 120000
    )
    out <- defaults

    # store additional params
    for (name in setdiff(names(out.args), "configs")) {
      out[[name]] <- out.args[[name]]
    }

    # update the defaults with the specified parameters and keep redundant parameters from configs
    for (name in names(configs)) {
        out[[name]] <- configs[[name]]
    }

    # update settings
    out[["early.stopping"]] <- ifelse(out[["early.stopping"]], out[['stopping.lag']], -1)
    if(is.null(out[['increase.size']]))  out[['increase.size']] <- out[['num.snps.batch']]/2

    # configure the temp file locations
    #   We will write some intermediate files to meta.dir and save.dir.
    #   those files will be deleted with snpnet::cleanUpIntermediateFiles() function.
    if (is.null(out[['results.dir']])) out[['results.dir']] <- tempdir(check = TRUE)
    dir.create(file.path(out[['results.dir']], out[["meta.dir"]]), showWarnings = FALSE, recursive = T)
    dir.create(file.path(out[['results.dir']], out[["save.dir"]]), showWarnings = FALSE, recursive = T)
    if(is.null(out[['gcount.full.prefix']])) out[['gcount.full.prefix']] <- file.path(
        out[['results.dir']], out[["meta.dir"]], out['gcount.basename.prefix']
    )

    out
}

updateConfigsWithFamily <- function(configs, family){
    out <- configs
    out[['family']] <- family
    out[['use.glmnetPlus']] <- checkGlmnetPlus(out[['use.glmnetPlus']], family)
    if (is.null(out[['metric']])) out[['metric']] <- setDefaultMetric(family)
    out
}

## logger functions

snpnetLogger <- function(message, log.time = NULL, indent=0, funcname='snpnet-ge'){
    if (is.null(log.time)) log.time <- Sys.time()
    cat('[', as.character(log.time), ' ', funcname, '] ', rep(' ', indent * 2), message, '\n', sep='')
}

timeDiff <- function(start.time, end.time = NULL) {
    if (is.null(end.time)) end.time <- Sys.time()
    paste(round(end.time-start.time, 4), units(end.time-start.time))
}

snpnetLoggerTimeDiff <- function(message, start.time, end.time = NULL, indent=0){
    if (is.null(end.time)) end.time <- Sys.time()
    snpnetLogger(paste(message, "Time elapsed:", timeDiff(start.time, end.time), sep=' '), log.time=end.time, indent=indent)
}


#' @importFrom methods as
as_dgCMatrix <- function(x) {
  as(as(x, "sparseMatrix"), "CsparseMatrix")
}

predict_sparsegl <- function(
    object, newx, s = NULL,
    type = c("link", "response", "coefficients", "nonzero", "class"),
    ...) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE)) {
      cli_abort(
        "You must supply a value for `newx` when `type` == '{type}'."
      )
    }
  }
  nbeta <- coef(object, s)
  if (type == "coefficients") {
    return(nbeta)
  }
  if (type == "nonzero") {
    return(nonzeroCoef(nbeta[-1, , drop = FALSE]))
  }
  if (inherits(newx, "sparseMatrix")) newx <- as_dgCMatrix(newx)
  dx <- dim(newx)
  p <- object$dim[1]
  if (is.null(dx)) newx <- matrix(newx, 1, byrow = TRUE)
  if (ncol(newx) != p) {
    cli_abort("The number of variables in `newx` must be {p}.")
  }
  #fit <- as.matrix(cbind2(1, newx) %*% nbeta)
  one_vec <- as.matrix(rep(1, nrow(newx)), nrow(newx), 1)
  fit <- as.matrix(newx %*% as.matrix(nbeta)[-1, , drop = FALSE] + one_vec %*% as.matrix(nbeta)[1, , drop = FALSE])
  fit
}


# Re-implemented the sparsegl:::calc_gamma function using Rcpp. 
# This modified version does not work for "ncols > 2" in the original code. 
Rcpp::cppFunction('
  NumericVector calc_gamma_Rcpp(NumericMatrix x, IntegerVector ix, IntegerVector iy, int bn) {
    // Helper function to calculate largest squared singular value
    auto maxeig2 = [](NumericMatrix& mat, int col1, int col2) {
      double mat11 = 0.0, mat12 = 0.0, mat22 = 0.0;

      // Cross-product of the matrix (mat\' * mat)
      for (int i = 0; i < mat.nrow(); ++i) {
        mat11 += mat(i, col1) * mat(i, col1);
        mat12 += mat(i, col1) * mat(i, col2);
        mat22 += mat(i, col2) * mat(i, col2);
      }

      // Trace and determinant
      double tr = mat11 + mat22;
      double det = mat11 * mat22 - mat12 * mat12;

      // Return the largest eigenvalue
      return (tr + sqrt(tr * tr - 4 * det)) / 2;
    };

    NumericVector gamma(bn);
    
    for (int g = 0; g < bn; ++g) {
      int start = ix[g] - 1; // Adjust 1-based index from R to 0-based for C++
      int end = iy[g] - 1;
      int ncols = end - start + 1;

      if (ncols == 2) {
        // Compute largest singular value directly without creating a submatrix
        gamma[g] = maxeig2(x, start, end);
      } else {
        // Sum of squares of the selected columns
        double sumsq = 0.0;
        for (int i = 0; i < x.nrow(); ++i) {
          for (int j = start; j <= end; ++j) {
            sumsq += x(i, j) * x(i, j);
          }
        }
        gamma[g] = sumsq;
      }
    }
    
    return gamma / x.nrow();
  }
')

integer_dc = dotCall64:::integer_dc
numeric_dc = dotCall64:::numeric_dc
getoutput = sparsegl:::getoutput

#' @importFrom dotCall64 integer_dc vector_dc numeric_dc
sgl_ls_modified <- function(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlam,
    flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
    lower_bnd, upper_bnd) {
  # call Fortran core
  is.sparse <- FALSE
  if (!is.numeric(y)) cli_abort("For family = 'gaussian', y must be numeric.")
  if (inherits(x, "sparseMatrix")) {
    is.sparse <- TRUE
    x <- as_dgCMatrix(x)
  }
  ym <- mean(y)
  if (intr) {
    y <- y - ym
    nulldev <- mean(y^2)
  } else {
    nulldev <- mean((y - ym)^2)
  }
  if (standardize) {
    sx <- sqrt(Matrix::colSums(x^2))
    sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
    xs <- 1 / sx
    x <- x %*% Matrix::Diagonal(x = xs)
  }
  if (is.sparse) {
    xidx <- as.integer(x@i + 1)
    xcptr <- as.integer(x@p + 1)
    xval <- as.double(x@x)
    nnz <- as.integer(utils::tail(x@p, 1))
  }

  gamma <- calc_gamma_Rcpp(x, ix, iy, bn)

  if (!is.sparse) {
    fit <- dotCall64::.C64(
      "sparse_four",
      SIGNATURE = c(
        "integer", "integer", "integer", "integer", "double",
        "integer", "integer", "double", "double", "double", "double",
        "integer", "integer", "integer", "double", "double",
        "double", "integer", "integer", "integer", "double", "double", "integer",
        "integer", "double", "integer", "integer", "double",
        "double", "double", "double"
      ),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(x), y = as.double(y), pf = pf,
      pfl1 = pfl1,
      # Read / write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin, ulam = ulam,
      eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam),
      activeGroup = integer_dc(pmax), nbeta = integer_dc(nlam),
      alam = numeric_dc(nlam), npass = integer_dc(1),
      jerr = integer_dc(1), mse = numeric_dc(nlam),
      # read only
      alsparse = asparse, lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 11), rep("rw", 8), rep("w", 9), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl"
    )
  } else { # sparse design matrix
    fit <- dotCall64::.C64(
      "spmat_four",
      SIGNATURE = c(
        "integer", "integer", "integer", "integer", "double",
        "integer", "integer", "double", "integer", "integer",
        "integer", "double", "double", "double", "integer", "integer",
        "integer", "double", "double", "double", "integer",
        "integer", "integer", "double", "double", "integer",
        "integer", "double", "integer", "integer", "double",
        "double", "double", "double"
      ),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(xval), xidx = xidx, xcptr = xcptr,
      nnz = nnz, y = as.double(y), pf = pf, pfl1 = pfl1,
      # Read write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin,
      ulam = ulam, eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam), activeGroup = integer_dc(pmax),
      nbeta = integer_dc(nlam), alam = numeric_dc(nlam),
      npass = integer_dc(1), jerr = integer_dc(1), mse = numeric_dc(nlam),
      # Read only
      alsparse = as.double(asparse), lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 14), rep("rw", 8), rep("w", 9), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl"
    )
  }

  # output
  outlist <- getoutput(x, group, fit, maxit, pmax, nvars, vnames, eps)
  if (standardize) outlist$beta <- outlist$beta * xs
  if (intr) {
    outlist$b0 <- outlist$b0 + ym
  } else {
    outlist$b0 <- rep(0, dim(outlist$beta)[2])
  }
  outlist$npasses <- fit$npass
  outlist$jerr <- fit$jerr
  outlist$group <- group
  outlist$mse <- fit$mse[seq(fit$nalam)]
  outlist$dev.ratio <- 1 - outlist$mse / nulldev
  outlist$nulldev <- nulldev
  class(outlist) <- c("lsspgl")
  outlist
}


#' @importFrom stats glm binomial gaussian
sgl_logit_modified <- function(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1,
    dfmax, pmax, nlam, flmin, ulam, eps,
    maxit, vnames, group, intr, asparse, standardize,
    lower_bnd, upper_bnd) {
  y <- as.factor(y)
  lev <- levels(y)
  ntab <- table(y)
  minclass <- min(ntab)
  if (minclass <= 1) {
    cli_abort("Binomial regression: one class has 1 or 0 observations; not supported")
  }
  if (length(ntab) != 2) {
    cli_abort("Binomial regression: more than two classes is not supported")
  }
  if (minclass < 8) {
    cli_warn("Binomial regression: one class has fewer than 8 observations; dangerous ground")
  }
  # TODO, enable prediction with class labels if factor is passed
  if (intr == 1L && flmin < 1) b0_first <- coef(glm(y ~ 1, family = binomial()))
  y <- 2 * (as.integer(y) - 1) - 1 # convert to -1 / 1

  is.sparse <- FALSE
  if (inherits(x, "sparseMatrix")) {
    is.sparse <- TRUE
    x <- as_dgCMatrix(x)
  }
  if (standardize) {
    sx <- sqrt(Matrix::colSums(x^2))
    sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
    xs <- 1 / sx
    x <- x %*% Matrix::Diagonal(x = xs)
  }
  if (is.sparse) {
    xidx <- as.integer(x@i + 1)
    xcptr <- as.integer(x@p + 1)
    xval <- as.double(x@x)
    nnz <- as.integer(utils::tail(x@p, 1))
  }

  gamma <- 0.25 * calc_gamma_Rcpp(x, ix, iy, bn)

  if (!is.sparse) {
    fit <- dotCall64::.C64(
      "log_sparse_four",
      SIGNATURE = c(
        "integer", "integer", "integer", "integer", "double",
        "integer", "integer", "double", "double", "double",
        "double", "integer", "integer", "integer", "double",
        "double", "double", "integer", "integer", "integer",
        "double", "double", "integer", "integer", "double",
        "integer", "integer", "double", "double", "double"
      ),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma,
      nobs = nobs, nvars = nvars, x = as.double(x), y = as.double(y), pf = pf,
      pfl1 = pfl1,
      # Read / write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin, ulam = ulam,
      eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam), activeGroup = integer_dc(pmax),
      nbeta = integer_dc(nlam), alam = numeric_dc(nlam), npass = integer_dc(1),
      jerr = integer_dc(1),
      # read only
      alsparse = asparse, lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 11), rep("rw", 8), rep("w", 8), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl"
    )
  } else {
    fit <- dotCall64::.C64(
      "log_spmat_four",
      SIGNATURE = c(
        "integer", "integer", "integer", "integer", "double",
        "integer", "integer", "double", "integer", "integer",
        "integer", "double", "double", "double", "integer",
        "integer", "integer", "double", "double", "double",
        "integer", "integer", "integer", "double", "double",
        "integer", "integer", "double", "integer", "integer",
        "double", "double", "double"
      ),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma,
      nobs = nobs, nvars = nvars, x = as.double(xval), xidx = xidx,
      xcptr = xcptr, nnz = nnz, y = as.double(y), pf = pf, pfl1 = pfl1,
      # Read / write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin,
      ulam = ulam, eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam), activeGroup = integer_dc(pmax),
      nbeta = integer_dc(nlam), alam = numeric_dc(nlam),
      npass = integer_dc(1), jerr = integer_dc(1),
      # Read only
      alsparse = as.double(asparse), lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 14), rep("rw", 8), rep("w", 8), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl"
    )
  }
  # output
  outlist <- getoutput(x, group, fit, maxit, pmax, nvars, vnames, eps)
  if (standardize) outlist$beta <- outlist$beta * xs

  outlist$b0 <- matrix(outlist$b0, nrow = 1)
  if (intr == 1L && flmin < 1) outlist$b0[1] <- b0_first
  outlist <- c(
    outlist,
    list(
      npasses = fit$npass,
      jerr = fit$jerr,
      group = group,
      classnames = lev
    )
  )
  class(outlist) <- c("logitspgl")
  outlist
}


#' Regularization paths for sparse group-lasso models
#'
#' @description
#' Fits regularization paths for sparse group-lasso penalized learning problems at a
#' sequence of regularization parameters `lambda`.
#' Note that the objective function for least squares is
#' \deqn{RSS/(2n) + \lambda penalty}
#' Users can also tweak the penalty by choosing a different penalty factor.
#'
#'
#' @param x Double. A matrix of predictors, of dimension
#'   \eqn{n \times p}{n * p}; each row
#'   is a vector of measurements and each column is a feature. Objects of class
#'   [`Matrix::sparseMatrix`] are supported.
#' @param y Double/Integer/Factor. The response variable.
#'   Quantitative for `family="gaussian"` and for other exponential families.
#'   If `family="binomial"` should be either a factor with two levels or
#'   a vector of integers taking 2 unique values. For a factor, the last level
#'   in alphabetical order is the target class.
#' @param group Integer. A vector of consecutive integers describing the
#'   grouping of the coefficients (see example below).
#' @param family Character or function. Specifies the generalized linear model
#'   to use. Valid options are:
#'   * `"gaussian"` - least squares loss (regression, the default),
#'   * `"binomial"` - logistic loss (classification)
#'
#'   For any other type, a valid [stats::family()] object may be passed. Note
#'   that these will generally be much slower to estimate than the built-in
#'   options passed as strings. So for example, `family = "gaussian"` and
#'   `family = gaussian()` will produce the same results, but the first
#'   will be much faster.
#' @param nlambda The number of \code{lambda} values - default is 100.
#' @param lambda.factor A multiplicative factor for the minimal lambda in the
#'   `lambda` sequence, where `min(lambda) = lambda.factor * max(lambda)`.
#'   `max(lambda)` is the smallest value of `lambda` for which all coefficients
#'   are zero. The default depends on the relationship between \eqn{n}
#'   (the number of rows in the matrix of predictors) and \eqn{p}
#'   (the number of predictors). If \eqn{n \geq p}, the
#'   default is `0.0001`.  If \eqn{n < p}, the default is `0.01`.
#'   A very small value of `lambda.factor` will lead to a
#'   saturated fit. This argument has no effect if there is user-defined
#'   `lambda` sequence.
#' @param lambda A user supplied `lambda` sequence. The default, `NULL`
#'   results in an automatic computation based on `nlambda`, the smallest value
#'   of `lambda` that would give the null model (all coefficient estimates equal
#'   to zero), and `lambda.factor`. Supplying a value of `lambda` overrides
#'   this behaviour. It is likely better to supply a
#'   decreasing sequence of `lambda` values than a single (small) value. If
#'   supplied, the user-defined `lambda` sequence is automatically sorted in
#'   decreasing order.
#' @param pf_group Penalty factor on the groups, a vector of the same
#'   length as the total number of groups. Separate penalty weights can be applied
#'   to each group of \eqn{\beta}{beta's}s to allow differential shrinkage.
#'   Can be 0 for some
#'   groups, which implies no shrinkage, and results in that group always being
#'   included in the model (depending on `pf_sparse`). Default value for each
#'   entry is the square-root of the corresponding size of each group.
#'   Because this default is typical, these penalties are not rescaled.
#' @param pf_sparse Penalty factor on l1-norm, a vector the same length as the
#'   total number of columns in `x`. Each value corresponds to one predictor
#'   Can be 0 for some predictors, which
#'   implies that predictor will be receive only the group penalty.
#'   Note that these are internally rescaled so that the sum is the same as
#'   the number of predictors.
#' @param dfmax Limit the maximum number of groups in the model. Default is
#'   no limit.
#' @param pmax Limit the maximum number of groups ever to be nonzero. For
#'   example once a group enters the model, no matter how many times it exits or
#'   re-enters model through the path, it will be counted only once.
#' @param eps Convergence termination tolerance. Defaults value is `1e-8`.
#' @param maxit Maximum number of outer-loop iterations allowed at fixed lambda
#'   value. Default is `3e8`. If models do not converge, consider increasing
#'   `maxit`.
#' @param intercept Whether to include intercept in the model. Default is TRUE.
#' @param asparse The relative weight to put on the \eqn{\ell_1}-norm in
#'   sparse group lasso. Default is `0.05` (resulting in `0.95` on the
#'   \eqn{\ell_2}-norm).
#' @param standardize Logical flag for variable standardization (scaling) prior
#'   to fitting the model. Default is TRUE.
#' @param lower_bnd Lower bound for coefficient values, a vector in length of 1
#'   or of length the number of groups. Must be non-positive numbers only.
#'   Default value for each entry is `-Inf`.
#' @param upper_bnd Upper for coefficient values, a vector in length of 1
#'   or of length the number of groups. Must be non-negative numbers only.
#'   Default value for each entry is `Inf`.
#' @param weights Double vector. Optional observation weights. These can
#'   only be used with a [stats::family()] object.
#' @param offset Double vector. Optional offset (constant predictor without a
#'   corresponding coefficient). These can only be used with a
#'   [stats::family()] object.
#' @param warm List created with [make_irls_warmup()]. These can only be used
#'   with a [stats::family()] object, and is not typically necessary even then.
#' @param trace_it Scalar integer. Larger values print more output during
#'   the irls loop. Typical values are `0` (no printing), `1` (some printing
#'   and a progress bar), and `2` (more detailed printing).
#'   These can only be used with a [stats::family()] object.
#'
#' @return An object with S3 class `"sparsegl"`. Among the list components:
#' * `call` The call that produced this object.
#' * `b0` Intercept sequence of length `length(lambda)`.
#' * `beta` A `p` x `length(lambda)` sparse matrix of coefficients.
#' * `df` The number of features with nonzero coefficients for each value of
#'     `lambda`.
#' * `dim` Dimension of coefficient matrix.
#' * `lambda` The actual sequence of `lambda` values used.
#' * `npasses` Total number of iterations summed over all `lambda` values.
#' * `jerr` Error flag, for warnings and errors, 0 if no error.
#' * `group` A vector of consecutive integers describing the grouping of the
#'     coefficients.
#' * `nobs` The number of observations used to estimate the model.
#'
#' If `sparsegl()` was called with a [stats::family()] method, this may also
#' contain information about the deviance and the family used in fitting.
#'
#'
#' @seealso [cv.sparsegl()] and the [`plot()`][plot.sparsegl()],
#'   [`predict()`][predict.sparsegl()], and [`coef()`][coef.sparsegl()]
#'   methods for `"sparsegl"` objects.
#'
#' @references Liang, X., Cohen, A., Sólon Heinsfeld, A., Pestilli, F., and
#'   McDonald, D.J. 2024.
#'   \emph{sparsegl: An `R` Package for Estimating Sparse Group Lasso.}
#'   Journal of Statistical Software, Vol. 110(6): 1–23.
#'   \doi{10.18637/jss.v110.i06}.
#'
#' @export
#'
#'
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit <- sparsegl(X, y, group = groups)
#'
#' yp <- rpois(n, abs(X %*% beta_star))
#' fit_pois <- sparsegl(X, yp, group = groups, family = poisson())
sparsegl_modified <- function(
    x, y, group = NULL, family = c("gaussian", "binomial"),
    nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
    lambda = NULL, pf_group = sqrt(bs), pf_sparse = rep(1, nvars),
    intercept = TRUE, asparse = 0.05, standardize = TRUE,
    lower_bnd = -Inf, upper_bnd = Inf,
    weights = NULL, offset = NULL, warm = NULL,
    trace_it = 0,
    dfmax = as.integer(max(group)) + 1L,
    pmax = min(dfmax * 1.2, as.integer(max(group))),
    eps = 1e-08, maxit = 3e+06) {
  this.call <- match.call()
  if (!is.matrix(x) && !inherits(x, "sparseMatrix")) {
    cli_abort("`x` must be a matrix.")
  }

  if (any(is.na(x))) cli_abort("Missing values in `x` are not supported.")

  y <- drop(y)
  if (!is.null(dim(y))) cli_abort("`y` must be a vector or 1-column matrix.")
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)

  if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")

  if (length(y) != nobs) {
    cli_abort("`x` has {nobs} rows while `y` has {length(y)}.")
  }

  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else {
    if (length(group) != nvars) {
      cli_abort(c(
        "The length of `group` is {length(group)}.",
        "It must match the number of columns in `x`: {nvars}"
      ))
    }
  }

  bn <- as.integer(max(group)) # number of groups
  bs <- as.integer(as.numeric(table(group))) # number of elements in each group

  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn))) {
    cli_abort("Groups must be consecutively numbered 1, 2, 3, ...")
  }

  if (asparse > 1) {
    cli_abort(c(
      "`asparse` must be less than or equal to 1.",
      i = "You may want {.fn glmnet::glmnet} instead."
    ))
  }

  if (asparse < 0) {
    asparse <- 0
    cli_warn("`asparse` must be in {.val [0, 1]}, running ordinary group lasso.")
  }
  if (any(pf_sparse < 0)) cli::cli_abort("`pf_sparse` must be non-negative.")
  if (any(is.infinite(pf_sparse))) {
    cli_abort(
      "`pf_sparse` may not be infinite. Simply remove the column from `x`."
    )
  }
  if (any(pf_group < 0)) cli_abort("`pf_group` must be non-negative.")
  if (any(is.infinite(pf_group))) {
    cli_abort(c(
      "`pf_group` must be finite.",
      i = "Simply remove the group from `x`."
    ))
  }
  if (all(pf_sparse == 0)) {
    if (asparse > 0) {
      cli_abort(
        "`pf_sparse` is identically 0 but `asparse` suggests some L1 penalty is desired."
      )
    } else {
      cli_warn("`pf_sparse` was set to 1 because `asparse` = {.val {0}}.")
      pf_sparse <- rep(1, nvars)
    }
  }

  ## Note: should add checks to see if any columns are completely unpenalized
  ## This is not currently expected.

  iy <- cumsum(bs) # last column of x in each group
  ix <- c(0, iy[-bn]) + 1 # first column of x in each group
  ix <- as.integer(ix)
  iy <- as.integer(iy)
  group <- as.integer(group)

  # parameter setup
  if (length(pf_group) != bn) {
    cli_abort(
      "The length of `pf_group` must be the same as the number of groups: {.val {bn}}."
    )
  }
  if (length(pf_sparse) != nvars) {
    cli_abort(
      "The length of `pf_sparse` must be equal to the number of predictors: {.val {nvars}}."
    )
  }

  pf_sparse <- pf_sparse / sum(pf_sparse) * nvars
  maxit <- as.integer(maxit)
  pf_group <- as.double(pf_group)
  pf_sparse <- as.double(pf_sparse)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)

  # lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) {
      cli::cli_abort("`lambda.factor` must be less than {.val {1}}.")
    }
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    # flmin = 1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) cli_abort("`lambda` must be non-negative.")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  intr <- as.integer(intercept)

  ### check on upper/lower bounds
  if (any(lower_bnd > 0)) cli_abort("`lower_bnd` must be non-positive.")
  if (any(upper_bnd < 0)) cli_abort("`upper_bnd` must be non-negative.")
  lower_bnd[lower_bnd == -Inf] <- -9.9e30
  upper_bnd[upper_bnd == Inf] <- 9.9e30
  if (length(lower_bnd) < bn) {
    if (length(lower_bnd) == 1) {
      lower_bnd <- rep(lower_bnd, bn)
    } else {
      cli_abort("`lower_bnd` must be length {.val {1}} or length {.val {bn}}.")
    }
  } else {
    lower_bnd <- lower_bnd[seq_len(bn)]
  }
  if (length(upper_bnd) < bn) {
    if (length(upper_bnd) == 1) {
      upper_bnd <- rep(upper_bnd, bn)
    } else {
      cli_abort("`upper_bnd` must be length {.val {1}} or length {.val {bn}}.")
    }
  } else {
    upper_bnd <- upper_bnd[seq_len(bn)]
  }
  storage.mode(upper_bnd) <- "double"
  storage.mode(lower_bnd) <- "double"

  # call R sub-function
  fam <- sparsegl:::validate_family(family)
  if (fam$check == "char") {
    family <- match.arg(family)
    if (!is.null(weights)) {
      cli_warn(c(
        "Currently, `weights` are only supported when `family` has class {.cls family}.",
        i = "Estimating unweighted sparse group lasso. See {.fn sparsegl::sparsegl}."
      ))
    }
    if (!is.null(offset)) {
      cli_warn(c(
        "Currently, `offset` is only supported when `family` has class {.cls family}.",
        i = "Estimating sparse group lasso without any offset. See {.fn sparsegl::sparsegl}."
      ))
    }
    fit <- switch(family,
      gaussian = sgl_ls_modified(
        bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
        dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
        as.double(asparse), standardize, lower_bnd, upper_bnd
      ),
      binomial = sgl_logit_modified(
        bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
        dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
        as.double(asparse), standardize, lower_bnd, upper_bnd
      )
    )
  }
  if (fam$check == "fam") {
    fit <- sgl_irwls(
      bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
      dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
      as.double(asparse), standardize, lower_bnd, upper_bnd, weights,
      offset, fam$family, trace_it, warm
    )
  }

  # output
  if (is.null(lambda)) fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  fit$asparse <- asparse
  fit$nobs <- nobs
  fit$pf_group <- pf_group
  fit$pf_sparse <- pf_sparse
  class(fit) <- c(class(fit), "sparsegl")
  fit
}