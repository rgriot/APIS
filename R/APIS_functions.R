#' APIS function that assigns with observed data
#'
#' This function performs the APIS procedure
#' @param off.genotype Offspring genotypes | Matrix (n*p)
#' where n = number of individuals
#' p = number of markers
#' rownames(offspring) = labels of offspring
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param sire.genotype Sire genotypes | Matrix (n*p)
#' where n = number of individuals
#' p = number of markers
#' rownames(sire) = labels of sires
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param dam.genotype Dam genotypes | Matrix (n*p)
#' where n = number of individuals
#' p = number of markers
#' rownames(dam) = labels of dams
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param error (default: NULL) The assignment error rate accepted by the user
#' @param exclusion.threshold (default: ncol(off.genotype)) Threshold for exclusion (number of mismatches allowed)
#' @param preselect.Parent (default: FALSE) Preselection of parents. Can be FALSE, an integer or a vector of two integers (number of sires, numbers of dams)
#' @keywords assignment APIS
#' @return pedigree
#' @return a log file
#' @examples data("genotype_APIS")
#'
#' result <- APIS(off.genotype = off_full,
#'                sire.genotype = sire_full,
#'                dam.genotype = dam_full,
#'                error = 0.05)
#' @export
APIS <- function(off.genotype, sire.genotype, dam.genotype, error = NULL, exclusion.threshold = ncol(off.genotype), preselect.Parent = F) {

  # Check inputs

  #	Check if all genotypes matrices have the same number of markers
  if (ncol(off.genotype) == ncol(sire.genotype) & ncol(off.genotype) == ncol(dam.genotype) & ncol(sire.genotype) == ncol(dam.genotype)) {
    cat("genotype matrices : OK")
    cat('\n')
  } else {
    stop("Your genotype matrices do not have the same number of markers")
  }

  #	Check if the number of mismatches allowed is lower than the number of markers and positive
  if ((0 <= exclusion.threshold) && (exclusion.threshold <= ncol(off.genotype))) {
    cat("exclusion threshold : OK")
    cat('\n')
  } else {
    stop("The exclusion threshold is greater than the number of markers")
  }

  # Check if the user-defined assignment error rate limit is a percentage
  if ((0 <= error) && (error <= 100)) {
    cat("assignment error rate : OK")
    cat('\n')
  } else {
    stop("The assignment error rate limit is NEGATIVE")
  }


  # Calculate the theoretical assignment power
  P <- assignmentPower(sire = sire.genotype, dam = dam.genotype)
  P2 <- substr(as.character(100 * P), 1, 6)
  cat("The assignment power of your marker set is ", P2, "%", sep = "")
  cat('\n')

  if (P >= 0.99) {
    cat("Theoretical assignment power : OK")
    cat('\n')
  } else {
    message("WARNING! Your marker set is not enough powerfull!")
  }

  # Assign with APIS
  assignResult 	<- assignmentFortran(offspring = off.genotype, sire = sire.genotype, dam = dam.genotype, thresh = exclusion.threshold)
  apisResult 		<- setThreshold(ped.log = assignResult$log.mendel, ped.exclu = assignResult$exclu, nb.mrk = assignResult$nb.mrk, error = error)

  pedigree 	<- apisResult$pedigree
  log 		  <- apisResult$log

  # Give recommendations according to the data set and the results
  cat('--------------------------------------', sep = '\n')
  cat('             APIS SUMMARY', sep = '\n')
  cat('--------------------------------------', sep = '\n')
  cat('Theoretical assignment power of the marker set : ', P2, "%", sep = "")
  cat('\n')
  cat('Assignment error rate accepted : ', error)
  cat('\n')
  assignmentRate <- length(pedigree$sire[which(pedigree$sire != '0')]) / nrow(pedigree)
  AR <- substr(as.character(100 * assignmentRate), 1, 6)
  cat('Assignment rate : ', AR, '%', sep = "")
  cat('\n')

  # Return outputs
  output <- list(pedigree = pedigree, log = log, error = error)
}

# ----------------------------------------------------------------------------------------------------------------
#' Assignment function to obtain the average Mendelian transmission probabilities
#'
#' This function calculates the average Mendelian transmission probabilities
#' @param offspring Offspring genotypes | Matrix (n*p) where n = number of individuals, p = number of markers
#' rownames(offspring) = labels of offspring
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param sire Sire genotypes | Matrix (n*p) where n = number of individuals, p = number of markers
#' rownames(sire) = labels of sires
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param dam Dam genotypes | Matrix (n*p) where n = number of individuals, p = number of markers
#' rownames(dam) = labels of dams
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param thresh (default: ncol(offspring) Threshold for exclusion (number of mismatches allowed)
#' @param preselect.Parent (default: FALSE) Preselection of parents. Can be FALSE, an integer or a vector of two integers (number of sires, numbers of dams)
#' @keywords assignment
#' @return intermidiate pedigree
#' @return log file for Mendelian transmission probabilities
#' @return log file for exclusion
#' @export
assignment <- function(offspring, sire, dam, thresh = ncol(offspring), preselect.Parent = F) {
  # DESCRIPTION
  # Function to calculate average Mendelian transmission probabilities

  # Stop if different number of markers are provided
  if (ncol(offspring)!=ncol(sire)&ncol(offspring)!=ncol(dam))
    stop('Genotypes must have the same number of markers')

  # Create the results matrix
  ped.log <- as.data.frame(matrix(NA,nrow = nrow(offspring), ncol = 16))
  colnames(ped.log) <- c('off', 'mrk_genotype','sire1', 'dam1', 'miss1', 'mendel1',
                         'sire2', 'dam2', 'miss2', 'mendel2', 'delta_Pmendel12',
                         'sire3', 'dam3', 'miss3', 'mendel3', 'delta_Pmendel23')
  ped.log[,1] <- rownames(offspring)
  ped.log[,2] <- as.numeric(apply(offspring,1, function(X) {length(X[X!='NA/NA'])}))

  ped.exclu <- as.data.frame(matrix(NA, nrow = nrow(offspring), ncol = 11))
  colnames(ped.exclu) <- c('off', 'mrk_genotype','sire1', 'dam1', 'miss1',
                           'sire2', 'dam2', 'miss2',
                           'sire3', 'dam3', 'miss3')

  # Parameters
  e <- 0.01 # Genotyping error of 1%

  # Estimate allele frequencies
  cat('Estimation of allele frequencies')
  cat('\n')
  Freq <- allFreq(offspring)

  x.col <- ncol(offspring)

  # Exclusion tables
  homo_exclu <- matrix(c(0,0,1,0,
                         0,0,1,0,
                         1,1,2,1,
                         0,0,1,0), nrow = 4, ncol = 4)
  rownames(homo_exclu) <- c('AA', 'AC', 'CC', 'miss')
  colnames(homo_exclu) <- c('AA', 'AC', 'CC', 'miss')

  hetero_exclu <- matrix(c(1,0,0,1,0,1,0,
                           0,0,0,0,0,1,0,
                           0,0,1,0,1,1,0,
                           1,0,0,1,0,1,0,
                           0,0,1,0,1,1,0,
                           1,1,1,1,1,2,1,
                           0,0,0,0,0,1,0), nrow = 7, ncol = 7)
  rownames(hetero_exclu) <- c('AA', 'AB', 'BB', 'AC', 'BC', 'CC', 'miss')
  colnames(hetero_exclu) <- c('AA', 'AB', 'BB', 'AC', 'BC', 'CC', 'miss')

  # Variables for parent selection
  parent.genotype <- rbind(sire, dam)

  parent.sex <- c(rep("M", times = nrow(sire)), rep("F", times = nrow(dam)))

  # Assignment Process
  cat('Assignment')
  cat('\n')

  # Set up the cluster for parallel iteration
  cl <- makeCluster(parallel::detectCores()-1)
  registerDoSNOW(cl)

  iterations <- nrow(offspring)
  pb <- txtProgressBar(min = 0, max = iterations, char = "><(((°> ", style = 3)

  progress <- function(n) {setTxtProgressBar(pb, n)}

  opts <- list(progress = progress)

  # Start
  A <- foreach(off=1:iterations, .multicombine = T,
               .packages = c('foreach', 'doParallel', 'doSNOW'), .options.snow = opts) %dopar% { # For each offspring

                 tmp <- offspring[off,, drop = F]


                 if (preselect.Parent == F) {
                   potential.sire <- rownames(sire)
                   potential.dam <- rownames(dam)
                 } else {
                   potential.parents <- selectParents(tmp, parent.genotype = parent.genotype,
                                                      parent.sex = parent.sex, n.Parent = preselect.Parent)

                   potential.sire <- potential.parents$sire_toKeep
                   potential.dam <- potential.parents$dam_toKeep
                 }


                 # Create temporary results
                 res <- matrix(NA, nrow = (length(potential.sire)*length(potential.dam)), ncol = 4)
                 colnames(res) <- c('sire', 'dam', 'score_exclu', 'P_mendel')
                 res[,1] <- rep(potential.sire, each = length(potential.dam))
                 res[,2] <- rep(potential.dam, times = length(potential.sire))

                 # Lielihood tables
                 table_like <- vector('list', x.col)
                 off_geno <- strsplit(tmp, split = "/", fixed = T)
                 for (i in 1:length(tmp)) {
                   if (off_geno[[i]][1]==off_geno[[i]][2]) {
                     fa <- Freq[i,which(colnames(Freq)==paste0('Freq_',off_geno[[i]][1]))]
                     fa <- ifelse(test = (fa==0|is.na(fa)), yes = Freq[i,which(colnames(Freq)==paste0('Freq_NA'))], no = fa)
                     table_like[[i]] <- matrix(c(1,0.5,e,fa,0.5,0.25,e,0.5*fa,e,e,e,e,fa,0.5*fa,e,fa*fa), nrow = 4, ncol = 4)
                   } else {
                     fa <- Freq[i,which(colnames(Freq)==paste0('Freq_',off_geno[[i]][1]))]
                     fa <- ifelse(test = (fa==0|is.na(fa)), yes = Freq[i,which(colnames(Freq)==paste0('Freq_NA'))], no = fa)
                     fb <- Freq[i,which(colnames(Freq)==paste0('Freq_',off_geno[[i]][2]))]
                     fb <- ifelse(test = (fb==0|is.na(fb)), yes = Freq[i,which(colnames(Freq)==paste0('Freq_NA'))], no = fb)

                     table_like[[i]] <- matrix(c(e,0.5,1,e,0.5,e,fb,
                                                 0.5,0.5,0.5,0.25,0.5,e,0.5*(fa+fb),
                                                 1,0.5,e,0.5,e,e,fa,
                                                 e,0.25,0.5,e,0.25,e,0.5*fb,
                                                 0.5,0.25,e,0.25,e,e,0.5*fb,
                                                 e,e,e,e,e,e,e,
                                                 fb,0.5*(fa+fb),fa,0.5*fb,0.5*fa,e,2*fa*fb), nrow = 7, ncol = 7)
                   }
                 }

                 t <- foreach(n = 1:nrow(res), .combine = rbind,
                              .multicombine = T, .packages = c('foreach', 'doParallel', 'doSNOW')) %dopar% { # For each parents pair

                                #p = sire // m = dam
                                p <- sire[which(rownames(sire)==res[n,1]),, drop = F]
                                m <- dam[which(rownames(dam)==res[n,2]),, drop = F]

                                # Keep likelihood and missmatches
                                sc_exclu <- vector(mode = 'numeric', length = x.col)
                                sc_like <- vector(mode = 'numeric', length = x.col)

                                off_split <- strsplit(tmp, split = "/", fixed = T)
                                p_split <- strsplit(p, split = "/", fixed = T)
                                m_split <- strsplit(m, split = "/", fixed = T)

                                mrk <- 1

                                while (mrk <= x.col & sum(sc_exclu) <= thresh) { # For each marker
                                  off_mrk <- off_split[[mrk]]
                                  p_mrk <- p_split[[mrk]]
                                  m_mrk <- m_split[[mrk]]

                                  if (off_mrk[1]=='NA'&off_mrk[2]=='NA') { # if offspring genotype is missing
                                    sc_exclu[mrk] <- 0
                                    sc_like[mrk] <- 1
                                  } else if (off_mrk[1]==off_mrk[2]&off_mrk[1]!='NA') { #If offspring is HOMOZYGOUS
                                    # SIRE ID
                                    if (p_mrk[1]==p_mrk[2]) { # If the sire is HOMOZYGOUS
                                      if (p_mrk[1]=='NA') {
                                        p_id <- 4 # sire NA/NA
                                      }else if (p_mrk[1]==off_mrk[1]) {
                                        p_id <- 1 # sire A/A
                                      } else {
                                        p_id <- 3 # sire "C/C"
                                      }
                                    } else { # If the sire is HETEROZYGOUS
                                      if ((p_mrk[1]!=off_mrk[1]&p_mrk[2]!=off_mrk[1])&(p_mrk[1]!=off_mrk[2]&p_mrk[2]!=off_mrk[2])) {
                                        p_id <- 3 # sire C/C
                                      } else {
                                        p_id <- 2 # sire A/C
                                      }
                                    }
                                    #DAM ID
                                    if (m_mrk[1]==m_mrk[2]) { #If the dam is HOMOZYGOUS
                                      if (m_mrk[1]=='NA') {
                                        m_id <- 4 # dam NA/NA
                                      }else if (m_mrk[1]==off_mrk[1]) {
                                        m_id <- 1 # dam A/A
                                      } else {
                                        m_id <- 3 # dam C/C
                                      }
                                    } else { # If the dam is HETEROZYGOUS
                                      if ((m_mrk[1]!=off_mrk[1]&m_mrk[2]!=off_mrk[1])&(m_mrk[1]!=off_mrk[2]&m_mrk[2]!=off_mrk[2])) {
                                        m_id <- 3 # dam C/C
                                      } else {
                                        m_id <- 2 # dam A/C
                                      }
                                    }

                                    sc_exclu[mrk] <- homo_exclu[m_id,p_id]
                                    sc_like[mrk] <- table_like[[mrk]][m_id, p_id]

                                  } else { #If the offspring is HETEROZYGOUS
                                    # SIRE ID
                                    if (p_mrk[1]==p_mrk[2]) { # If the sire is HOMOZYGOUS
                                      if (p_mrk[1]=='NA') {
                                        p_id <- 7 # sire NA/NA
                                      }else if (p_mrk[1]==off_mrk[1]) {
                                        p_id <- 1 # sire A/A
                                      } else if (p_mrk[1]==off_mrk[2]) {
                                        p_id <- 3 # sire B/B
                                      } else {
                                        p_id <- 6 # sire C/C
                                      }
                                    } else { # If the sire is HETEROZYGOUS
                                      if ((p_mrk[1]==off_mrk[1]|p_mrk[2]==off_mrk[1])&(p_mrk[1]==off_mrk[2]|p_mrk[2]==off_mrk[2])) {
                                        p_id <- 2 # sire A/B
                                      } else if ((p_mrk[1]==off_mrk[1]|p_mrk[2]==off_mrk[1])&(p_mrk[1]!=off_mrk[2]|p_mrk[2]!=off_mrk[2])){
                                        p_id <- 4 # sire A/C
                                      } else if ((p_mrk[1]!=off_mrk[1]|p_mrk[2]!=off_mrk[1])&(p_mrk[1]==off_mrk[2]|p_mrk[2]==off_mrk[2])) {
                                        p_id <- 5 # sire B/C
                                      } else {
                                        p_id <- 6 # sire C/C
                                      }
                                    }

                                    # DAM ID
                                    if (m_mrk[1]==m_mrk[2]) { # If the dam is HOMOZYGOUS
                                      if (m_mrk[1]=='NA') {
                                        m_id <- 7 # dam NA/NA
                                      }else if (m_mrk[1]==off_mrk[1]) {
                                        m_id <- 1 # dam A/A
                                      } else if (m_mrk[1]==off_mrk[2]) {
                                        m_id <- 3 # dam B/B
                                      } else {
                                        m_id <- 6 # dam C/C
                                      }
                                    } else { # If the dam is HETEROZYGOUS
                                      if ((m_mrk[1]==off_mrk[1]|m_mrk[2]==off_mrk[1])&(m_mrk[1]==off_mrk[2]|m_mrk[2]==off_mrk[2])) {
                                        m_id <- 2 # dam A/B
                                      } else if ((m_mrk[1]==off_mrk[1]|m_mrk[2]==off_mrk[1])&(m_mrk[1]!=off_mrk[2]|m_mrk[2]!=off_mrk[2])){
                                        m_id <- 4 # dam A/C
                                      } else if ((m_mrk[1]!=off_mrk[1]|m_mrk[2]!=off_mrk[1])&(m_mrk[1]==off_mrk[2]|m_mrk[2]==off_mrk[2])) {
                                        m_id <- 5 # dam B/C
                                      } else {
                                        m_id <- 6 # dam C/C
                                      }
                                    }

                                    # Get the score for tested marker
                                    sc_exclu[mrk] <- hetero_exclu[m_id,p_id]
                                    sc_like[mrk] <- table_like[[mrk]][m_id, p_id]
                                  }

                                  mrk <- mrk + 1

                                }

                                # Create the parallel loop output
                                r <- c(NA,NA)
                                r[1] <- sum(sc_exclu) # Number of missmatch
                                r[2] <- exp(sum(log(sc_like))/ped.log[off,2]) # Mendelian transimission probability

                                return(r)
                              }

                 # Working on the results
                 res[,3:4] <- t
                 res <- as.data.frame(res)
                 res$sire <- as.character(res$sire)
                 res$dam <- as.character(res$dam)
                 res$score_exclu <- as.numeric(as.character(res$score_exclu))
                 res$P_mendel <- as.numeric(as.character(res$P_mendel))

                 # Order by Mendelian transmission probabilities
                 res2 <- res[order(res[,4], res[,3], decreasing = T),]

                 delta_P12 <- res2[1,4]-res2[2,4]
                 delta_P23 <- res2[2,4]-res2[3,4]

                 p_fin <- res2[1,1]
                 m_fin <- res2[1,2]

                 ped.out <- c(rownames(tmp), NA, NA)
                 if (is.na(p_fin)&is.na(m_fin)) {
                   ped.out[2:3] <- c(NA, NA)
                 } else {
                   ped.out[2:3] <- c(p_fin,m_fin)
                 }

                 log.out <- unlist(c(ped.out[1], ped.log[off,2], ped.out[2:3], res2[1, 3:4],
                                     res2[2,1:4], delta_P12,
                                     res2[3,1:4], delta_P23))

                 # Order by mismatches
                 res2 <- res[order(res[,3], -res[,4], decreasing = F),]


                 exclu.out <- unlist(c(log.out[1],ped.log[off,2], res2[1,1:3], res2[2,1:3], res2[3,1:3]))

                 a <- list(ped.out, log.out, exclu.out)
               }

  # Stop the parallel cluster
  close(pb)
  stopCluster(cl)

  ped <- as.data.frame(t(as.data.frame(lapply(A, function(X) {t <- X[[1]]}))))
  colnames(ped) <- c('off', 'sire', 'dam')
  rownames(ped) <- c(1:nrow(ped))

  # Working on the data
  ped$off <- as.character(ped$off)
  ped$sire <- as.character(ped$sire)
  ped$dam <- as.character(ped$dam)

  # Create a log
  ped.log <- as.data.frame(t(as.data.frame(lapply(A, function(X) {t <- X[[2]]}))))
  colnames(ped.log) <- c('off', 'mrk_genotype','sire1', 'dam1', 'miss1', 'mendel1',
                         'sire2', 'dam2', 'miss2', 'mendel2', 'delta_Pmendel12',
                         'sire3', 'dam3', 'miss3', 'mendel3', 'delta_Pmendel23')
  rownames(ped.log) <- c(1:nrow(ped.log))

  ped.log[,] <- sapply(ped.log[,c(1:ncol(ped.log))], as.character)
  ped.log[,c(2, 5:6, 9:11, 14:16)] <- sapply(ped.log[,c(2, 5:6, 9:11, 14:16)], as.numeric)

  # Create a data frame from results by exclusion
  ped.exclu <- as.data.frame(t(as.data.frame(lapply(A, function(X) {t <- X[[3]]}))))
  colnames(ped.exclu) <- c('off', 'mrk_genotype','sire1', 'dam1', 'miss1',
                           'sire2', 'dam2', 'miss2',
                           'sire3', 'dam3', 'miss3')
  rownames(ped.exclu) <- c(1:nrow(ped.exclu))
  ped.exclu[,c(2,5,8,11)] <- sapply(ped.exclu[,c(2,5,8,11)], as.numeric)

  # Return the output
  return(list(pedigree = ped, log.mendel = ped.log, log.exclu = ped.exclu, nb.mrk = ncol(offspring)))
}

# ----------------------------------------------------------------------------------------------------------------
#' Assignment function to obtain the average Mendelian transmission probabilities using a Fortran library
#'
#' This function calculates the average Mendelian transmission probabilities
#' @param offspring Offspring genotypes | Matrix (n*p) where n = number of individuals, p = number of markers
#' rownames(offspring) = labels of offspring
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param sire Sire genotypes | Matrix (n*p) where n = number of individuals, p = number of markers
#' rownames(sire) = labels of sires
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param dam Dam genotypes | Matrix (n*p) where n = number of individuals, p = number of markers
#' rownames(dam) = labels of dams
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param thresh (default: ncol(offspring) Threshold for exclusion (number of mismatches allowed)
#' @param preselect.Parent (default: FALSE) Preselection of parents. Can be FALSE, an integer or a vector of two integers (number of sires, numbers of dams)
#' @keywords assignment
#' @return intermidiate pedigree
#' @return log file for Mendelian transmission probabilities
#' @return log file for exclusion
#' @export
assignmentFortran <- function(offspring, sire, dam, thresh = ncol(offspring)) {

  # DESCRIPTION
  # Function to calculate average Mendelian transmission probabilities

  # Stop if different number of markers are provided
  if (ncol(offspring)!=ncol(sire)&ncol(offspring)!=ncol(dam))
    stop('Genotypes must have the same number of markers')

  # Create the results matrix
  ped.log <- as.data.frame(matrix(NA,nrow = nrow(offspring), ncol = 16))
  colnames(ped.log) <- c('off', 'mrk_genotype','sire1', 'dam1', 'miss1', 'mendel1',
                         'sire2', 'dam2', 'miss2', 'mendel2', 'delta_Pmendel12',
                         'sire3', 'dam3', 'miss3', 'mendel3', 'delta_Pmendel23')
  ped.log[,1] <- rownames(offspring)
  ped.log[,2] <- as.numeric(apply(offspring,1, function(X) {length(X[X!='NA/NA'])}))

  ped.exclu <- as.data.frame(matrix(NA, nrow = nrow(offspring), ncol = 11))
  colnames(ped.exclu) <- c('off', 'mrk_genotype','sire1', 'dam1', 'miss1',
                           'sire2', 'dam2', 'miss2',
                           'sire3', 'dam3', 'miss3')

  # Parameters
  e <- 0.01 # Genotyping error of 1%

  # Estimate allele frequencies
  cat('Estimation of allele frequencies')
  cat('\n')
  Freq <- allFreq(offspring)

  x.col <- ncol(offspring)
  iterations <- nrow(offspring)

  cat('Recoding datasets')
  cat('\n')

  # Recode for Fortran subroutine
  recodeFortran <- function(mrk, list.mrk) {
    tmp <- unlist(strsplit(mrk, '/'))
    all1 <- list.mrk[which(list.mrk[,1] == tmp[1]),2]
    all2 <- list.mrk[which(list.mrk[,1] == tmp[2]),2]

    return(c(all1, all2))
  }

  recodeFreq <- function(name.col, list.mrk) {
    tmp <- unlist(strsplit(name.col, '_'))[2]
    return(list.mrk[which(list.mrk[,1] == tmp), 2])
  }

  variant <- unique(unlist(strsplit(as.vector(offspring), '/')))
  variant <- variant[-which(variant == "NA")]

  variant.corres <- data.frame(variant = as.character(variant),
                               recode = c(1:length(variant)))
  variant.corres$variant <- as.character(variant.corres$variant)
  variant.corres <- rbind(variant.corres, c(as.character("NA"), 0))

  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)


  recode.off <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
    tmp <- as.numeric(as.vector(sapply(offspring[i,, drop = F], recodeFortran, list.mrk = variant.corres)))
  }

  recode.sire <- foreach(i = 1:nrow(sire), .combine = rbind) %dopar% {
    tmp <- as.numeric(as.vector(sapply(sire[i,, drop = F], recodeFortran, list.mrk = variant.corres)))
  }

  recode.dam <- foreach(i = 1:nrow(dam), .combine = rbind) %dopar% {
    tmp <- as.numeric(as.vector(sapply(dam[i,, drop = F], recodeFortran, list.mrk = variant.corres)))
  }

  stopCluster(cl)

  rownames(recode.off) <- rownames(offspring)
  rownames(recode.sire) <- rownames(sire)
  rownames(recode.dam) <- rownames(dam)

  Freq <- Freq[,-which(colnames(Freq) == "Freq_NA")]
  Freq <- Freq[, c((floor(ncol(Freq)/2)+2):ncol(Freq))]

  colnames(Freq) <- sapply(colnames(Freq), recodeFreq, list.mrk = variant.corres)
  Freq <- rbind(colnames(Freq), Freq)
  Freq <- apply(Freq, 2, as.numeric)
  Freq <- Freq[,order(Freq[1,])]
  Freq <- Freq[-1,]

  # Assignment Process
  cat('Assignment')
  cat('\n')

  pb <- txtProgressBar(min = 0, max = iterations, char = "><(((°> ", style = 3)

  for (off in 1:nrow(offspring)) { #Pour chaque descendant
    tmp <- recode.off[off,]

    #Resultats
    res <- matrix(NA, nrow = (nrow(sire)*nrow(dam)), ncol = 4)
    colnames(res) <- c('sire', 'dam', 'score_exclu', 'P_mendel')
    res[,1] <- rep(rownames(sire), each = nrow(dam))
    res[,2] <- rep(rownames(dam), times = nrow(sire))
    res <- as.data.frame(res)

    nMrk         = as.integer(x.col)
    nSires       = as.integer(nrow(recode.sire))
    nDams        = as.integer(nrow(recode.dam))
    nVariant     = as.integer(ncol(Freq))
    output_sires = vector(mode = 'integer', length = nSires*nDams)
    output_dams  = vector(mode = 'integer', length = nSires*nDams)
    output_score = vector(mode = 'numeric', length = nSires*nDams)
    output_miss  = vector(mode = 'numeric', length = nSires*nDams)

    outputFortran <- .Fortran("likelihoodcalculation", as.integer(tmp), as.integer(recode.sire),
                              as.integer(recode.dam), as.integer(nMrk), as.integer(nVariant), as.integer(nSires), as.integer(nDams),
                              as.double(Freq), output_sires, output_dams, as.double(output_score), as.integer(output_miss))


    res[,3] <- as.integer(outputFortran[[12]])
    res[,4] <- exp(as.numeric(as.numeric(outputFortran[[11]]))/ped.log[off,2])

    #Working on the results
    res$sire <- as.character(res$sire)
    res$dam <- as.character(res$dam)
    res$score_exclu <- as.numeric(as.character(res$score_exclu))
    res$P_mendel <- as.numeric(as.character(res$P_mendel))

    # Order by Mendelian transmission probabilities
    res2 <- res[order(res[,4], res[,3], decreasing = T),]

    delta_P12 <- res2[1,4]-res2[2,4]
    delta_P23 <- res2[2,4]-res2[3,4]

    p_fin <- res2[1,1]
    m_fin <- res2[1,2]

    ped.log[off,3:ncol(ped.log)] <- unlist(c(p_fin, m_fin, res2[1, 3:4], res2[2,1:4], delta_P12, res2[3,1:4], delta_P23))

    #Order by mismatches
    res2 <- res[order(res[,3], -res[,4], decreasing = F),]

    ped.exclu[off,] <- unlist(c(ped.log[off,c(1,2)], res2[1,1:3], res2[2,1:3], res2[3,1:3]))

    setTxtProgressBar(pb, off)
  }

  cat('\n')

  # Create the outputs
  ped <- data.frame(off = ped.log$off,
                    sire = ped.log$sire1,
                    dam = ped.log$dam1)
  colnames(ped) <- c('off', 'sire', 'dam')
  rownames(ped) <- c(1:nrow(ped))

  ped$off <- as.character(ped$off)
  ped$sire <- as.character(ped$sire)
  ped$dam <- as.character(ped$dam)

  # Create a log
  ped.log[,] <- sapply(ped.log[,c(1:ncol(ped.log))], as.character)
  ped.log[,c(2, 5:6, 9:11, 14:16)] <- sapply(ped.log[,c(2, 5:6, 9:11, 14:16)], as.numeric)

  # Create a data frame from results by exclusion
  ped.exclu[,c(2,5,8,11)] <- sapply(ped.exclu[,c(2,5,8,11)], as.numeric)

  # Return the output
  return(list(pedigree = ped, log.mendel = ped.log, log.exclu = ped.exclu, nb.mrk = ncol(offspring)))
}

# -------------------------------------------------------------------------------------------------------------------
#' Set the APIS threshold
#'
#' This function calculates the threshold for APIS
#' @param ped.log log.like for assignment function
#' @param ped.exclu log.exclu for assignment function
#' @param nb.mrk Number of markers
#' @param error (default: NULL) The assignment error rate accepted by the user
#' @keywords assignment
#' @return pedigree
#' @return log file
#' @export
setThreshold <- function(ped.log, ped.exclu, nb.mrk, error = NULL) {

  cat('===================================================', sep = '\n')
  cat('           ___   _____   _   _____  ', sep = '\n')
  cat('          /   | |  _  \\ | | /  ___/ ', sep = '\n')
  cat('         / /| | | |_| | | | | |___  ', sep = '\n')
  cat('        / / | | |  ___/ | | \\ __  \\ ', sep = '\n')
  cat('       / /  | | | |     | |  ___| | ', sep = '\n')
  cat('      /_/   |_| |_|     |_| /_____/ ', sep = '\n')
  cat('\n')
  cat('---------------------------------------------------', sep = '\n')
  cat('AUTO-ADAPTIVE PARENTAGE INFERENCE SOFTWARE', sep = '\n')
  cat('---------------------------------------------------', sep = '\n')

  # Create the pedigree output
  ped <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(ped.log)))
  colnames(ped) <- c('off', 'sire', 'dam')
  ped[,1] <- ped.log[,1]

  # Plot of Mendelian transmission probability distributions
  mendel <- rep(NA, times = 2*nrow(ped.log))
  mendel[seq(1, length(mendel), 2)] <- ped.log$mendel1
  mendel[seq(2, length(mendel), 2)] <- ped.log$mendel2

  data.mendel <- data.frame(mendel = mendel,
                            P = rep(c("P1", "P2"), times = nrow(ped.log)))

  delta <- rep(NA, times = 2*nrow(ped.log))
  delta[seq(1, length(delta), 2)] <- ped.log$delta_Pmendel12
  delta[seq(2, length(delta), 2)] <- ped.log$delta_Pmendel23

  data.delta <- data.frame(delta = delta,
                            P = rep(c("delta1", "delta2"), times = nrow(ped.log)))

  miss <- rep(NA, times = 2*nrow(ped.log))
  miss[seq(1, length(miss), 2)] <- ped.log$miss1
  miss[seq(2, length(miss), 2)] <- ped.log$miss2

  data.miss <- data.frame(miss = miss,
                           P = rep(c("miss1", "miss2"), times = nrow(ped.log)))

  if (is.null(error)) {
    error <- as.numeric(readline(prompt = 'What assignment error rate do you accept : '))
  } else {
    error <- error
  }

  # Calculate the median of P2m(:)
  med_0 <- median(ped.log$mendel2)

  mendel2_o <- sort(ped.log$mendel2)
  mendel1_o <- sort(ped.log$mendel1)

  N1_0 <- length(ped.log$mendel1[which(ped.log$mendel1<=med_0)])
  N0 <- length(ped.log$mendel1)

  while(TRUE) {
    seuil <- min(ped.log$mendel2)

    N2_l <- round(length(ped.log$mendel2[which(ped.log$mendel2<=seuil)]) -
                    (length(ped.log$mendel3[which(ped.log$mendel3<=seuil)]) * ((2*N1_0)/N0)))

    vu <- c()

    cpt <- 1
    while(N2_l<((N0-(2*N1_0))/2) & cpt<=nrow(ped.log)) {
      seuil <- mendel2_o[cpt]

      N2_l <- round(length(ped.log$mendel2[which(ped.log$mendel2<=seuil)]) -
                      (length(ped.log$mendel3[which(ped.log$mendel3<=seuil)]) * ((2*N1_0)/N0)))

      cpt <- cpt+1
    }
    med_1 <- seuil
    N1_0 <- length(ped.log$mendel1[which(ped.log$mendel1<=med_0)])
    N1_1 <- length(ped.log$mendel1[which(ped.log$mendel1<=med_1)])

    diffN <- N1_1 - N1_0
    if (diffN<=1) {
      break
    } else {
      med_0 <- med_1
    }
  }

  N1_1min <- length(which(ped.log$mendel1<=median(ped.log$mendel2)))
  N1_1 <- ifelse(test = N1_1>round(nrow(ped.log)/2), yes = round(nrow(ped.log)/2), no = N1_1)

  cat('Estimated number of offspring with at least one missing parent : between',2*N1_1min,'and',2*N1_1)
  cat('\n')

  #####-----------------------------------------------------
  ##### THRESHOLD
  #####-----------------------------------------------------

  # If the number of offspring with at least one missing parent is LOWER than the user-defined error
  if ((2*N1_1)<=round(error*nrow(ped.log))) {
    cat('--------------------------------------', sep = '\n')
    cat('     BEST MENDELIAN PROBABILITY', sep = '\n')
    cat('--------------------------------------', sep = '\n')

    ped[,2:3] <- ped.log[,3:4]
    ped$assign <- 'assign'

    thresh.mendel <- min(ped.log$delta_Pmendel12)

  } else {
    # If the number of offspring with at least one missing parent is GREATER than the user-defined error
    cat('--------------------------------------', sep = '\n')
    cat('    DELTA OF MENDELIAN PROBABILITY', sep = '\n')
    cat('--------------------------------------', sep = '\n')

    s.delta23 <- sort(ped.log$delta_Pmendel23, decreasing = T)

    thresh.mendel <- quantile(s.delta23[1:(nrow(ped.log) - 2*N1_1min)], probs = (1-error), type = 5, na.rm = T)
    cat('Threshold for delta :', thresh.mendel)
    cat('\n')

    ped[,2:3] <- ped.log[,c(3:4)]
    ped$assign <- ifelse(test = ped.log$delta_Pmendel12 >= thresh.mendel, yes = 'assign', no = 'no.assign')
  }

  # Reformate outputs
  ped.final <- ped
  ped.final[which(ped.final$assign == 'no.assign'), 2:3] <- c(NA, NA)
  ped.final <- ped.final[,-4]

  # Plot the distributions
  plot_mendel <- ggplot(data = data.mendel, aes(x = mendel, fill = P)) +
    geom_histogram(data = subset(data.mendel, P == 'P2'), bins = 30) +
    geom_histogram(data = subset(data.mendel, P == 'P1'), alpha = 0.8, bins = 30) +
    xlab(label = "average Mendelian tranmission probability") +
    ylab(label = "number of individuals") +
    theme(axis.title.x = element_text(margin = margin(20, 0, 0, 0))) +
    theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0))) +
    guides(fill = guide_legend(title = "Mendelian probability"))

  plot_delta <- ggplot(data = data.delta, aes(x = delta, fill = P)) +
    geom_histogram(data = subset(data.delta, P == 'delta2'), bins = 30) +
    geom_histogram(data = subset(data.delta, P == 'delta1'), alpha = 0.8, bins = 30) +
    geom_vline(xintercept = thresh.mendel) +
    xlab(label = "delta") +
    ylab(label = "number of individuals") +
    theme(axis.title.x = element_text(margin = margin(20, 0, 0, 0))) +
    theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0))) +
    guides(fill = guide_legend(title = "Delta"))

  plot_miss <- ggplot(data = data.miss, aes(x = miss, fill = P)) +
    geom_histogram(data = subset(data.miss, P == 'miss2'), bins = 30) +
    geom_histogram(data = subset(data.miss, P == 'miss1'), alpha = 0.8, bins = 30) +
    xlab(label = "number of mismatches") +
    ylab(label = "number of individuals") +
    theme(axis.title.x = element_text(margin = margin(20, 0, 0, 0))) +
    theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0))) +
    guides(fill = guide_legend(title = "missmatches"))

  grid.arrange(plot_delta, plot_mendel, plot_miss, nrow = 3, ncol = 1)

  # Return the output
  return(list(pedigree = ped.final, log = ped.log, error = error))
}

# -------------------------------------------------------------------------------------------------------------------
#' Estimate the allele frequencies
#'
#' This function estimates allele frequencies
#' @param genotype A matrix of genotypes (n*p)
#' n = number of individuals
#' p = number of markers (coded as "All1/All2", ex: "A/A" or "NA/NA" for missing genotype)
#' @keywords allele frequencies
#' @examples data("genotype_APIS")
#' freq <- allFreq(off)
#' @return allele frequencies
#' @export
allFreq <- function(genotype) {
  # DESCRIPTION
  # Estimate allele frequencies based on genotype matrix


  # Create the genotype matrix for new coding genotypes (2 columns)
  mat.geno <- matrix(NA, nrow = nrow(genotype), ncol = 2*ncol(genotype))

  imp <- seq(1,ncol(mat.geno),2)

  # Divide each genotype (coded A/A) into 2 columns
  for (i in c(1:ncol(genotype))) {
    tmp <- strsplit(genotype[,i], split = '/', fixed = T)
    M <- t(mapply(FUN = function(X) {X}, tmp))
    mat.geno[,(imp[i]:(imp[i]+1))] <- M
  }

  # List of the different alleles
  variant <- sort(unique(unlist(as.list(apply(mat.geno,2,unique)))))

  # Create the results matrix
  mat.res <- matrix(0, nrow = ncol(genotype), ncol = length(variant))
  rownames(mat.res) <- colnames(genotype)
  colnames(mat.res) <- variant

  for (n in 1:nrow(mat.res)) {
    tmp <- table(mat.geno[,(imp[n]:(imp[n]+1))])
    mat.res[n,match(names(tmp), colnames(mat.res))] <- tmp
  }

  # Calculte the allele frequencies
  mat.freq <- mat.res/(rowSums(mat.res[,which(colnames(mat.res)!='NA')]))
  colnames(mat.freq) <- paste0('Freq_',colnames(mat.res))

  # Merge the results
  res <- cbind(mat.res, tot = rowSums(mat.res), mat.freq)

  # Return the result
  return(res)
}

# ------------------------------------------------------------------------------------------------------------------
#' calculte the theoretical assignment power
#'
#' This function calculates the theoretical assignment power of the marker set
#' @param sire Sire genotypes | Matrix (n*p) where n = number of individuals, p = number of markers
#' rownames(sire) = labels of sires
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @param dam Dam genotypes | Matrix (n*p) where n = number of individuals, p = number of markers
#' rownames(dam) = labels of dams
#' marker coding = "All1/All2" example: "A/A", "A/B", "NA/NA" (for missing genotype)
#' @return Theoretical assignment power of the marker set
#' @examples data("genotype_APIS")
#' assignmentPower(sire_full, dam_full)
#' @keywords assignment exclusion power
#' @export
assignmentPower <- function(sire, dam) {
  # DESCRIPTION
  # This function calculates the theoretical assignment power as proposed in Vandeputte, M (2012)

  pop <- rbind(sire, dam)

  # Importe the allFreq function and calculate the allele frequencies
  freq <- allFreq(as.matrix(pop))

  col <- which(colnames(freq)=='tot')
  freq.calc <- as.data.frame(freq[,((col+1):ncol(freq))])
  freq.calc <- freq.calc[,-which(colnames(freq.calc) == "Freq_NA")]

  mcol <- ncol(freq.calc)

  # Calculate Q1 and Q3 for each marker
  freq.calc$Q1i <- 1 - 2*rowSums(freq.calc[,1:mcol]^2) +
    rowSums(freq.calc[,1:mcol]^3) + 2*rowSums(freq.calc[,1:mcol]^4) -
    2*rowSums(freq.calc[,1:mcol]^2)^2 - 3*rowSums(freq.calc[,1:mcol]^5) +
    3*rowSums(freq.calc[,1:mcol]^3)*rowSums(freq.calc[,1:mcol]^2)

  freq.calc$Q3i <- 1 + 4*rowSums(freq.calc[,1:mcol]^4) -
    4*rowSums(freq.calc[,1:mcol]^5) - 3*rowSums(freq.calc[,1:mcol]^6) -
    8*rowSums(freq.calc[,1:mcol]^2)^2 + 2*rowSums(freq.calc[,1:mcol]^3)^2 +
    8*rowSums(freq.calc[,1:mcol]^3)*rowSums(freq.calc[,1:mcol]^2)

  # Calculate the global Q1 and Q3
  Q1 <- 1 - prod(1-freq.calc$Q1i)
  Q3 <- 1 - prod(1-freq.calc$Q3i)

  # Calculate the assignment power
  Pu <- Q1^(nrow(dam)+nrow(sire)-2)*Q3^((nrow(dam)-1)*(nrow(sire)-1))

  # Return the result
  return(Pu)
}

# ------------------------------------------------------------------------------------------------------------------
#' Establish personal threshold
#'
#' This function allows the user to set up his own threshold
#' @param APIS.result APIS function output
#' @param method the method for the new threshold | 'delta' for deltas,
#' 'Pmendel' for Mendelian porbabilities, 'exclusion' for mismatches
#' @param threshold personal threshold | default values are implemented
#' @return new pedigree from the new threshold
#' @examples data("genotype_APIS")
#'
#' result <- APIS(off.genotype = off_full,
#'                sire.genotype = sire_full,
#'                dam.genotype = dam_full,
#'                error = 0.05)
#'
#' new.result <- personalThreshold(result, method = 'Pmendel')
#' @keywords assignment APIS threshold
#' @export
personalThreshold <- function(APIS.result, method, threshold = NULL) {
  # Get the result
  pedigree <- APIS.result$pedigree
  log <- APIS.result$log
  error <- APIS.result$error

  # Check the method
  if (method == "delta") {
    # DO the delta threshold
    col_toKeep <- 14
    if (is.null(threshold)) {
      return(APIS.result)
    } else {
      threshold <- threshold
    }
  } else if (method == "Pmendel") {
    # DO the Pmendel threshold
    col_toKeep <- 7
    if (is.null(threshold)) {
      threshold <- quantile(x = log[,12], probs = (1 - error), type = 5)
    } else {
      threshold <- threshold
    }
  } else if (method == "exclusion") {
    # DO the exclusion threshold
    col_toKeep <- 5
    if (is.null(threshold)) {
      threshold <- ceiling(0.05 * max(log$mrk_genotype))
    } else {
      threshold <- threshold
    }
  } else {
    stop("Invalid method")
  }

  # Re-assgin the offspring
  tmp <- log[,c(1,3:4)]

  if (method == 'exclusion') {
    tmp$assign <- ifelse(test = log[,col_toKeep] <= threshold, yes = 'assign', no = 'no.assign')
  } else {
    tmp$assign <- ifelse(test = log[,col_toKeep] >= threshold, yes = 'assign', no = 'no.assign')
  }

  tmp[which(tmp$assign == 'no.assign'), 2:3] <- c(NA, NA)
  pedigree <- tmp[,1:3]

  # Plot results
  mendel <- rep(NA, times = 2*nrow(log))
  mendel[seq(1, length(mendel), 2)] <- log$mendel1
  mendel[seq(2, length(mendel), 2)] <- log$mendel2

  data.mendel <- data.frame(mendel = mendel,
                            P = rep(c("P1", "P2"), times = nrow(log)))

  delta <- rep(NA, times = 2*nrow(log))
  delta[seq(1, length(delta), 2)] <- log$delta_Pmendel12
  delta[seq(2, length(delta), 2)] <- log$delta_Pmendel23

  data.delta <- data.frame(delta = delta,
                           P = rep(c("delta1", "delta2"), times = nrow(log)))

  miss <- rep(NA, times = 2*nrow(log))
  miss[seq(1, length(miss), 2)] <- log$miss1
  miss[seq(2, length(miss), 2)] <- log$miss2

  data.miss <- data.frame(miss = miss,
                          P = rep(c("miss1", "miss2"), times = nrow(log)))

  plot_mendel <- ggplot(data = data.mendel, aes(x = mendel, fill = P)) +
    geom_histogram(data = subset(data.mendel, P == 'P2'), bins = 30) +
    geom_histogram(data = subset(data.mendel, P == 'P1'), alpha = 0.8, bins = 30) +
    xlab(label = "average Mendelian tranmission probability") +
    ylab(label = "number of individuals") +
    theme(axis.title.x = element_text(margin = margin(20, 0, 0, 0))) +
    theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0))) +
    guides(fill = guide_legend(title = "Mendelian probability"))

  plot_delta <- ggplot(data = data.delta, aes(x = delta, fill = P)) +
    geom_histogram(data = subset(data.delta, P == 'delta2'), bins = 30) +
    geom_histogram(data = subset(data.delta, P == 'delta1'), alpha = 0.8, bins = 30) +
    xlab(label = "delta") +
    ylab(label = "number of individuals") +
    theme(axis.title.x = element_text(margin = margin(20, 0, 0, 0))) +
    theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0))) +
    guides(fill = guide_legend(title = "Delta"))

  plot_miss <- ggplot(data = data.miss, aes(x = miss, fill = P)) +
    geom_histogram(data = subset(data.miss, P == 'miss2'), bins = 30) +
    geom_histogram(data = subset(data.miss, P == 'miss1'), alpha = 0.8, bins = 30) +
    xlab(label = "number of mismatches") +
    ylab(label = "number of individuals") +
    theme(axis.title.x = element_text(margin = margin(20, 0, 0, 0))) +
    theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0))) +
    guides(fill = guide_legend(title = "missmatches"))

  if (method == 'delta') {
    # Add delta threshold

    plot_mendel <- plot_mendel +
      geom_vline(xintercept = threshold)

  } else if (method == 'Pmendel') {
    # Add Pmendel threshold

    plot_mendel <- plot_mendel +
      geom_vline(xintercept = threshold)

  } else {
    # Add mismatches threshold

    plot_miss <- plot_miss +
      geom_vline(xintercept = threshold)

  }

  grid.arrange(plot_delta, plot_mendel, plot_miss, nrow = 3, ncol = 1)

  return(list(pedigree = pedigree, log = log, error = error, threshold = threshold))
}

# ------------------------------------------------------------------------------------------------------------------
#' Select most likely parents for potent parent pairs tests
#'
#' This function allows the selection of the most likely parents for assignment, reducing computation time
#' @param off.genotype genotype of one offspring
#' @param parent.genotype genotype matrix of parent genotypes
#' @param parent.sex vector of parents sex
#' @param n.Parent vector of number of sires and dams to select
#' @return list of potential sires and dams
#' @keywords assignment APIS threshold
#' @export
selectParents <- function(off.genotype, parent.genotype, parent.sex, n.Parent) {
  # Initialize the variables
  if (length(n.Parent) == 1) {
    sire.keep <- n.Parent
    dam.keep <- n.Parent
  } else {
    sire.keep <- n.Parent[1]
    dam.keep <- n.Parent[2]
  }


  i <- 1

  output <- data.frame(parent = rownames(parent.genotype),
                       sex = parent.sex,
                       mismatch = NA)
  output$parent <- as.character(output$parent)

  # Create the probability tables
  off.geno <- strsplit(off.genotype, split = '/')
  probability_table <- vector('list', length(off.genotype))

  for (m in c(1:length(off.genotype))) {
    # If the offspring is homozygous
    if (off.geno[[m]][1] == off.geno[[m]][2]) {
      probability_table[[m]] <- c(0, 0, 1, NA)
    } else {
      probability_table[[m]] <- c(0, 0, 0, 0, 0, 1, NA)
    }
  }

  # Loop over all the parents
  for (p in output$parent) {
    p.geno <- parent.genotype[which(rownames(parent.genotype) %in% p),]
    p.geno <- strsplit(p.geno, split = '/')
    # Probability
    parent_probability <- rep(NA, length(p.geno))

    # Loop over all the markers
    for (m in c(1:length(off.genotype))) {

      off.mrk <- off.geno[[m]]
      p.mrk <- p.geno[[m]]

      if (off.mrk[1] == off.mrk[2] & off.mrk[1] == 'NA') { # If offspring is NA/NA
        parent_probability[m] <- 1
      } else if (off.mrk[1] == off.mrk[2] & off.mrk[1] != 'NA') { # If offspring is homozygous
        # Check parent genotype
        if (p.mrk[1] == p.mrk[2] & p.mrk[1] == 'NA') { # parent NA/NA
          parent_probability[m] <- probability_table[[m]][4]
        } else if (p.mrk[1] == p.mrk[2] & p.mrk[1] != 'NA') {
          if (p.mrk[1] == off.mrk[1]) { # parent A/A
            parent_probability[m] <- probability_table[[m]][1]
          } else { # parent C/C
            parent_probability[m] <- probability_table[[m]][3]
          }
        } else { # parent heterozygous
          if ((p.mrk[1] != off.mrk[1] & p.mrk[2] != off.mrk[1]) & (p.mrk[1] != off.mrk[2] & p.mrk[2] != off.mrk[2])) {
            parent_probability[m] <- probability_table[[m]][3]
          } else {
            parent_probability[m] <- probability_table[[m]][2]
          }
        }
      } else { # If offspring is heterozygous
        if (p.mrk[1] == p.mrk[2]) { # parent homzygous
          if (p.mrk[1] == 'NA') { # parent NA/NA
            parent_probability[m] <- probability_table[[m]][7]
          } else if (p.mrk[1] == off.mrk[1]) { # parent A/A
            parent_probability[m] <- probability_table[[m]][1]
          } else if (p.mrk[1] == off.mrk[2]) { # parent B/B
            parent_probability[m] <- probability_table[[m]][3]
          } else { # parent C/C
            parent_probability[m] <- probability_table[[m]][6]
          }
        } else { # parent heterozygous
          if ((p.mrk[1] == off.mrk[1] | p.mrk[2] == off.mrk[1]) & (p.mrk[1] == off.mrk[2] | p.mrk[2] == off.mrk[2])) {
            # parent A/B
            parent_probability[m] <- probability_table[[m]][2]
          } else if ((p.mrk[1] == off.mrk[1] | p.mrk[2] == off.mrk[1]) & (p.mrk[1] != off.mrk[2] | p.mrk[2] != off.mrk[2])) {
            # parent A/C
            parent_probability[m] <- probability_table[[m]][4]
          } else if ((p.mrk[1] != off.mrk[1] | p.mrk[2] != off.mrk[1]) & (p.mrk[1] == off.mrk[2] | p.mrk[2] == off.mrk[2])) {
            # parent B/C
            parent_probability[m] <- probability_table[[m]][5]
          } else {
            # parent C/C
            parent_probability[m] <- probability_table[[m]][6]
          }
        }
      }
    }

    # Calculte the average probability for the parent
    # Write the result
    output[i, 3] <- sum(parent_probability, na.rm = T)
    i <- i + 1
  }

  output <- output[order(output$mismatch),]


  sire_toKeep <- output[which(output$sex == 'M'), ]
  s.keep <- sire_toKeep[which(sire_toKeep$mismatch >= min(output$mismatch) &  sire_toKeep$mismatch <= (min(output$mismatch) + 2)),]
  select.sire <- ifelse(test = sire.keep > nrow(sire_toKeep), yes = nrow(sire_toKeep), no = sire.keep)

  if (nrow(s.keep) < sire.keep) {
    s.keep <- rbind(s.keep, sire_toKeep[c((nrow(s.keep) + 1):select.sire),])
  } else {
    s.keep <- s.keep
  }
  sire_toKeep <- s.keep$parent

  dam_toKeep <- output[which(output$sex == 'F'), ]
  d.keep <- dam_toKeep[which(dam_toKeep$mismatch >= min(output$mismatch) & dam_toKeep$mismatch <= (min(output$mismatch) + 2)),]
  select.dam <- ifelse(test = dam.keep > nrow(dam_toKeep), yes = nrow(dam_toKeep), no = dam.keep)


  if (nrow(d.keep) < dam.keep) {
    d.keep <- rbind(d.keep, dam_toKeep[c((nrow(d.keep) + 1):select.dam),])
  } else {
    d.keep <- d.keep
  }
  dam_toKeep <- d.keep$parent

  # Return the most likely parents
  return(list(sire_toKeep = sire_toKeep, dam_toKeep = dam_toKeep))
}

