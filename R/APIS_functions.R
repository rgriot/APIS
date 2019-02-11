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
#' @keywords assignment APIS
#' @return pedigree
#' @return a log file
#' @examples data("genotype_APIS")
#'
#' result <- APIS(off.genotype = off,
#'                sire.genotype = sire,
#'                dam.genotype = dam,
#'                error = 0.05)
#' @export
APIS <- function(off.genotype, sire.genotype, dam.genotype, error = NULL, exclusion.threshold = ncol(off.genotype)) {

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
  cat("The assignment power of your marker set is ", P)
  cat('\n')

  if (P >= 0.99) {
    cat("Theoretical assignment power : OK")
    cat('\n')
  } else {
    message("WARNING! Your marker set is not enough powerfull!")
  }

  # Assign with APIS
  assignResult 	<- assignment(off = off.genotype, sire = sire.genotype, dam = dam.genotype, thresh = exclusion.threshold)
  apisResult 		<- setThreshold(ped.log = assignResult$log.mendel, ped.exclu = assignResult$exclu, nb.mrk = assignResult$nb.mrk, error = error)

  pedigree 	<- apisResult$pedigree
  log 		<- apisResult$log

  # Give recommendations according to the data set and the results


  # Return outputs
  output <- list(pedigree = pedigree, log = log)
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
#' @keywords assignment
#' @return intermidiate pedigree
#' @return log file for Mendelian transmission probabilities
#' @return log file for exclusion
#' @export
assignment <- function(offspring, sire, dam, thresh = ncol(offspring)) {
  # DESCRIPTION
  # Function to calculate average Mendelian transmission probabilities

  # Stop if different number of markers are provided
  if (ncol(offspring)!=ncol(sire)&ncol(offspring)!=ncol(dam))
    stop('Genotypes must have the same number of markers')

  # Create the results matrix
  ped <- matrix(NA, nrow = nrow(offspring), ncol = 7)
  colnames(ped) <- c('off', 'sire1', 'dam1', 'nb_exclu1', 'log_like1', 'P_mendel1', 'assign')
  ped[,1] <- rownames(offspring)

  ped.log <- as.data.frame(matrix(NA,nrow = nrow(offspring), ncol = 21))
  colnames(ped.log) <- c('off', 'mrk_genotype','sire1', 'dam1', 'miss1', 'like1', 'mendel1',
                         'sire2', 'dam2', 'miss2', 'like2', 'mendel2', 'delta_like12', 'delta_Pmendel12',
                         'sire3', 'dam3', 'miss3', 'like3', 'mendel3', 'delta_like23', 'delta_Pmendel23')
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

  cat('Assignment')
  cat('\n')

  # Set up the cluster for parallel iteration
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)

  # Start
  A <- foreach(off=1:nrow(offspring), .multicombine = T,
               .packages = c('foreach', 'doParallel')) %dopar% { # For each offspring
                 tmp <- offspring[off,, drop = F]

                 # Create temporary results
                 res <- matrix(NA, nrow = (nrow(sire)*nrow(dam)), ncol = 5)
                 colnames(res) <- c('sire', 'dam', 'score_exclu', 'score_like', 'P_mendel')
                 res[,1] <- rep(rownames(sire), each = nrow(dam))
                 res[,2] <- rep(rownames(dam), times = nrow(sire))

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
                              .multicombine = T, .packages = c('foreach', 'doParallel')) %dopar% { # For each parents pair

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
                                r <- c(NA,NA,NA)
                                r[1] <- sum(sc_exclu) # Number of missmatch
                                r[2] <- sum(log(sc_like)) # Log likelihood
                                r[3] <- exp(sum(log(sc_like))/ped.log[off,2]) # Mendelian transimission probability

                                return(r)
                              }

                 # Working of the results
                 res[,3:5] <- t
                 res <- as.data.frame(res)
                 res$sire <- as.character(res$sire)
                 res$dam <- as.character(res$dam)
                 res$score_exclu <- as.numeric(as.character(res$score_exclu))
                 res$score_like <- as.numeric(as.character(res$score_like))
                 res$P_mendel <- as.numeric(as.character(res$P_mendel))

                 # Order by Mendelian transmission probabilities
                 res2 <- res[order(res[,5], res[,3], decreasing = T),]
                 res2 <- res2[which(res2$P_mendel!=-Inf),]

                 delta_log12 <- res2[1,4]-res2[2,4]
                 delta_log23 <- res2[2,4]-res2[3,4]

                 delta_P12 <- res2[1,5]-res2[2,5]
                 delta_P23 <- res2[2,5]-res2[3,5]

                 p_fin <- res2[1,1]
                 m_fin <- res2[1,2]

                 ped.out <- c(rownames(tmp),NA,NA,NA,NA,NA)
                 if (is.na(p_fin)&is.na(m_fin)) {
                   ped.out[2:3] <- c('no_match', 'no_match')
                   ped.out[4:6] <- c(NA, NA, NA)
                 } else {
                   ped.out[2:3] <- c(p_fin,m_fin)
                   nb_exclu <- res2[which(res2$sire==p_fin&res2$dam==m_fin),3]
                   ped.out[4:6] <- c(nb_exclu, res2[1,4], res2[1,5])
                 }

                 log.out <- unlist(c(ped.out[1], ped.log[off,2],ped.out[2:6], res2[2,1:5], delta_log12, delta_P12, res2[3,1:5], delta_log23, delta_P23))

                 cat("\r", off, "on", nrow(offspring))

                 # Order by mismatches
                 res2 <- res[order(res[,3], -res[,4], decreasing = F),]
                 res2 <- res2[which(res2$score_like!=-Inf),]


                 exclu.out <- unlist(c(log.out[1],ped.log[off,2], res2[1,1:3], res2[2,1:3], res2[3,1:3]))

                 a <- list(ped.out, log.out, exclu.out)
               }

  # Stop the parallel cluster
  stopCluster(cl)

  ped <- as.data.frame(t(as.data.frame(lapply(A, function(X) {t <- X[[1]]}))))
  colnames(ped) <- c('off', 'sire1', 'dam1', 'nb_exclu1', 'log_like1', 'P_mendel1')
  rownames(ped) <- c(1:nrow(ped))

  # Working on the data
  ped$off <- as.character(ped$off)
  ped$sire1 <- as.character(ped$sire1)
  ped$dam1 <- as.character(ped$dam1)
  ped$nb_exclu1 <- as.numeric(as.character(ped$nb_exclu1))
  ped$log_like1 <- as.numeric(as.character(ped$log_like1))
  ped$P_mendel1 <- as.numeric(as.character(ped$P_mendel1))

  # Create a log
  ped.log <- as.data.frame(t(as.data.frame(lapply(A, function(X) {t <- X[[2]]}))))
  colnames(ped.log) <- c('off', 'mrk_genotype','sire1', 'dam1', 'miss1', 'like1', 'mendel1',
                         'sire2', 'dam2', 'miss2', 'like2', 'mendel2', 'delta_like12', 'delta_Pmendel12',
                         'sire3', 'dam3', 'miss3', 'like3', 'mendel3', 'delta_like23', 'delta_Pmendel23')
  rownames(ped.log) <- c(1:nrow(ped.log))
  ped.log[,] <- sapply(ped.log[,c(1:ncol(ped.log))], as.character)
  ped.log[,c(2,5:7,10:14, 17:21)] <- sapply(ped.log[,c(2,5:7,10:14, 17:21)], as.numeric)

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
  ped <- as.data.frame(matrix(NA, ncol = 6, nrow = nrow(ped.log)))
  colnames(ped) <- c('off', 'sire', 'dam', 'miss', 'like', 'P_mendel')
  ped[,1] <- ped.log[,1]

  # Plot of Mendelian transmission probability distributions
  log.mendel1 <- sort(ped.log$mendel1)
  log.mendel2 <- sort(ped.log$mendel2)
  max.l1 <- max(max(log.mendel1), max(log.mendel2))
  min.l2 <- min(min(log.mendel1), min(log.mendel2))

  res.graph <- matrix(NA, nrow = (nrow(ped.log)+1), ncol = 3)
  colnames(res.graph) <- c('pos', 'P_mendel1', 'P_mendel2')

  slide.win <- (max.l1 - min.l2)/20
  it <- (max.l1 - min.l2)/nrow(ped.log)
  pos.in <- min.l2
  i <- 1

  while(pos.in < max.l1) {
    l.win <- pos.in + slide.win
    cpt.log1 <- length(log.mendel1[which(log.mendel1>pos.in&log.mendel1<=l.win)])
    cpt.log2 <- length(log.mendel2[which(log.mendel2>pos.in&log.mendel2<=l.win)])

    pos <- (l.win + pos.in)/2
    res.graph[i,] <- c(pos, cpt.log1, cpt.log2)

    i <- i+1
    pos.in <- pos.in + it
  }

  res.graph <- as.data.frame(na.omit(res.graph))
  res.graph$sum <- res.graph$P_mendel1 + res.graph$P_mendel2

  cat('\n')

  plot(res.graph$pos, res.graph$P_mendel1, col = 'red', type = 'l', lty = 3,
       xlim = sort(c(min.l2, max.l1)), ylim = c(0,max(res.graph$P_mendel2,res.graph$P_mendel1, res.graph$sum)),
       xlab = 'Mendelian probability', ylab = 'Frequency')
  points(res.graph$pos, res.graph$P_mendel2, col = 'blue', type = 'l', lty = 3)
  legend('topleft', legend = c('P1', 'P2'), lty = c(3,3),
         col = c('red', 'blue'))

  if (is.null(error)) {
    error <- as.numeric(readline(prompt = 'What assignment error rate do you accept : '))
  } else {

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

    ped[,2:6] <- ped.log[,3:7]
    ped$assign <- 'assigne'

  } else {
    # If the number of offspring with at least one missing parent is GREATER than the user-defined error
    cat('--------------------------------------', sep = '\n')
    cat('    DELTA OF MENDELIAN PROBABILITY', sep = '\n')
    cat('--------------------------------------', sep = '\n')

    s.delta23 <- sort(ped.log$delta_Pmendel23, decreasing = T)

    thresh.mendel <- quantile(s.delta23[1:(nrow(ped.log) - 2*N1_1min)], probs = (1-error), type = 5, na.rm = T)
    cat('Threshold for delta :', thresh.mendel)
    cat('\n')

    ped[,2:6] <- ped.log[,c(3:7)]
    ped$assign <- ifelse(test = ped.log$delta_Pmendel12 >= thresh.mendel, yes = 'assigne', no = 'non.assigne')
  }


  # Return the output
  return(list(pedigree = ped, log = ped.log, error = error))
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
#' assignmentPower(sire, dam)
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
