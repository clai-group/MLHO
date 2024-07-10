#' Benchmark for rolling events with history/past/last timeframes vs one-hot encoding for rolling events (patients with several labels at different timestamps)
#'
#' @param MLHO.dat your dbmart goes here
#' @param labels your labeldt goes here
#' @param binarize if you want the outcome to be binary
#' @param sparsity if you want to apply sparsity
#' @param jmi if you want to do jmi
#' @param topn the number of features to be selected
#' @param patients vector of patients
#' @param multicore if you want to parallelize the jmi
#' @param encounterLevel set to true if you have multiple date assigned labels per patient
#' @param valuesToMerge set to true if you have a "value" column for phenx and want to unite these columns into one
#' @param timeBuffer a vector containing three numbers (h,p,l) describing the time buffer that should be added to the dates used to calculate the last, past and history event groups (only used when encounterLevel=TRUE).
#'
#' @return
#' @export
#'

benchmark.revolving <- function(MLHO.dat,
                      labels,
                      binarize=FALSE,
                      sparsity=NA,
                      jmi=TRUE,
                      topn=200,
                      patients,
                      multicore=FALSE,
                      valuesToMerge=FALSE,
                      timeBufffer=c(h=0,p=0,l=0)){

  require("dplyr")
  require("DT")
  require('tidyr')
  require('foreach')

  labels <- subset(labels, labels$patient_num %in% unique(MLHO.dat$patient_num))
  dbmart.atemporal <- foreach(i=1:nrow(labels), .combine=rbind) %do% {
    MLHO.dat %>% dplyr::filter(start_date <= labels[i,]$start_date - timeBufffer[3] & patient_num == labels[i,]$patient_num) %>%
      dplyr::mutate(patient_num=paste(patient_num,labels[i,]$start_date,sep='_'))
  }
  labels_atemporal <- labels %>% dplyr::mutate(patient_num=paste(patient_num,start_date,sep='_')) %>% dplyr::select(-start_date)



  dbmart.atemporal.wide <- MSMR.lite(dbmart.atemporal, labels_atemporal, binarize, sparsity, jmi, topn, unique(dbmart.atemporal$patient_num), multicore, encounterLevel=FALSE, valuesToMerge, timeBufffer)

  dbmart.temporal.wide <- MSMR.lite(MLHO.dat, labels, binarize, sparsity, jmi, topn, patients, multicore, encounterLevel=TRUE, valuesToMerge, timeBufffer)

  out <- list()
  out$atemporal.wide <- dbmart.atemporal.wide
  out$temporal.wide<- dbmart.temporal.wide

  return (out)

}
