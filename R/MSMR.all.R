#' MSMR all
#'
#' @param MLHO.dat your dbmart goes here
#' @param labels your labeldt goes here
#' @param multiExposure if you have multiple exposures/treatments
#' @param benchmark description
#' @param binarize if you want the outcome to be binary
#' @param sparsity if you want to apply sparsity
#' @param jmi if you want to do jmi
#' @param topn the number of features to be selected
#' @param patients vector of patients
#' @param multicore if you want to parallelize the jmi
#' @param encounterLevel set to true if you have multiple date assigned labels per patient
#' @param valuesToMerge set to true if you have a "value" column for phenx and want to unite these columns into one
#'
#' @return
#' @export
#'

MSMR.all <- function(MLHO.dat = dat.train,
                      labels = test_label,
                      multiExposure = TRUE,
                      benchmark = F, 
                      binarize=FALSE,
                      sparsity=NA,
                      jmi=TRUE,
                      topn=200,
                      patients = uniqpats.train,
                      multicore=FALSE,
                      encounterLevel=FALSE,
                      valuesToMerge=FALSE,
                      timeBufffer=c(h=0,p=0,l=0,o=0)){
  if (multiExposure ==T){
    if (benchmark==T){
      dbmart.wide = benchmark.revolving(MLHO.dat, 
                                      labels, 
                                      binarize, 
                                      sparsity, 
                                      jmi, 
                                      topn, 
                                      uniqpats.train, 
                                      multicore,  
                                      valuesToMerge, 
                                      timeBufffer)
    }else{
      dbmart.wide <- MSMR.me(MLHO.dat, 
                            labels, 
                            binarize, 
                            sparsity, 
                            jmi, 
                            topn, 
                            uniqpats.train, 
                            multicore, 
                            encounterLevel=T, 
                            valuesToMerge, 
                            timeBufffer)
    }
  }else{
    dbmart.wide <- MSMR.lite(MLHO.dat, 
                             labels, 
                             binarize, 
                             sparsity, 
                             jmi, 
                             topn, 
                             uniqpats.train, 
                             multicore)
  }
  
  return (dbmart.wide)
}



