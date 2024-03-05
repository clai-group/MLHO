#' MSMR lite
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
#'
#' @return
#' @export
#'

MSMSR.lite <- function(MLHO.dat,
                       labels,
                       binarize=FALSE,
                       sparsity=NA,
                       jmi=TRUE,
                       topn=200,
                       patients,
                       multicore=FALSE,
                       encounterLevel=FALSE,
                       valuesToMerge=FALSE){

  require("dplyr")
  require("DT")
  require('tidyr')

  if(valuesToMerge){
    MLHO.dat <- MLHO.dat %>%
      tidyr::unite("phenx",phenx,value, na.rm=TRUE, remove=TRUE)
  }
  MLHO.dat.table <- MLHO.dat %>%
    select(patient_num,phenx) %>%
    unique() %>%
    pull(phenx) %>%
    table()

  MLHO.dat.agg <- data.frame(phenx=names(MLHO.dat.table),distinct_patients=as.vector(MLHO.dat.table))

  if(!is.na(sparsity)){
    print("step - 1: sparsity screening!")
    ##remove low-prevalence features
    #avrs <- c(as.character(subset(MLHO.dat.agg$phenx,MLHO.dat.agg$distinct_patients > round(length(patients)*sparsity))))
    avrs <- MLHO.dat.table[MLHO.dat.table > round(length(patients) * sparsity)] %>%
      names()
    MLHO.dat <- subset(MLHO.dat,MLHO.dat$phenx %in% avrs)
  }


  ##if(MLHO.dat has Value column) -> merge with phenX
  MLHO.encounter.data = data.frame(patient_num = character(),
                                   phenx = character(),
                                   start_date=as.Date(character()))
  MLHO.encounter.data$start_date <- as.Date(MLHO.encounter.data$start_date)


  if(encounterLevel){

    print("Applying encounter based transformations!")
    MLHO.dat$start_date<-as.Date(MLHO.dat$start_date)
    labels$start_date<-as.Date(labels$start_date)

    for(patient in patients){
      #get data and labels for current patient
      patient.data <- MLHO.dat %>%
        dplyr::filter(patient_num == patient) %>%
        dplyr::select(everything())

      encounters <-labels %>%
        dplyr::filter(patient_num == patient) %>%
        dplyr::select(patient_num, start_date)

      #summarize all events before the first encounter as history
      first.encounter <- min(encounters$start_date)
      patient_labels <-labels %>%
        dplyr::filter(patient_num == patient)
      history.data <- patient.data %>%
        dplyr::filter(start_date < first.encounter) %>%
        dplyr::select(patient_num, phenx, start_date) %>%
        dplyr::mutate(phenx = paste0(phenx,'_history', sep=""))



      past.encounter <- NULL

      for(i in seq(length(encounters$start_date))){
        encounter = encounters$start_date[i]
        last.encounter <- encounter
        if(is.null(past.encounter)){
          encounter.data <- patient.data %>%
            dplyr::filter (start_date >= first.encounter & start_date <= last.encounter) %>%
            dplyr::select(patient_num, phenx, start_date) %>%
            dplyr::mutate(phenx = paste0(phenx,'_last', sep=""))
        } else{
          encounter.data <- patient.data %>%
            dplyr::filter (start_date > past.encounter & start_date <= last.encounter) %>%
            dplyr::select(patient_num, phenx, start_date) %>%
            dplyr::mutate(phenx = paste0(phenx,'_last', sep=""))
        }
        if(!is.null(past.encounter)){
          encounter.data <- patient.data %>%
            dplyr::filter (start_date <= past.encounter & start_date>=first.encounter) %>%
            dplyr::select(patient_num, phenx, start_date) %>%
            dplyr::mutate(phenx = paste0(phenx,'_past', sep="")) %>%
            rbind(encounter.data)
        }

        #set new past encounter
        past.encounter <- last.encounter
        #append history, past and last merge patient_num with encounter_date

        MLHO.encounter.data$start_date <- as.Date(MLHO.encounter.data$start_date)
        encounter.data$start_date <-as.Date(encounter.data$start_date)

        MLHO.encounter.data <- rbind(history.data, encounter.data) %>%
          dplyr::mutate(patient_num = paste0(patient,"_",encounter)) %>%
          rbind(MLHO.encounter.data)
      }
    }
    #wide table
    MLHO.dat.wide <- MLHO.encounter.data %>%
      dplyr::group_by(patient_num,phenx) %>%
      dplyr::summarise(n=n(),.groups = "drop") %>%
      tidyr::pivot_wider(id_cols = "patient_num",names_from = "phenx",values_from = n,values_fill = 0)


    labels <- labels %>%
      dplyr::mutate(patient_num = paste0(patient_num,"_" ,start_date)) %>%
      dplyr::select(-start_date)
  }else{ #old code
    MLHO.dat.wide <- MLHO.dat %>%
      dplyr::group_by(patient_num,phenx) %>%
      dplyr::summarise(n=n(),.groups = "drop") %>%
      tidyr::pivot_wider(id_cols = "patient_num",names_from = "phenx",values_from = n,values_fill = 0)
    MLHO.dat.wide <- MLHO.dat.wide[, !(names(MLHO.dat.wide) %in% c("NA"))]
  }

  AVR <- MLHO.dat.wide

  if(jmi == TRUE){
    if(multicore==TRUE){
      ###setup parallel backend for case 4 loops
      cores<-detectCores()
      cl <- makeCluster(cores[1]-2)
      registerDoParallel(cl)
    }

    print("step 2: JMI dimensionality reduction!")
    # AVR <- merge(AVR,labels,by="patient_num")
    AVR <- dplyr::left_join(AVR,labels,by="patient_num")
    require("praznik")
    # AVR <- as.data.frame(AVR)
    ##calculating JOINT mutual information
    JMIs.AVR <- as.data.frame(JMI(AVR[,!(names(AVR) %in% names(labels))],
                                  as.factor(AVR$label),k = dim(AVR)[2]-ncol(labels), threads = 0)$score)

    JMIs.AVR$features <- rownames(JMIs.AVR)
    rownames(JMIs.AVR) <- NULL
    colnames(JMIs.AVR) <- c("JMI.score","phenx")

    # JMIs.AVR <- merge(JMIs.AVR,MLHO.dat.agg,by="phenx",all.x = TRUE)
    JMIs.AVR <- dplyr::left_join(JMIs.AVR,MLHO.dat.agg,by="phenx")
    JMIs.AVR$JMI.score <- round(JMIs.AVR$JMI.score,3)
    JMIs.AVR$rank <- rank(order(order(JMIs.AVR$JMI.score,-JMIs.AVR$distinct_patients)))


    topfeatures.AVR <- subset(JMIs.AVR,JMIs.AVR$rank <= topn)
    topfeatures.AVR <- c(as.character(unique(topfeatures.AVR$phenx)))
    AVR <- AVR[, names(AVR) %in% topfeatures.AVR | names(AVR) %in% names(labels)]

  }

  if (jmi==FALSE){
    # AVR <- merge(AVR,labels,by="patient_num")
    AVR <- dplyr::left_join(AVR,labels,by="patient_num")
  }

  #binarize?
  if(binarize==TRUE){
    AVR[,2:dim(AVR)[2]] <- +(AVR[,2:dim(AVR)[2]] > 0)
  }

  AVR <- as.data.frame(AVR)

  return(AVR)

}
