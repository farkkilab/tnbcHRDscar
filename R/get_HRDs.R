#' New HRD scar score

get.tnbcHRDscars <- function(seg, chrominfo = "grch38", LOH_windos=c(10,30), LST_segSize=5e6, LST_mindistance=2e6){
  seg <- preparing.input(seg)
  
  if (chrominfo == "grch38"){
    chrom = chrominfo_grch38
  }
  if (chrominfo == "grch37"){
    chrom = chrominfo_grch37
  }
  
  #Calculating HRD-LOH
  HRD_LOHs <- features.LOH(seg, MbSizes=LOH_windos)
  HRD_LOHs <- HRD_LOHs[,1]

  #Calculating LSTs
  LSTs <- LSTs(seg, chrominfo=chrom, segsizes=LST_segSize, mindistance=LST_mindistance)

  #Calculating nTAIs
  res_ai<- calc.ai_new(seg, chrom)
  nTAIs <- res_ai[,1]

  tnbcHRDscar <- HRD_LOHs + LSTs + nTAIs
  
  HR.status <- ifelse(tnbcHRDscar >= 53, "HRD","HRP")

  #Concatenating results
  HRDresulst <- cbind(HRD_LOHs, LSTs, nTAIs, tnbcHRDscar, HR.status)
  colnames(HRDresulst) <- cbind("nLOH","LSTs","nTAIs", "tnbcHRDscar","HR.status")
  assign("HRDresulst",as.data.frame(HRDresulst),envir = .GlobalEnv)
  return(HRDresulst)
}