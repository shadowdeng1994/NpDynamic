#' Reconstruct under variant mutation rates
#'

#' @export
GetNpDynamic_VarMutRate <- function(FFFly,TTTreeAnn,TTTipInfo,FFFilterOrgan){
  tmp.CorrectedLTT <- list()

  # Analysis with different mutation rate.
  for(MMMutRate in seq(0.5,1.5,0.25)){
    tmp.step <- 1
    tmp.window <- GetWindow(TTTipInfo,tmp.step)*MMMutRate
    for(ooo in FFFilterOrgan$Organ %>% unique){
      tmp.LTT <- GetAllLTT(TTTreeAnn,TTTipInfo,tmp.step,tmp.window,ooo)
      tmp.Delta_n <- (sss+1):max(tmp.LTT$MutBin) %>%
        lapply(function(ttt){
          GetDelta_n(TTTreeAnn,tmp.LTT,tmp.step,ttt)
        }) %>% bind_rows
      tmp.CorrectedLTT[[paste0(ooo,"_",MMMutRate)]] <-
        GetCorrectedLTT(FFFly,ooo,tmp.LTT,tmp.Delta_n) %>%
        mutate(MutRate=MMMutRate)
    }
  }

  tmp.out <- tmp.CorrectedLTT %>% bind_rows %>% mutate(Fly=FFFly)
  return(tmp.out)
}



#' @export
RunNpVarMutRate <- function(){
  message("## 1/4 ## Load in OrganInfor.ann")
  load("PackagedData.RData")
  OrganInfo.ann <- LoadInOrganInfo()
  tmp.FilterOrgan <- GetFilterOrgan(var.TipInfo.L5,var.TipInfo.L6)

  message("## 2/4 ## Bootstraping")
  var.NpDynamic_VarMutRate <-
    rbind(
      GetNpDynamic_VarMutRate("L5",var.TreeAnn.L5,var.TipInfo.L5,tmp.FilterOrgan),
      GetNpDynamic_VarMutRate("L6",var.TreeAnn.L6,var.TipInfo.L6,tmp.FilterOrgan)
    )

  message("## 3/4 ## Cume Np dynamic")
  var.cumeNp_VarMutRate <- GetCumeNp_VarMutRate(var.NpDynamic_VarMutRate)

  message("## 4/4 ## Output to NpVarMutRate.RData")
  save(var.NpDynamic_VarMutRate,var.cumeNp_VarMutRate,file="NpVarMutRate.RData")
}
