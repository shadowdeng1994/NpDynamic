#' Bootrapping with 90% Sampling
#'

#' @export
# down-sampling
GetNpDynamic_90Sampling <- function(FFFly,TTTreeAnn,TTTipInfo,FFFilterOrgan){
  tmp.CorrectedLTT <- list()

  tmp.step <- 1
  tmp.window <- GetWindow(TTTipInfo,tmp.step)
  for(rrr in 1:100){
    for(ooo in FFFilterOrgan$Organ %>% unique){
      tmp.LTT <- GetAllLTT_90Sampling(TTTreeAnn,TTTipInfo,tmp.step,tmp.window,ooo)
      tmp.Delta_n <- (sss+1):max(tmp.LTT$MutBin) %>%
        lapply(function(ttt){
          GetDelta_n(TTTreeAnn,tmp.LTT,tmp.step,ttt)
        }) %>% bind_rows
      tmp.CorrectedLTT[[paste0(ooo,"_",rrr)]] <-
        GetCorrectedLTT(FFFly,ooo,tmp.LTT,tmp.Delta_n) %>%
        mutate(RID=rrr)
    }
  }

  tmp.out <- tmp.CorrectedLTT %>% bind_rows %>% mutate(Fly=FFFly)
  return(tmp.out)
}



#' @export
# Proceed down-sampling for many times
RunNpBootrap <- function(){
  message("## 1/4 ## Load in OrganInfor.ann")
  load("PackagedData.RData")
  OrganInfo.ann <- LoadInOrganInfo()
  tmp.FilterOrgan <- GetFilterOrgan(var.TipInfo.L5,var.TipInfo.L6)

  message("## 2/4 ## Bootstraping")
  var.NpDynamic_90Sampling <-
    rbind(
      GetNpDynamic_90Sampling("L5",var.TreeAnn.L5,var.TipInfo.L5,tmp.FilterOrgan),
      GetNpDynamic_90Sampling("L6",var.TreeAnn.L6,var.TipInfo.L6,tmp.FilterOrgan)
      )

  message("## 3/4 ## Cume Np dynamic")
  var.cumeNp_90Sampling <- GetCumeNp_90Sampling(var.NpDynamic_90Sampling)

  message("## 4/4 ## Output to RunNpBoot_90.RData")
  save(var.NpDynamic_90Sampling,var.cumeNp_90Sampling,file="RunNpBoot_90.RData")
}
