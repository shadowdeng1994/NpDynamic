#' Assess sampling sufficiency for each organ
#'

#' @export
# Down-sampling for each organ of each individual
GetNpDynamic_DownSampling <- function(FFFly,TTTreeAnn,TTTipInfo,FFFilterOrgan){
  tmp.CorrectedLTT <- list()

  tmp.step <- 1
  tmp.window <- GetWindow(TTTipInfo,tmp.step)
  for(nnn in seq(10,1500,2)){
    for(ooo in FFFilterOrgan$Organ %>% unique){
      tmp.LTT <- GetAllLTT_DownSampling(TTTreeAnn,TTTipInfo,tmp.step,tmp.window,ooo,nnn)
      tmp.Delta_n <- (sss+1):max(tmp.LTT$MutBin) %>%
        lapply(function(ttt){
          GetDelta_n(TTTreeAnn,tmp.LTT,tmp.step,ttt)
        }) %>% bind_rows
      tmp.CorrectedLTT[[paste0(ooo,"_",nnn)]] <-
        GetCorrectedLTT(FFFly,ooo,tmp.LTT,tmp.Delta_n) %>%
        mutate(SampleSize=nnn)
    }
  }

  tmp.out <- tmp.CorrectedLTT %>% bind_rows %>% mutate(Fly=FFFly)
  return(tmp.out)
}



#' @export
# Estimate Np dynamic
RunDownSamplingNp <- function(){
  message("## 1/6 ## Load in OrganInfor.ann")
  load("PackagedData.RData")
  OrganInfo.ann <- LoadInOrganInfo()
  tmp.MostInfoOrgan <- GetMostInfoOrgan(var.TipInfo.L5,var.TipInfo.L6)

  message("## 2/6 ## Down-sampling")
  tmp.NpDynamic_DownSampling <-
    rbind(
      GetNpDynamic_DownSampling("L5",var.TreeAnn.L5,var.TipInfo.L5,tmp.FilterOrgan),
      GetNpDynamic_DownSampling("L6",var.TreeAnn.L6,var.TipInfo.L6,tmp.FilterOrgan)
    )

  message("## 3/6 ## Cume Np dynamic")
  tmp.cumeNp_DownSampling <- GetCumeNp_DownSampling(tmp.NpDynamic_DownSampling)

  message("## 4/6 ## K-S test")
  tmp.ks_test <- RunKSTest(tmp.cumeNp_DownSampling)

  message("## 5/6 ## Extract library size")
  tmp.LibSize <- GetLibSize(var.TipInfo.L5,var.TipInfo.L6)

  message("## 6/6 ## Output to Resample.RData")
  var.Resample <- list()
  var.Resample[["LTT"]] <- tmp.NpDynamic_DownSampling
  var.Resample[["CumeNp"]] <- tmp.cumeNp_DownSampling
  var.Resample[["ks_test"]] <- tmp.ks_test
  var.Resample[["LibSize"]] <- tmp.LibSize

  save(var.Resample,file="Resample.RData")
}
