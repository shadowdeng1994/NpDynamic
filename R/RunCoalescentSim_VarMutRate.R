#' Simulation of coalescent method under different mutation rates
#'




#' @export
GetNpDynamic_Sim_VarMutRate <- function(DDData,MMMutRate){
  message("## 1/2 ## Load in SimTree")
  load(paste0(DDData,".RData"))

  tmp.step <- 1
  tmp.window <- tmp.step*MMMutRate

  message("## 2/2 ## Extract subtree and reconstruct Np dynamic")
  c(1/1000,5/1000,10/1000,50/1000) %>%
  lapply(function(nnn){
    tmp.Subtree <- GetSubTreeWithNode2Height(tmp.Tree,tmp.Node2Height)

    tmp.SubNode2Height <- GetNode2Height(tmp.Subtree)
    tmp.SubNode2Tip <- GetNode2Tips(tmp.Subtree)

    tmp.SubLTT <- GetAllLTT_Sim(tmp.SubNode2Height,tmp.step,tmp.window)
    tmp.Delta_n <- (sss+1):max(tmp.SubLTT$MutBin) %>%
      lapply(function(ttt){
        GetDelta_n_Sim(tmp.SubNode2Tip,tmp.SubLTT,tmp.step,ttt)
        }) %>% bind_rows
    tmp.CorrectedLTT <- GetCorrectedLTT_Sim(tmp.SubLTT,tmp.Delta_n)

    tmp.out <- tmp.CorrectedLTT %>% mutate(Coverage=nnn)
    return(tmp.out)
  }) %>% bind_rows %>% mutate(Data=DDData,Window=tmp.window)
}


RunCoalescentSim_VarMutRate <- function(){
  tmp.Dynamic <- 2^c(0:14,rep(15,5),14:5)

  var.Sim_VarMutRate <-
  rbind(
    GetNpDynamic_Sim_VarMutRate("Data_1",1),
    GetNpDynamic_Sim_VarMutRate("Data_2",0.75),
    GetNpDynamic_Sim_VarMutRate("Data_3",0.5),
    GetNpDynamic_Sim_VarMutRate("Data_4",0.25)
  )

  save(tmp.Dynamic,var.Sim_VarMutRate,file="Sim_VarMutRate.RData")
}







