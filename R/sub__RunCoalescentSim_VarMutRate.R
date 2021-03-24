#' Subfunctions for RunCoalescentSim_VarMutRate.R
#'

GetSubTreeWithNode2Height <- function(TTTree,NNNode2Height,RRRatio){
  NNNode2Height %>%
    select(Node,NodeHeight,Tip) %>% filter(!is.na(Tip)) %>%
    .$Tip %>% sample(round((TTTree %>% .$tip.label %>% length)*RRRatio)) %>%
    get_subtree_with_tips(TTTree,.,force_keep_root = T) %>% .$subtree
}



GetAllLTT_Sim <- function(SSSubNode2Height,SSStep,WWWindow){
  SSSubNode2Height %>% fun.GetLTT(.,WWWindow-1e-5,SSSubNode2Height$NodeHeight %>% max-1e-5,WWWindow) %>%
    mutate(MutBin=dense_rank(Height))
}

GetDelta_n_Sim <- function(SSSubNode2Tip,SSSubLTT,SSStep,ttt){
  tmp.out <- SSSubLTT %>% filter(MutBin==ttt-SSStep|MutBin==ttt) %>%
    left_join(SSSubNode2Tip,by="Node") %>%
    group_by(Tip) %>% filter(n()==2) %>%
    group_by(TTT=c("AncLTT","LTT")[dense_rank(MutBin)]) %>% summarise(LTT=length(Node %>% unique)) %>%
    spread(TTT,LTT) %>% summarise(MutBin=ttt,LTT=LTT,Delta_n=LTT-AncLTT,Step=SSStep)
  if(tmp.out$Delta_n==0){
    while(tmp.out$Delta_n==0){
      SSStep <- SSStep+1
      if((ttt-SSStep)<=0){ break }
      tmp.out <- SSSubLTT %>% filter(MutBin==ttt-SSStep|MutBin==ttt) %>%
        left_join(SSSubNode2Tip,by="Node") %>%
        group_by(Tip) %>% filter(n()==2) %>%
        group_by(TTT=c("AncLTT","LTT")[dense_rank(MutBin)]) %>% summarise(LTT=length(Node %>% unique)) %>%
        spread(TTT,LTT) %>% summarise(MutBin=ttt,LTT=LTT,Delta_n=LTT-AncLTT,Step=SSStep)
    }
  }
  return(tmp.out)
}

GetCorrectedLTT_Sim <- function(LLLTTT,DDDelta_n){
  LLLTTT %>% group_by(MutBin,Height) %>% summarise %>% group_by %>%
    left_join(DDDelta_n) %>%
    filter(!is.na(LTT)) %>%
    mutate(Score=2*Delta_n/LTT/(LTT-1),NorScore=Score/Step*sss,Np=(1/NorScore+1)/2)
}


















