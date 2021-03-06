#' Subfunctions for GetNpDynamic.R
#'

# Get window size with per-generation mutation rate
GetWindow <- function(TTTipInfo,SSStep){
  1/TTTipInfo$Height2MutNum$coefficients[[2]]/SSStep
}

# Get LTT with tree height
GetLTT <- function(NNNode2Height,From,To,By){
  seq(From,To,By) %>%
    lapply(function(hhh){
      NNNode2Height %>% filter(FromHeight<hhh,NodeHeight>=hhh) %>% group_by(Node) %>% summarise(Height=hhh)
    }) %>% bind_rows
}

# Get all LTT from Q0 to Q95
GetAllLTT <- function(TTTreeAnn,TTTipInfo,sss,www,ooo){
  TTTreeAnn$Node2Organ %>%
    filter(Organ==ooo) %>%
    group_by(Node) %>% summarise %>%
    left_join(TTTreeAnn$Node2Height %>%
                mutate(Add=0,Add=replace(Add,Node<=(TTTreeAnn$Tree %>% .$tip.label %>% length),www*1.01)) %>%
                mutate(NodeHeight=NodeHeight+Add)) %>%
    # from 0 to Q95
    GetLTT(.,0,TTTipInfo$Tip2Height %>% filter(Organ==ooo) %>% .$NodeHeight %>% quantile(0.95),www) %>%
    mutate(MutBin=dense_rank(Height))
}

# down-sampling
GetAllLTT_90Sampling <- function(TTTreeAnn,TTTipInfo,sss,www,ooo){
  TTTipInfo$Tip2Height %>%
    filter(Organ==ooo) %>%
    slice(sample(1:n(),round(n()*0.9))) %>%
    left_join(TTTreeAnn$Node2Tips) %>%
    group_by(Node) %>% summarise %>%
    left_join(TTTreeAnn$Node2Height %>%
                mutate(Add=0,Add=replace(Add,Node<=(TTTreeAnn$Tree %>% .$tip.label %>% length),www*1.01)) %>%
                mutate(NodeHeight=NodeHeight+Add)) %>%
    # from 0 to Q95
    GetLTT(.,0,TTTipInfo$Tip2Height %>% filter(Organ==ooo) %>% .$NodeHeight %>% quantile(0.95),www) %>%
    mutate(MutBin=dense_rank(Height))
}

# Get all LTT after down-sampling
GetAllLTT_DownSampling <- function(TTTreeAnn,TTTipInfo,sss,www,ooo,nnn){
  TTTipInfo$Tip2Height %>%
    filter(Organ==ooo) %>%
    slice(sample(1:n(),nnn)) %>%
    left_join(TTTreeAnn$Node2Tips) %>%
    group_by(Node) %>% summarise %>%
    left_join(TTTreeAnn$Node2Height %>%
                mutate(Add=0,Add=replace(Add,Node<=(TTTreeAnn$Tree %>% .$tip.label %>% length),www*1.01)) %>%
                mutate(NodeHeight=NodeHeight+Add)) %>%
    # from 0 to Q95
    GetLTT(.,0,TTTipInfo$Tip2Height %>% filter(Organ==ooo) %>% .$NodeHeight %>% quantile(0.95),www) %>%
    mutate(MutBin=dense_rank(Height))
}

# Calcuating minimum LTT decrease on each generation
GetDelta_n <- function(TTTreeAnn,LLLTT,sss,ttt){
  tmp.step <- sss
  # LTT decrease between g and g-1
  tmp.out <- LLLTT %>% filter(MutBin==ttt-tmp.step|MutBin==ttt) %>%
    left_join(TTTreeAnn$Node2Tips,by="Node") %>%
    group_by(Tip) %>% filter(n()==2) %>%
    group_by(TTT=c("AncLTT","LTT")[dense_rank(MutBin)]) %>% summarise(LTT=length(Node %>% unique)) %>%
    spread(TTT,LTT) %>% summarise(MutBin=ttt,LTT=LTT,Delta_n=LTT-AncLTT,Step=tmp.step)
  # if D(g,g-1) equal to 0, move to g-2. Keep going until D(g,g-j) over 0 or g-j less 0.
  while(tmp.out$Delta_n==0){
    tmp.step <- tmp.step+1
    if((ttt-tmp.step)<=0){ break }
    tmp.out <- LLLTT %>% filter(MutBin==ttt-tmp.step|MutBin==ttt) %>%
      left_join(TTTreeAnn$Node2Tips,by="Node") %>%
      group_by(Tip) %>% filter(n()==2) %>%
      group_by(TTT=c("AncLTT","LTT")[dense_rank(MutBin)]) %>% summarise(LTT=length(Node %>% unique)) %>%
      spread(TTT,LTT) %>% summarise(MutBin=ttt,LTT=LTT,Delta_n=LTT-AncLTT,Step=tmp.step)
  }
  return(tmp.out)
}

# Transform LTT to Np.
GetCorrectedLTT <- function(FFFly,OOOrgan,LLLTTT,DDDelta_n){
  LLLTTT %>% group_by(MutBin,Height) %>% summarise %>% group_by %>%
    left_join(DDDelta_n) %>%
    filter(!is.na(LTT)) %>%
    mutate(Score=2*Delta_n/LTT/(LTT-1),NorScore=Score/Step*sss,Np=(1/NorScore+1)/2) %>%
    mutate(Fly=FFFly,Organ=OOOrgan)
}

# Calculating cumulative sum of Np
GetCumeNp <- function(NNNpDynamic_L5,NNNpDynamic_L6){
  rbind(NNNpDynamic_L5$LTT,NNNpDynamic_L6$LTT) %>%
    filter(MutBin>2) %>%
    filter(is.finite(Np)) %>%
    arrange(Fly,Organ,MutBin) %>%
    group_by(Fly,Organ) %>% mutate(CumeNp=cumsum(Np)) %>%
    group_by
}

# down-sampling
GetCumeNp_90Sampling <- function(NNNpDynamic_90Sampling){
  NNNpDynamic_90Sampling %>%
    filter(MutBin>2) %>%
    filter(is.finite(Np)) %>%
    arrange(Fly,Organ,RID,MutBin) %>%
    group_by(Fly,Organ,RID) %>% mutate(CumeNp=cumsum(Np)) %>%
    group_by
}

# donw-sampling
GetCumeNp_DownSampling <- function(NNNpDynamic_DownSampling){
  NNNpDynamic_DownSampling %>%
    filter(MutBin>2) %>%
    filter(is.finite(Np)) %>%
    arrange(Fly,Organ,SampleSize,MutBin) %>%
    group_by(Fly,Organ,SampleSize) %>% mutate(CumeNp=cumsum(Np)) %>%
    group_by
}

# Calculating under variant mutation rate.
GetCumeNp_VarMutRate <- function(NNNpDynamic_VarMutRate){
  NNNpDynamic_VarMutRate %>%
    filter(MutBin>2) %>%
    filter(is.finite(Np)) %>%
    arrange(Fly,Organ,MutRate,MutBin) %>%
    group_by(Fly,Organ,MutRate) %>% mutate(CumeNp=cumsum(Np)) %>%
    group_by
}




