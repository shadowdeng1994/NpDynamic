#' Subfunctions for Build Tree with Np dynamic
#'


# Get Np dynamic table
GetNpDynamicTable <- function(DDDynamic){
  DDDynamic %>% tbl_df %>%
    rename(NpSize=value) %>%
    mutate(ApproSize=round(DDDynamic),Generation=1:n())
}

# Get mutation pool for mutation assignment
GetMutPool <- function(RRRef){
  RRRef %>%
    mutate(Alt=NA,Alt=replace(Alt,Base=="C","T"),Alt=replace(Alt,Base=="G","A")) %>%
    filter(!is.na(Alt)) %>%
    unite(Base,Posi,Alt,col="MutInfo",sep="_") %>%
    mutate(MutID=1:n())
}

# Simulation cell division and mutation accumulation
AssignDivisionAndMutations <- function(DDDynamicTable,PPPool,MMMutRate){
  tmp.Division <- list()
  tmp.Division[[1]] <- data.frame("CellID"=1:DDDynamicTable$NpSize[1]) %>% tbl_df %>%
    mutate(Generation=1,Layer=1:n(),Layer=Layer-mean(Layer),
           Parent=0,ParentGeneration=0,ParentLayer=0,AllParent=0 %>% as.character) %>%
    mutate(MutNew=sample(1:nrow(PPPool),1,replace = T) %>% as.character,MutNew=replace(MutNew,runif(n())>MMMutRate,NA),MutAnn=MutNew)
  for(iii in 2:nrow(DDDynamicTable)){
    tmp.Division[[iii]] <-
      tmp.Division[[iii-1]] %>%
      select(Parent=CellID,ParentGeneration=Generation,ParentLayer=Layer,ParentMut=MutAnn,AllParent) %>%
      slice(sample(1:n(),DDDynamicTable$NeSize[iii-1])) %>%
      mutate(MutNew=sample(1:nrow(PPPool),n(),replace = T) %>% as.character,MutNew=replace(MutNew,runif(n())>MMMutRate,NA),HAHA=MutNew) %>%
      unite(ParentMut,HAHA,col="MutAnn",sep=":") %>%
      mutate(AllParent=paste(AllParent,Parent,sep=":")) %>%
      rbind(.,.) %>% arrange(Parent) %>%
      mutate(CellID=1:n(),CellID=CellID+nrow(tmp.Division %>% bind_rows),Generation=iii,Layer=1:n(),Layer=Layer-mean(Layer))
  }
  tmp.out <- tmp.Division %>% bind_rows

  return(tmp.out)
}

# Simulate cell division
GetDivisionTable <- function(DDDivision){
  DDDivision %>% arrange(Parent) %>%
    mutate(Tip=paste0("C",CellID,"T:0")) %>%
    mutate(Inter=paste0("C",Parent,"T:0")) %>%
    group_by(Inter,Generation,MutNew) %>% summarise(Newick=paste(Tip,collapse = ",")) %>%
    group_by %>% mutate(LLL=1,LLL=replace(LLL,is.na(MutNew),0),Newick=paste0("(",Newick,"):",LLL)) %>%
    arrange(Generation)
}

# Transfer into newick format
GetNewickString <- function(DDDivisionTable){
  tmp.NewickString <- "(C0T:0);"
  for(iii in 1:nrow(DDDivisionTable)){
    tmp.NewickString <- gsub(DDDivisionTable$Inter[iii],DDDivisionTable$Newick[iii],tmp.NewickString)
  }
  return(tmp.NewickString)
}

# Merging subpopulation
GetDivisionTable_Structured <- function(DDDivision,PPPop){
  DDDivision %>% arrange(Parent) %>%
    mutate(Tip=paste0("C",CellID,"T_",PPPop,":0")) %>%
    mutate(Inter=paste0("C",Parent,"T_",PPPop,":0")) %>%
    group_by(Inter,Generation,MutNew) %>% summarise(Newick=paste(Tip,collapse = ",")) %>%
    group_by %>% mutate(LLL=1,LLL=replace(LLL,is.na(MutNew),0),Newick=paste0("(",Newick,"):",LLL)) %>%
    arrange(Generation)
}

# Transfer into newick format
GetNewickString_Structured <- function(DDDivisionTable,PPPop){
  tmp.NewickString <- paste0("(C0T_",PPPop,":0);")
  for(iii in 1:nrow(DDDivisionTable)){
    tmp.NewickString <- gsub(DDDivisionTable$Inter[iii],DDDivisionTable$Newick[iii],tmp.NewickString)
  }
  return(tmp.NewickString)
}








