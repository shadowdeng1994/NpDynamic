#' Subfunctions for GetTreeAnn.R
#'
#'

# load in tree file
LoadInTree <- function(FFFly){
  system.file("TreeData", paste0(FFFly,".iqtree.nwk"), package = "NpDynamic") %>%
    read.tree
}

# load in node's information
LoadInTreeAnn <- function(FFFly){
  system.file("TreeData", paste0(FFFly,".iqtree.tsv"), package = "NpDynamic") %>%
    read_tsv(col_names = F) %>% tbl_df %>%
    mutate(AAA=X2) %>% separate(AAA,LETTERS[1:3],sep="-") %>%
    select(Tip=X1,Organ=A,TID=X2,Status=X3,Filter=X4)
}

# get tree height for each internal node
GetNode2Height <- function(TTTree){
  cbind(TTTree$edge,TTTree %>% nodeHeights) %>% as.data.frame %>% tbl_df %>%
    rename(From=V1,Node=V2,FromHeight=V3,NodeHeight=V4)
}

# get tree height for each terminal node
GetNode2Tips <- function(TTTree){
  rbind(
    TTTree$tip.label %>% tbl_df %>% mutate(Node=1:n()) %>% select(Node,Tip=value),
    ((1:TTTree$Nnode)+length(TTTree$tip.label)) %>%
      lapply(function(nnn){
        extract.clade(TTTree,nnn)$tip.label %>% tbl_df %>% mutate(Node=nnn)
      }) %>% bind_rows %>% select(Node,Tip=value)
  ) %>% tbl_df
}

# get organ component for each internal node
GetNode2Organ <- function(NNNode2Tip,TTTip){
  NNNode2Tip %>% left_join(TTTip) %>%
    group_by(Node,Organ) %>% summarise(Count=n()) %>%
    group_by
}
