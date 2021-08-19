#' Subfunctions for GetPackagedData.R
#'

# load in organ information
LoadInOrganInfo <- function(){
  tmp.Mark <-
    c("Legs",
      "Appendages",
      "Nervous system",
      "Digestive system",
      "Reproductive system",
      "Others",
      "")
  system.file("TreeData", "OrganAnnotation.tsv", package = "NpDynamic") %>%
    read.delim(col.names=c("Stage","FullName","Chinese","Organ","Germlayer","Mark")) %>%
    tbl_df %>%
    mutate(Mark=factor(Mark,tmp.Mark))
}

GetMostInfoOrgan <- function(TTTipInfo_L5,TTTipInfo_L6){
  rbind(
    TTTipInfo_L5$Tip2Height %>%
      group_by(Organ) %>% summarise(Fly="L5",Count=n()),
    TTTipInfo_L6$Tip2Height %>%
      group_by(Organ) %>% summarise(Fly="L6",Count=n())
  ) %>%
    # Select the most informatic organs
    right_join(
      rbind(TTTipInfo_L5$Tip2Height %>% mutate(Fly="L5"),TTTipInfo_L6$Tip2Height %>% mutate(Fly="L6")) %>%
        group_by(Fly,Organ) %>% summarise(Count=n()) %>%
        group_by(Organ) %>% filter(n()==2,all(Count>100)) %>% summarise())
}

GetFilterOrgan <- function(TTTipInfo_L5,TTTipInfo_L6){
  GetMostInfoOrgan(TTTipInfo_L5,TTTipInfo_L6) %>%
    # exclude sampling-insufficient organs
    filter(Organ!="Mg",Organ!="Vc")
}


