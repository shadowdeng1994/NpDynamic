#' Subfunctions for GetGst.R
#'


GetPopComponent <- function(TTTreeAnn){
  TTTreeAnn$Tree %>% nodeHeights %>%
    as.data.frame %>% .$V1 %>% unique %>%
    lapply(function(hhh){
      tmp1 <- TTTreeAnn$Tree %>% split_tree_at_height(HHH) %>% .$subtrees
      tmp2 <-
        1:length(tmp1) %>% lapply(function(III){
          tmp1[[III]] %>% .$tree %>% .$tip.label %>% tbl_df %>%
            rename(Tip=value) %>%
            left_join(TTTreeAnn$TipInfo %>% select(Tip,Organ)) %>%
            group_by(Organ) %>% summarise(CladeID=III,Count=n())
        }) %>% bind_rows %>% mutate(CutHeight=HHH)
      }) %>% bind_rows
}

GetPopComponent_Shuffle <- function(TTTreeAnn){
  tmp.shffule <- TTTreeAnn$TipInfo %>% mutate(Organ=sample(Organ))
  TTTreeAnn$Tree %>% nodeHeights %>%
    as.data.frame %>% .$V1 %>% unique %>%
    lapply(function(hhh){
      tmp1 <- TTTreeAnn$Tree %>% split_tree_at_height(HHH) %>% .$subtrees
      tmp2 <-
        1:length(tmp1) %>% lapply(function(III){
          tmp1[[III]] %>% .$tree %>% .$tip.label %>% tbl_df %>%
            rename(Tip=value) %>%
            left_join(tmp.shffule %>% select(Tip,Organ)) %>%
            group_by(Organ) %>% summarise(CladeID=III,Count=n())
        }) %>% bind_rows %>% mutate(CutHeight=HHH)
    }) %>% bind_rows
}


GetGstWithPopComp <- function(FFFly,PPPopcomp){
  PPPopcomp[[FFFly]] %>% .$CutHeight %>% unique %>%
    lapply(function(hhh){
      tmp <-
        PPPopcomp %>% .[[FFFly]] %>%
        filter(CutHeight==hhh) %>%
        spread(Organ,Count,0) %>% gather(Organ,Count,-CladeID,-CutHeight) %>%
        arrange(CutHeight,CladeID) %>%
        group_by(CutHeight) %>% mutate(TotalSize=sum(Count)) %>%
        group_by(CutHeight,CladeID) %>% mutate(W_i=sum(Count)/TotalSize,P_k_in_i=Count/sum(Count))

      left_join(
        group_by(tmp,CutHeight,CladeID,W_i) %>% summarise(J_i=sum(P_k_in_i^2)) %>%
          group_by(CutHeight) %>% summarise(J_sub=sum(W_i*J_i)),
        group_by(tmp,CutHeight,Organ) %>% summarise(P_k_in_all_2=sum(W_i*P_k_in_i)^2) %>%
          group_by(CutHeight) %>% summarise(J_total=sum(P_k_in_all_2))
      ) %>%
        mutate(Gst=(J_sub-J_total)/(1-J_total)) %>%
        mutate(Fly=FFFly,CutHeight=hhh)
    }) %>% bind_rows
}


GetGstWithPopComp_Organ <- function(FFFly,PPPopcomp){
  tmp2 <- list()
  for(ooo in PPPopcomp[[FFFly]] %>% .$Organ %>% unique){
    tmp1 <- PPPopcomp[[FFFly]] %>%
      mutate(OOODet=Organ==ooo,Organ=c("Others",ooo)[OOODet+1]) %>%
      group_by(Organ,CladeID,CutHeight) %>% summarise(Count=sum(Count)) %>%
      arrange(CutHeight,CladeID)

    tmp.out[[ooo]] <-
      tmp1 %>% .$CutHeight %>% unique %>%
      lapply(function(hhh){
        tmp <-
          tmp1 %>%
          filter(CutHeight==hhh) %>%
          spread(Organ,Count,0) %>% gather(Organ,Count,-CladeID,-CutHeight) %>%
          arrange(CutHeight,CladeID) %>%
          group_by(CutHeight) %>% mutate(TotalSize=sum(Count)) %>%
          group_by(CutHeight,CladeID) %>% mutate(W_i=sum(Count)/TotalSize,P_k_in_i=Count/sum(Count))

        left_join(
          group_by(tmp,CutHeight,CladeID,W_i) %>% summarise(J_i=sum(P_k_in_i^2)) %>%
            group_by(CutHeight) %>% summarise(J_sub=sum(W_i*J_i)),
          group_by(tmp,CutHeight,Organ) %>% summarise(P_k_in_all_2=sum(W_i*P_k_in_i)^2) %>%
            group_by(CutHeight) %>% summarise(J_total=sum(P_k_in_all_2))
        ) %>%
          mutate(Gst=(J_sub-J_total)/(1-J_total)) %>%
          mutate(Fly=FFFly,Organ=ooo,CutHeight=hhh)
      }) %>% bind_rows
  }

  tmp.out <- tmp2 %>% bind_rows
  return(tmp.out)
}
















