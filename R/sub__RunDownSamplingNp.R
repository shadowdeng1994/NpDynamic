#' Subfunctions for RunDownSamplingNp.R
#'

RunKSTest <- function(CCcumeNp_DownSampling){
  CCcumeNp_DownSampling %>%
    select(Fly,Organ,SampleSize,MutBin,CumeNp) %>%
    group_by(Fly,Organ) %>%
    left_join(filter(.,dense_rank(-SampleSize)==1),by=c("Fly","Organ","MutBin")) %>%
    group_by(Fly,Organ,SampleSize.x,SampleSize.y) %>% filter(n()>20) %>% summarise(pval=ks.test(CumeNe.x,CumeNe.y)$p.value) %>%
    group_by(Organ) %>% mutate(padj=p.adjust(pval,method="fdr")) %>% arrange(padj) %>%
    group_by %>% mutate(Sig=(padj<0.05))
}


GetLibSize <- function(NNNpDynamic_L5,NNNpDynamic_L6){
  rbind(NNNpDynamic_L5$Tip2Height %>% mutate(Fly="L5"),
        NNNpDynamic_L6$Tip2Height %>% mutate(Fly="L6")) %>%
    group_by(Fly,Organ) %>% summarise(Count=n()) %>%
    group_by(Organ) %>% filter(all(Count>100),n()==2)
}
