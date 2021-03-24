#' Subfunctions for FitWithSShapeCurve.R
#'

FitGrowthCurve <- function(FFFly,FFFilterOrgan){
  FFFilterOrgan %>% .$Organ %>% unique %>%
    lapply(function(ooo){
      tmp1 <- var.cumeNp %>% filter(Fly==FFFly,Organ==ooo)
      tmp2 <- growthcurver::SummarizeGrowth(tmp1$MutBin,tmp1$CumeNp,bg_correct="none")
      rbind(
        tmp2$vals %>%
          unlist %>% as.data.frame %>% rename(Value='.') %>%
          rownames_to_column("ColName") %>%
          mutate(Group="Parameter"),
        data.frame(
          ColName=tmp2$data$t,
          Value=tmp2$model %>% predict %>% as.numeric,
          Group="Predict"
        )
      ) %>% mutate(Fly=FFFly,Organ=ooo)
    }) %>% bind_rows
}
