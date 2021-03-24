#' Build Tree with Np dynamic
#'

#' @export
GetTreeWithNpDynamic <- function(FFFile,DDDynamic,MMMutRate){
  message("## 1/6 ## Load in Ref")
  var.Ref <- LoadInRef()
  tmp.Pool <- GetMutPool()

  message("## 2/6 ## Get Np dynamic")
  tmp.Dynamic <- DDDynamic
  tmp.DynamicTable <- GetNpDynamicTable(tmp.Dynamic)

  message("## 3/6 ## Get cell division evenets and mutations")
  tmp.Division <- AssignDivisionAndMutations(tmp.DynamicTable,tmp.Pool,MMMutRate)
  tmp.DivisionTable <- GetDivisionTable(tmp.Division)

  message("## 4/6 ## Get *.newick")
  tmp.NewickString <- GetNewickString(tmp.DivisionTable)
  tmp.Tree <- read.newick(text = tmp.NewickString)

  message("## 5/6 ## Get Node2Height")
  tmp.Node2Height <- GetNode2Height(tmp.Tree) %>%
    left_join(tmp.Tree$tip.label %>% tbl_df %>% rename(Tip=value) %>% mutate(Node=1:n()))

  message("## 6/6 ## Output")
  save(tmp.Dynamic,tmp.NewickString,tmp.Tree,tmp.Node2Height,
       file=paste0(FFFile,".RData"))
}

#' @export
GetSimTree_VarMutRate <- function(){
  tmp.Dynamic <- 2^c(0:14,rep(15,5),14:5)
  GetTreeWithNpDynamic("Data_1",tmp.Dynamic,1)
  GetTreeWithNpDynamic("Data_2",tmp.Dynamic,0.75)
  GetTreeWithNpDynamic("Data_3",tmp.Dynamic,0.5)
  GetTreeWithNpDynamic("Data_4",tmp.Dynamic,0.25)
}



