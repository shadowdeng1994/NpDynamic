#' Simulation of coalescent method with different population structure
#'


#' @export
GetSimTree_Panmictic <- function(FFFile,DDDynamic_1,DDDynamic_2,MMMutRate){
  message("## 1/6 ## Load in Ref")
  var.Ref <- LoadInRef()
  tmp.Pool <- GetMutPool()

  message("## 2/6 ## Get Np dynamic")
  tmp.Dynamic_1 <- DDDynamic_1
  tmp.Dynamic_2 <- DDDynamic_2
  tmp.Dynamic <- tmp.Dynamic_1+tmp.Dynamic_2
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
GetSimTree_Structured <- function(FFFile,DDDynamic_1,DDDynamic_2,MMMutRate){
  message("## 1/6 ## Load in Ref")
  var.Ref <- LoadInRef()
  tmp.Pool <- GetMutPool()

  message("## 2/6 ## Get Np dynamic")
  tmp.Dynamic_1 <- DDDynamic_1
  tmp.Dynamic_2 <- DDDynamic_2
  tmp.DynamicTable_1 <- GetNpDynamicTable(tmp.Dynamic_1)
  tmp.DynamicTable_2 <- GetNpDynamicTable(tmp.Dynamic_2)

  message("## 3/6 ## Get cell division evenets and mutations")
  tmp.Division_1 <- AssignDivisionAndMutations(tmp.DynamicTable_1,tmp.Pool,MMMutRate)
  tmp.DivisionTable_1 <- GetDivisionTable_Structured(tmp.Division_1,"p1")
  tmp.Division_2 <- AssignDivisionAndMutations(tmp.DynamicTable_2,tmp.Pool,MMMutRate)
  tmp.DivisionTable_2 <- GetDivisionTable_Structured(tmp.Division_2,"p2")

  message("## 4/6 ## Get *.newick")
  tmp.NewickString_1 <- GetNewickString_Structured(tmp.DivisionTable_1,"p1")
  tmp.NewickString_2 <- GetNewickString_Structured(tmp.DivisionTable_2,"p2")
  tmp.Tree_1 <- read.newick(text = tmp.NewickString_1)
  tmp.Tree_2 <- read.newick(text = tmp.NewickString_2)
  tmp.Tree <- bind.tree(tmp.Tree_1,tmp.Tree_2)

  message("## 5/6 ## Get Node2Height")
  tmp.Node2Height <- GetNode2Height(tmp.Tree) %>%
    left_join(tmp.Tree$tip.label %>% tbl_df %>% rename(Tip=value) %>% mutate(Node=1:n()))

  message("## 6/6 ## Output")
  save(tmp.Dynamic_1,tmp.Dynamic_2,
       tmp.NewickString_1,tmp.NewickString_2,
       tmp.Tree,
       tmp.Node2Height,
       file=paste0(FFFile,".RData"))
}

#' @export
RunCoalescentSim_PopStructure <- function(){
  tmp.NpDynamic.p1 <- 2^c(1:10,rep(10,15),11:15)
  tmp.NpDynamic.p2 <- 2^c(rep(1,10),1:15,rep(15,5))

  GetSimTree_Panmictic("Panmictic",tmp.NpDynamic.p1,tmp.NpDynamic.p2,1)
  GetSimTree_Structured("Structured",tmp.NpDynamic.p1,tmp.NpDynamic.p2,1)
}
