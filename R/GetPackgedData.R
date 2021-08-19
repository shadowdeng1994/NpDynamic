#' extract Np dynamic from *.newick and save as PackagedData.RData
#'


#' @export
# Subfunction for preprocessing the result of iqtree
GetTreeAnn <- function(FFFly){
  tmp.tree <- LoadInTree(FFFly)
  tmp.ann <- LoadInTreeAnn(FFFly)
  tmp.tip <- tmp.ann %>% slice(grep("^T",Tip))
  tmp.node <- tmp.ann %>% slice(grep("^I",Tip))

  tmp.Node2Height <- GetNode2Height(tmp.tree)
  tmp.Node2Tips <- GetNode2Tips(tmp.tree)
  tmp.Node2Organ <- GetNode2Organ(tmp.Node2Tips,tmp.tip)

  tmp.TreeAnn <- list()
  tmp.TreeAnn[["Fly"]] <- FFFly
  tmp.TreeAnn[["Tree"]] <- tmp.tree
  tmp.TreeAnn[["Node2Height"]] <- tmp.Node2Height
  tmp.TreeAnn[["Node2Tips"]] <- tmp.Node2Tips
  tmp.TreeAnn[["Node2Organ"]] <- tmp.Node2Organ
  tmp.TreeAnn[["TipInfo"]] <- tmp.tip
  tmp.TreeAnn[["NodeInfo"]] <- tmp.node

  return(tmp.TreeAnn)
}

#' @export
# Subfunction for preprocessing the ccs data
GetRawReads <- function(FFFly,TTTreeAnn){
  tmp.RawReads <- list()
  tmp.RawReads[["Fly"]] <- FFFly
  tmp.RawReads[["Tip2Seq"]] <-
    LoadInCentroid(FFFly) %>%
    right_join(TTTreeAnn$TipInfo) %>% select(Tip,Seq)

  return(tmp.RawReads)
}

# Subfunction for summarising the tips' information
GetTipInfo <- function(FFFly,TTTreeAnn,RRRawReads,RRRef){
  tmp.Tip2AllMut <- GetTip2AllMut(FFFly,RRRawReads,RRRef)
  tmp.Tip2SNP <- GetTip2SNP(tmp.Tip2AllMut,TTTreeAnn)
  tmp.Tip2Height <- GetTip2Height(TTTreeAnn)
  tmp.Height2MutNum <- GetHeight2MutNum(tmp.Tip2SNP,tmp.Tip2Height)

  tmp.TipInfo <- list()
  tmp.TipInfo[["Fly"]] <- FFFly
  tmp.TipInfo[["Tip2AllMut"]] <- tmp.Tip2AllMut
  tmp.TipInfo[["Tip2SNP"]] <- tmp.Tip2SNP
  tmp.TipInfo[["Tip2Height"]] <- tmp.Tip2Height
  tmp.TipInfo[["Height2MutNum"]] <- tmp.Height2MutNum

  return(tmp.TipInfo)
}

#' @export
# Inferring Np dynamic
GetNpDynamic <- function(FFFly,TTTreeAnn,TTTipInfo,FFFilterOrgan){
  # Split tree height into cell division according to per-generation mutation rate
  tmp.step <- 1
  tmp.window <- GetWindow(TTTipInfo,tmp.step)

  # Calculate corrected Np
  tmp.CorrectedLTT <- list()
  for(ooo in FFFilterOrgan$Organ %>% unique){
    tmp.LTT <- GetAllLTT(TTTreeAnn,TTTipInfo,tmp.step,tmp.window,ooo)
    tmp.Delta_n <- (sss+1):max(tmp.LTT$MutBin) %>%
      lapply(function(ttt){
        GetDelta_n(TTTreeAnn,tmp.LTT,tmp.step,ttt)
        }) %>% bind_rows
    tmp.CorrectedLTT[[ooo]] <- GetCorrectedLTT(FFFly,ooo,tmp.LTT,tmp.Delta_n)
  }

  tmp.NpDynamic <- list()
  tmp.NpDynamic[["Fly"]] <- FFFly
  tmp.NpDynamic[["MinWindow"]] <- tmp.window
  tmp.NpDynamic[["LTT"]] <- tmp.CorrectedLTT %>% bind_rows

  return(tmp.NpDynamic)
}

#' @export
GetPackagedData <- function(){
  message("## 1/6 ## Load in OrganInfor.ann")
  OrganInfo.ann <- LoadInOrganInfo()

  message("## 2/6 ## Get var.TreeAnn.*")
  var.TreeAnn.L5 <- GetTreeAnn("L5")
  var.TreeAnn.L6 <- GetTreeAnn("L6")

  message("## 3/6 ## Get var.TipInfo.*")
  var.Ref <- LoadInRef()
  var.RawReads.L5 <- GetRawReads("L5",var.TreeAnn.L5)
  var.RawReads.L6 <- GetRawReads("L6",var.TreeAnn.L6)
  var.TipInfo.L5 <- GetTipInfo("L5",var.TreeAnn.L5,var.RawReads.L5,var.Ref)
  var.TipInfo.L6 <- GetTipInfo("L6",var.TreeAnn.L6,var.RawReads.L6,var.Ref)

  message("## 4/6 ## Get NpDynamic.*")
  tmp.FilterOrgan <- GetFilterOrgan(var.TipInfo.L5,var.TipInfo.L6)
  var.NpDynamic.L5 <- GetNpDynamic("L5",var.TreeAnn.L5,var.TipInfo.L5,tmp.FilterOrgan)
  var.NpDynamic.L6 <- GetNpDynamic("L6",var.TreeAnn.L6,var.TipInfo.L6,tmp.FilterOrgan)

  message("## 5/6 ## Get CumeNpDynamic")
  var.cumeNp <- GetCumeNp(var.NpDynamic.L5,var.NpDynamic.L6)

  message("## 6/6 ## Output to PackagedData.RData")
  save(
    var.TreeAnn.L5,
    var.TreeAnn.L6,
    var.TipInfo.L5,
    var.TipInfo.L6,
    var.Ref,
    var.NpDynamic.L5,
    var.NpDynamic.L6,
    var.cumeNp,
    file="PackagedData.RData")
}

