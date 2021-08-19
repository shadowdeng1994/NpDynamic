#' Estimate tissue similarity with Gst
#'


#' @export
GetGst <- function(){
  # load in preprocessed tree.
  load("PackagedData.RData")

  # Extract the organs' componence for each node
  var.PopCompOnTree <- list()
  var.PopCompOnTree[["L5"]] <- GetPopComponent(var.TreeAnn.L5)
  var.PopCompOnTree[["L6"]] <- GetPopComponent(var.TreeAnn.L6)
  
  # Extract the organs' componence for each node after shuffling
  var.PopCompOnTree_Shuffle <- list()
  var.PopCompOnTree_Shuffle[["L5"]] <- GetPopComponent_Shuffle(var.TreeAnn.L5)
  var.PopCompOnTree_Shuffle[["L6"]] <- GetPopComponent_Shuffle(var.TreeAnn.L6)

  # Calculate Gst
  var.Gst <-
    rbind(
      GetGstWithPopComp("L5",var.PopCompOnTree),
      GetGstWithPopComp("L6",var.PopCompOnTree)
    )

  # Calculate Gst after shuffling
  var.Gst_Shuffle <-
    rbind(
      GetGstWithPopComp("L5",var.PopCompOnTree_Shuffle),
      GetGstWithPopComp("L6",var.PopCompOnTree_Shuffle)
    )
  
  # Calcualte one-to-others Gst
  var.Gst_Organ <-
    rbind(
      GetGstWithPopComp_Organ("L5",var.PopCompOnTree),
      GetGstWithPopComp_Organ("L6",var.PopCompOnTree)
    )

  save(var.Gst,var.Gst_Shuffle,var.Gst_Organ,file="Gst.RData")
}

