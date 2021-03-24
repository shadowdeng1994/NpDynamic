#' Fitting Np dynamic with S-shape curve
#'


#' @export
FitWithSShapeCurve <- function(){
  message("## 1/4 ## Load in OrganInfor.ann")
  load("PackagedData.RData")
  tmp.FilterOrgan <- GetFilterOrgan(var.TipInfo.L5,var.TipInfo.L6)

  var.SShapeFitting <-
    rbind(
      FitGrowthCurve("L5",tmp.FilterOrgan),
      FitGrowthCurve("L6",tmp.FilterOrgan)
    )

  save(var.SShapeFitting,file="SShapeFitting.RData")
}
