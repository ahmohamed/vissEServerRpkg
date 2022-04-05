#' @export
runExample <- function(method) {
  ds = c(ora="ora_example", gsea="gsea_example", visium="visium_example", sc="sc_example")[[method]]
  getdata(ds)
}
