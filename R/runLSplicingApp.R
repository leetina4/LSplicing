# runLSplicingApp.R

#' \code{runLSplicingApp} run the shiny app \code{./inst/shiny-scripts/lsplicingApp/}.
#' @return NULL (invisible) Invoked for its side-effect of launching a shiny app.
#' @examples
#' # No runnable example applies - but BiocCheck() requires one. So ...
#' NULL
#'\dontrun{
#' runLSplicingApp()
#'}
#' @export

runLSplicingApp <- function() {
  appDir <- system.file("shiny-scripts", "lsplicingApp", package = "LSplicing")
  shiny::runApp(appDir, display.mode = "normal")
  return(invisible(NULL))
}

# [END]
