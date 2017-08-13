shinyInfl <- function() {
  appDir <- system.file("shiny", package = "reverseR")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `reverseR`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}