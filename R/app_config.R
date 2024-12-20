#' Retrieve a System File Path for the App
#'
#' Provides the system file path for resources within the `zAMPExplorer` package.
#'
#' @param ... Additional arguments passed to `system.file()` to specify file paths.
#' @return A character string containing the requested file path.
#' @noRd
app_sys <- function(...) {
  system.file(..., package = "zAMPExplorer")
}

#' Read Application Configuration
#'
#' Retrieves configuration values from the `golem-config.yml` file.
#'
#' @param value The configuration key to retrieve.
#' @param config The active configuration profile. Defaults to the value of the
#'   `GOLEM_CONFIG_ACTIVE` environment variable. If unset, uses `R_CONFIG_ACTIVE`, or `"default"`.
#' @param use_parent Logical. If `TRUE`, searches parent directories for the configuration file.
#'   Default is `TRUE`.
#' @param file The path to the configuration file. Defaults to `golem-config.yml` in the package.
#' @return The requested configuration value.
#' @noRd
get_golem_config <- function(
    value,
    config = Sys.getenv(
      "GOLEM_CONFIG_ACTIVE",
      Sys.getenv(
        "R_CONFIG_ACTIVE",
        "default"
      )
    ),
    use_parent = TRUE,
    file = app_sys("golem-config.yml")) {
  config::get(
    value = value,
    config = config,
    file = file,
    use_parent = use_parent
  )
}
