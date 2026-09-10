#' Deprecated alias for EXNEX fitting
#'
#' `exnex_surv()` is deprecated in favour of [pooling_surv()]. It emits a
#' one-line deprecation message and forces `pooling = "exnex"`, so existing
#' calls continue to fit the EXNEX hierarchy unchanged.
#'
#' @inheritParams pooling_surv
#' @return An object of class `pooling_surv`, identical to
#'   `pooling_surv(..., pooling = "exnex")`.
#' @examples
#' \dontrun{
#' fit <- exnex_surv(survival::Surv(time, event) ~ group, data = d)
#' }
#' @export
exnex_surv <- function(x, ...) {
  message(
    "exnex_surv() is deprecated; use pooling_surv() instead. ",
    "Fitting with pooling = 'exnex'."
  )
  pooling_surv(x, ..., pooling = "exnex")
}
