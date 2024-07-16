#' Downsample a signal by an integer factor
#'
#' @inheritParams gsignal::downsample
#' @return
#' @export
#' @importFrom stats is.ts is.mts deltat start
#' @importFrom gsignal downsample
#'
#' @examples
#' x=1:10
#' x
#' downsample(x, 2)
#' x3=matrix(1:30, ncol=3, dimnames=list(NULL, c("x", "y", "z")))
#' x3
#' downsample(x3, 2)
downsample <- function(x, n, phase = 0) UseMethod('downsample')

#' @rdname downsample
#' @export
downsample.default <- function(x, n, phase = 0) {
  gsignal::downsample(x, n=n, phase = phase)
}

#' @rdname downsample
#' @export
#' @examples
#' xts=stats::ts(1:10, start=0, deltat=0.1)
#' xts
#' downsample(xts, 2)
#'
#' x3ts=stats::ts(matrix(1:30, ncol=3, dimnames=list(NULL, c("x", "y", "z"))), start=0, deltat=0.1)
#' x3ts
#' downsample(x3ts, 2)
downsample.ts <- function(x, n, phase = 0) {
  stopifnot(is.ts(x))
  dx=if(is.mts(x)) as.matrix(x) else as.vector(x)
  ts(gsignal::downsample(dx, n=n, phase = phase), deltat = deltat(x)*n, start = start(x))
}
