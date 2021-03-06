#' Compute (total) weighted hypervolume given a set of rectangles
#' 
#' The function `whv_rect()` calculates the hypervolume weighted by a set of rectangles (with zero weight outside the rectangles). The function `total_whv_rect()` calculates the total weighted hypervolume as `hypervolume() + scalefactor * abs(prod(reference - ideal)) * whv_rect()`. The details of the computation are given by \citet{DiaLop2020ejor}.
#' 
#' @template arg_data
#'
#' @param rectangles (`matrix()`) weighted rectangles that will bias the
#'   computation of the hypervolume. Maybe generated by [eafdiff()] with
#'   `rectangles=TRUE` or by [choose_eafdiff()].
#'
#' @template arg_refpoint
#' 
#' @template arg_maximise
#'
#' @details
#'   TODO
#'
#' @return  A single numerical value.
#' 
#' @seealso    [read_datasets()], [eafdiff()], [choose_eafdiff()], [whv_hype()]
#' 
#' @examples
#'
#'
#' rectangles <- as.matrix(read.table(header=FALSE, text='
#'  1.0  3.0  2.0  Inf    1
#'  2.0  3.5  2.5  Inf    2
#'  2.0  3.0  3.0  3.5    3
#' '))
#' whv_rect (matrix(2, ncol=2), rectangles, reference = 6)
#' whv_rect (matrix(c(2, 1), ncol=2), rectangles, reference = 6)
#' whv_rect (matrix(c(1, 2), ncol=2), rectangles, reference = 6)
#'
#' @references
#' \insertAllCited{}
#'
#'@export
#'@md
whv_rect <- function(data, rectangles, reference, maximise = FALSE)
{
  data <- check_dataset(data)
  nobjs <- ncol(data)
  npoints <- nrow(data)
  if (is.null(reference)) stop("reference cannot be NULL")
  if (length(reference) == 1) reference <- rep_len(reference, nobjs)
  # FIXME: This is wrong for maximisation
  stopifnot(maximise == FALSE)
  # FIXME: Do this in C code!
  rectangles_a <- rectangles[,c(1,3), drop=FALSE]
  rectangles_a[rectangles_a > reference[1]] <- reference[1]
  rectangles_b <- rectangles[,c(2,4), drop=FALSE]
  rectangles_b[rectangles_b > reference[2]] <- reference[2]
  rectangles[,c(1,3)] <- rectangles_a
  rectangles[,c(2,4)] <- rectangles_b
  # Remove empty rectangles maybe created above.
  rectangles <- rectangles[ (rectangles[,1] != rectangles[,3]) & (rectangles[,2] != rectangles[,4]),
                         , drop = FALSE]
  rectangles_nrows <- nrow(rectangles)

  if (nobjs != 2) stop("sorry: only 2 objectives supported")
    
  if (ncol(rectangles) != 5) stop("rectangles: invalid number of columns")
  
  if (any(maximise)) {
    if (length(maximise) == 1) {
      data <- -data
      reference <- -reference
      rectangles[,1:4] <- -rectangles[,1:4, drop = FALSE]
      
    } else if (length(maximise) != nobjs) {
      stop("length of maximise must be either 1 or ncol(data)")
    }
    data[,maximise] <- -data[,maximise]
    reference[maximise] <- -reference[maximise]
    pos <- as.vector(matrix(1:4, nrow=2)[,maximise])
    rectangles[,pos] <- -rectangles[,pos]
  }
  return(.Call(rect_weighted_hv2d_C,
               as.double(t(data)),
               as.integer(npoints),
               as.double(t(rectangles)),
               as.integer(rectangles_nrows)))
}


#' @template arg_ideal
#' 
#' @param scalefactor (`numeric(1)`) real value within \eqn{(0,1]} that scales
#'   the overall weight of the differences. This is parameter psi (\eqn{\psi}) in \citet{DiaLop2020ejor}.
#'
#' @examples
#' total_whv_rect (matrix(2, ncol=2), rectangles, reference = 6, ideal = c(1,1))
#' total_whv_rect (matrix(c(2, 1), ncol=2), rectangles, reference = 6, ideal = c(1,1))
#' total_whv_rect (matrix(c(1, 2), ncol=2), rectangles, reference = 6, ideal = c(1,1))
#'
#'@rdname whv_rect
#'
#' @export
#'@md
total_whv_rect <- function(data, rectangles, reference, maximise = FALSE, ideal = NULL, scalefactor = 0.1)
{
  nobjs <- ncol(data) 
  maximise <- as.logical(rep_len(maximise, nobjs))
  if (nobjs != 2) stop("sorry: only 2 objectives supported")
  if (ncol(rectangles) != 5) stop("rectangles: invalid number of columns")
  if (scalefactor <= 0 || scalefactor > 1) stop("scalefactor must be within (0,1]")

  hv <- hypervolume(data, reference, maximise = maximise)
  whv <- whv_rect(data, rectangles, reference, maximise = maximise)
  if (is.null(ideal)) {
    # FIXME: Should we include the range of the rectangles here?
    minmax <- apply(data, 2, range)
    lower <- minmax[1,]
    upper <- minmax[2,]
    ideal <- ifelse(maximise, upper, lower)
  }
  if (length(ideal) != nobjs) {
    stop("ideal should have same length as nobjs")
  }
  beta <- scalefactor * abs(prod(reference - ideal))
  #cat("beta: ", beta, "\n")
  return (hv + beta * whv)
}

#' Approximation of the (weighted) hypervolume by Monte-Carlo sampling
#' 
#' Return an estimation of the hypervolume of the space dominated by the input
#' data following the procedure described by \citet{AugBadBroZit2009gecco}. A
#' weight distribution describing user preferences may be specified.
#'
#' @template arg_data
#'
#' @template arg_refpoint
#' 
#' @template arg_maximise
#'
#' @template arg_ideal
#'
#' @param nsamples (`integer(1)`) number of samples for Monte-Carlo sampling.
#'
#' @param dist (`list()`) weight distribution. See Details.
#'
#' @details
#' A weight distribution  \citep{AugBadBroZit2009gecco} can be provided via the `dist` argument. The ones currently supported are:
#'  * `type="point"` describes a goal in the objective space, where `mu` gives the coordinates of the goal. The resulting weight distribution is a multivariate normal distribution centred at the goal. 
#' * `type="exponential"` describes an exponential distribution with rate parameter `1/mu`, i.e., \eqn{\lambda = \frac{1}{\mu}}.
#'
#' @return A single numerical value.
#'
#' @references
#' \insertAllCited{}
#' 
#' @seealso    [read_datasets()], [eafdiff()], [whv_rect()]
#' 
#' @examples
#' 
#' whv_hype (matrix(2, ncol=2), reference = 4, ideal = 1)
#'
#' whv_hype (matrix(c(3,1), ncol=2), reference = 4, ideal = 1)
#'
#' whv_hype (matrix(2, ncol=2), reference = 4, ideal = 1,
#'           dist = list(type="exponential", mu=0.2))
#'
#' whv_hype (matrix(c(3,1), ncol=2), reference = 4, ideal = 1,
#'           dist = list(type="exponential", mu=0.2))
#'
#' whv_hype (matrix(2, ncol=2), reference = 4, ideal = 1,
#'           dist = list(type="point", mu=c(1,1)))
#'
#' whv_hype (matrix(c(3,1), ncol=2), reference = 4, ideal = 1,
#'           dist = list(type="point", mu=c(1,1)))
#'
#'@export
#'@md
whv_hype <- function(data, reference, ideal, maximise = FALSE,
                     dist = list(type = "uniform"), nsamples = 1e5L)
{
  data <- check_dataset(data)
  nobjs <- ncol(data)
  npoints <- nrow(data)
  if (is.null(reference)) stop("reference cannot be NULL")
  if (length(reference) == 1) reference <- rep_len(reference, nobjs)
  if (is.null(ideal)) stop("ideal cannot be NULL")
  if (length(ideal) == 1) ideal <- rep_len(ideal, nobjs)
  if (nobjs != 2) {
    stop("sorry: only 2 objectives supported")
  }

  if (any(maximise)) {
    if (length(maximise) == 1) {
      data <- -data
      reference <- -reference
      ideal <- -ideal
    } else if (length(maximise) != nobjs) {
      stop("length of maximise must be either 1 or ncol(data)")
    }
    data[,maximise] <- -data[,maximise]
    reference[maximise] <- -reference[maximise]
    ideal[maximise] <- -ideal[maximise]
  }
  seed <- get_seed()
  return(.Call(whv_hype_C,
               as.double(t(data)),
               as.integer(npoints),
               as.double(ideal),
               as.double(reference),
               dist,
               as.integer(seed),
               as.integer(nsamples)))
}

get_seed <- function() sample.int(.Machine$integer.max, 1)
