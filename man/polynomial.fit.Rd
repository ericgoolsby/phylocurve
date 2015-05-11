\name{polynomial.fit}
\alias{polynomial.fit}
\title{
Fit polynomial curves to species data using stepwise regression
}
\description{
Transforms raw data into regression curves using stepwise (AIC or BIC) polynomial regression. The maximum number of parameters (nterms), response data can be constrained between minima and maxima (for example, the default sets any negative predicted y value to 0).
}
\usage{
polynomial.fit(data, x_variable, y_variable, method = "BIC",
  nterms = 2, min_x = -Inf, max_x = Inf, min_y = 0, max_y = Inf, eval_length = 30)
}
\arguments{
  \item{data}{
Data frame with function-valued species data along with a column named "species"
}
  \item{x_variable}{
The name of the column in \code{data} corresponding to the predictor (x) variable
}
  \item{y_variable}{
The name of the column in \code{data} corresponding to the response (y) variable
}
  \item{method}{
Metric for stepwise polynomial regression. Either "AIC" or "BIC" (the default).
}
  \item{nterms}{
Maximum number of terms to allow in polynomial regression.
}
  \item{min_x}{
Lower range to remove from raw x data
}
  \item{max_x}{
Upper range to remove from raw x data
}
  \item{min_y}{
Lower range to remove from predicted response (for example, for response data in which only positive values are allowed, the min_y can be set to 0 (the default)).
}
  \item{max_y}{
Upper range to remove from predicted response.
}
  \item{eval_length}{
How many points to evaluate along the range of x.
}
}
\value{
\item{X }{X (formatted for phylocurve.generalized)}
\item{Y }{Y (formatted for phylocurve.generalized)}
}
\examples{
x_length <- 50
nspecies <- 25
sim_data <- sim.curves(nspecies = nspecies,x_length = x_length)

# Fit polynomial curves to raw data
poly_data <- polynomial.fit(data = sim_data$data,x_variable = "x",y_variable = "y",min_y = 0)

# Align curve data using dynamic time warping
# NOT RUN (takes several minutes)
# P <- phylocurve.generalized(tree = sim_data$tree,X = poly_data$X,Y = poly_data$Y)

# Perform comparative analyses on aligned curves
# rate.mult.fitted <- rate.mult(tree = P$tree,Y = P$aligned_data)
# K.mult(rate.mult.fitted) 
}