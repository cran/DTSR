% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DTSR.R
\name{DTSR}
\alias{DTSR}
\title{Distributed Trimmed Scores Regression (DTSR) for Handling Missing Data}
\usage{
DTSR(data0, data.sample, data.copy, mr, km, D)
}
\arguments{
\item{data0}{The original dataset containing the response variable and features.}

\item{data.sample}{The dataset used for sampling, which may contain missing values.}

\item{data.copy}{A copy of the original dataset, used for comparison or validation.}

\item{mr}{Indices of the rows with missing values that need to be predicted.}

\item{km}{The number of clusters for k-means clustering.}

\item{D}{The number of blocks to divide the data into.}
}
\value{
A list containing:
\item{XDTSR}{The imputed dataset.}
\item{RMSEDTSR}{The Root Mean Squared Error.}
\item{MAEDTSR}{The Mean Absolute Error.}
\item{REDTSR}{The Relative Eelative Error.}
\item{GCVDTSR}{The DTSR for Generalized Cross-Validation.}
\item{timeDTSR}{The DTSR algorithm execution time.}
}
\description{
This function performs DTSR to handle missing data by dividing the dataset into D blocks,
applying the Trimmed Scores Regression (TSR) method to each block, and then combining
the results. It calculates various evaluation metrics including RMSE, MMAE, RRE, and
Consistency Proportion Index (CPP) using different hierarchical clustering methods.
}
\examples{
# Create a sample matrix with random values and introduce missing values
set.seed(123)
n <- 100
p <- 10
D <- 2
data.sample <- matrix(rnorm(n * p), nrow = n)
data.sample[sample(1:(n-10), (p-2))] <- NA
data.copy <- data.sample
data0 <- data.frame(data.sample, response = rnorm(n))
mr <- sample(1:n, 10)  # Sample rows for evaluation
km <- 3  # Number of clusters
# Perform DTSR imputation
result <- DTSR(data0, data.sample, data.copy, mr, km,D)
# Print the results
print(result$XDTSR)
}
\seealso{
\code{\link{TSR}} for the original TSR function.
}
\keyword{DTSR}
\keyword{PCA}
\keyword{SVD}
\keyword{TSR}
\keyword{imputation}
