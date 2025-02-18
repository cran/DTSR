% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DRPCA.R
\name{DRPCA}
\alias{DRPCA}
\title{Distributed Robust Principal Component Analysis (DRPCA) for Handling Missing Data}
\usage{
DRPCA(data0, data.sample, data.copy, mr, km, D)
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
\item{XDRPCA}{The imputed dataset.}
\item{RMSEDRPCA}{The Root Mean Squared Error.}
\item{MAEDRPCA}{The Mean Absolute Error.}
\item{REDRPCA}{The Relative Eelative Error.}
\item{GCVDRPCA}{Distributed DRPCA Imputation for Generalized Cross-Validation.}
\item{timeDRPCA}{The DRPCA algorithm execution time.}
}
\description{
This function performs DRPCA to handle missing data by dividing the dataset into D blocks,
applying the Robust Principal Component Analysis (RPCA) method to each block, and then combining
the results. It calculates various evaluation metrics including RMSE, MMAE, RRE, and
Generalized Cross-Validation (GCV) using different hierarchical clustering methods.
}
\examples{
# Create a sample dataset with missing values
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
result <- DRPCA(data0, data.sample, data.copy, mr, km, D)
#Print the results
print(result$XDRPCA)
}
\seealso{
\code{\link{RPCA}} for the original RPCA function.
}
\keyword{DRPCA}
\keyword{PCA}
\keyword{RPCA}
\keyword{SVD}
\keyword{imputation}
