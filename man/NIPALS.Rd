% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NIPALS.R
\name{NIPALS}
\alias{NIPALS}
\title{NIPALS Algorithm with RPCA and Clustering}
\usage{
NIPALS(data0, data.sample, data.copy, mr, km)
}
\arguments{
\item{data0}{The original dataset containing the response variable and features.}

\item{data.sample}{The dataset used for sampling, which may contain missing values.}

\item{data.copy}{A copy of the original dataset, used for comparison or validation.}

\item{mr}{Indices of the rows with missing values that need to be predicted.}

\item{km}{The number of clusters for k-means clustering.}
}
\value{
A list containing:
\item{Xnew}{The imputed dataset.}
\item{RMSE}{The Root Mean Squared Error.}
\item{CPP1}{The K-means clustering Consistency Proportion Index.}
\item{CPP2}{The Hierarchical Clustering Complete Linkage Consistency Proportion Index.}
\item{CPP3}{The Hierarchical Clustering Single Linkage Consistency Proportion Index.}
\item{CPP4}{The Hierarchical Clustering Average Linkage Consistency Proportion Index.}
\item{CPP5}{The Hierarchical Clustering Centroid linkage Consistency Proportion Index.}
\item{CPP6}{The Hierarchical Clustering Median Linkage Consistency Proportion Index.}
\item{CPP7}{The Hierarchical Clustering Ward's Method Consistency Proportion Index.}
\item{timeNIPALS}{The NIPALS algorithm execution time.}
}
\description{
This function performs the NIPALS (Nonlinear Iterative Partial Least Squares) algorithm to handle missing data by imputing the missing values based on the correlation structure within the data. It also calculates the RMSE and Consistency Proportion Index (CPP) using different hierarchical clustering methods.
}
\examples{
# Create a sample matrix with random values and introduce missing values
set.seed(123)
n <- 100
p <- 5
data.sample <- matrix(rnorm(n * p), nrow = n)
data.sample[sample(1:(n*p), 20)] <- NA
data.copy <- data.sample
data0 <- data.frame(data.sample, response = rnorm(n))
mr <- sample(1:n, 10)  # Sample rows for evaluation
km <- 3  # Number of clusters
# Perform NIPALS imputation
result <- NIPALS(data0, data.sample, data.copy, mr, km)
# Print the results
print(result$RMSE)
print(result$CPP1)
print(result$Xnew)

}
\seealso{
\code{\link{princomp}} and \code{\link{svd}} for more information on PCA and SVD.
}
\keyword{NIPALS}
\keyword{PCA}
\keyword{RPCA}
\keyword{SVD}
\keyword{imputation}
