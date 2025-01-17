% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterMotifs.R
\name{clusterMotifs}
\alias{clusterMotifs}
\title{Assign each motif to a cluster}
\usage{
clusterMotifs(
  motifs,
  type = c("PPM", "ICM"),
  method = c("PCC", "EUCL", "SW", "KL", "ALLR", "BHAT", "HELL", "SEUCL", "MAN",
    "ALLR_LL", "WEUCL", "WPCC"),
  power = 1,
  agglom = "complete",
  thresh = 0.2,
  return_d = FALSE,
  plot = FALSE,
  labels = FALSE,
  cex = 1,
  linecol = "red",
  ...
)
}
\arguments{
\item{motifs}{A list of universalmotifs or a list of PWMs}

\item{type}{Can be ICM or PPM}

\item{method}{The method to be used for determining similarity/distances}

\item{power}{Raise correlation matrices to this power before converting to a
distance matrix. Only applied if method is either "PCC" or "WPCC"}

\item{agglom}{Method to be used for agglomeration by \link[stats]{hclust}}

\item{thresh}{Tree heights below which motifs are formed into a cluster}

\item{return_d}{logical(1) Return the distance matrices for each cluster}

\item{plot}{Show tree produced by \link[stats]{hclust}. If requested the
value set by thresh will be shown as a horizontal line}

\item{labels, cex}{Passed to \link[stats]{plot.hclust}}

\item{linecol}{Passed to \link[graphics]{abline} as the argument \code{col}}

\item{...}{passed to \link[universalmotif]{compare_motifs}}
}
\value{
Named vector with numeric values representing which cluster each
motif has been assigned to.

If setting return_d = TRUE, a named list will be returned with the clusters
as the element \code{cl} and a list with distance matrices for each cluster as the
element \code{d}
}
\description{
Cluster related motifs for testing as a group
}
\details{
This builds on \link[universalmotif]{compare_motifs}, enabling the
assignment of each PWM to a cluster, and subsequent testing of motifs as
a cluster, rather than returning individual results.

Internally all matrices are converted to distance matrices and
\link[stats]{hclust} is used to form clusters. By default, options such as
"EUCL", "MAN" produce distances, whilst similarity matrices are produced
when choosing "PCC" and other correlation based metrics. In these cases, the
distance matrix is obtained by taking 1 - similarity.

By default PWM labels are hidden (labels = FALSE), however these can be shown
using labels = NULL as explained in \link[stats]{plot.hclust}.

Raising the threshold will lead to fewer, larger clusters whilst leaving
this value low will return a more conservative approach, with more smaller
clusters.
The final decision as the best clustering strategy is highly subjective and
left to the user.
Manual inspection of motifs within a cluster can be performed using
\link[universalmotif]{view_motifs}, as shown in the vignette.
}
\examples{
# Load the example motifs
data("ex_pfm")

# Return a vector with each motif assigned a cluster
# The default uses Pearson's Correlation Coefficient
clusterMotifs(ex_pfm)

# Preview the settings noting that showing labels can clutter the plot
# with large numbers of motifs. The defaults for Euclidean distance
# show the threshold may need raising
clusterMotifs(ex_pfm, plot = TRUE, labels = NULL, method = "EUCL")

}
