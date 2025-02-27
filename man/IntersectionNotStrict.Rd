% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tileCount.R
\name{IntersectionNotStrict}
\alias{IntersectionNotStrict}
\title{Count reads overlapping genomic ranges}
\usage{
IntersectionNotStrict(
  features,
  reads,
  ignore.strand = TRUE,
  inter.feature = FALSE
)
}
\arguments{
\item{features}{A object of \link[GenomicRanges:GRanges-class]{GRanges}
representing the
feature regions to be counted.}

\item{reads}{An object that represents the data to be counted. See
\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}.
If reads are more than 1 bam files,
it should be a vector of character with full path, otherwise current
working directory
is the default directory. For paired end reads,}

\item{ignore.strand}{logical(1). ignore strand?}

\item{inter.feature}{not used. This parameter is required by
\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}.}
}
\value{
return a summarized experiment object with chromosome-level depth
information for each input sample as metadata.
}
\description{
Count reads overlapping a set of genimc features represented as
genomic ranges. This function does not work for parallel.
}
