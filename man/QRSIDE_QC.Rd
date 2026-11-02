\name{QRSIDE_QC}
\alias{QRSIDE_QC}
\title{Filter spatial counts and align spot coordinates}
\usage{QRSIDE_QC(geneexpr, pos, threshold)}
\arguments{
\item{geneexpr}{Nonnegative integer count matrix with named spots in rows and genes in columns.}
\item{pos}{Two numeric coordinate columns with the same unique spot names.}
\item{threshold}{Keep rows and columns whose total counts exceed this nonnegative threshold.}
}
\description{Filters spots and genes and returns their matching spatial coordinates in the original spot order. Matrix dimensions are preserved when one row or column remains.}
\value{A list containing the filtered count matrix and aligned coordinate data frame.}
