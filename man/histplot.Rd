% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/histplot.R
\name{histplot}
\alias{histplot}
\title{histplot}
\usage{
histplot(
  x,
  color = "#99d6e1",
  rug_color = "#f79999",
  freq = FALSE,
  lwd = 2,
  rug_lwd = 1,
  main = "",
  xlab = "",
  ylab = "",
  cex.main = 1.5,
  cex.lab = 1.5,
  cex.axis = 1.5,
  breaks = "Sturges",
  ylim = NULL,
  xpos = 0.03,
  ypos = 0,
  cex.text = 1.2
)
}
\arguments{
\item{x}{Input phenotype data.}

\item{color}{The color of histogram.}

\item{rug_color}{The color of rug under the histogram.}

\item{freq}{If FALSE, the histogram graphic is a representation of frequencies; if TRUE, the histogram graphic is a representation of probability densitie. (Default: FALSE)}

\item{lwd}{The line width of histogram. (Default: 2)}

\item{rug_lwd}{The line width of rug under the histogram. (Default: 1)}

\item{main}{The title of plot.}

\item{xlab}{The X axis labels.}

\item{ylab}{The Y axis labels}

\item{cex.main}{The magnification to be used for title. (Default: 1.5)}

\item{cex.lab}{The magnification to be used for axis labels. (Default: 1.5)}

\item{cex.axis}{The magnification to be used for axis annotation. (Default: 1.5)}

\item{breaks}{The number of bars in the histogram.}

\item{ylim}{Y axis ranges.}

\item{xpos}{The horizontal position of the pvalue label. (Default: 0.03)}

\item{ypos}{The vertical position of the pvalue label. (Default: 0)}

\item{cex.text}{The magnification to be used for pvalue labels. (Default: 1.2)}
}
\value{
Histogram and p-value of Shapiro-Wilk Normality Test.
}
\description{
Histogram drawing
}
\examples{
data("wheatds")
inlier <- outlier(wheatds, sample = "Line", loc = "Env", rep = "Rep", phe = "DS", mode = "blup")
stat_out <- stat(x = inlier, sample = "Sample", phe = "inlier")
histplot(x = stat_out$mean)
}
\author{
Peng Zhao <pengzhao@nwafu.edu.cn>
}
