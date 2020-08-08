#' Histogram drawing
#'
#' @title histplot
#' @param x Input phenotype data.
#' @param color The color of histogram.
#' @param rug_color The color of rug under the histogram.
#' @param freq If FALSE, the histogram graphic is a representation of frequencies; if TRUE, the histogram graphic is a representation of probability densitie. (Default: FALSE)
#' @param lwd The line width of histogram. (Default: 2)
#' @param rug_lwd The line width of rug under the histogram. (Default: 1)
#' @param main The title of plot.
#' @param xlab The X axis labels.
#' @param ylab The Y axis labels
#' @param cex.main The magnification to be used for title. (Default: 1.5)
#' @param cex.lab The magnification to be used for axis labels. (Default: 1.5)
#' @param cex.axis The magnification to be used for axis annotation. (Default: 1.5)
#' @param breaks The number of bars in the histogram.
#' @param ylim Y axis ranges.
#' @param xpos The horizontal position of the pvalue label. (Default: 0.03)
#' @param ypos The vertical position of the pvalue label. (Default: 0)
#' @param cex.text The magnification to be used for pvalue labels. (Default: 1.2)
#'
#' @return Histogram and p-value of Shapiro-Wilk Normality Test.
#' @export
#' @importFrom graphics hist lines rug text
#' @importFrom stats dnorm shapiro.test
#' @examples
#' data("wheatds")
#' inlier <- outlier(wheatds, sample = "Line", loc = "Env", rep = "Rep", phe = "DS", mode = "blup")
#' stat_out <- stat(x = inlier, sample = "Sample", phe = "inlier")
#' histplot(x = stat_out$mean)
#' @author Peng Zhao <pengzhao@nwafu.edu.cn>

histplot <- function(x, color = "#99d6e1", rug_color = "#f79999", freq = FALSE, lwd = 2, rug_lwd = 1, main = "", xlab = "", ylab = "", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, breaks="Sturges", ylim = NULL, xpos = 0.03, ypos = 0, cex.text = 1.2){
  x <- na.omit(x)
  normal_test <- shapiro.test(x)
  dp <- options(digits=3)
  on.exit(options(dp))
  pvalue <- format(normal_test$p.value, scientific=TRUE)
  p_plot <- paste("p-value = ", pvalue, sep = "\n")
  h <- hist(x, col = color, freq = freq, lwd = lwd, main = main, xlab = xlab, ylab = ylab, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, breaks = breaks, ylim = ylim)
  rug(jitter(x), side=1, col=rug_color, lwd=rug_lwd)
  xfit <- seq(min(x), max(x), length=40)
  meannum <- mean(x, na.rm=TRUE)
  sdnum <- sd(x, na.rm=TRUE)
  yfit <- dnorm(xfit, mean = meannum, sd=sdnum)
  yfit <- yfit*diff(h$mids[1:2]*length(x))
  lines(xfit, yfit, col="black", lwd=2)
  text((max(h$breaks) - xpos*max(h$breaks)), (max(h$counts) + ypos * max(h$counts)), p_plot, cex = cex.text)
}
