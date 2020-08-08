#' Calculate statistical indicators of phenotypic data
#'
#' @title stat
#' @param x Input phenotype data file.
#' @param sample The column name of the sample name in phenotypic data. (Default: NULL)
#' @param phe The column name of the phenotypic value in data. (Default: NULL)
#'
#' @return Mean, median, standard deviation, standard error of phenotypic data for each sample.
#' @export
#'
#' @importFrom stats median na.omit quantile sd
#' @examples
#' data("wheatds")
#' inlier <- outlier(wheatds, sample = "Line", loc = "Env", rep = "Rep", phe = "DS", mode = "blup")
#' stat_out <- stat(x = inlier, sample = "Sample", phe = "inlier")
#' @author Peng Zhao <pengzhao@nwafu.edu.cn>
stat <- function(x, sample = NULL, phe = NULL) {
	x.Sample <- NULL
	names(x)[names(x) == sample] <- "Sample"
  names(x)[names(x) == phe] <- "Phe"
  phe <- x$Phe
  input <- data.frame(x$Sample, phe)
  unique_sample <- unique(input$x.Sample)
  output <- data.frame()
  for (i in seq_along(unique_sample)) {
    value <- (subset(input,x.Sample==unique_sample[i]))$phe
    mean_value <- mean(value)
    median_value <- median(value)
    sd_value <- sd(value)
    se_value <- sd(value) / sqrt(length(value))
    n_value <- length(value)
    output_value <- data.frame(Sample = unique_sample[i], mean = mean_value, median = median_value, sd = sd_value, se = se_value, n = n_value)
    output <- rbind(output, output_value)
  }
  return(output)
}
