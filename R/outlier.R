#' Remove outliers from phenotypic data
#'
#' @title outlier
#' @param x Input phenotype data file.
#' @param sample The column name of the sample name in phenotypic data. (Default: NULL)
#' @param year The column name of the year in phenotypic data. (Default: NULL)
#' @param loc The column name of the location in phenotypic data. (Default: NULL)
#' @param rep The column name of the replication in phenotypic data. (Default: NULL)
#' @param phe The column name of the phenotypic value in data. (Default: NULL)
#' @param fold Fold before inter-quartile range. (Default: 1.5)
#' @param mode Type of input phenotypic data. "normal" means normal data, "blup" means data containing year/location/repetition. (Default: "normal")
#'
#' @return phenotypic data with outliers removed.
#' @export
#'
#' @importFrom stats quantile
#' @importFrom tidyr separate
#' @examples
#' data("wheatds")
#' inlier <- outlier(wheatds, sample = "Line", loc = "Env", rep = "Rep", phe = "DS", mode = "blup")
#' @author Peng Zhao <pengzhao@nwafu.edu.cn>

outlier <- function(x, sample = NULL, year = NULL, loc = NULL, rep = NULL, phe = NULL, fold = 1.5, mode = "normal") {
  x.Sample <- NULL
  Sample <- NULL
  names(x)[names(x) == sample] <- "Sample"
  names(x)[names(x) == phe] <- "Phe"
  if(mode == "normal") {
	  phe <- x$Phe
	  input <- data.frame(x$Sample, phe)
	  unique_sample <- unique(input$x.Sample)
	  output <- data.frame()
	  for (i in seq_along(unique_sample)) {
	    phe_value <- (subset(input,x.Sample==unique_sample[i]))$phe
	    q <- quantile(phe_value, na.rm = T)
	    iqr <- q[4] - q[2]
	    min_value <- q[2] - fold*iqr
	    max_value <- q[4] + fold*iqr
	    inlier.pos <- which(phe_value >= min_value & phe_value <= max_value)
	    value <- phe_value[inlier.pos]
	    if(length(value) > 0) {
	      output_value <- data.frame(Sample=unique_sample[i], inlier=value)
	      output <- rbind(output, output_value)
	    }
	  }
		return(output)
	} else if (mode == "blup") {
		## Null number
	  null_num <- 0
	  no_null_list <- list(sample, year, loc, rep, phe)
	  for(i in 1:length(no_null_list)){
			if(is.null(no_null_list[[i]]) == FALSE){
				null_num = null_num + 1
		  }
		}
	  ## Year/Loc/Rep
	  if (null_num == 5) {
	  	names(x)[names(x) == sample] <- "Sample"
	  	names(x)[names(x) == year] <- "Year"
	  	names(x)[names(x) == loc] <- "Loc"
	  	names(x)[names(x) == rep] <- "Rep"
	  	names(x)[names(x) == phe] <- "Phe"
	    phe <- x$Phe
	    ## Merge sample
	    sample_name <- paste(x$Sample,x$Year,x$Loc,sep="&")
	    input <- data.frame(sample_name,phe)
	    unique_sample <- unique(sample_name)
	    output <- data.frame()
	    ## outlier filter
	    for (i in seq_along(unique_sample)) {
	      phe_value <- (subset(input,sample_name==unique_sample[i]))$phe
	      q <- quantile(phe_value,na.rm = T)
	      iqr <- q[4] - q[2]
	      min_value <- q[2] - fold*iqr
	      max_value <- q[4] + fold*iqr
	      inlier.pos <- which(phe_value >= min_value & phe_value <= max_value)
	      value <- phe_value[inlier.pos]
	      if(length(value) > 0) {
	      	rep_value <- c()
					for (r in 1:length(value)) {
						r_num <- paste("rep",r,sep = "")
						rep_value <- c(rep_value, r_num)
					}
	        output_value <- data.frame(Sample=unique_sample[i],Rep=rep_value,inlier=value)
	        output <- rbind(output,output_value)
	      }
	    }
	    out_split <- separate(data=output,col=Sample,into=c("Sample","Year","Loc"),sep="&")
	    return(out_split)
	  } else if (null_num == 4) {
	  	names(x)[names(x) == sample] <- "Sample"
	  	names(x)[names(x) == rep] <- "Rep"
	  	names(x)[names(x) == phe] <- "Phe"
	  	## Loc/Rep
	  	if (is.null(year) == TRUE) {
		  	names(x)[names(x) == loc] <- "Loc"
		    phe <- x$Phe
		    sample_name <- paste(x$Sample,x$Loc,sep="&")
		    input <- data.frame(sample_name,phe)
		    unique_sample <- unique(sample_name)
		    output <- data.frame()
		    ## outlier filter
		    for (i in seq_along(unique_sample)) {
		      phe_value <- (subset(input,sample_name==unique_sample[i]))$phe
		      q <- quantile(phe_value,na.rm = T)
		      iqr <- q[4] - q[2]
		      min_value <- q[2] - fold*iqr
		      max_value <- q[4] + fold*iqr
		      inlier.pos <- which(phe_value >= min_value & phe_value <= max_value)
		      value <- phe_value[inlier.pos]
		      if(length(value) > 0) {
		      	rep_value <- c()
						for (r in 1:length(value)) {
							r_num <- paste("rep",r,sep = "")
							rep_value <- c(rep_value, r_num)
						}
		        output_value <- data.frame(Sample=unique_sample[i],Rep=rep_value,inlier=value)
		        output <- rbind(output,output_value)
		      }
		    }
		    out_split <- separate(data=output,col=Sample,into=c("Sample","Loc"),sep="&")
		    return(out_split)
		  } else if (is.null(loc) == TRUE) {
		  	names(x)[names(x) == year] <- "Year"
		    phe <- x$Phe
		    sample_name <- paste(x$Sample,x$Year,sep="&")
		    input <- data.frame(sample_name,phe)
		    unique_sample <- unique(sample_name)
		    output <- data.frame()
		    ## outlier filter
		    for (i in seq_along(unique_sample)) {
		      phe_value <- (subset(input,sample_name==unique_sample[i]))$phe
		      q <- quantile(phe_value,na.rm = T)
		      iqr <- q[4] - q[2]
		      min_value <- q[2] - fold*iqr
		      max_value <- q[4] + fold*iqr
		      inlier.pos <- which(phe_value >= min_value & phe_value <= max_value)
		      value <- phe_value[inlier.pos]
		      if(length(value) > 0) {
						rep_value <- c()
						for (r in 1:length(value)) {
							r_num <- paste("rep",r,sep = "")
							rep_value <- c(rep_value, r_num)
						}
		        output_value <- data.frame(Sample=unique_sample[i],Rep=rep_value,inlier=value)
		        output <- rbind(output,output_value)
		      }
		    }
		    out_split <- separate(data=output,col=Sample,into=c("Sample","Year"),sep="&")
		    return(out_split)
		  }
	  } else if (null_num == 3) {
	  	names(x)[names(x) == sample] <- "Sample"
	  	names(x)[names(x) == phe] <- "Phe"
	  	## Loc
	  	if (is.null(year) == TRUE) {
			  names(x)[names(x) == loc] <- "Loc"
			  phe <- x$Phe
		    sample_name <- paste(x$Sample,x$Loc,sep="&")
		    input <- data.frame(sample_name,phe)
		    unique_sample <- unique(sample_name)
		    output <- data.frame()
		    ## outlier filter
		    for (i in seq_along(unique_sample)) {
		      phe_value <- (subset(input,sample_name==unique_sample[i]))$phe
		      q <- quantile(phe_value,na.rm = T)
		      iqr <- q[4] - q[2]
		      min_value <- q[2] - fold*iqr
		      max_value <- q[4] + fold*iqr
		      inlier.pos <- which(phe_value >= min_value & phe_value <= max_value)
		      value <- phe_value[inlier.pos]
		      if(length(value) > 0) {
		        output_value <- data.frame(Sample=unique_sample[i],inlier=value)
		        output <- rbind(output,output_value)
		      }
		    }
		    out_split <- separate(data=output,col=Sample,into=c("Sample","Loc"),sep="&")
		    return(out_split)
		  } else if (is.null(loc) == TRUE) {
		  	names(x)[names(x) == year] <- "Year"
			  phe <- x$Phe
		    sample_name <- paste(x$Sample,x$Year,sep="&")
		    input <- data.frame(sample_name,phe)
		    unique_sample <- unique(sample_name)
		    output <- data.frame()
		    ## outlier filter
		    for (i in seq_along(unique_sample)) {
		      phe_value <- (subset(input,sample_name==unique_sample[i]))$phe
		      q <- quantile(phe_value,na.rm = T)
		      iqr <- q[4] - q[2]
		      min_value <- q[2] - fold*iqr
		      max_value <- q[4] + fold*iqr
		      inlier.pos <- which(phe_value >= min_value & phe_value <= max_value)
		      value <- phe_value[inlier.pos]
		      if(length(value) > 0) {
		        output_value <- data.frame(Sample=unique_sample[i],inlier=value)
		        output <- rbind(output,output_value)
		      }
		    }
		    out_split <- separate(data=output,col=Sample,into=c("Sample","Year"),sep="&")
		    return(out_split)
		  }
	  }
	}
}
