#' Performs the Best Linear Unbiased Prediction (BLUP)
#'
#' @title blup
#' @param x Input phenotype data file.
#' @param sample The column name of the sample name in phenotypic data. (Default: NULL)
#' @param year The column name of the year in phenotypic data. (Default: NULL)
#' @param loc The column name of the location in phenotypic data. (Default: NULL)
#' @param rep The column name of the replication in phenotypic data. (Default: NULL)
#' @param phe The column name of the phenotypic value in data. (Default: NULL)
#' @param fold Fold before inter-quartile range. (Default: 1.5)
#'
#' @return Estimate BLUPs for a phenotypic data with outliers removed on a per sample basis.
#' @export
#' @importFrom lme4 lmer
#' @importFrom lme4 ranef
#' @importFrom tidyr separate
#'
#' @examples
#' data("wheatds")
#' blup_out <- blup(wheatds, sample = "Line", loc = "Env", rep = "Rep", phe = "DS")
#' @author Peng Zhao <pengzhao@nwafu.edu.cn>

blup <- function(x, sample = NULL, year = NULL, loc = NULL, rep = NULL, phe = NULL, fold = 1.5) {
	## To prevent an error
	no_error <- options(lmerControl=list(check.nobs.vs.rankZ = "warning",check.nobs.vs.nlev = "warning",check.nobs.vs.nRE = "warning",check.nlev.gtreq.5 = "warning",check.nlev.gtr.1 = "warning"))
  on.exit(options(no_error))
  Sample <- NULL
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
    ## BLUP
    out_split <- separate(data=output,col=Sample,into=c("Sample","Year","Loc"),sep="&")
    SAMPLE <- as.factor(out_split$Sample)
    YEAR <- as.factor(out_split$Year)
    LOC <- as.factor(out_split$Loc)
    REP <- as.factor(out_split$Rep)
    TRAIT <- as.numeric(out_split$inlier)
    model <- lmer(TRAIT ~ (1|SAMPLE) + (1|LOC) + (1|YEAR) + (1|REP%in%LOC:YEAR) + (1|SAMPLE:LOC) + (1|SAMPLE:YEAR))
    traitsblup <- ranef(model)
    traitslineblup <- traitsblup$SAMPLE
    ## BLUP result summary
    sample <- row.names(traitsblup$SAMPLE)
    blup_value <- traitsblup$SAMPLE$`(Intercept)`
    blup_out <- data.frame(sample,blup_value)
    sample_phe <- data.frame(SAMPLE,TRAIT)
    unique_sample <- unique(SAMPLE)
    mean_out <- data.frame()
    for (i in seq_along(unique_sample)) {
      phe_value <- (subset(sample_phe,SAMPLE==unique_sample[i]))$TRAIT
      mean_value <- mean(phe_value)
      output_mean <- data.frame(sample=unique_sample[i],mean=mean_value)
      mean_out <- rbind(mean_out,output_mean)
    }
    blup_mean <- merge(x = blup_out, y = mean_out, by = "sample", all = TRUE)
    mean_adj <- blup_mean$blup_value + blup_mean$mean
    blup_output <- data.frame(Sample=blup_mean$sample,Mean=blup_mean$mean,Blup=blup_mean$blup_value,Adj_mean=mean_adj)
    ## Calculate heritability
    sum_model <- summary(model)
    var_temp <- data.frame(sum_model$varcor)
    Vg <- var_temp$vcov[3]
    Vge <- var_temp$vcov[2]
    Vgy <- var_temp$vcov[1]
    Ve <- var_temp$vcov[7]
    Year_num <- sum_model[["ngrps"]][["YEAR"]]
    Loc_num <- sum_model[["ngrps"]][["LOC"]]
    Rep_num <- length(summary(REP))
    h <- Vg / (Vg + Vge/Loc_num + Vgy/Year_num + Ve/(Year_num*Rep_num*Loc_num))
    h_output <- paste("Heritability = ", h, sep = "")
    cat(h_output,"\n")
    return(blup_output)
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
	    ## BLUP
	    out_split <- separate(data=output,col=Sample,into=c("Sample","Loc"),sep="&")
	    SAMPLE <- as.factor(out_split$Sample)
	    LOC <- as.factor(out_split$Loc)
	    REP <- as.factor(out_split$Rep)
	    TRAIT <- as.numeric(out_split$inlier)
	    model <- lmer(TRAIT ~ (1|SAMPLE) + (1|LOC) + (1|REP%in%LOC) + (1|SAMPLE:LOC))
	    traitsblup <- ranef(model)
	    traitslineblup <- traitsblup$SAMPLE
	    ## BLUP result summary
	    sample <- row.names(traitsblup$SAMPLE)
	    blup_value <- traitsblup$SAMPLE$`(Intercept)`
	    blup_out <- data.frame(sample,blup_value)
	    sample_phe <- data.frame(SAMPLE,TRAIT)
	    unique_sample <- unique(SAMPLE)
	    mean_out <- data.frame()
	    for (i in seq_along(unique_sample)) {
	      phe_value <- (subset(sample_phe,SAMPLE==unique_sample[i]))$TRAIT
	      mean_value <- mean(phe_value)
	      output_mean <- data.frame(sample=unique_sample[i],mean=mean_value)
	      mean_out <- rbind(mean_out,output_mean)
	    }
	    blup_mean <- merge(x = blup_out, y = mean_out, by = "sample", all = TRUE)
	    mean_adj <- blup_mean$blup_value + blup_mean$mean
	    blup_output <- data.frame(Sample=blup_mean$sample,Mean=blup_mean$mean,Blup=blup_mean$blup_value,Adj_mean=mean_adj)
	    ## Calculate heritability
	    sum_model <- summary(model)
	    var_temp <- data.frame(sum_model$varcor)
	    Vg <- var_temp$vcov[2]
	    Vge <- var_temp$vcov[1]
	    Ve <- var_temp$vcov[5]
	    Loc_num <- sum_model[["ngrps"]][["LOC"]]
	    Rep_num <- length(summary(REP))
	    h <- Vg / (Vg + Vge/Loc_num + Ve/(Rep_num*Loc_num))
	    h_output <- paste("Heritability = ", h, sep = "")
	    cat(h_output,"\n")
	    return(blup_output)
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
	    ## BLUP
	    out_split <- separate(data=output,col=Sample,into=c("Sample","Year"),sep="&")
	    SAMPLE <- as.factor(out_split$Sample)
	    YEAR <- as.factor(out_split$Year)
	    REP <- as.factor(out_split$Rep)
	    TRAIT <- as.numeric(out_split$inlier)
	    model <- lmer(TRAIT ~ (1|SAMPLE) + (1|YEAR) + (1|REP%in%YEAR) + (1|SAMPLE:YEAR))
	    traitsblup <- ranef(model)
	    traitslineblup <- traitsblup$SAMPLE
	    ## BLUP result summary
	    sample <- row.names(traitsblup$SAMPLE)
	    blup_value <- traitsblup$SAMPLE$`(Intercept)`
	    blup_out <- data.frame(sample,blup_value)
	    sample_phe <- data.frame(SAMPLE,TRAIT)
	    unique_sample <- unique(SAMPLE)
	    mean_out <- data.frame()
	    for (i in seq_along(unique_sample)) {
	      phe_value <- (subset(sample_phe,SAMPLE==unique_sample[i]))$TRAIT
	      mean_value <- mean(phe_value)
	      output_mean <- data.frame(sample=unique_sample[i],mean=mean_value)
	      mean_out <- rbind(mean_out,output_mean)
	    }
	    blup_mean <- merge(x = blup_out, y = mean_out, by = "sample", all = TRUE)
	    mean_adj <- blup_mean$blup_value + blup_mean$mean
	    blup_output <- data.frame(Sample=blup_mean$sample,Mean=blup_mean$mean,Blup=blup_mean$blup_value,Adj_mean=mean_adj)
	    ## Calculate heritability
	    sum_model <- summary(model)
	    var_temp <- data.frame(sum_model$varcor)
	    Vg <- var_temp$vcov[2]
	    Vge <- var_temp$vcov[1]
	    Ve <- var_temp$vcov[5]
	    Year_num <- sum_model[["ngrps"]][["YEAR"]]
	    Rep_num <- length(summary(REP))
	    h <- Vg / (Vg + Vge/Year_num + Ve/(Rep_num*Year_num))
	    h_output <- paste("Heritability = ", h, sep = "")
	    cat(h_output,"\n")
	    return(blup_output)
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
	    ## BLUP
	    out_split <- separate(data=output,col=Sample,into=c("Sample","Loc"),sep="&")
	    SAMPLE <- as.factor(out_split$Sample)
	    LOC <- as.factor(out_split$Loc)
	    TRAIT <- as.numeric(out_split$inlier)
	    model <- lmer(TRAIT ~ (1|SAMPLE) + (1|LOC))
	    traitsblup <- ranef(model)
	    traitslineblup <- traitsblup$SAMPLE
	    ## BLUP result summary
	    sample <- row.names(traitsblup$SAMPLE)
	    blup_value <- traitsblup$SAMPLE$`(Intercept)`
	    blup_out <- data.frame(sample,blup_value)
	    sample_phe <- data.frame(SAMPLE,TRAIT)
	    unique_sample <- unique(SAMPLE)
	    mean_out <- data.frame()
	    for (i in seq_along(unique_sample)) {
	      phe_value <- (subset(sample_phe,SAMPLE==unique_sample[i]))$TRAIT
	      mean_value <- mean(phe_value)
	      output_mean <- data.frame(sample=unique_sample[i],mean=mean_value)
	      mean_out <- rbind(mean_out,output_mean)
	    }
	    blup_mean <- merge(x = blup_out, y = mean_out, by = "sample", all = TRUE)
	    mean_adj <- blup_mean$blup_value + blup_mean$mean
	    blup_output <- data.frame(Sample=blup_mean$sample,Mean=blup_mean$mean,Blup=blup_mean$blup_value,Adj_mean=mean_adj)
	    ## Calculate heritability
	    sum_model <- summary(model)
	    var_temp <- data.frame(sum_model$varcor)
	    Vg <- var_temp$vcov[1]
	    Ve <- var_temp$vcov[3]
	    Loc_num <- sum_model[["ngrps"]][["LOC"]]
	    h <- Vg / (Vg + Ve/Loc_num)
	    h_output <- paste("Heritability = ", h, sep = "")
	    cat(h_output,"\n")
	    return(blup_output)
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
	    ## BLUP
	    out_split <- separate(data=output,col=Sample,into=c("Sample","Year"),sep="&")
	    SAMPLE <- as.factor(out_split$Sample)
	    YEAR <- as.factor(out_split$Year)
	    TRAIT <- as.numeric(out_split$inlier)
	    model <- lmer(TRAIT ~ (1|SAMPLE) + (1|YEAR))
	    traitsblup <- ranef(model)
	    traitslineblup <- traitsblup$SAMPLE
	    ## BLUP result summary
	    sample <- row.names(traitsblup$SAMPLE)
	    blup_value <- traitsblup$SAMPLE$`(Intercept)`
	    blup_out <- data.frame(sample,blup_value)
	    sample_phe <- data.frame(SAMPLE,TRAIT)
	    unique_sample <- unique(SAMPLE)
	    mean_out <- data.frame()
	    for (i in seq_along(unique_sample)) {
	      phe_value <- (subset(sample_phe,SAMPLE==unique_sample[i]))$TRAIT
	      mean_value <- mean(phe_value)
	      output_mean <- data.frame(sample=unique_sample[i],mean=mean_value)
	      mean_out <- rbind(mean_out,output_mean)
	    }
	    blup_mean <- merge(x = blup_out, y = mean_out, by = "sample", all = TRUE)
	    mean_adj <- blup_mean$blup_value + blup_mean$mean
	    blup_output <- data.frame(Sample=blup_mean$sample,Mean=blup_mean$mean,Blup=blup_mean$blup_value,Adj_mean=mean_adj)
	    ## Calculate heritability
	    sum_model <- summary(model)
	    var_temp <- data.frame(sum_model$varcor)
	    Vg <- var_temp$vcov[1]
	    Ve <- var_temp$vcov[3]
	    Year_num <- sum_model[["ngrps"]][["YEAR"]]
	    h <- Vg / (Vg + Ve/Year_num)
	    h_output <- paste("Heritability = ", h, sep = "")
	    cat(h_output,"\n")
	    return(blup_output)
	  }
  }
}
