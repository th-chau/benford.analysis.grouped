#### Package benford.analysis.grouped ####

##' Stremlining Benford Analysis for data validation and forensic analytics by group.
##'
##' Modified to enable group analysis functionalities; original functions are mostly kept untouched, with slight modifications in efficiency of the code and to reduce dependencies.
##' Goal: new functions that enable efficient analysis and forewarning of irregularities that are clustered in certain groups
##' 
##' @references Alexander, J. (2009). Remarks on the use of Benford's Law. Working Paper, Case Western Reserve University, Department of Mathematics and Cognitive Science.
##' \cr\cr Berger, A. and Hill, T. (2011). A basic theory of Benford's Law. Probability Surveys, 8, 1-126.
##' \cr\cr Hill, T. (1995). A statistical derivation of the significant-digit law. Statistical Science, 10(4), 354-363.
##' \cr\cr Nigrini, M. J. (2012). Benford's Law: Application for Forensic Accounting, Auditing and Fraud Detection. Wiley and Sons: New Jersey.
##' \cr\cr Nigrini, M. J. (2011). Forensic Analyticis: Methods and Techniques for Forensic Accounting Investigations.Wiley and Sons: New Jersey.
##' \cr\cr 
##' 
##' @docType package
##' @aliases benford.analysis.grouped-package
##' @name benford.analysis.grouped
##' 
NULL

##' @title Benford Analysis of a dataset
##' @description This function validates a dataset using Benford's Law.
##' Its main purposes are to find out where the dataset deviates from Benford's Law and 
##' to identify suspicious data that need further verification. 
##' 
##' For a more complete example, see the package help at \link{benford.analysis}.
##' @usage
##' benford(data, number.of.digits = 2, sign = "positive", discrete=TRUE, round=3)
##' @param data a numeric vector.
##' @param number.of.digits how many first digits to analyse .
##' @param sign  The default value for sign is "positive" and it analyzes only data greater than zero. 
##' There are also the options "negative" and "both" that will analyze only negative values or both positive and negative values of the data,
##' respectively. For large datasets with both positive and negative numbers, 
##' it is usually recommended to perform a separate analysis for each group,
##' for the incentives to manipulate the numbers are usually different.
##' @param discrete most real data - like population numbers or accounting data - are discrete, so 
##' the default is TRUE. This paramater sets rounding to the differences of the ordered data to avoid floating point number
##' errors in the second order distribution, that usually occurs when data is discrete
##' and the ordered numbers are very close to each other. If your data is continuous
##' (like a simulated lognormal) you should run with discrete = FALSE. 
##' @param round it defines the number of digits that the rounding will use if discrete = TRUE.
##' @return An object of class Benford containing the results of the analysis. It is a list of 
##' eight objects, namely:
##' 
##' \item{info}{general information, including \itemize{
##' \item data.name: the name of the data used.
##' \item n: the number of observations used.
##' \item n.second.order: the number of observations used for second order analysis.
##' \item number.of.digits: the number of first digits analysed.
##' }}
##' 
##' \item{data}{a data frame with: \itemize{
##' \item lines.used: the original lines of the dataset.
##' \item data.used: the data used.
##' \item data.mantissa: the log data's mantissa.
##' \item data.digits: the first digits of the data.
##' }}
##' 
##' \item{s.o.data}{a data frame with: \itemize{
##' \item data.second.order: the differences of the ordered data.
##' \item data.second.order.digits: the first digits of the second order analysis.
##' }}
##' 
##' \item{bfd}{a data frame with: \itemize{
##' \item digits: the groups of digits analysed.
##' \item data.dist: the distribution of the first digits of the data.
##' \item data.second.order.dist: the distribution of the first digits of the second order analysis.
##' \item benford.dist: the theoretical benford distribution.
##' \item data.second.order.dist.freq: the frequency distribution of the first digits of the second order analysis.
##' \item data.dist.freq: the frequency distribution of the first digits of the data.
##' \item benford.dist.freq: the theoretical benford frequency distribution.
##' \item benford.so.dist.freq: the theoretical benford frequency distribution of the second order analysis.
##' \item data.summation: the summation of the data values grouped by first digits.
##' \item abs.excess.summation: the absolute excess summation of the data values grouped by first digits.
##' \item difference: the difference between the data and benford frequencies.
##' \item squared.diff: the chi-squared difference between data and benford frequencies.
##' \item absolute.diff: the absolute difference between data and benford frequencies.
##' \item ks
##' }}
##' 
##' \item{mantissa}{a data frame with: \itemize{
##' \item mean.mantissa: the mean of the mantissa.
##' \item var.mantissa: the variance of the mantissa.
##' \item ek.mantissa: the excess kurtosis of the mantissa.
##' \item sk.mantissa: the skewness of the mantissa.
##' }}
##' 
##' \item{MAD}{the mean absolute deviation.}
##' \item{distortion.factor}{the distortion factor} 
##' 
##' \item{stats}{list of "htest" class statistics: \itemize{
##' \item chisq: Pearson's Chi-squared test.
##' \item mantissa.arc.test: Mantissa Arc Test.
##' }}
##' @examples 
##' data(corporate.payment) #loads data
##' bfd.cp <- benford(corporate.payment$Amount) #generates benford object
##' bfd.cp #prints
##' plot(bfd.cp) #plots
##' 
##' @export

benford <- function(data, number.of.digits = 2, sign = "positive", discrete=TRUE, round=3){
  
  data.name <- as.character(deparse(substitute(data)))
  
  benford.digits <- generate.benford.digits(number.of.digits)
    
  benford.dist <- generate.benford.distribution(benford.digits)
  
  empirical.distribution <- generate.empirical.distribution(data, number.of.digits,sign, second.order = FALSE, benford.digits)
  
  n <- length(empirical.distribution$data)
  
  second.order <- generate.empirical.distribution(data, number.of.digits,sign, second.order = TRUE, benford.digits, discrete = discrete, round = round)
  
  n.second.order <- length(second.order$data)
  
  benford.dist.freq <- benford.dist*n
  
  ## calculating useful summaries and differences
  difference <- empirical.distribution$dist.freq - benford.dist.freq
  
  squared.diff <- ((empirical.distribution$dist.freq - benford.dist.freq)^2)/benford.dist.freq
  
  absolute.diff <- abs(empirical.distribution$dist.freq - benford.dist.freq)
  
  ### chi-squared test
  chisq.bfd <- chisq.test.bfd(squared.diff, data.name)
  
  ### MAD
  mean.abs.dev <- sum(abs(empirical.distribution$dist - benford.dist)/(length(benford.dist)))
  

  ### Summation
  summation <- generate.summation(benford.digits,empirical.distribution$data, empirical.distribution$data.digits)
  abs.excess.summation <- abs(summation - mean(summation))
  
  ### Mantissa
  mantissa <- extract.mantissa(empirical.distribution$data)
  mean.mantissa <- mean(mantissa)
  var.mantissa <- var(mantissa)
  ek.mantissa <- excess.kurtosis(mantissa)
  sk.mantissa <- skewness(mantissa)
  
  ### Mantissa Arc Test
  mat.bfd <- mantissa.arc.test(mantissa, data.name)
  
  ### Distortion Factor
  distortion.factor <- DF(empirical.distribution$data)  
  
  ## recovering the lines of the numbers
  if (sign == "positive") lines <- which(data > 0 & !is.na(data))
  if (sign == "negative") lines <- which(data < 0 & !is.na(data))
  if (sign == "both")     lines <- which(data != 0 & !is.na(data))
  #lines <- which(data %in% empirical.distribution$data)
  
  ## output
  output <- list(info = list(data.name = data.name,
                             n = n,
                             n.second.order = n.second.order,
                             number.of.digits = number.of.digits),
                 
                 data = data.table(lines.used = lines,
                                   data.used = empirical.distribution$data,
                                   data.mantissa = mantissa,
                                   data.digits = empirical.distribution$data.digits),
                 
                 s.o.data = data.table(second.order = second.order$data,
                                       data.second.order.digits = second.order$data.digits),
                 
                 bfd = data.table(digits = benford.digits,
                                  data.dist = empirical.distribution$dist,
                                  data.second.order.dist = second.order$dist,
                                  benford.dist = benford.dist,
                                  data.second.order.dist.freq = second.order$dist.freq,
                                  data.dist.freq = empirical.distribution$dist.freq,
                                  benford.dist.freq = benford.dist.freq,
                                  benford.so.dist.freq = benford.dist*n.second.order,
                                  data.summation = summation,
                                  abs.excess.summation = abs.excess.summation,
                                  difference = difference,
                                  squared.diff = squared.diff,
                                  absolute.diff = absolute.diff),
                 
                 mantissa = data.table(statistic = c("Mean Mantissa", 
                                                     "Var Mantissa", 
                                                     "Ex. Kurtosis Mantissa",
                                                     "Skewness Mantissa"),
                                       values = c(mean.mantissa = mean.mantissa,
                                                  var.mantissa = var.mantissa,
                                                  ek.mantissa = ek.mantissa,
                                                  sk.mantissa = sk.mantissa)),
                 MAD = mean.abs.dev,
                 
                 distortion.factor = distortion.factor,
                 
                 stats = list(chisq = chisq.bfd,
                              mantissa.arc.test = mat.bfd)
  )
  
  class(output) <- "Benford"
  
  return(output)
  
}

##' @title Plot method for Benford Analysis
##' @description The \code{plot} method for "Benford" objects.
##' @usage 
##' 
##' \method{plot}{Benford}(x,except=c("chi square","abs diff"), multiple=TRUE, ...) 
##' @param  x a "Benford" object
##' @param except it specifies which plots are not going to be plotted.
##' Currently, you can choose from 7 plots: "digits", "second order", "summation",
##' "mantissa", "chi square", "abs diff", "ex summation". If you want to plot all, just
##' put except = "none". The default is not to plot the "mantissa" and "abs diff".
##' @param multiple if TRUE, all plots are grouped in the same window.
##' @param ... arguments to be passed to generic plot functions,
##' @return Plots the Benford object.
##' @export
##' @importFrom graphics abline axis barplot legend lines par plot
##' @importFrom stats pchisq var
##' @importFrom utils head
##' @importFrom stats setNames
plot.Benford<- function(x, except=c("chi square", "abs diff"), multiple=TRUE ,...){
  
  old.par <- par(no.readonly=TRUE)
  
  if(class(x)!="Benford") stop("Class(x) must be 'Benford'")
  
  if(!any(except %in% c("digits", "second order", "summation",
                     "mantissa", "chi squared", "abs diff","none"))) {stop("Invalid except name. Type ?plot.Benford for help.")}
  
  if(multiple){
    
    nGraphics <- 8 - length(except)
    
    if(nGraphics<4){rows=1; cols=nGraphics}
    if(nGraphics >=4 & nGraphics<=6){rows=2; cols=3}
    if(nGraphics>6){rows=2; cols=4}
    
    par(mfrow=c(rows,cols))
  }
  
  if(all(except!="digits")){
  plotting.data.vs.benford(x, ...)
  }
  
  if(all(except!="second order")){
  plotting.second.order(x, ...)
  }
  
  if(all(except!="summation")){
  plotting.summation(x, ...)
  }
  
  if(all(except!="mantissa")){
  plotting.ordered.mantissa(x, ...)
  }
  
  if(all(except!="chi squared")){
  plotting.chi_squared(x, ...)
  }
  
  if(all(except!="abs diff")){
  plotting.abs.diff(x, ...)
  }
  
  if(all(except!="ex summation")){
  plotting.ex.summation(x, ...)
  }
  plotting.legend(x)
  
  par(old.par)
  
}

##' @title Print method for Benford Analysis
##' @description The \code{print} method for "Benford" objects.
##' @usage 
##' 
##' \method{print}{Benford}(x, how.many=5, ...)
##' @param  x a "Benford" object.
##' @param how.many a number that defines how many of the biggest absolute differences to show.
##' @param ... arguments to be passed to generic print functions.
##' @return Prints the Benford object.
##' @export
print.Benford <- function(x,how.many=5,...){
  
  if(class(x)!="Benford") stop("Class(x) must be 'Benford'")
  cat("\nBenford object:\n",
      "\nData:", x[["info"]]$data.name,
      "\nNumber of observations used =", x[["info"]]$n,
      "\nNumber of obs. for second order =", x[["info"]]$n.second.order,
      "\nFirst digits analysed =", x[["info"]]$number.of.digits)    
  
  cat("\n\nMantissa: \n")
  cat("\n")
  pretty_print <- x$mantissa
  pretty_print$statistic <- gsub("Mantissa|\\s", "", pretty_print$statistic)
  pretty_print <- setNames(pretty_print, c("Statistic", "Value"))
  print.data.frame(pretty_print, row.names = FALSE, digits = 2)
  cat("\n")
  
  cat("\nThe", how.many, "largest deviations: \n")
  cat("\n")
  print.data.frame(round(head(x[["bfd"]][order(absolute.diff, decreasing=TRUE)][,list(digits, absolute.diff)], how.many), 2))
  cat("\n")
  
  cat("Stats:\n")
  print(x[["stats"]]$chisq)
  print(x[["stats"]]$mantissa.arc.test)
  cat("Mean Absolute Deviation:",x[["MAD"]])
  cat("\nDistortion Factor:", x[["distortion.factor"]])
  cat("\n\nRemember: Real data will never conform perfectly to Benford's Law. You should not focus on p-values!")
  values <- NULL
  absolute.diff <- NULL
  digits <- NULL
}

##' @title Enhanced Benford function by group
##' @description This runs \code{benford} by group, and allows selective reporting of results. This returns a list of Benford objects by group
##' For a more complete example, see the package help at \link{benford.analysis}.
##'
##' @usage
##' benford.by.group(data, groupnamevect, number.of.digits = 2, sign = "positive", discrete=TRUE, round=3, minlength = 100)
##' @param data a numeric vector.
##' @param groupnamevect a character vector that allows the numeric vector to be grouped into different groups
##' @param number.of.digits how many first digits to analyse .
##' @param sign  The default value for sign is "positive" and it analyzes only data greater than zero. 
##' There are also the options "negative" and "both" that will analyze only negative values or both positive and negative values of the data,
##' respectively. For large datasets with both positive and negative numbers, 
##' it is usually recommended to perform a separate analysis for each group,
##' for the incentives to manipulate the numbers are usually different.
##' @param discrete most real data - like population numbers or accounting data - are discrete, so 
##' the default is TRUE. This paramater sets rounding to the differences of the ordered data to avoid floating point number
##' errors in the second order distribution, that usually occurs when data is discrete
##' and the ordered numbers are very close to each other. If your data is continuous
##' (like a simulated lognormal) you should run with discrete = FALSE. 
##' @param round it defines the number of digits that the rounding will use if discrete = TRUE.
##' @param minlength it defines the minimum elements within a group in order for a benford object to be returned
##' @return A list containing objects of class Benford containing the results of the analysis. It is a list of 
##' Benford objects grouped by unique values of groupnamevect, and in each of the Benford object within the list it contains:
##' 
##' \item{info}{general information, including \itemize{
##' \item data.name: the name of the data used.
##' \item n: the number of observations used.
##' \item n.second.order: the number of observations used for second order analysis.
##' \item number.of.digits: the number of first digits analysed.
##' }}
##' 
##' \item{data}{a data frame with: \itemize{
##' \item lines.used: the original lines of the dataset.
##' \item data.used: the data used.
##' \item data.mantissa: the log data's mantissa.
##' \item data.digits: the first digits of the data.
##' }}
##' 
##' \item{s.o.data}{a data frame with: \itemize{
##' \item data.second.order: the differences of the ordered data.
##' \item data.second.order.digits: the first digits of the second order analysis.
##' }}
##' 
##' \item{bfd}{a data frame with: \itemize{
##' \item digits: the groups of digits analysed.
##' \item data.dist: the distribution of the first digits of the data.
##' \item data.second.order.dist: the distribution of the first digits of the second order analysis.
##' \item benford.dist: the theoretical benford distribution.
##' \item data.second.order.dist.freq: the frequency distribution of the first digits of the second order analysis.
##' \item data.dist.freq: the frequency distribution of the first digits of the data.
##' \item benford.dist.freq: the theoretical benford frequency distribution.
##' \item benford.so.dist.freq: the theoretical benford frequency distribution of the second order analysis.
##' \item data.summation: the summation of the data values grouped by first digits.
##' \item abs.excess.summation: the absolute excess summation of the data values grouped by first digits.
##' \item difference: the difference between the data and benford frequencies.
##' \item squared.diff: the chi-squared difference between data and benford frequencies.
##' \item absolute.diff: the absolute difference between data and benford frequencies.
##' \item ks
##' }}
##' 
##' \item{mantissa}{a data frame with: \itemize{
##' \item mean.mantissa: the mean of the mantissa.
##' \item var.mantissa: the variance of the mantissa.
##' \item ek.mantissa: the excess kurtosis of the mantissa.
##' \item sk.mantissa: the skewness of the mantissa.
##' }}
##' 
##' \item{MAD}{the mean absolute deviation.}
##' \item{distortion.factor}{the distortion factor} 
##' 
##' \item{stats}{list of "htest" class statistics: \itemize{
##' \item chisq: Pearson's Chi-squared test.
##' \item mantissa.arc.test: Mantissa Arc Test.
##' }}
##' @examples 
##' data(corporate.payment) #loads data
##' bfd.cp <- benford.by.group(corporate.payment$Amount, corporate.payment$Date) #generates list of benford objects
##' bfd.cp[1] #prints
##' plot(bfd.cp[1]) #plots
##' 
##' @export
benford.by.group <- function(data, groupnamevect, number.of.digits = 2, sign = "positive", discrete=TRUE, round=3, minlength = 100){
  
  output <- list()

  # Select only non-NA values
  groupnamevect <- groupnamevect[!is.na(data)]
  data <- data[!is.na(data)]

  for(i in 1:length(unique(groupnamevect))){
      ElementsToSelect <- groupnamevect == unique(groupnamevect)[[i]]
      DataSubSet <- data[ElementsToSelect]
      DataSubSet <- DataSubSet[!is.na(DataSubSet) & !DataSubSet == 0]
        if(!is.na(unique(groupnamevect)[[i]]) & length(DataSubSet) > minlength){
          names(DataSubSet) <- paste0(unique(groupnamevect)[[i]])
          temp <- benford(DataSubSet, number.of.digits = 2, sign = "positive", discrete=TRUE, round=3)
          output[[unique(groupnamevect)[[i]]]] <- temp
    }
  }

  return(output)

}


##' @title Subset Benford objects based on anomalous results
##' @description This takes a list of \code{benford} objects and returns a list of \code{benford} objects that pass the specified criteria.
##'
##' @usage
##' benford.select.group(ListOfBenfordObjects, AreaThreshold = 0.1, MantissaThreshold = 0.05, LogicalRelationship = "AND")
##' @param ListOfBenfordObjects a list of benford objects.
##' @param AreaThreshold this sets the *minimum* value for an object to be flagged.
##' This is the summed deviation of the discrete empirical distribution curve from the discrete Benford curve; in other words, this is the total area in whcih the two curves are not in agreement. 
##' @param This sets the *maximum* p-value of the Mantissa arc test for an object to be flagged.
##' @param LogicalRelationship this takes a string - "AND", "OR" - to specify whether an object should be flagged when it satisfies one or both criteria
##' @return A list containing objects of class Benford that satisfy one or both of the criteria.
##' 
##' @examples 
##' Sublist <- benford.select.group(ListOfBenfordObjects, AreaThreshold = 0.1, MantissaThreshold = 0.05, LogicalRelationship = "AND")
##' plot(ListOfBenfordObjects[[1]])
##' 
##' @export

benford.select.group <- function(ListOfBenfordObjects, AreaThreshold = 0.1, MantissaThreshold = 0.05, LogicalRelationship = "AND"){

if(!(LogicalRelationship == "AND" | LogicalRelationship == "and" | LogicalRelationship == "or" | LogicalRelationship == "OR")){
  stop("Logical relationship between the two conditions must be a character vector 'AND' / 'OR'")
}

SubsetOutput <- list()

  for(i in 1:length(names(ListOfBenfordObjects))){

    AreaDeviation <- sum(abs(ListOfBenfordObjects[[i]][["bfd"]][["data.dist"]] - ListOfBenfordObjects[[i]][["bfd"]][["benford.dist"]]))
    Mantissa.p.Value <- ListOfBenfordObjects[[i]][["stats"]][["mantissa.arc.test"]][["p.value"]]

    if(LogicalRelationship == "AND" | LogicalRelationship == "and"){
      if(AreaDeviation > AreaThreshold & Mantissa.p.Value < MantissaThreshold){
            SubsetOutput[[names(ListOfBenfordObjects)[[i]]]] <- ListOfBenfordObjects[[i]]
      }
    }

    if(LogicalRelationship == "Or" |LogicalRelationship == "or"){
      if(AreaDeviation > AreaThreshold | Mantissa.p.Value < MantissaThreshold){
            SubsetOutput[[names(ListOfBenfordObjects)[[i]]]] <- ListOfBenfordObjects[[i]]
      }
    }

  }

return(SubsetOutput)

}


##' @title Summary table of a list of Benford objects
##' @description This takes a list of \code{benford} objects and returns a summary table of the list, containing the following for each object: name, observations, number of digits used, Mantissa statistic, p-value of MantissaStat arc test and the total deviation of the empirical distribution from the Benford distribution.
##'
##' @usage
##' benford.sum.table(ListOfBenfordObjects)
##' @param ListOfBenfordObjects a list of benford objects.
##' 
##' @export

benford.sum.table <- function(ListOfBenfordObjects){

    for(i in 1:length(names(ListOfBenfordObjects))){

    AreaDeviation <- sum(abs(ListOfBenfordObjects[[i]][["bfd"]][["data.dist"]] - ListOfBenfordObjects[[i]][["bfd"]][["benford.dist"]]))
    Mantissa.p.Value <- ListOfBenfordObjects[[i]][["stats"]][["mantissa.arc.test"]][["p.value"]]

    if(i == 1){

      OutputDF <- data.frame(Name = names(ListOfBenfordObjects)[[i]], Observations = ListOfBenfordObjects[[i]][["info"]][["n"]], NoOfDigits = ListOfBenfordObjects[[i]][["info"]][["number.of.digits"]], MantissaStat = ListOfBenfordObjects[[i]][["stats"]][["mantissa.arc.test"]][["statistic"]], MantissaPValue = Mantissa.p.Value, AreaDeviation = AreaDeviation)

    }else{

    Temp <- data.frame(Name = names(ListOfBenfordObjects)[[i]], Observations = ListOfBenfordObjects[[i]][["info"]][["n"]], NoOfDigits = ListOfBenfordObjects[[i]][["info"]][["number.of.digits"]], MantissaStat = ListOfBenfordObjects[[i]][["stats"]][["mantissa.arc.test"]][["statistic"]], MantissaPValue = Mantissa.p.Value, AreaDeviation = AreaDeviation)

      OutputDF <- rbind(OutputDF, Temp)
    }

  }

  return(OutputDF)

}

##' @title Plots a list of Benford objects
##' @description This takes a list of \code{benford} objects and generates plots for each of the objects, while providing additional information compared to the default plot.Benford function. 
##' @usage
##' Benford.list.plot(ListOfBenfordObjects, except=c("summation", "abs diff"), multiple=TRUE ,...))
##' @param ListOfBenfordObjects a list of benford objects.
##' @param except a character vector containing what plots not to plot. For details use ?plot.Benford
##' @param multiple plot multiple plots on the same page
##' @export

Benford.list.plot <- function(List, except=c("summation", "abs diff"), multiple=TRUE ,...){
  

for(i in 1:length(names(List))){


  x <- List[[i]]

  old.par <- par(no.readonly=TRUE)
  
  if(class(x)!="Benford") stop("Class(x) must be 'Benford'")
  
  if(!any(except %in% c("digits", "second order", "summation",
                     "mantissa", "chi squared", "abs diff","none"))) {stop("Invalid except name. Type ?plot.Benford for help.")}
  
  if(multiple){
    
    nGraphics <- 8 - length(except)
    
    if(nGraphics<4){rows=1; cols=nGraphics}
    if(nGraphics >=4 & nGraphics<=6){rows=2; cols=3}
    if(nGraphics>6){rows=2; cols=4}
    
    par(mfrow=c(rows,cols))
  }
  
  if(all(except!="digits")){
  plotting.data.vs.benford.v2(x, ...)
  }
  
  if(all(except!="second order")){
  plotting.second.order(x, ...)
  }
  
  if(all(except!="summation")){
  plotting.summation(x, ...)
  }
  
  if(all(except!="mantissa")){
  plotting.ordered.mantissa(x, ...)
  }
  
  if(all(except!="chi squared")){
  plotting.chi_squared.v2(x, ...)
  }
  
  if(all(except!="abs diff")){
  plotting.abs.diff(x, ...)
  }
  
  if(all(except!="ex summation")){
  plotting.ex.summation(x, ...)
  }
  
  plotting.legend.list(x, names(List)[[i]])
  
  par(old.par)
  
}}