#' Perform differential splicing analysis
#'
#' Bayesian inference followed by differential analysis of
#' posterior distributions with respect to PSI/PSU/PIR.
#' With replicate data, joint posterior distributions for a sample are
#' estimated from emperical posterior distributions of the
#' replicates using maximum-likelihood (MLE) fitting.
#'
#' @details
#' Originally designed to accept PSI values generated from
#' \href{https://github.com/vastgroup/vast-tools}{VAST-TOOLS}.
#'
#' Other formats are planned.
#'
#' \strong{Statistics Options}
#'
#' Probably the most important extra options to consider are -r PROB (--prob), -m MINDIFF (--minDiff) and -e MINREADS (--minReads) These represent the stringency criterion for filtering of visual output and textual data sent to file.
#'
#' The -r flag represents the minimal probability of acceptance that is required to consider a comparison to be 'believable'. By default this is 0.95, but it can be altered depending on stringency requirements.
#'
#' The -m flag represents the minimum difference between psi1 and psi2 that you will accept, such that we are are sure with at least probability -r that there is a difference of at least -m. -m does not currently alter the output sent to STDOUT, but does filter what is plotted to PDF and printed to file.
#'
#' The -e flag specifies the minimum number of reads for a sample/event to be compared. In cases where the prior distribution has been methodically calculated and/or is believable beyond an uninformative prior (like the uniform default), this may not be necessary, however it is still highly recommended. The default value for -e is 10, though this could easily be higher.
#'
#' Additionally, diff allows you to alter the parameters of the conjugate beta prior distribution, which is set as a uniform beta with shape parameters --alpha and --beta as 1 and 1 respectively. Beta shape parameters greater than one alter this probability distribution, and may be more or less applicable to certain uses, see: \href{http://www.wolframalpha.com/input/?i=beta+distribution}{beta distribution} NOTE: When considering differential analysis of event types like intron retention it may be more appropriate to use a custom prior model that is able to more accurately reflect the lower expectation of inclusion levels.
#'
#' In the case that you have paired samples, where NormalA is dependent on PerturbationA, it is appropriate to use the --paired=TRUE flag. For example when considering NormalA and NormalB, to compare to PerturbationA and PerturbationB, the probability that P( joint_psi1 - joint_psi2 > -m ) is calculated such that NormalA is only compared to PerturbationA, and then NormalB is compared to PerturbationB. No MLE fitting is used in this case.
#'
#' In all multireplicate cases where --paired=FALSE, the posterior distributions of the individual replicates are used to estimate a 'best fit joint posterior' distribution over psi for each sample.
#'
#' \strong{Performance Options}
#'
#' The -s flag can be used to specify the -s SIZE of the emperical posterior distribution to sample, lower numbers decrease accuracy but increase performance.
#'
#' The diff command is also able to run in parallel.., specify the number of cores to use with -c INT Obviously more cores will increase the speed of diff, though it may increase the RAM usage as well..
#'
#' Using the -n flag to specify the number of lines to read/process at a time, will set a max threshold to the RAM used by parallel processing with the -c flag. A lower number means that diff will use significantly less memory, however by decreasing -n you have increased the number of times that the mclapply function must calculate the parallel processing overhead. The default is 100, which works well.
#'
#' @param input Exact or Partial match to PSI table in output directory
#' @param replicate_A Comma-separated string of sample names for group A
#' @param replicate_B Comma-separate string of sample names for group B
#' @param name_A Name of the first replicate set A
#' @param name_B Name of the second replicate set B
#' @param num_lines Number of lines to read/process in parallel at a time... lower number = less memory = greater overhead
#' @param paired Samples are paired, -a pairOneA,pairTwoA,.. -b pairOneB,pairTwoB,..
#' @param filter Filter output for differential events only
#' @param pdf Plot visual output (pdf) for differential events into FILE
#' @param output Output directory
#' @param prob Probability threshold for P( (psi1 - psi2) > x ) > threshold
#' @param min_diff Threshold for min diff where P( (psi1 - psi2) > threshold ) > --prob
#' @param min_reads Threshold for min reads in a sample (use this flag unless you believe the prior)
#' @param alpha First shape parameter for the Beta prior distribution P(psi), Uniform by default
#' @param beta Second shape parameter for the Beta prior distribution P(psi), Uniform by default
#' @param size Size of the posterior emperical distribution over psi, lower = faster...
#' @param cores Number of cores to use for plot processing
#' @param seed Seed the RNG for a deterministic result
#' @param verbose Enable verbose. Default \code{TRUE}
#'
#' @export
#' @import grid
#' @import MASS
#' @import reshape2
#' @import ggplot2
#' @import parallel
#' @examples
#' \dontrun{
#'  diffr("INCLUSION_LEVELS.tab", "Rep1,Rep2", "Rep3,Rep4")
#' }
diffr <- function(input, replicate_A, replicate_B,
                  name_A = NULL, name_B = NULL,
                  num_lines = 10000, paired = FALSE,
                  filter = TRUE, pdf = "input.DIFF_plots",
                  output = NULL, prob = 0.95, min_diff = 0.1, min_reads = 10,
                  alpha = 1, beta = 1, size = 500,
                  cores = 1, seed = 10, verbose = TRUE) {
  ## seed RNG
  set.seed(seed)

  ## try and find the input file if they aren't exact
  if(!file.exists(input)) {
    potentialFiles <- Sys.glob( paste(c("*",input,"*"), collapse="") )
    if( length( potentialFiles ) >= 1) {
      # now sort and take the one with the 'biggest' number of samples
      potentialFiles_sort <- rev( sort( potentialFiles ) )
      input <- potentialFiles_sort[1]
    } else {
      # Still can't find input after searching...
      stop("[vast diff error]: No input file given!")
    }
  }

  ## Setting input files.
  inputFile <- file( input, 'r' )


  firstRepSet <- unlist(strsplit( as.character(replicate_A) , "," ))
  secondRepSet <- unlist(strsplit( as.character(replicate_B), "," ))

  if( length(firstRepSet) <= 0 ||
      length(secondRepSet) <= 0) {
    stop("[vast diff error]: No replicate sample names given!!! -a sampA,sampB -b sampC,sampD")
  }

  # Set number of replicates
  firstRepN <- length(firstRepSet)
  secondRepN <- length(secondRepSet)

  # Make sure there are sample names
  if(is.null(name_A ) ) {
    name_A <- firstRepSet[1]
  }
  if(is.null( name_B) ) {
    name_B<- secondRepSet[1]
  }
  # Set output sample names for plot
  sampOneName <- substr(name_A, 1, 9)
  sampTwoName <- substr(name_B, 1, 9)

  ## INITIALIZE LISTS ##
  shapeFirst <- vector("list", firstRepN)
  shapeSecond <- vector("list", secondRepN)

  psiFirst <- vector("list", firstRepN)
  psiSecond <- vector("list", secondRepN)

  # Get header
  head <- readLines( inputFile, n=1 )
  head_n <- unlist( strsplit( head, "\t" ) )

  #DEPRECATED -TSW 03/26/2015
  # if we are to filter to stdout, then print header
  if( filter ) {
    #  writeLines(head, stdout())
    writeLines(sprintf("GENE\tEVENT\t%s\t%s\tE[dPsi]\tMV[dPsi]_at_%s", name_A, name_B, prob), stdout())
  }

  # check if header is correct..  TODO

  # Indexes of samples of interest
  repAind <- which( head_n %in% firstRepSet  )
  repBind <- which( head_n %in% secondRepSet )

  if(length(repAind) == 0 ||
     length(repBind) == 0) {
    stop("[vast diff error]: Incorrect sampleNames given, One or more do not exist!!!\n")
  }

  # Indexes of Quals
  repA.qualInd <- repAind + 1
  repB.qualInd <- repBind + 1

  # make sure this succeeded  TODO

  # CONST
  alphaList <- seq(0,1,0.01)

  ### TMP OUT
  if(pdf == "input.DIFF_plots") {
    pdfname <- sub("\\.[^.]*(\\.gz)?$", ".DIFF_plots.pdf", basename(input))
    signame <- sub("\\.[^.]*(\\.gz)?$", ".DIFF_sig.txt", basename(input))
  } else {
    pdfname <- paste(c(pdf, ".pdf"), collapse="")
    signame <- paste(c(pdf, ".txt"), collapse="")
  }

  sighandle <- file(signame, "w")

  pdf(pdfname, width=7, height=3.5, family="sans", compress=FALSE)
  write(head, sighandle)

  ### BEGIN READ INPUT ###
  # Iterate through input, 'num_lines' at a time to reduce overhead/memory
  while(length( lines <- readLines(inputFile, n=num_lines) ) > 0) {

    # use parallel computing to store plots in plotListed
    # then print them to the pdf afterwards before next chunk of num_lines from file.
    plotListed <- vector("list", length(lines))
    eventTitleListed <- vector("list", length(lines))

    plotListed <- mclapply(1:length(lines), function(i) {

      tabLine <- unlist( strsplit( lines[i], "\t" ) )
      #writeLines(paste(tabLine[repA.qualInd], collapse="\t"), stderr());

      # Posterior parameters... Prior given from command line --alpha, --beta
      shapeFirst <- lapply( tabLine[repA.qualInd], function(x) {
        parseQual(x, alpha, beta)
      } )
      shapeSecond <- lapply( tabLine[repB.qualInd], function(x) {
        parseQual(x, alpha, beta)
      } )

      totalFirst <- unlist(lapply( shapeFirst, function(x) { x[1] + x[2] }))
      totalSecond <- unlist(lapply( shapeSecond, function(x) { x[1] + x[2] }))

      # if no data, next;
      if( all(totalFirst < (min_reads + alpha + beta)) ||
          all(totalSecond < (min_reads + alpha + beta)) ) {
        return(NULL)
      }

      firstShapeMat <- do.call(rbind, shapeFirst)
      secondShapeMat <- do.call(rbind, shapeSecond)

      firstShapeAve <- c( mean(firstShapeMat[,1]), mean(firstShapeMat[,2]) )
      secondShapeAve <- c( mean(secondShapeMat[,1]), mean(secondShapeMat[,2]) )

      # Sample Posterior Distributions
      psiFirst <- lapply( shapeFirst, function(x) {
        #sample here from rbeta(N, alpha, beta) if > -e
        if(x[1]+x[2] < min_reads) { return(NULL) }
        rbeta(size, shape1=x[1], shape2=x[2])
      })

      psiSecond <- lapply( shapeSecond, function(x) {
        #sample here from rbeta(N, alpha, beta)
        if(x[1]+x[2] < min_reads) { return(NULL) }
        rbeta(size, shape1=x[1], shape2=x[2])
      })

      # calculate expected value of psi for each replicate
      expFirst <- unlist(lapply(shapeFirst, function(x) {
        if(x[1]+x[2] < min_reads) { return(NULL) }
        x[1] / (x[1] + x[2])
      }))
      expSecond <- unlist(lapply(shapeSecond, function(x) {
        if(x[1]+x[2] < min_reads) { return(NULL) }
        x[1] / (x[1] + x[2])
      }))

      if(paired) { #make sure both samples have a non-NULL replicate
        for(lstInd in 1:length(psiFirst)) {
          if(is.null(psiFirst[[lstInd]]) || is.null(psiSecond[[lstInd]])) {
            psiFirst[[lstInd]] <- NULL
            psiSecond[[lstInd]] <- NULL
          }
        }
      }

      # Create non-parametric Joint Distributions
      psiFirstComb <- do.call(c, psiFirst)
      psiSecondComb <- do.call(c, psiSecond)

      if( length(psiFirstComb) <= 0 || length(psiSecondComb) <= 0 ) { return(NULL) }

      #    print(length(psiFirstComb))

      # if they aren't paired, then shuffle the joint distributions...
      if( !paired ) {
        paramFirst <- try (suppressWarnings(
          fitdistr(psiFirstComb,
                   "beta",
                   list( shape1=firstShapeAve[1], shape2=firstShapeAve[2])
          )$estimate ), TRUE )
        paramSecond <- try (suppressWarnings(
          fitdistr(psiSecondComb,
                   "beta",
                   list( shape1=secondShapeAve[1], shape2=secondShapeAve[2])
          )$estimate ), TRUE )
        # if optimization fails its because the distribution is too narrow
        # in which case our starting shapes should already be good enough
        if(class(paramFirst) != "try-error") {
          psiFirstComb <- rbeta(size, shape1=paramFirst[1], shape2=paramFirst[2])
        }
        if(class(paramSecond) != "try-error") {
          psiSecondComb <- rbeta(size, shape1=paramSecond[1], shape2=paramSecond[2])
        }
      }

      # get emperical posterior median of psi
      medOne <- median(psiFirstComb)
      medTwo <- median(psiSecondComb)

      # look for a max difference given prob cutoff...
      if(medOne > medTwo) {
        max <- maxDiff(psiFirstComb, psiSecondComb, prob)
      } else {
        max <- maxDiff(psiSecondComb, psiFirstComb, prob)
      }
      #    writeLines(lines[i], stderr()) ### DEBUGGING

      # SIGNIFICANT from here on out:

      if( filter ) {
        filtOut <- sprintf("%s\t%s\t%f\t%f\t%f\t%s", tabLine[1], tabLine[2], medOne, medTwo, medOne - medTwo, round(max,2))
      } else {
        filtOut <- NULL
      }

      # check for significant difference
      if(max < min_diff) {
        # non-sig, return null plots and text output
        return(list(NULL, NULL, NULL, NULL, filtOut))
      } else {
        sigInd <- i

        eventTitle <- paste(c("Gene: ", tabLine[1], "  Event: ", tabLine[2]), collapse="")
        eventCoord <- paste(c("Coordinates: ", tabLine[3]), collapse="")
        #    eventTitleListed[[i]] <- paste(c("Gene: ", tabLine[1], "     ", "Event: ", tabLine[2]), collapse="")

        # Print visual output to pdf;
        if( medOne > medTwo ) {
          retPlot <- plotDiff(psiFirstComb, psiSecondComb, expFirst, expSecond, max, medOne, medTwo, sampOneName, sampTwoName , FALSE, alphaList)
        } else {
          retPlot <- plotDiff(psiSecondComb, psiFirstComb, expFirst, expSecond, max, medTwo, medOne, sampTwoName, sampOneName , TRUE, alphaList)
        }
        # sig event return
        return(list(retPlot, eventTitle, eventCoord, sigInd, filtOut))  #return of mclapply function
      }
    }, mc.cores=cores, mc.preschedule=TRUE, mc.cleanup=TRUE) #End For

    for(it in 1:length(lines)) {
      if(is.null(plotListed[[it]])) { next; }
      # PRINT MAIN OUTPUT
      if(!is.null(plotListed[[it]][[5]])) {
        write( plotListed[[it]][[5]], stdout() )
      }

      # PRINT SIG OUTPUT
      if(!is.null(plotListed[[it]][[4]])) {
        writeLines( lines[ plotListed[[it]][[4]] ], sighandle )
      }

      # PRINT LIST OF PLOTS.
      if(is.null(plotListed[[it]][[1]])) { next; }
      plotPrint(plotListed[[it]][[2]], plotListed[[it]][[3]], plotListed[[it]][[1]])
    }

  } #End While

  garbage <- dev.off()

  flush(sighandle)
  close(sighandle)
  close(inputFile)

}
