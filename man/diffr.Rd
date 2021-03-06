% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/diffr.R
\name{diffr}
\alias{diffr}
\title{Perform differential splicing analysis}
\usage{
diffr(input, replicate_A, replicate_B, name_A = NULL, name_B = NULL,
  num_lines = 10000, paired = FALSE, filter = TRUE,
  pdf = "input.DIFF_plots", output = NULL, prob = 0.95, min_diff = 0.1,
  min_reads = 10, alpha = 1, beta = 1, size = 500, cores = 1,
  seed = 10, verbose = TRUE)
}
\arguments{
\item{input}{Exact or Partial match to PSI table in output directory}

\item{replicate_A}{Comma-separated string of sample names for group A}

\item{replicate_B}{Comma-separate string of sample names for group B}

\item{name_A}{Name of the first replicate set A}

\item{name_B}{Name of the second replicate set B}

\item{num_lines}{Number of lines to read/process in parallel at a time... lower number = less memory = greater overhead}

\item{paired}{Samples are paired, -a pairOneA,pairTwoA,.. -b pairOneB,pairTwoB,..}

\item{filter}{Filter output for differential events only}

\item{pdf}{Plot visual output (pdf) for differential events into FILE}

\item{output}{Output directory}

\item{prob}{Probability threshold for P( (psi1 - psi2) > x ) > threshold}

\item{min_diff}{Threshold for min diff where P( (psi1 - psi2) > threshold ) > --prob}

\item{min_reads}{Threshold for min reads in a sample (use this flag unless you believe the prior)}

\item{alpha}{First shape parameter for the Beta prior distribution P(psi), Uniform by default}

\item{beta}{Second shape parameter for the Beta prior distribution P(psi), Uniform by default}

\item{size}{Size of the posterior emperical distribution over psi, lower = faster...}

\item{cores}{Number of cores to use for plot processing}

\item{seed}{Seed the RNG for a deterministic result}

\item{verbose}{Enable verbose. Default \code{TRUE}}
}
\description{
Bayesian inference followed by differential analysis of
posterior distributions with respect to PSI/PSU/PIR.
With replicate data, joint posterior distributions for a sample are
estimated from emperical posterior distributions of the
replicates using maximum-likelihood (MLE) fitting.
}
\details{
Originally designed to accept PSI values generated from
\href{https://github.com/vastgroup/vast-tools}{VAST-TOOLS}.

Other formats are planned.

\strong{Statistics Options}

Probably the most important extra options to consider are -r PROB (--prob), -m MINDIFF (--minDiff) and -e MINREADS (--minReads) These represent the stringency criterion for filtering of visual output and textual data sent to file.

The -r flag represents the minimal probability of acceptance that is required to consider a comparison to be 'believable'. By default this is 0.95, but it can be altered depending on stringency requirements.

The -m flag represents the minimum difference between psi1 and psi2 that you will accept, such that we are are sure with at least probability -r that there is a difference of at least -m. -m does not currently alter the output sent to STDOUT, but does filter what is plotted to PDF and printed to file.

The -e flag specifies the minimum number of reads for a sample/event to be compared. In cases where the prior distribution has been methodically calculated and/or is believable beyond an uninformative prior (like the uniform default), this may not be necessary, however it is still highly recommended. The default value for -e is 10, though this could easily be higher.

Additionally, diff allows you to alter the parameters of the conjugate beta prior distribution, which is set as a uniform beta with shape parameters --alpha and --beta as 1 and 1 respectively. Beta shape parameters greater than one alter this probability distribution, and may be more or less applicable to certain uses, see: \href{http://www.wolframalpha.com/input/?i=beta+distribution}{beta distribution} NOTE: When considering differential analysis of event types like intron retention it may be more appropriate to use a custom prior model that is able to more accurately reflect the lower expectation of inclusion levels.

In the case that you have paired samples, where NormalA is dependent on PerturbationA, it is appropriate to use the --paired=TRUE flag. For example when considering NormalA and NormalB, to compare to PerturbationA and PerturbationB, the probability that P( joint_psi1 - joint_psi2 > -m ) is calculated such that NormalA is only compared to PerturbationA, and then NormalB is compared to PerturbationB. No MLE fitting is used in this case.

In all multireplicate cases where --paired=FALSE, the posterior distributions of the individual replicates are used to estimate a 'best fit joint posterior' distribution over psi for each sample.

\strong{Performance Options}

The -s flag can be used to specify the -s SIZE of the emperical posterior distribution to sample, lower numbers decrease accuracy but increase performance.

The diff command is also able to run in parallel.., specify the number of cores to use with -c INT Obviously more cores will increase the speed of diff, though it may increase the RAM usage as well..

Using the -n flag to specify the number of lines to read/process at a time, will set a max threshold to the RAM used by parallel processing with the -c flag. A lower number means that diff will use significantly less memory, however by decreasing -n you have increased the number of times that the mclapply function must calculate the parallel processing overhead. The default is 100, which works well.
}
\examples{
\dontrun{
 diffr("INCLUSION_LEVELS.tab", "Rep1,Rep2", "Rep3,Rep4")
}
}

