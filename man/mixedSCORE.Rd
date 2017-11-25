\name{mixedSCORE}
\alias{mixedSCORE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Mixed Membership Estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{mixedSCORE(A,K,...)} estimates the membership matrix in a network by the mixed SCORE (Spectral Clustering On Ratio-of-Eigenvector) algorithm. 
}
\usage{
mixedSCORE(A, K, L_theo = FALSE, mode = c("collapse", "each", "mutual") )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{an igraph object that represents an undirected network. }
  \item{K}{an integer indicating the number of communities in the network}
  \item{L_theo}{logical. If \code{TRUE}, use the theoretical approach to estimate the tuning parameter L, which is the number of cluster to be used in the k-means step.}
  \item{mode}{an optional character constant that defines the conversion algorithm if the input igraph object A is directed}
}

\details{
The input graph should be undirected. If the graph is directed, \code{mixedSCORE} will use \code{as.undirected} in the \code{igraph} package for conversion. By default, \code{mixedSCORE} uses the collapse mode. For details on the conversion algorithm, see \code{\link[igraph]{as.undirected}}.

The leading eigenvector of the adjacency matrix of A should have strictly positive components. Otherwise, \code{mixedSCORE} will extract the giant component of the network, and positive components of the leading eigenvector are guaranteed by the Perron-Frobenius theorem. 
}

%%  ~~ If necessary, more details than the description above ~~
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\code{mixedSCORE} returns an object of class \code{mSCORE}, which is a list containing the following components:
\item{direct}{a character indicating the conversion algorithm if the graph is converted. If \code{NULL}, the input graph is undirected and does not need conversion. }
\item{gc}{a list of character that contains the names of nodes that are not in the giant component, if the giant component is extracted from the network. If \code{NULL}, the giant component is not extracted.}
\item{member}{the estimate of the n by K mixed membership matrix, where n is the number of nodes and K is the number of communities. The (i,j) entry corresponds to the probability that the i-th node belongs to the j-th community. }
\item{l}{the optimal tuning parameter, which is an integer greater than K}
\item{R}{the n by K-1 matrix whose columns represent the ratio of eigenvectors of the adjacency matrix.}
\item{vertice}{the matrix whose columns are the estimated vertices in the vertice hunting step }
\item{P}{the estimate of the K by K matrix P in the Degree Corrected Block model.}
\item{B}{the estimate of the K by K matrix B in the Degree Corrected Block model. }
}

\references{
Jin, Ke and Luo (2017) "Estimating Network Membership By Simplex Vertices Hunting", Manuscript.
%% ~put references to the literature/web site here ~
}
\author{
Tracy Ke, Lijia Zhou and Qi Zhu. 

Maintainer: Lijia Zhou <zlj@uchicago.edu>, Qi Zhu <qizhu@uchicago.edu>.

%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[igraph]{as.undirected}}, \code{\link{SCORE}}, and \code{\link{topicSCORE}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
# football contains the network of American football games 
# between Division IA colleges during regular season 
# Fall 2000, as compiled by M. Girvan and M. Newman.

# M. Girvan and M. E. J. Newman, Community structure in 
# social and biological networks, Proc. Natl. Acad. Sci. 
# USA 99, 7821-7826 (2002).

data(football)

fit = mixedSCORE(football,4) 
summary(fit)

fit$member

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
