\name{SCORE}
\alias{SCORE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fast Community Detection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{SCORE(A,K,...)} detects communities in a network by the SCORE (Spectral Clustering On Ratio-of-Eigenvector) algorithm.
}
\usage{
SCORE(A,K,mode="collapse")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{an igraph object that represents an undirected network. }
  \item{K}{an integer indicating the number of communities in the network}
  \item{mode}{an optional character constant that defines the conversion algorithm if the input igraph object A is directed}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
The input graph should be undirected. If the graph is directed, \code{mixedSCORE} will use \code{as.undirected} in the \code{igraph} package for conversion. By default, \code{mixedSCORE} uses the collapse mode. For details on the conversion algorithm, see \code{\link[igraph]{as.undirected}}.

The leading eigenvector of the adjacency matrix of A should have strictly positive components. Otherwise, \code{mixedSCORE} will extract the giant component of the network, and positive components of the leading eigenvector are guaranteed by the Perron-Frobenius theorem. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\code{SCORE} returns an object of class \code{SCORE}, which is a list containing the following components:
\item{community}{the n-dimensional community vector, where n is the number of nodes. The i-th component of the vector indicate which community the i-th node belongs to. }
\item{direct}{a character indicating the conversion algorithm if the graph is converted. If \code{NULL}, the input graph is undirected and does not need conversion. }
\item{gc}{a list of character that contains the names of nodes that are not in the giant component, if the giant component is extracted from the network. If \code{NULL}, the giant component is not extracted.}
\item{R}{the n by K-1 matrix whose columns represent the ratio of eigenvectors of the adjacency matrix.}
}
\references{
%% ~put references to the literature/web site here ~
Jiashun Jin (2012) "Fast Community Detection By SCORE", Manuscript.
}

\author{
Tracy Ke, Lijia Zhou and Qi Zhu. 

Maintainer: Lijia Zhou <zlj@uchicago.edu>, Qi Zhu <qizhu@uchicago.edu>.

%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[igraph]{as.undirected}}, \code{\link{mixedSCORE}}, and \code{\link{topicSCORE}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
# karate contains the network of friendships between the 34
# members of a karate club at a US university, as described 
# by Wayne Zachary in 1977.

# W. W. Zachary, An information flow model for conflict 
# and fission in small groups, Journal of Anthropological
# Research 33, 452-473 (1977).

data(karate)

fit = SCORE(karate,2)
summary(fit)

fit$community

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
