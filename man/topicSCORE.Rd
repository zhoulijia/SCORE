\name{topicSCORE}
\alias{topicSCORE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal Topic Estimation 
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{topicSCORE} uses an SVD approach to perform an optimal topic estimation.

}
\usage{
topicSCORE(D,K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{the p by n text corpus matrix, where p is the number of common words and n is the number of documents. The (i,j) entry corresponds to the observed fraction of word i in document j.}
  \item{K}{an integer indicating the number of topics}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

It returns an object of class \code{tSCORE}, which is a list containing the following components:
\item{topic}{the estimate of the p by K topic matrix, whose columns correspond to the expected frequencies of words in a document that discuss a certain topic. In particular, the (i,j) entry is the expected frequency of word i in a document that discuss topic j. }
\item{weight}{the p by K weight matrix Pi, whose rows correspond to the weight vector of a certain word. }
\item{R}{the n by K-1 matrix whose columns represent the ratio of singular vectors of the text corpus matrix.}
\item{vertice}{the matrix whose columns are the estimated vertices in the vertice hunting step }
}

\references{
%% ~put references to the literature/web site here ~
Ke and Wang (2017) "A New SVD Approach To Optimal Topic Estimation".
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# AP (Harman, 1993) consists of 2246 news articles with a vocabulary of 
# 10473 words. After preprocessing, approximately 8000 words are kept.
# See reference for more details.

data(AP)           # the AP file is stored in .mat format
obj = readMat(AP)
D = obj$D          # text corpus matrix
vocab = obj$volc   # list of vocabulary

fit = topicSCORE(D,3)
summary(fit)

fit$topic
fit$weight

word1 = fit$topic[,1]
idx = c()
w = tail(sort(word1),15)
for(i in 15:1){
  idx = c(idx,which(word1 == w[i]))
}

vocab[idx]

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
