#' NNS Term Matrix
#'
#' Generates a term matrix for text classification use in \link{NNS.reg}.
#'
#' @param x Text A two column dataset should be used.  Concatenate text from original sources to comply with format.  Also note the possiblity of factors in \code{"DV"}, so \code{"as.numeric(as.character(...))"} is used to avoid issues.
#' @param oos Out-of-sample text dataset to be classified.
#' @param names Column names for \code{"IV"} and \code{"oos"}.  Defaults to FALSE.
#' @return Returns the text as independent variables \code{"IV"} and the classification as the dependent variable \code{"DV"}.  Out-of-sample independent variables are returned with \code{"OOS"}.
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995}
#' @examples
#' x <- data.frame(cbind(c("sunny", "rainy"), c(1, -1)))
#' NNS.term.matrix(x)
#'
#' ### Concatenate Text with space seperator, cbind with "DV"
#' x <- data.frame(cbind(c("sunny", "rainy"), c("windy", "cloudy"), c(1, -1)))
#' x <- data.frame(cbind(paste(x[ , 1], x[ , 2], sep = " "), as.numeric(as.character(x[ , 3]))))
#' NNS.term.matrix(x)
#'
#'
#' ### NYT Example
#' \dontrun{
#' require(RTextTools)
#' data(NYTimes)
#'
#' ### Concatenate Columns 3 and 4 containing text, with column 5 as DV
#' NYT=data.frame(cbind(paste(NYTimes[ , 3], NYTimes[ , 4], sep = " "),
#'                      as.numeric(as.character(NYTimes[ , 5]))))
#' NNS.term.matrix(NYT)}
#' @export


NNS.term.matrix <- function(x, oos = NULL, names = FALSE){

  p <- length(oos)

  x <- t(t(x))
  n <- nrow(x)

  #Remove commas, etc.
  mgsub <- function(pattern, x, ...) {
      result <- x
      for (i in 1:length(pattern)) {
          result <- gsub(pattern[i], "", result, ...)
      }
      result
  }

  #Use all lowercase to simplify instances
  x[ , 1] <- tolower(mgsub(c(",", ";", ":", "'s", " . "), x[ , 1]))

  unique.vocab <- unique(unlist(strsplit(as.character(x[ , 1]), " ", fixed = TRUE)))

  #Sub with a longer .csv to be called to reduce IVs
  prepositions <- c("a", "in", "of", "our", "the", "is", "for", "with", "we", "this", "it", "but", "was",
                   "at", "to", "on", "aboard", "aside", "by", "means", "spite", "about", "as", "concerning",
                   "instead", "above", "at", "considering", "into", "according", "atop", "despite", "view",
                   "across", "because", "during", "near", "like", "across", "after", "against", "ahead", "along",
                   "alongside", "amid", "among", "apart", "around", "out", "outside", "over", "owing", "past", "prior",
                   "before", "behind", "below", "beneath", "beside", "besides", "between", "beyond", "regarding",
                   "round", "since", "through", "throughout", "till", "down", "except", "from", "addition",
                   "back", "front", "place", "regard", "inside", "together", "toward", "under", "underneath", "until",
                   "nearby", "next", "off", "account", "onto", "top", "opposite", "out", "unto", "up", "within", "without", "what")



  #Remove prepositions
  preps <- unique.vocab %in% c(prepositions)

  if(sum(preps) > 0){
      unique.vocab <- unique.vocab[!preps]
  }

  if(!is.null(oos)){
      oos.preps <- oos %in% c(prepositions)
      if(sum(oos.preps) > 0){
          oos <- oos[!oos.preps]
      }
  }

  NNS.TM <- (t(sapply(1 : length(x[ , 1]), function(i) stringr::str_count(x[i, 1], unique.vocab))))

  if(names){
      colnames(NNS.TM) <- c(unique.vocab)
  }

  if(!is.null(oos)){
      OOS.TM <- (t(sapply(1 : length(oos), function(i) stringr::str_count(oos[i], unique.vocab))))
  if(names){
      colnames(OOS.TM) <- c(unique.vocab)
  }

  return(list("IV" = NNS.TM,
              "DV" = as.numeric(as.character(x[ , 2])),
              "OOS" = OOS.TM))
  } else {
      return(list("IV" = NNS.TM,
                "DV" = as.numeric(as.character(x[ , 2]))))
  }

}
