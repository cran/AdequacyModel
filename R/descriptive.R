descriptive <- function(x)
{
    
  statistics   <- list(mean = round(mean(x), 5),median = round(median(x), 5),
                       mode = round(searchMode(x)$modalValue, 5), variance = round(var(x), 5),
                       Skewness = round(skewness(x), 5), Kurtosis = round(kurtosis(x), 5), minimum = round(min(x), 5),
                       maximum = round(max(x), 5), n = length(x))
  return(statistics)
}