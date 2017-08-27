km.means.quantiles <- function(km.object){
    print("When assuming career span equals 20")
    print(km.object, print.rmean=TRUE, rmean=20)
    print("When assumming career span equals 35")
    print(km.object, print.rmean=TRUE, rmean=35) # Projecting to a tradition career length
    print(paste("25 Percent of the sample (or this category) lasts",quantile(km.object,.25)$quantile,"or fewer years"))
    print(paste("10 Percent of the sample (or this category) lasts",quantile(km.object,.10)$quantile,"or fewer years"))
}
