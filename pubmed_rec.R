get_pubmed_ids<-function (pubmed_query_string, api_key = NULL)
{
    old_warn <- options()$warn
    options(warn = -1)
    t_0 <- Sys.time()
    myQuery <- as.character(pubmed_query_string)
    myQuery <- gsub(" ", "+", myQuery, fixed = TRUE)
    myPubmedURL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
        "db=pubmed&term=", myQuery, "&usehistory=y", sep = "")
    if (!is.null(api_key)) {
        myPubmedURL <- paste(myPubmedURL, "&api_key=", api_key,
            sep = "")
    }
    idXML <- NULL
    try_num <- 1
    while (is.null(idXML)) {
        if (try_num > 1)
            Sys.sleep(time = 2)
        t_1 <- Sys.time()
        if (as.numeric(difftime(t_1, t_0, units = "mins")) >
            2) {
            message("Killing the request! Something is not working. Please, try again later")
            return()
        }
        idXML <- tryCatch({
            IDconnect <- suppressWarnings(url(myPubmedURL, open = "rb",
                encoding = "UTF8"))
            idXML <- suppressWarnings(readLines(IDconnect, warn = FALSE,
                encoding = "UTF8"))
            close(IDconnect)
            idXML <- paste(idXML, collapse = "")
            if (grepl("<ERROR>", substr(idXML, 1, 250))) {
                NULL
            }
            else {
                idXML
            }
        }, error = function(e) {
            NULL
        })
        myIDlist <- NULL
        if (!is.null(idXML)) {
            tryCatch({
                myIDlist <- list()
                my_tags <- c("Count", "RetMax", "RetStart", "QueryKey",
                  "WebEnv", "IdList", "TranslationSet", "QueryTranslation")
                for (j in 1:length(my_tags)) {
                  ttag <- my_tags[j]
                  xx <- custom_grep(idXML, tag = ttag, "char")
                  myIDlist[[ttag]] <- xx[1]
                }
                nutag <- "Id"
                xx <- myIDlist[["IdList"]]
                xx <- custom_grep(xx, "Id", format = "list")
                names(xx) <- rep("Id", length(xx))
                myIDlist[["IdList"]] <- xx
                xx <- myIDlist[["TranslationSet"]]
                myIDlist[["TranslationSet"]] <- list()
                nutag <- c("From", "To")
                for (z in nutag) {
                  yy <- custom_grep(xx, z, format = "char")
                  myIDlist[["TranslationSet"]][[z]] <- yy[1]
                }
            }, error = function(e) {
                idXML <- NULL
            })
        }
        if (!is.list(myIDlist)) {
            idXML <- NULL
        }
        try_num <- try_num + 1
    }
    myIDlist[["OriginalQuery"]] <- myQuery
    myIDlist[["APIkey"]] <- api_key
    options(warn = old_warn)
    return(myIDlist)
}