ancestorcount <- function(x, ancestor) {length(ancestor[[x]])}



xx <- as.list(GOBPANCESTOR)
gos <- names(xx)
#bpancestors <- lapply(X = gos, FUN = ancestorcount, ancestor = GOBPANCESTOR)
bpanc <- data.frame(go = gos,  bp = unlist(bpancestors),row.names = gos)

xx <- as.list(GOCCANCESTOR)
gos <- names(xx)
#ccancestors <- lapply(X = gos, FUN = ancestorcount, ancestor = GOCCANCESTOR)
ccanc <- data.frame(go = gos,  bp = unlist(ccancestors),row.names = gos)

xx <- as.list(GOMFANCESTOR)
gos <- names(xx)
#mfancestors <- lapply(X = gos, FUN = ancestorcount, ancestor = GOMFANCESTOR)
mfanc <- data.frame(go = gos,  bp = unlist(mfancestors),row.names = gos)


#save(list = c("bpanc","ccanc","mfanc"), file = "~/Documents/R/packages/GONetwork/R/sysdata.rda",compress = TRUE)

