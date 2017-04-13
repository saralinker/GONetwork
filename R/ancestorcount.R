# ancestorcount <- function(x, ancestor) {length(ancestor[[x]])}
# getTerm <- function(g){
#  xx[[g]]@Term
#}

# 
# require(GO.db)
# 
# 
# xx <- as.list(GOBPANCESTOR)
# gos <- names(xx)
# #bpancestors <- lapply(X = gos, FUN = ancestorcount, ancestor = GOBPANCESTOR)
# bpanc <- data.frame(go = gos,  bp = unlist(bpancestors),row.names = gos)
# 
# xx <- as.list(GOCCANCESTOR)
# gos <- names(xx)
# #ccancestors <- lapply(X = gos, FUN = ancestorcount, ancestor = GOCCANCESTOR)
# ccanc <- data.frame(go = gos,  bp = unlist(ccancestors),row.names = gos)
# 
# xx <- as.list(GOMFANCESTOR)
# gos <- names(xx)
# #mfancestors <- lapply(X = gos, FUN = ancestorcount, ancestor = GOMFANCESTOR)
# mfanc <- data.frame(go = gos,  bp = unlist(mfancestors),row.names = gos)


 
# GO_space <- data.frame(id = c(as.character(bpanc$go), as.character(ccanc$go), as.character(mfanc$go)), 
#                        space = c(rep("biological_process", nrow(bpanc)),
#                                      rep("cellular_component", nrow(ccanc)),
#                                          rep("molecular_function", nrow(mfanc))
#                                      ),
#                        parents = c(bpanc$bp, ccanc$bp, mfanc$bp)
# )

 #xx <- as.list(GOTERM)
 #GO_space$term <- unlist(lapply(X = GO_space$id,FUN = getTerm))
 
#save(list = c("mouse_gaf","human_gaf","GO_space"), file = "~/Documents/R/packages/GONetwork/R/sysdata.rda",compress = TRUE)
# 
