#Corylus avellana (European hazelnut)

coav_genome <-  readDNAStringSet("data/GCA_901000735.2_CavTom2PMs-1.0_genomic.fna")
fexc_genome <- fexc_genome[1:11,]

coav_CACTA <- packSearch(DNAString("CACTA"), coav_genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
coav_CACTC <- packSearch(DNAString("CACTC"), coav_genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
coav_CACTG <- packSearch(DNAString("CACTG"), coav_genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
coav_CACTT <- packSearch(DNAString("CACTT"), coav_genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)

coav_AGAGT <- packSearch(DNAString("AGAGT"), coav_genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
coav_ATAGG <- packSearch(DNAString("ATAGG"), coav_genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
coav_ATCCT <- packSearch(DNAString("ATCCT"), coav_genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)

coav_5nt <- data.frame("species"="Cave",
                       "seq" = c("CACTA", "CACTC", "CACTG", "CACTT",
                                 "AGAGT", "ATAGG", "ATCCT"),
                       "freq" = c(as.integer(dim(coav_CACTA)[1]),
                                  as.integer(dim(coav_CACTC)[1]),
                                  as.integer(dim(coav_CACTG)[1]),
                                  as.integer(dim(coav_CACTT)[1]),
                                  as.integer(dim(coav_AGAGT)[1]),
                                  as.integer(dim(coav_ATAGG)[1]),
                                  as.integer(dim(coav_ATCCT)[1])))

cave_tsd1 <- packSearch(DNAString("CACTT"), coav_genome[1],
                   elementLength = c(300, 3500), tsdLength = 1)
cave_tsd2 <- packSearch(DNAString("CACTT"), coav_genome[1],
                   elementLength = c(300, 3500), tsdLength = 2)
cave_tsd4 <- packSearch(DNAString("CACTT"), coav_genome[1],
                   elementLength = c(300, 3500), tsdLength = 4)
cave_tsd5 <- packSearch(DNAString("CACTT"),coav_genome[1],
                   elementLength = c(300, 3500), tsdLength = 5)
cave_tsd <- data.frame("spec" = as.factor("Cave"),
                       "freq" = c(as.integer(dim(cave_tsd1)[1]),
                                  as.integer(dim(cave_tsd2)[1]),
                                  as.integer(dim(coav_CACTT)[1]),
                                  as.integer(dim(cave_tsd4)[1]),
                                  as.integer(dim(cave_tsd5)[1])),
                       "tsd" = c(1,2,3,4,5))



