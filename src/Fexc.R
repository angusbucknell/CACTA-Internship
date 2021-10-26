fexc_genome <-  readDNAStringSet("data/GCA_019097785.1_FRAX_001_PL_genomic.fna")
fexc_genome <- fexc_genome[1:23,]

fexc_CACTA <- packSearch(DNAString("CACTA"), fexc_genome[1],
                         elementLength = c(300, 3500), tsdLength = 3)
fexc_CACTC <- packSearch(DNAString("CACTC"), fexc_genome[1],
                         elementLength = c(300, 3500), tsdLength = 3)
fexc_CACTG <- packSearch(DNAString("CACTG"), fexc_genome[1],
                         elementLength = c(300, 3500), tsdLength = 3)
fexc_CACTT <- packSearch(DNAString("CACTT"), fexc_genome[1],
                         elementLength = c(300, 3500), tsdLength = 3)

fexc_AGAGT <- packSearch(DNAString("AGAGT"), fexc_genome[1],
                         elementLength = c(300, 3500), tsdLength = 3)
fexc_ATAGG <- packSearch(DNAString("ATAGG"), fexc_genome[1],
                         elementLength = c(300, 3500), tsdLength = 3)
fexc_ATCCT <- packSearch(DNAString("ATCCT"), fexc_genome[1],
                         elementLength = c(300, 3500), tsdLength = 3)

fexc_5nt <- data.frame("species"= "Fexc",
                       "seq" = c("CACTA", "CACTC", "CACTG", "CACTT",
                                 "AGAGT", "ATAGG", "ATCCT"),
                       "freq" = c(as.integer(dim(fexc_CACTA)[1]),
                                  as.integer(dim(fexc_CACTC)[1]),
                                  as.integer(dim(fexc_CACTG)[1]),
                                  as.integer(dim(fexc_CACTT)[1]),
                                  as.integer(dim(fexc_AGAGT)[1]),
                                  as.integer(dim(fexc_ATAGG)[1]),
                                  as.integer(dim(fexc_ATCCT)[1])))


fexc_tsd1 <- packSearch(DNAString("CACTT"), fexc_genome[1],
                        elementLength = c(300, 3500), tsdLength = 1)
fexc_tsd2 <- packSearch(DNAString("CACTT"), fexc_genome[1],
                        elementLength = c(300, 3500), tsdLength = 2)
fexc_tsd4 <- packSearch(DNAString("CACTT"), fexc_genome[1],
                        elementLength = c(300, 3500), tsdLength = 4)
fexc_tsd5 <- packSearch(DNAString("CACTT"),fexc_genome[1],
                        elementLength = c(300, 3500), tsdLength = 5)
fexc_tsd <- data.frame("seq" = as.factor("Fexc"),
                       "freq" = c(as.integer(dim(fexc_tsd1)[1]),
                                  as.integer(dim(fexc_tsd2)[1]),
                                  as.integer(dim(fexc_CACTT)[1]),
                                  as.integer(dim(fexc_tsd4)[1]),
                                  as.integer(dim(fexc_tsd5)[1])),
                       "tsd" = c(1,2,3,4,5))
