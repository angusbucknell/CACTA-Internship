library(Biostrings)
library(packFinder)
library(ggplot2)
library(rtracklayer)


genome <- readDNAStringSet("C:/Users/pbuck/OneDrive/Documents/Qrob_PM1N.fa")
TE <- readGFFAsGRanges("C:/Users/pbuck/OneDrive/Documents/Qrob_PM1N_refTEs.gff")


#### 5nt packFinder ####
CACTA <- packSearch(DNAString("CACTA"), genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
CACTC <- packSearch(DNAString("CACTC"), genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
CACTG <- packSearch(DNAString("CACTG"), genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
CACTT <- packSearch(DNAString("CACTT"), genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)

AGAGT <- packSearch(DNAString("AGAGT"), genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
ATAGG <- packSearch(DNAString("ATAGG"), genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)
ATCCT <- packSearch(DNAString("ATCCT"), genome[1],
                    elementLength = c(300, 3500), tsdLength = 3)

dims_5nt <- data.frame("seq" = c("CACTA", "CACTC", "CACTG", "CACTT",
                                 "AGAGT", "ATAGG", "ATCCT"),
                       "freq" = c(as.integer(dim(CACTA)[1]),
                                  as.integer(dim(CACTC)[1]),
                                  as.integer(dim(CACTG)[1]),
                                  as.integer(dim(CACTT)[1]),
                                  as.integer(dim(AGAGT)[1]),
                                  as.integer(dim(ATAGG)[1]),
                                  as.integer(dim(ATCCT)[1])))
sem_5nt <- sd(dims_5nt$freq)/sqrt(dim(dims_5nt)[1])
barplot(dims_5nt$freq)

##### 5nt frequencies on chromosome 1 GFX #####
dims_5nt_gfx <- ggplot(dims_5nt, aes(seq, freq)) + geom_bar(stat = "identity") +
  xlab("Search sequence") + ylab("Frequency of matches") +
  geom_errorbar(aes(ymin = freq - sem_5nt, ymax = freq + sem_5nt), width=.2,
                position = position_dodge(.9)) + ylim(0, 6000)


#### CACTT TSD analysis ####
tsd1 <- packSearch(DNAString("CACTT"), genome[1],
                   elementLength = c(300, 3500), tsdLength = 1)
tsd2 <- packSearch(DNAString("CACTT"), genome[1],
                   elementLength = c(300, 3500), tsdLength = 2)
tsd3 <- packSearch(DNAString("CACTT"), genome[1],
                   elementLength = c(300, 3500), tsdLength = 3)
tsd4 <- packSearch(DNAString("CACTT"), genome[1],
                   elementLength = c(300, 3500), tsdLength = 4)
tsd5 <- packSearch(DNAString("CACTT"), genome[1],
                   elementLength = c(300, 3500), tsdLength = 5)

tsd_dims <- data.frame("seq" = as.factor("CACTT"),
                       "freq" = c(as.integer(dim(tsd1)[1]),
                                  as.integer(dim(tsd2)[1]),
                                  as.integer(dim(tsd3)[1]),
                                  as.integer(dim(tsd4)[1]),
                                  as.integer(dim(tsd5)[1])),
                       "tsd" = c(1,2,3,4,5))
plot(tsd_dims$freq ~ tsd_dims$tsd)

##### Control TSD analysis #####

AGAGT_tsd1 <- packSearch(DNAString("AGAGT"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 1)
AGAGT_tsd2 <- packSearch(DNAString("AGAGT"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 2)
AGAGT_tsd4 <- packSearch(DNAString("AGAGT"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 4)
AGAGT_tsd5 <- packSearch(DNAString("AGAGT"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 5)

ATAGG_tsd1 <- packSearch(DNAString("ATAGG"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 1)
ATAGG_tsd2 <- packSearch(DNAString("ATAGG"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 2)
ATAGG_tsd4 <- packSearch(DNAString("ATAGG"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 4)
ATAGG_tsd5 <- packSearch(DNAString("ATAGG"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 5)

ATCCT_tsd1 <- packSearch(DNAString("ATCCT"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 1)
ATCCT_tsd2 <- packSearch(DNAString("ATCCT"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 2)
ATCCT_tsd4 <- packSearch(DNAString("ATCCT"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 4)
ATCCT_tsd5 <- packSearch(DNAString("ATCCT"), genome[1],
                         elementLength = c(300, 3500), tsdLength = 5)


seq_vector <- c("AGAGT","AGAGT","AGAGT","AGAGT","AGAGT",
                "ATAGG","ATAGG","ATAGG","ATAGG","ATAGG",
                "ATCCT","ATCCT","ATCCT","ATCCT","ATCCT")
freq_vector <- c(dim(AGAGT_tsd1)[1],dim(AGAGT_tsd2)[1],dim(AGAGT)[1],dim(AGAGT_tsd4)[1],dim(AGAGT_tsd5)[1],
                 dim(ATAGG_tsd1)[1],dim(ATAGG_tsd2)[1],dim(ATAGG)[1],dim(ATAGG_tsd4)[1],dim(ATAGG_tsd5)[1],
                 dim(ATCCT_tsd1)[1],dim(ATCCT_tsd2)[1],dim(ATCCT)[1],dim(ATCCT_tsd4)[1],dim(ATCCT_tsd5)[1])
tsd_dims <- rbind(tsd_dims, data.frame("seq" = as.factor(seq_vector),
                                       "freq" = freq_vector,
                                       "tsd" = c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)))

##### 5nt TSD GFX #####
tsd_gfx <- ggplot(tsd_dims, aes(x=tsd, y=freq, group=seq, color=seq)) + 
  geom_line(size=1) + geom_point(size=2) +
  xlab("Terminal Site Duplication (TSD) length") + ylab("Frequency of matches")

#### CACTT genome analysis ####

CACTT_genome <- packSearch(DNAString("CACTT"), genome,
                           elementLength = c(300, 3500), tsdLength = 3)
CACTT_genome <- CACTT_genome[startsWith(CACTT_genome$seqnames, "Qrob_Chr"),]
write.csv(CACTT_genome,"C:/Users/pbuck/OneDrive/Documents/CACTT.csv")

genome_freq <- as.data.frame(table(as.factor(CACTT_genome[startsWith(CACTT_genome$seqnames, "Qrob_Chr"),]$seqnames))/(width(genome[1:12])/1000000))
sem_chr <- sd(genome_freq$Freq)/sqrt(dim(genome_freq)[1])

##### Match frequencies across chromosomes GFX #####
chr_gfx <- ggplot(genome_freq, aes(Var1, Freq)) + geom_bar(stat="identity") +
  xlab("Chromosome") + ylab("Matches per 1 Mbp of chromosome length") +
  ylim(0, 100) + scale_x_discrete(labels= 1:12)

##### Long CACTT analysis #####
CACTT_long <- packSearch(DNAString("CACTT"), genome,
                           elementLength = c(3000, 5500), tsdLength = 3)
CACTT_long <- makeGRangesFromDataFrame(CACTT_long)
CACTT_long <- CACTT_long[startsWith(as.character(seqnames(CACTT_long)), "Qrob_Chr")]
long_seq <- genome[CACTT_long]



#### vmatch genome & GFF ####

vmatch <- vmatchPattern("CACTT", genome)
vmatch <- unlist(vmatch[startsWith(names(vmatch), "Qrob_Chr")])

vmatch_TE <- vmatchPattern("CACTT", genome[TE])
vmatch_TE <- vmatch_TE[startsWith(names(vmatch_TE), "Qrob_Chr")]
vmatch_TE <- vmatch_TE[elementNROWS(vmatch_TE) > 0]
vmatch_TE <- unique(as.data.frame(unlist(vmatch_TE)))
vmatch_TE <- IRanges(start = vmatch_TE$start, end = vmatch_TE$end,
                     width = vmatch_TE$width, names = vmatch_TE$names)
