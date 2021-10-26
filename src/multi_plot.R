all_5nt <- rbind(rbind(dims_5nt, coav_5nt), fexc_5nt)
all_5nt$species <- as.factor(all_5nt$species)
cave_sem <- sd(coav_5nt$freq)/sqrt(dim(coav_5nt)[1])
fexc_sem <- sd(fexc_5nt$freq)/sqrt(dim(fexc_5nt)[1])
all_5nt$sem[1:7] <- sem_5nt
all_5nt$sem[8:14] <- cave_sem 
all_5nt$sem[15:21] <- fexc_sem
matches_5nt_gfx <- ggplot(all_5nt, aes(fill=species, x=seq, y=freq)) + 
  geom_bar(stat="identity", position = "dodge", color="black") +
  labs(x="Search sequence", y="Frequency of matches", fill = "Species") +
  scale_fill_discrete(labels=c("C. avellana", "F. excelsior", "Q. robur")) +
  ylim(c(0, 6300)) +
  geom_errorbar(aes(ymin=freq-sem, ymax=freq+sem), width=.3,position=position_dodge(.9))


all_tsd <- rbind(as.data.frame(table(CACTT$TSD)), as.data.frame(table(coav_CACTT$TSD)))
all_tsd <- rbind(all_tsd, as.data.frame(table(fexc_CACTT$TSD)))
all_tsd$spec <- "Qrob"
all_tsd$spec[63:126] <- "Cave"
all_tsd$spec[127:190] <- "Fexc"
all_tsd_gfx <- ggplot(all_tsd, aes(fill=spec, x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  scale_fill_discrete(labels=c("C. avellana", "F. excelsior", "Q. robur")) +
  labs(x="Search sequence", y="Frequency of matches", fill = "Species")

tsd_line <- rbind(cave_tsd, fexc_ts)



#https://www.ncbi.nlm.nih.gov/assembly/GCF_000002775.4/