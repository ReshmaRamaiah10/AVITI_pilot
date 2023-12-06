library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly=TRUE)

dat <- read.csv(args[1])
tab3 <- dat[dat$Rat>log(2),]
tab3 <- tab3[order(tab3$Rat, decreasing=T),]

png(file=args[2], width=800, height=800)
ggplot(dat, aes(x=AVITI, y=Illumina, color=Diff)) +
    geom_point(size=2) +
    scale_color_manual(values=c("False"="coral1", "True"="azure4")) +
    xlab("ln TPM AVITI") +
    ylab("ln TPM Illumina") +
    theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=20), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)), axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_text_repel(aes(x=AVITI, y=Illumina, label=gene), data=tab3[1:20,], color="black", max.overlaps=20, size=5) +
    annotate("text", x=1.5, y=10, label=paste0("Pearson's r = ", args[3]), size=10)
dev.off()
