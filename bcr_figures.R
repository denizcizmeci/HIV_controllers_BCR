# Extended somatic hypermutation and clonal evolution in B-cells of HIV controllers with broadly neutralizing antibodies
# Deniz Cizmeci, 07/03/2020
# Ragon Institute

# Load data ---------------------------------------------------------------
#setwd("~/Documents/projects/bcr/hiv")

load("./data/bcr.RData")

col.tn <- "#fdc086"
col.nn <- "#a6cee3"

library(ggpubr)
library(ggsci)
library(ggrepel)
library(patchwork)
library(paletteer)

library(dplyr)
library(tidyr)

bcr$gr <- as.vector(bcr$group)
bcr$gr[which(bcr$gr == "Non-Neutralizer")] <- "NN"
bcr$gr[which(bcr$gr == "Mid-Neutralizer")] <- "MN"
bcr$gr[which(bcr$gr == "Top-Neutralizer")] <- "TN"

bcr$gr <- factor(bcr$gr, levels = c("NN", "MN", "TN"))

hiv <- droplevels(subset(bcr, bcr$group != "Mid-Neutralizer"))

subjects <- unique(hiv[,c("donor_identifier", "gr", "n_totSeq")])

# cdr3 length -------------------------------------------------------------

hiv$h_cdr3_length <- nchar(as.vector(hiv$h_cdr3))

table(hiv$h_cdr3_length, hiv$gr)

cdr3 <- as.data.frame(table(hiv$h_cdr3_length, hiv$gr))

nn.cdr3 <- subset(cdr3, cdr3$Var2 == "NN")
nn.cdr3$perc <- (nn.cdr3$Freq/2707)*100
tn.cdr3 <- subset(cdr3, cdr3$Var2 == "TN")
tn.cdr3$perc <- (tn.cdr3$Freq/5771)*100

cdr3 <- rbind(nn.cdr3, tn.cdr3)


pdf("cdr3length.pdf", width= 7, height =3)
ggplot(cdr3, aes(Var1, perc)) +
  geom_bar(stat = "identity", aes(fill = Var2), color = "black", size = 0.15,
           position = position_dodge(0.8), width = 0.7, show.legend = FALSE) +
  #coord_flip()+scale_x_discrete(limits = rev(levels(fam$h_cdr3_length)))+
  scale_fill_manual(values = c(col.nn, col.tn))+ xlab("")+
  theme_classic()+ ylab(bquote('Percent of all sequences (%)'))+
  xlab("CDRH3 length (aa)")+
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size =10),
        axis.text.y = element_text(face = "bold", size = 12))
dev.off()


# Figure2 number of mutations ---------------------------------------------

mut <- as.data.frame(table(hiv$h_mut_total, hiv$gr))
table(hiv$gr)
nn.mut <- subset(mut, mut$Var2 == "NN")
nn.mut$perc <- (nn.mut$Freq/2707)*100
tn.mut <- subset(mut, mut$Var2 == "TN")
tn.mut$perc <- (tn.mut$Freq/5771)*100

mut <- rbind(nn.mut, tn.mut)

mut$count <- as.numeric(mut$Var1)

pdf("hm.pdf", width= 8, height =3)
ggplot(mut, aes(count, perc)) +
  geom_bar(stat = "identity", aes(fill = Var2), color = "black", size = 0.15,
           position = position_dodge(0.8), width = 0.7, show.legend = FALSE) +
  #coord_flip()+scale_x_discrete(limits = rev(levels(fam$h_cdr3_length)))+
  scale_fill_manual(values = c(col.nn, col.tn))+ xlab("")+
  theme_classic()+ ylab(bquote('Percent of\nall sequences (%)'))+
  xlab("Number of mutations")+
  #scale_x_continuous(breaks=seq(0, 150, by = 50))
  theme(axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size =12, angle = 0), 
        axis.text.y = element_text( size = 12))
dev.off()

library(plyr)

bcr$Neutralization.breadth <- paste0(round_any(bcr$breadth, 10, f = ceiling), "%")
bcr$Neutralization.breadth <- factor(bcr$Neutralization.breadth, 
                                     levels = c("0%","10%","20%","30%", "40%",
                                                "50%","60%","70%","80%","90%","100%"))

col.breadth <- c("royalblue4","royalblue1", "#778BD9", "#87d5f5",  "#c7e6f2", "#edf1f2", 
                 "#f7f1d7", "lemonchiffon1", "#FAE389","#F5A463", "#AD5E3E")

pdf("IgH_mutations.pdf", width= 8, height =3)
p <- ggplot(bcr, aes(donor_identifier, h_mut_total))
p + geom_violin(scale = "width", aes(fill = Neutralization.breadth),alpha = 0.8)+ geom_jitter(height = 0, width = 0.1, size = 0.1 , 
                                                                                              color = "gray32")+ 
  geom_boxplot(outlier.shape = NA, width = 0.2, aes(fill = Neutralization.breadth), alpha = 0.8,size = 0.25, show.legend = FALSE) + 
  #facet_wrap(. ~ group,scales = "free_x") +
  scale_x_discrete(limits = rev(levels(bcr$donor_identifier)))+
  scale_fill_manual(drop = FALSE,values= col.breadth,
                    name="Neutralization\nbreadth",
                    breaks=levels(bcr$Neutralization.breadth),
                    labels=c("0%","10%","20%","30%", "40%",
                             "50%","60%","70%","80%","90%","100%")) +
  theme_classic() + 
  theme(axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(color = "black",size = 12, angle = 45, hjust = 1), 
        axis.text.y = element_text( size = 12),
        plot.title = element_text(size = 6, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.position="right",
        strip.text.x = element_text(size = 12),
        legend.justification="left",
        legend.text=element_text(size=12),
        legend.key.size = unit(0.4, "cm"))+
  #guides(fill=guide_legend(ncol=4)) +
  guides(fill = guide_legend(reverse=T))+
  xlab("Study Participant ID") + ylab("Number of mutations")
dev.off()

#### Light chain

mut <- as.data.frame(table(hiv$l_mut_total, hiv$gr))
table(hiv$gr)
nn.mut <- subset(mut, mut$Var2 == "NN")
nn.mut$perc <- (nn.mut$Freq/2707)*100
tn.mut <- subset(mut, mut$Var2 == "TN")
tn.mut$perc <- (tn.mut$Freq/5771)*100

mut <- rbind(nn.mut, tn.mut)

mut$count <- as.numeric(mut$Var1)

pdf("lm.pdf", width= 8, height =3)
ggplot(mut, aes(count, perc)) +
  geom_bar(stat = "identity", aes(fill = Var2), color = "black", size = 0.15,
           position = position_dodge(0.8), width = 0.7, show.legend = FALSE) +
  #coord_flip()+scale_x_discrete(limits = rev(levels(fam$h_cdr3_length)))+
  scale_fill_manual(values = c(col.nn, col.tn))+ xlab("")+
  theme_classic()+ ylab(bquote('Percent of\nall sequences (%)'))+
  xlab("Number of mutations")+
  #scale_x_continuous(breaks=seq(0, 150, by = 50))
  theme(axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size =12, angle = 0), 
        axis.text.y = element_text( size = 12))
dev.off()

pdf("IgL_mutations.pdf", width= 8, height =3)
p <- ggplot(bcr, aes(donor_identifier, l_mut_total))
p + geom_violin(scale = "width", aes(fill = Neutralization.breadth),alpha = 0.8)+ geom_jitter(height = 0, width = 0.1, size = 0.1 , 
                                                                                              color = "gray32")+ 
  geom_boxplot(outlier.shape = NA, width = 0.2, aes(fill = Neutralization.breadth), alpha = 0.8,size = 0.25, show.legend = FALSE) + 
  #facet_wrap(. ~ group,scales = "free_x") +
  scale_fill_manual(drop = FALSE,values= col.breadth,
                    name="Neutralization\nbreadth",
                    breaks=levels(bcr$Neutralization.breadth),
                    labels=c("0%","10%","20%","30%", "40%",
                             "50%","60%","70%","80%","90%","100%")) +
  theme_classic() + 
  scale_x_discrete(limits = rev(levels(bcr$donor_identifier)))+
  theme(axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(color = "black",size = 12, angle = 45, hjust = 1), 
        axis.text.y = element_text( size = 12),
        plot.title = element_text(size = 6, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.position="right",
        strip.text.x = element_text(size = 12),
        legend.justification="left",
        legend.text=element_text(size=12),
        legend.key.size = unit(0.4, "cm"))+
  #guides(fill=guide_legend(ncol=4)) +
  guides(fill = guide_legend(reverse=T))+
  xlab("Study Participant ID") + ylab("Number of mutations")
dev.off()
