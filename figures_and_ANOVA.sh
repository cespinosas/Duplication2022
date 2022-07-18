#Essential code for the figures in the article. You will need to arrange the data as indicated.library(ggplot2)require(reshape2)library(ggExtra)#FIG 1#####################################################################################Fig 1a# from a file:# In "data1a[,1]" mutational robustness of GRNs before adding a duplicate gene.# In "data1a[,2]" mutational robustness of GRNs after adding a duplicate gene.# In "data1a[,3]" you can consider 1 if the GRN preserves the original phenotype after adding a duplicate gene, and 0 if don't.data1a <- read.csv(file = "filename.txt")setEPS()postscript("/Fig1a.eps", width = 2.3, height = 2.3, colormodel = "rgb", paper = "special")p <- ggplot(data1a, aes(x = data1a[,1], y = data1a[,2], color = data1a[,3])) + geom_point(size=0.3) + labs(x=expression(~italic("R\u00b5")*" before duplication"), y= expression(~italic("R\u00b5")*" after duplication")) + theme(text = element_text(size=8)) + geom_abline(intercept = 0, size = 0.2) + theme(legend.position = "none")ggMarginal(p, type = "histogram")dev.off()#Fig 1b# In the file you need to order similarity in:# In data1b[,1] the average similarity of all the GRNs before adding a gene.# In data1b[,2] the average similarity of the duplication-resistant GRNs before duplication.# In data1b[,3] the average similarity of the duplication-resistant GRNs after duplication.# In data1b[,2] the average similarity of the addition-resistant GRNs before adding a gene.# In data1b[,3] the average similarity of the addition-resistant GRNs after adding a gene.data1b <- read.csv(file = "filename.txt")setEPS()postscript("/Fig1b.eps", width = 2.3, height = 2.25, colormodel = "rgb", paper = "special")ggplot(melt(data1b), aes(variable, value)) + geom_boxplot(fill="gray", outlier.size = 0.1, size = 0.2) + labs(x=" ", y=expression(italic("S\u00b5"))) + stat_summary(fun=mean, geom="point", shape=20, size=1, color="black", fill="black") + scale_x_discrete(labels=c("Esp" = "Original\nnetworks", "BeforeDup" = "Before", "AfterDup" = "After", "BeforeCI" = "Before", "AfterCI" = "After")) + theme(text = element_text(size=8))dev.off()#####################################################################################FIG 2#Fig 2a# In the file you need to order the data in 3 columns:# In data2a[,1] goes the "Regulator" of the interaction; which can be D for duplicate and N for non-duplicate.# In data2a[,2] goes the "Target" of the interaction; which can be D for duplicate, N for non-duplicate and S for structural gene.# In data2a[,3] goes the "Value" of the similarity after delete that interaction.data2a <- read.csv(file ="filename.txt")setEPS()postscript("/Fig2a.eps", width = 2, height = 2.3, colormodel = "rgb", paper = "special")SoilSciGuylabs <- c("D", "N", "S")ggplot(data2a, aes(x=Target, y=Value, fill=Regulator, color=Regulator)) + geom_boxplot(color = "black", outlier.size = 0.1, size = 0.2)  + stat_summary(fun="mean", geom="point", size=1, position=position_dodge(width=0.75), color="black") + labs(x="Target", y=expression(~italic("S\u00b5")*" (deleting interactions)"))+ xlab("Target") + scale_x_discrete(labels= SoilSciGuylabs) + theme(legend.position = "none") +  theme(text = element_text(size=8))dev.off()# In the file you need to order the data in 3 columns:# In data2b[,1] goes the "Regulator" of the interaction; which can be D for duplicate and N for non-duplicate.# In data2b[,2] goes the "Target" of the interaction; which can be D for duplicate, N for non-duplicate and S for structural gene.# In data2b[,3] goes the "Value" of the average similarity after adding that interaction.data2b <- read.csv(file ="filename.txt")setEPS()postscript("/Fig2b.eps", width = 2.7, height = 2.3, colormodel = "rgb", paper = "special")SoilSciGuylabs <- c("D", "N", "S")ggplot(data2b, aes(x=Target, y=Value, fill=Regulator, color=Regulator)) + geom_boxplot(color = "black", outlier.size = 0.1, size = 0.2) + stat_summary(fun="mean", geom="point", size=1, position=position_dodge(width=0.75), color="black") + labs(x="Target", y=expression(~italic("S\u00b5")*" (adding interactions)")) + scale_x_discrete(labels= SoilSciGuylabs) +  theme(text = element_text(size=8))dev.off()#####################################################################################FIG 3#Fig 3a# In data3a: data3a[,1] have the number of incoming interactions of the duplicate gene for duplication-resistant GRNs, and data3a[,2] the number of incoming interactions of the duplicate gene for duplication-susceptible GRNs.data3a <- read.csv(file ="filename.txt")setEPS()postscript("/Fig3a.eps", width = 2.3, height = 2.3, colormodel = "rgb", paper = "special")ggplot(melt(data3a), aes(variable, value)) + geom_boxplot(fill="gray", outlier.size = 0.1, size = 0.2) + xlab(" ") + ylab("Incoming interactions") + stat_summary(fun=mean, geom="point", shape=20, size=1, color="black", fill="black") + theme(text = element_text(size=8)) + scale_x_discrete(labels=c("Resistant" = "Duplication\nresistant", "Susceptible" = "Duplication\nsusceptible"))dev.off()# In data3b: data3b[,1] have the number of outgoing interactions of the duplicate gene for duplication-resistant GRNs, and data3b[,2] the number of outgoing interactions of the duplicate gene for duplication-susceptible GRNs.data3b <- read.csv(file ="filename.txt")setEPS()postscript("/Fig3b.eps", width = 2.3, height = 2.3, colormodel = "rgb", paper = "special")ggplot(melt(data3b), aes(variable, value)) + geom_boxplot(fill="gray", outlier.size = 0.1, size = 0.2) + xlab(" ") + ylab("Outgoing interactions") + stat_summary(fun=mean, geom="point", shape=20, size=1, color="black", fill="black") + theme(text = element_text(size=8)) + scale_x_discrete(labels=c("Resistant" = "Duplication\nresistant", "Susceptible" = "Duplication\nsusceptible"))dev.off()#####################################################################################FIG 4# In data4 must be 2 columns: in data4[,1] with the same similarity than in data1b[,1], and in data4[,4] the number of accessible phenotypes before adding a gene.data4 <- read.csv(file ="filename.txt")setEPS()postscript("/Fig4.eps", width = 2.3, height = 2.3, colormodel = "rgb", paper = "special")p <- ggplot(a, aes(x = a[,1], y = a[,2], color = "gray")) + geom_point(size=0.3) + theme(legend.position = "none") + geom_abline(intercept = 0, size = 0.2) + labs(x=expression("GNRs' "~italic("S\u00b5")), y="New phenotypes after mutations") + theme(text = element_text(size=8))ggMarginal(p, type = "histogram")dev.off()#####################################################################################FIG 5#Fig 5alibrary(eulerr)setEPS()postscript("/Users/yuridiaposadas/Documents/Doctorado/Resultados/paper/Figuras/Fig5a.eps", width = 2, height = 2, colormodel = "rgb", paper = "special")plot(euler(c("Before"=5,"After"=4,"Before&After"=2)),quantities = FALSE,fills=c("#66C2A5","#FC8D62", "#8DA0CB"))dev.off()#Fig 5b# The data in data5b must be in 4 columns:# In data5b[,1] the number ("Num") in the order in which they should appear. In this case, "Total" is 1, "Resistant" 2 and "Susceptible" 3.# In data5b[,2] the colum "Via" order the data as "Total" (that are resistant and susceptible GRNs), "Resistant" (if they are duplication-resistant GRNs) and "Susceptible" (if they are duplication-susceptible GRNs). So, the data is doubled in this file.# In data5b[,3] is the kind of accessible phenotype: "Lost", "Recurrent" or "Newcomer".# In data5b[,3] is the number of accessible phenotypes.data5b <- read.csv(file ="filename.txt")setEPS()postscript("/Users/yuridiaposadas/Documents/Doctorado/Resultados/paper/Figuras/Fig5b.eps", width = 2.8, height = 2.3, colormodel = "rgb", paper = "special")SoilSciGuylabs <- c("Total", "Resistant", "Susceptible")ggplot(data5b, aes(x=Num, y=Datos, fill=Accessible.phenotypes)) + geom_boxplot(outlier.size = 0.1, size = 0.2) + stat_summary(fun="mean", geom="point", size=1, position=position_dodge(width=0.75), color="black") + aes(stringr::str_wrap(Num, 3), Datos) + ylab("Accessible phenotypes") + xlab(" ") + scale_x_discrete(labels= SoilSciGuylabs) + labs(fill='Accessible\nphenotypes') + scale_fill_brewer(palette="Set2") + theme(text = element_text(size=8))dev.off()#Fig 5c# In data5c must be in 4 columns:# In data5c[,1] the number ("Num") in the order in which they should appear. In this case, "Total" is 1, "Resistant" 2 and "Susceptible" 3.# In data5c[,2] the colum "Via" order the data as "Total" (that are resistant and susceptible GRNs), "Resistant" (if they are duplication-resistant GRNs) and "Susceptible" (if they are duplication-susceptible GRNs). So, the data is doubled in this file.# In data5c[,3] is the kind of accessible phenotype: "Lost" or "Recurrent".# In data5c[,3] is the number of mutation leading to the accessible phenotypes.data5c <- read.csv(file ="filename.txt")setEPS()postscript("/Users/yuridiaposadas/Documents/Doctorado/Resultados/paper/Figuras/Fig5c.eps", width = 2.3, height = 2.3, colormodel = "rgb", paper = "special")SoilSciGuylabs <- c("Total", "Resistant", "Susceptible")ggplot(a, aes(x=Num, y=Datos, fill=Accessible.phenotypes)) + geom_boxplot(outlier.size = 0.1, size = 0.2) + stat_summary(fun="mean", geom="point", size=1, position=position_dodge(width=0.75), color="black") + aes(stringr::str_wrap(Num, 3), Datos) + ylab("Mutations leading to accessible\nphenotypes per phenotype") + xlab(" ") + scale_x_discrete(labels= SoilSciGuylabs) + labs(fill='Accessible\nphenotypes') + scale_fill_manual(values=c("#66C2A5", "#8DA0CB")) + theme(text = element_text(size=8)) + theme(legend.position = "none")dev.off()#Fig 5d# In data5d must be in 3 columns:# In data5d[,1] the average similarity between recurrently accessible phenotypes and the original phenotype.# In data5d[,2] the average similarity between lost accessible phenotypes and the original phenotype.# In data5d[,3] the average similarity between newcomer accessible phenotypes and the original phenotype.data5d <- read.csv(file ="filename.txt")setEPS()postscript("/Users/yuridiaposadas/Documents/Doctorado/Resultados/paper/Figuras/Fig5d.eps", width = 2.3, height = 2.3, colormodel = "rgb", paper = "special")SoilSciGuylabs <- c("Recurrent", "Lost", "Newcomer")ggplot(melt(df), aes(variable, value)) + geom_boxplot(fill=c("#8DA0CB", "#66C2A5", "#FC8D62"), outlier.size = 0.1, size = 0.2) + xlab(" ") + ylab("Phenotypes' S\u00b5") + stat_summary(fun=mean, geom="point", shape=20, size=1, color="black", fill="black") + theme(text = element_text(size=8)) + ylim(0, 1) + xlab("Accessible phenotypes") + scale_x_discrete(labels= SoilSciGuylabs)dev.off()########################################################Type II ANOVA######################################################library(car)#Type II ANOVA for deletion#From Fig 2a:data2a <- read.csv(file ="filename.txt")time.lm <- lm(formula = Value ~ Regulator * Target, data = data2a)time.II.aov <- car::Anova(time.lm, type = 2)#Effect size of independent on dependentlibrary(lsr)etaSquared(time.lm, type = 2, anova = TRUE)#Type II ANOVA for addition#From Fig 2b:data2b <- read.csv(file ="filename.txt")time.lm <- lm(formula = Value ~ Regulator * Target, data = data2a)time.II.aov <- car::Anova(time.lm, type = 2)#Effect size of independent on dependentetaSquared(time.lm, type = 2, anova = TRUE)