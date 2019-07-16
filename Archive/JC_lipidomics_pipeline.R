# Read in file
lipidomics <- read.csv(file.choose("~/Deskstop/documents/Felipe/Lipidomics/20171127_JC_SeipinKO/20171127_JC_SeipinKO.csv"), header = T)

# Count A, B, C and D's
lipidomics$A <- apply(lipidomics[,30:47], 1, function(x) length(which(x=="A")))
lipidomics$B <- apply(lipidomics[,30:47], 1, function(x) length(which(x=="B")))
lipidomics$C <- apply(lipidomics[,30:47], 1, function(x) length(which(x=="C")))
lipidomics$D <- apply(lipidomics[,30:47], 1, function(x) length(which(x=="D")))
lipidomics$No_grade <- apply(lipidomics[,30:47], 1, function(x) length(which(x=="-")))

# Count APValue<=0.001
lipidomics$APvalue_001 <- apply(lipidomics[,11:28], 1, function(x) length(which(x<=0.001)))

# Filter the dataset based on your criteria
filtered_lipidomics <- subset(lipidomics, 
                              Rej == 0 & 
                                (A >= 3 | B >= 3) &
                                (APvalue_001 >=3))

# Check how many lipids per class that passed the filtering
summary(filtered_lipidomics$Class)

#bb <- filtered_lipidomics %>% rowwise %>% filter(Class=="PS")

# QC PLOT 1 - Abundance vs. retention time [requires ggplot2]
ggplot(data=filtered_lipidomics,aes(x = log10(MainArea.s1.), y = BaseRt)) +
  geom_point() +
  theme_bw() +
  facet_grid(.~Class) + ggtitle("there is his")


ggplot(data=filtered_lipidomics,aes(x = log10(), y = BaseRt)) +
  geom_point() +
  theme_bw() +
  facet_grid(.~Class) + ggtitle("there is his")

# Preparation for pair-wise correlations
inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(inf2NA(x), inf2NA(y), use = "pairwise.complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(is.na(p<0.01)) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

# QC PLOT 2 - Pair-wise correlation between replicates
pairs(log10(filtered_lipidomics[,c(49, 59:60)]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor)
pairs(log10(filtered_lipidomics[,61:63]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor)
pairs(log10(filtered_lipidomics[,64:66]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor)
pairs(log10(filtered_lipidomics[,50:52]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor)
pairs(log10(filtered_lipidomics[,53:55]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor)
pairs(log10(filtered_lipidomics[,56:58]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor)

# Formatting the table for PCA
filtered_lipidomics_PCA <- as.data.frame(t(subset(filtered_lipidomics[,49:66])))
colnames(filtered_lipidomics_PCA) <- filtered_lipidomics$LipidMolec
filtered_lipidomics_PCA$Group <- c("Wt", "Wt_24h", "Wt_24h", "Wt_24h", "Wt_48h", "Wt_48h", "Wt_48h", "Wt_A", "Wt_A", "Wt_A", "Wt", "Wt", "D", "D", "D", "Q", "Q", "Q")

# Perform PCA [requires FactoMineR]
res.pca = PCA(filtered_lipidomics_PCA, scale.unit=TRUE, ncp=5, quali.sup=629, graph=T)
concat = cbind.data.frame(filtered_lipidomics_PCA[,629],res.pca$ind$coord)
ellipse.coord = coord.ellipse(concat,bary=T)
plot.PCA(res.pca,habillage=629,ellipse=ellipse.coord,cex=0.8, label="all")

# Output the filtered lipidome
write.csv(filtered_lipidomics, "filtered_lipidomics.csv")
