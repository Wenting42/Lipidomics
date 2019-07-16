volc = ggplot(input, aes(logFC, -log10(adj.P.Val))) + #volcanoplot with log2Foldchange versus pvalue
  geom_rect(aes(xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf),
            fill = "grey90")+
  geom_vline(xintercept=0, linetype="solid", colour="grey60", size=0.4) +
  #Add a horizontal line for P-value cut-off
  geom_hline(yintercept=-log10(0.05), linetype="solid", colour="grey60", size=0.4)+
  geom_hline(yintercept=-log10(0.01), linetype="solid", colour="grey60", size=0.4)+
  #scale_y_continuous(breaks=seq(,by=0.5))+  # Ticks from 0-10, every .25
  geom_point(aes(col=sig, size=AveExpr)) + #add points colored by significance, point size by average value
  scale_size_continuous(range = c(0.25,2.5))+ #add size limit
  scale_color_manual(values=c("dodgerblue3", "darkorange3", 'chartreuse4', 'grey45')) + 
  #ggtitle ('Total Lipid species') +
  #xlim(c(-2.7, 2.7)) + ylim(c(0, 2)) + #set the x,y axis limits.
  xlab("Log2 fold change (KO / WT)") + ylab("-Log10(q value)") + 
  theme_bw()+ # remove background
  theme(line=element_blank())+
  geom_text_repel(data=subset(input,adj.P.Val<0.05 & abs(mutateddf$logFC) > 5), aes(label=gene), size=2.5)#adding text for the FDR points.
#Add a vertical line for fold change cut-offs



volc 
