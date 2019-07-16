message("which class you want to visualize for barplot?")
pick.class <- as.character(readline("Please input the class for barplot visualization: "))


lapply(sample.list, function(x) individuleClassPlot(x, filtered.lipidomics))


individuleClassPlot <- function(pick.class, sample.list, filtered.lipidomics){
filtered.class <- filtered.lipidomics %>% 
  select(Class, LipidMolec, sample.list) %>% 
  filter(grepl(pick.class, Class))

class.data <- aggregate(.~Class+LipidMolec, data=filtered.class, FUN=sum) %>% select(-Class)
data.wide <- data.frame(t(class.data[,-1]), group.repeats)
molec.names <- class.data[,1] %>% unlist() %>% str_extract_all(., "\\(\\d*.*\\)") %>% str_remove_all(., "[\\(\\)]") 

names(data.wide) <- c(molec.names, "experiment.group")


class.sd <- aggregate(.~experiment.group, data=data.wide, function(x) sd(x)) 
class.sd <- class.sd %>% gather( Aceyl, sd, -experiment.group)


class.long <- data.wide %>% group_by(experiment.group) %>% summarise_all(funs(sum)) 

class.data.long <- class.long %>% gather(Aceyl, AreaValue, -experiment.group)

class.data.long <- bind_cols(class.data.long, sd=class.sd$sd)



ggplot(class.data.long, aes(x=Aceyl, fill=experiment.group, y=AreaValue)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=AreaValue, ymax=AreaValue+sd), width=.8,
                position=position_dodge(.9)) +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),      ###### base of 10
                labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_continuous(labels = scales::scientific)+                      ###### scientific notation
  labs(ylab="Main Area", fill="Groups", xlab="AT") +
  ggtitle(pick.class) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Aceyl composition") +
  coord_cartesian(ylim=c(1e5,1e10))


message("Please input the name you want to store for the graph. e.g. TG.pdf")
class.type <- readline("The lipid class plot name: ")
ggsave(filename = class.type, path = 'plot/', device = "pdf")
}
