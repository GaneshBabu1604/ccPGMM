library(ggplot2)
library(agricolae)
library(patchwork)
library(reshape2) 
#Loading the simulation results
#Scenario 1
Clean_Data = read.csv('wellseparatedcluster_consolidated_results.csv')
#Scenario 2
Clean_Data = read.csv('heavyoverlapcluster_consolidated_results.csv')
#Scenario 3
Clean_Data = read.csv('mildoverlapcluster_consolidated_results.csv')
#Scenario 4
Clean_Data = read.csv('syntheticcereal_consolidated_results.csv')


#Adjusted Rand index plot
Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='d = 10',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y1 = ggplot(Clean_Data_10_0, aes(x = M, y = ARI, color = Model)) +
  geom_boxplot() +
  facet_wrap(~d) + theme_bw()+ ylim(-0.1,1)+ theme_bw()+ xlab('M')+ylab('Adjusted Rand index')+theme(axis.text.y = element_blank(), 
                                                                                                     axis.ticks.y = element_blank(), 
                                                                                                     axis.title.y = element_blank())
Y1
Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='d = 20',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y2 = ggplot(Clean_Data_10_0, aes(x = M, y = ARI, color = Model)) +
  geom_boxplot() +
  facet_wrap(~d) + theme_bw()+ ylim(-0.1,1)+ theme_bw()+ xlab('M')+ylab('Adjusted Rand index')+theme(axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                                                                                                     axis.title.y = element_blank())
Y2
Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='PGMM',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y3 = ggplot(Clean_Data_10_0, aes(x = M, y = ARI)) + ylab('Adjusted Rand index') +
  geom_boxplot(color = 'blue')  + ylim(-0.1,1) + theme_bw()+ theme(axis.text.x = element_blank(), 
                                                                   axis.ticks.x = element_blank(), 
                                                                   axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                                                                   axis.title.y = element_blank()
  )+ facet_wrap(~d) 
Y3
Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='DBSCAN',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y4 = ggplot(Clean_Data_10_0, aes(x = M, y = ARI)) +
  geom_boxplot(color = 'blue')  + ylim(-0.1,1) + theme_bw()+ theme(axis.text.x = element_blank(), 
                                                                   axis.ticks.x = element_blank(), 
                                                                   axis.title.x = element_blank()
  )+ facet_wrap(~d) +ylab('Adjusted Rand index')
Y4

Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='GMM',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y5 = ggplot(Clean_Data_10_0, aes(x = M, y = ARI)) +
  geom_boxplot(color = 'blue')  + ylim(-0.1,1) + theme_bw()+ theme(axis.text.x = element_blank(), 
                                                                   axis.ticks.x = element_blank(), 
                                                                   axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                                                                   axis.title.y = element_blank()
  )+ facet_wrap(~d) +ylab('Adjusted Rand index')
Y5
y3 =  Y4+ plot_spacer() + Y5+ plot_spacer()  +Y3 + plot_spacer() + Y1 + plot_spacer() + Y2 + plot_layout(widths = c(1.0,-0.3,1.0,-0.3,1.0,-0.3, 5,-0.3,5),guides = 'collect')
pdf(file = 'heavyoverlap_sim_equal.pdf',10,5)
y3
dev.off()


#Time Plot
Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='d = 10',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y1 = ggplot(Clean_Data_10_0, aes(x = M, y = Parallel.Time, color = Model)) +
  geom_boxplot() +
  facet_wrap(~d) + theme_bw()+ ylim(0,25)+ theme_bw()+ xlab('M')+ylab('Time (Minutes)')+theme(axis.text.y = element_blank(), 
                                                                                              axis.ticks.y = element_blank(), 
                                                                                              axis.title.y = element_blank())
Y1
Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='d = 20',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y2 = ggplot(Clean_Data_10_0, aes(x = M, y = Parallel.Time, color = Model)) +
  geom_boxplot() +
  facet_wrap(~d) + theme_bw()+ ylim(0,25)+ theme_bw()+ xlab('M')+ylab('Time (Minutes)')+theme(axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                                                                                              axis.title.y = element_blank())
Y2
Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='PGMM',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y3 = ggplot(Clean_Data_10_0, aes(x = M, y = Parallel.Time)) + ylab('Time (Minutes)') +
  geom_boxplot(color = 'blue')  + ylim(0,25) + theme_bw()+ theme(axis.text.x = element_blank(), 
                                                                 axis.ticks.x = element_blank(), 
                                                                 axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                                                                 axis.title.y = element_blank()
  )+ facet_wrap(~d) 
Y3
Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='DBSCAN',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y4 = ggplot(Clean_Data_10_0, aes(x = M, y = Parallel.Time)) +
  geom_boxplot(color = 'blue')  + ylim(0,25) + theme_bw()+ theme(axis.text.x = element_blank(), 
                                                                 axis.ticks.x = element_blank(), 
                                                                 axis.title.x = element_blank()
  )+ facet_wrap(~d) +ylab('Time (Minutes)')
Y4

Clean_Data_10_0 = Clean_Data[Clean_Data$Entropy == 0 & Clean_Data$d=='GMM',]
Clean_Data_10_0$R = factor(Clean_Data_10_0$r,levels = c(10,25,50,100))
Y5 = ggplot(Clean_Data_10_0, aes(x = M, y = Parallel.Time)) +
  geom_boxplot(color = 'blue')  + ylim(0,25) + theme_bw()+ theme(axis.text.x = element_blank(), 
                                                                 axis.ticks.x = element_blank(), 
                                                                 axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                                                                 axis.title.y = element_blank()
  )+ facet_wrap(~d) +ylab('Time (Minutes)')
Y5
y3 = Y4+ plot_spacer()+ Y5+ plot_spacer()  +Y3 + plot_spacer() + Y1 + plot_spacer() + Y2 + plot_layout(widths = c(1.0,-0.3,1.0,-0.3,1.0,-0.3, 5,-0.3,5),guides = 'collect')
pdf(file = 'heavyoverlap_sim_equal_time.pdf',10,5)
y3
dev.off()