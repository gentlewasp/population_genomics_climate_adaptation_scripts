library(tidyverse)
chosen_start=7915976
chosen_end=8390282
between_ld <- read_tsv("tp.thin100.maf5.Chr15.tp_chr5_within.windows.ld.gz")%>% filter(win1 != win2)
within_ld <- read_tsv("tp.thin100.maf5.Chr15.tp_chr15_within.windows.ld.gz") %>%
             rename(win1 = win2, win2 = win1) %>%
             filter(win1 != win2)
ld_plot<-  rbind(between_ld,within_ld) %>% 
      filter(n > 5)
str(ld_plot)
max_ld <- max(ld_plot[,7])
max_ld
data <- ld_plot %>%
      ggplot(.,aes()) + 
      geom_tile(aes(x=win1/1000000,y=win2/1000000,fill=max_2_r2)) +
      annotate("segment",y=0,x=chosen_start/1000000,
               yend=0,xend=chosen_end/1000000,
               color="#9B2374",size=10) +
      annotate("segment",y=chosen_start/1000000,x=0,
               yend=chosen_end/1000000,xend=0,
               color="#9B2374",size=10) +
      scale_fill_viridis_c(limits=c(0,max_ld),name="LD") +
      theme_linedraw() +ylab("Mbp") + xlab("Mbp") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            panel.border = element_blank(),
            legend.position = "bottom",
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            plot.margin = margin(40, 40, 40, 40, "pt")) +
      scale_x_continuous(limits = c(0,max(between_ld$win2)/1000000), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0,max(between_ld$win2)/1000000), expand = c(0, 0)) +
      coord_cartesian(clip = 'on') +
      theme(plot.tag.position = c(0.1, 0.5),
            plot.tag=element_text(size=2,color="red"))
pdf("tp_chr15_win50k_max_2_r2.pdf",width=16,height=16) 
plot(data)
dev.off()

