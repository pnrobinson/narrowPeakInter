library(ggplot2)
library(ggpubr)

dat <- read.table("../h3k27ac-enhancers.txt")
colnames(dat) <- c('group', 'H3K27ac')


g <- ggplot(dat, aes(x = group, y = H3K27ac)) +  geom_boxplot()
g <- g + theme_classic()
#g <- g + scale_fill_brewer(palette = "Dark2")
g <- g + scale_y_continuous(trans='log10')
g <- g + stat_compare_means()
g <- g + theme(axis.title.x = element_blank(),
       	 axis.title.y = element_text(size=24),
       	 		    axis.text = element_text(size=18))
g <- g + scale_fill_brewer(palette = "Dark2")

ggsave("h3k27ac.pdf")