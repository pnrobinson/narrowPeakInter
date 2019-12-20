library(ggplot2)
library(ggpubr)





dat <- read.table("../h3k27ac-enhancer.txt")
colnames(dat) <- c('group', 'H3K27ac')


g <- ggplot(dat, aes(x = group, y = H3K27ac)) +  geom_boxplot(notch = TRUE,aes(fill=group))
g <- g + theme_classic()
g <- g + scale_x_discrete(c("CGI", "Non-CGI"))
#g <- g + scale_fill_brewer(palette = "Dark2")
g <- g + scale_y_continuous(trans='log10')
g <- g + stat_compare_means(label.x=1.3, label.y=4, size=6, legend = "none")
g <- g + theme(axis.title.x = element_blank(),
               axis.title.y = element_text(size=24),
               axis.text.x = element_text(size=28),
               axis.text.y = element_text(size=22),
               legend.position = "none")

ggsave("h3k27ac-enhancer.pdf")


dat2 <-  read.table("../h3k27ac-promoter.txt")
colnames(dat2) <- c('group', 'H3K27ac')


g <- ggplot(dat2, aes(x = group, y = H3K27ac)) +  geom_boxplot(notch = TRUE,aes(fill=group))
g <- g + theme_classic()
g <- g + scale_x_discrete(c("CGI", "Non-CGI"))
#g <- g + scale_fill_brewer(palette = "Dark2")
g <- g + scale_y_continuous(trans='log10')
g <- g + stat_compare_means(label.x=1.3, label.y=4, size=6, legend = "none")
g <- g + theme(axis.title.x = element_blank(),
               axis.title.y = element_text(size=24),
               axis.text.x = element_text(size=28),
               axis.text.y = element_text(size=22),
                    legend.position = "none")


ggsave("h3k27ac-promoter.pdf")

