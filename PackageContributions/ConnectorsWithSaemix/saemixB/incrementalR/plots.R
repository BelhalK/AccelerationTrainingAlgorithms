seplot <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  graf <- ggplot(df)+geom_line(aes(iterations,value,by=value,colour = df$algo),show.legend = legend) +
  xlab("iterations") + ylab('value') + facet_wrap(~variable,scales = "free_y") + theme_bw() 
  grid.arrange(graf)
  # do.call("grid.arrange", c(graf, ncol=1, top=title))
}
