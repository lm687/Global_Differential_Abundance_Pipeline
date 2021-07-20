probs = c(0.2, 0.1, 0.1, 0.05, 0.45, 0.08, 0.02)
m1 = rmultinom(1, 300, probs)
m2 = rmultinom(1, 200, probs)

plot(unlist(outer(as.vector(m1), as.vector(m1), function(i,j) log(i/j) )),
     unlist(outer(as.vector(m2), as.vector(m2), function(i,j) log(i/j) )))

ggplot(cbind.data.frame(logR1 = as.vector(unlist(outer(as.vector(m1), as.vector(m1), function(i,j) log(i/j) ))),
                        logR2 = as.vector(unlist(outer(as.vector(m2), as.vector(m2), function(i,j) log(i/j) ))),
                        width=rep(m1/sum(m1), each=length(probs)),
                        height=rep(m1/sum(m1), length(probs))
                        ))+
  geom_point(aes(x=logR1, y=logR2))+
  geom_segment(aes(x=logR1-(width/2), xend=logR1+(width/2),
                   y=logR2-(height/2), yend=logR2+(height/2)))+theme_bw()

outer(as.vector(m1), as.vector(m1), function(i,j) log(i/j) )
outer(as.vector(m2), as.vector(m2), function(i,j) log(i/j) )