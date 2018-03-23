gg.gauge <- function(pos,breaks=c(0,16.7,33.4,50.1,66.8,83.5,100)) {
  require(ggplot2)
  get.poly <- function(a,b,r1=0.5,r2=1.0) {
    th.start <- pi*(1-a/100)
    th.end   <- pi*(1-b/100)
    th       <- seq(th.start,th.end,length=100)
    x        <- c(r1*cos(th),rev(r2*cos(th)))
    y        <- c(r1*sin(th),rev(r2*sin(th)))
    return(data.frame(x,y))
  }
  if(pos <= 16.7){
      rating <- "VL"
    }
  else if(pos > 16.7 && pos <= 33.4){
      rating <- "L"
    }
  else if(pos > 33.4 && pos <= 50.1){
      rating <- "M"
    }
  else if(pos > 50.1 && pos <= 66.8){
      rating <- "SG"
    }
  else if(pos > 66.8 && pos <= 83.5){
      rating <- "H"
    }
  else if(pos > 83.5){
      rating <- "SV"
    }
  else {
      rating <- "NA"
    }
  ggplot()+ 
    geom_polygon(data=get.poly(breaks[1],breaks[2]),aes(x,y),fill="green4")+
    geom_polygon(data=get.poly(breaks[2],breaks[3]),aes(x,y),fill="green2")+
    geom_polygon(data=get.poly(breaks[3],breaks[4]),aes(x,y),fill="yellow")+
    geom_polygon(data=get.poly(breaks[4],breaks[5]),aes(x,y),fill="orange")+
    geom_polygon(data=get.poly(breaks[5],breaks[6]),aes(x,y),fill="red2")+
    geom_polygon(data=get.poly(breaks[6],breaks[7]),aes(x,y),fill="red4")+
    geom_polygon(data=get.poly(pos-1,pos+1,0.2),aes(x,y))+
    annotate("text",x=0,y=0,label=rating,vjust=0,size=5,fontface="bold")+
    coord_fixed()+
    theme_bw()+
    theme(axis.text=element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank()) 
}