plotBars <- function(x, 
                     quantile=c(0.0275, 0.975),
                     at = at,
                     widthBar = 0.1,
                     quantile1=c(0.25, 0.75),
                     col="red", 
                     col1=NULL, 
                     alpha= 0.5,
                     alpha1= 0.8,
                     horiz=FALSE,
                     borderCol=NA,
                     lwd=1
                     ){

quantile95 <- quantile(x, prob=c(0.0275, 0.975))
if(is.null(col1)){
  col1 <-col
}
if(horiz==FALSE){
polygon(x = c(at- widthBar, at+ widthBar,
              at+ widthBar, at- widthBar ),
        y = c(quantile95[1], quantile95[1],
              quantile95[2], quantile95[2]), 
        col=adjustcolor(col, alpha),
        border= NA)
if(!is.null(quantile1)){
quantile50 <- quantile(x, prob=c(0.25, 0.75))

polygon(x = c(at-widthBar, at+widthBar,
              at+widthBar, at-widthBar ),
        y = c(quantile50[1], quantile50[1],
              quantile50[2], quantile50[2]), 
        col=adjustcolor(col1, alpha1),
        border= NA)

}
  
  if(!is.na(borderCol)){
    polygon(x = c(at- widthBar, at+ widthBar,
                  at+ widthBar, at- widthBar ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=NA,
            border= borderCol,lwd=lwd)
  }  
  
}else{
  polygon(x =c(quantile95[1], quantile95[1],
               quantile95[2], quantile95[2]) ,
          y = c(at- widthBar, at+ widthBar,
            at+ widthBar, at- widthBar ), 
          col=adjustcolor(col, alpha),
          border= NA)
  if(!is.null(quantile1)){
    quantile50 <- quantile(x, prob=c(0.25, 0.75))
    
    polygon(x = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]),
            y = c(at-widthBar, at+widthBar,
                  at+widthBar, at-widthBar ), 
            col=adjustcolor(col1, alpha1),
            border= NA)
    
  }
  
  if(!is.na(borderCol)){
    polygon(x =c(quantile95[1], quantile95[1],
                 quantile95[2], quantile95[2]) ,
            y = c(at- widthBar, at+ widthBar,
                  at+ widthBar, at- widthBar ), 
            col=NA,
            border= borderCol,lwd=lwd)
  }
}

}