getBIN2<-
  function(x, BIN_IN=c(0,10,100),divid=1,type=1){
   
    sapply(x,function(x,BIN= BIN_IN){
      
      up  <- uptxt<- BIN[x<BIN][1]
      if(is.na(up)){
        up    <- rev(BIN[x>=BIN])[1] + mean(BIN[-1]-BIN[1:(length(BIN)-1)])
        uptxt <- paste0("Plus","+")
      }
        
      dwn <- dwntxt <- rev(BIN[x>=BIN])[1]
      if(is.na(dwn)){
        dwn <- BIN[x<BIN][1]- mean(BIN[-1]-BIN[1:(length(BIN)-1)])
        dwntxt <- paste0((BIN)[1],"-")
      }
       
      mid   <- mean(c(up,dwn),na.rm=T)
      range <- paste0("[",dwntxt,",",uptxt,")")
      range <- paste0("L>=",dwntxt,"<",uptxt)
      if(type==1)
        return(mid=mid/divid)
      if(type==2)
        return(range)
      if(type==3)
        return(list(mid=mid/divid,range=range))
    })
  }
# used in makePropB:
#-------------------------------
getBIN<-function(x, BigBININ=BigBIN,divideby = 10){
  sapply(x,function(x,bins= BigBININ){
    up  <- BigBININ[x<BigBININ][1]
    if(is.na(up)) 
      up <- BigBININ[x<=BigBININ][1]
    dwn <- rev(BigBININ[x>BigBININ])[1]
    if(is.na(dwn)) 
      dwn <- rev(BigBININ[x>=BigBININ])[1]
    mid <- (up+dwn)/2
    return(mid=mid/divideby)
  })
}
# 
# assign_byBIN<-
#   function(x, BIN_IN=c(390,400),divid=1,type=1){
# 
#     bin_out <- x*0
#     
#     bin_out[x>=BIN_IN[1]&x<BIN_IN[2]]<- BIN_IN[1]
#     
#   
#     
#       up  <- uptxt<- BIN[x<BIN][1]
#       if(is.na(up)){
#         up    <- rev(BIN[x>=BIN])[1] + mean(BIN[-1]-BIN[1:(length(BIN)-1)])
#         uptxt <- paste0("Plus","+")
#       }
#       
#       dwn <- dwntxt <- rev(BIN[x>=BIN])[1]
#       if(is.na(dwn)){
#         dwn <- BIN[x<BIN][1]- mean(BIN[-1]-BIN[1:(length(BIN)-1)])
#         dwntxt <- paste0((BIN)[1],"-")
#       }
#       
#       mid <- mean(c(up,dwn),na.rm=T)
#       range <- paste0("L>=",dwntxt,"<",uptxt)
#       if(type==1)
#         return(mid=mid/divid)
#       if(type==2)
#         return(range)
#       
#     }
#   