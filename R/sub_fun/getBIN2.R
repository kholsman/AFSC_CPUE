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
       
      mid <- mean(c(up,dwn),na.rm=T)
      range <- paste0("L>=",dwntxt,"<",uptxt)
      if(type==1)
        return(mid=mid/divid)
      if(type==2)
        return(range)
      
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