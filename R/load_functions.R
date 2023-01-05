#load_functions.R


for(d in dir("R/sub_fun")) 
  source(paste0("R/sub_fun/",d))  # source(file.path(code.path,"R/AGG_PREY_FUN.R"))

# # load general functions
# source("R/funKir.R") 

# Minor functions:
  first<-function(x,na.rm=TRUE){
    if(na.rm==TRUE) x<-x[-which(is.na(x))]
    return(x[1])
  }

# used in get_prey_data() function:
  
  expandMat  <- function(dat=bywt2,prey_listIN=prey_list){
    bigmat   <-  data.frame(matrix(0,dim(dat)[1],length(prey_listIN)+1))
    colnames(bigmat)  <- c("ID",prey_listIN)
    for(i in 1:(dim(dat)[2]))
      bigmat[,which(colnames(bigmat)==colnames(dat)[i])]<-dat[,i]
    return(bigmat)
  }

  # used in makePropB:
  #-------------------------------
  getBIN  <- function(x, BigBININ=BigBIN){
    sapply(x,function(x,bins= BigBININ){
      up  <- BigBININ[x<BigBININ][1]
      if(is.na(up)) 
        up <- BigBININ[x<=BigBININ][1]
      dwn <- rev(BigBININ[x>BigBININ])[1]
      if(is.na(dwn)) 
        dwn <- rev(BigBININ[x>=BigBININ])[1]
      mid <- (up+dwn)/2
      return(mid=mid/10)
    })
  }
  
  replace_space <- function(x){
    sapply(prey, function(x){
      tmp <- unlist(strsplit(x,split=""))
      if(any(tmp==" "))
        tmp[tmp==" "]<-"."
      return(paste0(tmp,collapse=""))
    })}
  
  list_match <- function (x,listIN=splist){ 
    sapply(x, function(x,listin=listIN){
      nout <- NA
      for(i in 1:length(listin))
        if(is.element(x,listin[[i]]))
          nout <- i
        return(nout)
        
    })}
  
  # # used to preview shapefiles:
  # preview<-function(shp,epsg=3995){
  #   ggplot() + 
  #     theme_light()+
  #     geom_sf(data=st_transform(shp,crs=paste0("+init=epsg:",epsg)))
  # }
  # 
  ####################################################################
  ## SQL R Functions
  ####################################################################
  length.na     <-function(x){
    sub<-x[is.na(x)==FALSE]
    return(length(sub))
  }
  saveif        <-function(nm,race_only1=race_only){
    
    if(race_only1==1){
      txt<-file.path(path,paste(nm,"_noObs.Rdata",sep=""))
      eval(parse(text=paste("save(",nm,",file=txt)",sep="")))
    }else{
      txt<-file.path(path,paste(nm,".Rdata",sep=""))
      eval(parse(text=paste("save(",nm,",file=txt)",sep="")))
    }
  }
  save2DB       <-function(nm){
    tables2<-sqlTables(con,schema="HOLSMANK")$TABLE_NAME
    if(any(tables2==nm)) sqlDrop(con, nm)
    saveAndSub("LWA_ALL",path1=file.path(data.path, subDir))
    eval(parse(text=paste("sqlSave(con,", nm,",safer=FALSE)",sep="")))
  }
  summarizeTots <-function(predpreyfile=predprey_anne){
    yrs<-sort(unique(predpreyfile$YEAR))
    preds<-unique(predpreyfile$GOAPOLL_PRED)
    smry<-data.frame(matrix(0,length(yrs),length(preds)))
    colnames(smry)<-preds;rownames(smry)<-yrs
    
    for(y in 1:length(yrs)){
      
      for(p in 1:length(preds)){
        sub<-predpreyfile[predpreyfile$YEAR==yrs[y]&predpreyfile$GOAPOLL_PRED==preds[p]&predpreyfile$TWT>0,]
        smry[y,p]<-length(unique(sub$PPID))
        
      }
    }
    return(smry)
  }
  
  
  