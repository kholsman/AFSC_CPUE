#' Kirstin Holsman
#' kirstin.holsman@noaa.gov
#' July 2020
#' based on code from Steve Barbeaux, Dan Nichol

  for(l in 1:dim(srvys)[1]){
    for(s in 1:dim(species_lkup)[1]){
      cat("_____________________________________________");cat("\n")
      cat(paste("running query for",species_lkup$sp[s],"\n"))
      tt  <- Get_DATA(
        username,
        password,
        species = species_lkup$num[s],
        survey  = srv$num[l],
        bins    = sp_bins[[s]],
        LW_a    = LW_a_in[s],
        LW_b    = LW_b_in[s],
        plotit  = F)
      
      tmptx   <- paste0(species_lkup$sp[s],".cpue.",srvys$reg[l])
      eval(parse(text=paste0(tmptx, " <- tt")))
    
      fltxt   <- file.path(path,paste(tmptx,".Rdata",sep=""))
      eval(parse(text=paste("save(",tmptx,",file=fltxt)",sep="")))
      
      if(plotit){
        pdf(file=paste0(fltxt,"_Plot.pdf"))
        plot_dist_num(data= tt, Lbin=c(100,400), plotT=1,bynum=F)
        dev.off()
      }
      rm(list=list(fltxt,tmptx,tt))
    }
  }


