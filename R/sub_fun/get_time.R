#'
#'
#'
#'get_time.R
#'K Holsman
#'

get_time<- function(x, nm = "year",formatIN = "%Y-%m-%d %H:%M:%S"){
  a <- strptime(x, format = formatIN)
  if(is.null(nm)){
    out <- a
  }else{
    out <- eval(parse(text = paste0("a$",nm)))
  }
  return(out)
}
