readinteger <- function()
{ 
  n <- readline(prompt="Continue plotting: ")
  if(!grepl("^[0-9]+$",n))
  {
    return(readinteger())
  }
  
  return(as.integer(n))
}
