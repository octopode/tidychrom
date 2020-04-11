# utility function to convert JRW filename to sample code
filename2samp <- function(filename){
  prefix <- substr(filename, 0, 3)
  samp <- str_split(filename, "_") %>%
    .[[1]] %>% .[1] %>%
    substr(4, str_length(.)) %>%
    str_pad(., 4, "0", side = "left") %>%
    paste(prefix, ., sep="")
  return(samp)
}
