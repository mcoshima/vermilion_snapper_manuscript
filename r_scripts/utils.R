## Split recombine function

splt.recombine <- function(df1, df2, group_var, N){
  
  group_var <- enquo(group_var)
  col <- quo_name(group_var)
  
  names <- df1 %>%
    select(!!group_var) %>%
    unique() %>%
    pull() %>%
    sort()
  
  split.dat <- df1 %>%
    group_by(!!group_var) %>%
    group_split() %>%
    setNames(names)
  
  names.2 <- df2 %>%
    select(!!group_var) %>%
    unique() %>%
    pull() %>%
    sort()
  
  split.new <- df2 %>%
    #mutate_at(vars(!!col), function(x) x = as.numeric(paste(x))) %>%
    group_by(!!group_var) %>%
    group_split() %>%
    setNames(names.2)
  
  for(i in 1:N){
    
    x <- which(names(split.dat) == names(split.new)[i])
    
    split.dat[[x]] <- rbind((split.dat[x][[1]]),
                            (split.new[i][[1]]))
    
    colnames(split.dat[[x]]) <- colnames(df1)
    
  }
  
  return(do.call(rbind, split.dat))
  
}

### Function for running SS and specifying the file name based on the operating system
### for os specify "windows" for windows, "osx" for Mac, or "lin" for Linux
### for xtras you can add additional commands, such as -nohess or -mceval
run_ss <- function(dir., os = "osx", xtras=NULL){
  out <- tryCatch(
    expr = {
      
      message("Running SS now")
      if(os == "lin"){
        system(paste("cd", dir., "&& ./ss_linux", xtras, "> /dev/null 2>&1", sep = " "))
      }
      if(os == "windows"){
        shell(paste("cd/d", dir., "&& ss_win", xtras, " >NUL 2>&1", sep = " "))
        }
     if(os == "osx"){
       system(paste("cd", dir., "&& ./ss_osx", xtras, " > /dev/null 2>&1", sep = " "))
     }
      
      assign("error", FALSE, envir = globalenv())
      
    },
    warning = function(w){
      print(paste("There was a warning", w))
      assign("error", TRUE, envir = globalenv())
      return(error)
    },
    error = function(e){
      
      print(paste("The model didn't converge", e))
      assign("error", TRUE, envir = globalenv())
      return(error)
    },
    finally = {
      return(error)
    }
    
  )
  return(out)
}
