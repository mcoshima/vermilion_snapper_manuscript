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

