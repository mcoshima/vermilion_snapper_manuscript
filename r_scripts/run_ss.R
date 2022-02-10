#' This runs SS and returns TRUE if there was an error that prevented the model from converging.
#' @param dir. the directory of where to run ss from
#' @param os the type of operating system ss is running on ("win" for windows, "osx" for mac, "lin" for linux)
#' @param xtra_args any additional arguments to specify for running ss (ie -nohess, -mcmc, etc.)
#' @keywords SS
#' @export
#'

run_ss <- function(dir., os = "win", xtra_args = NULL){
  out <- tryCatch(
    expr = {

      message("Running SS now")
      if(os == "osx"){
        system(paste("cd", dir., "&& ./SS3", xtra_args, "> /dev/null 2>&1", sep = " "))
        }
         if(os == "lin"){
        system(paste("cd", dir., "&& ./SS3", xtra_args, "> /dev/null 2>&1", sep = " "))
        }
      if(os == "win"){
          shell(paste("cd/d", dir., "&& ss", xtra_args, ">NUL 2>&1", sep = " "))
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