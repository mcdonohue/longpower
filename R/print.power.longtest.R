print.power.longtest <-
function(x, ...)
{
    cat("\n    ", x$method, "\n\n")
    note <- x$note
    R <- x$R
    x[c("method","note","R")] <- NULL
    cat(paste(format(names(x), width= 15, justify = "right"),
	      format(x), sep= " = "), sep= "\n")
    if(!is.null(note))
	cat("\n", "NOTE:", note, "\n")
    else
	cat("\n")
    if(!is.null(R)){
	  cat("\n", "R:\n")
	  print(R)
	  cat("\n")
	}
    else
	cat("\n")
    invisible(x)
}