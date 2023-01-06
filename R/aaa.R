ignored_prefixes <- paste("^", c(
  'WARN',
  "Loading",
  'loading', 'downloading', 'retrieving', 'require', # AnnotationHub messages
  'snapshot',
  'Attaching',
  'package',
  'The following object',
  'Welcome to Bioconductor'
), sep = '', collapse = '|')

suppressively <- function(f) {
  function(...) suppressPackageStartupMessages(f(...))
}


with_handlers <- function(code, .logger=NULL, .console=TRUE) {
  messages <- character()
  mHandler <- function(m) {
    m = trimws(m$message)
    if (grepl(ignored_prefixes, m) ||
        m == "" || grepl('deprecate', m, ignore.case = T)) {
      return(invokeRestart("muffleMessage"))
    }
    if (!is.null(.logger)) .logger(m)
    fmt = futile.logger::flog.layout()(futile.logger::INFO, m)
    if (.console) message(fmt)
    messages <<- c(messages, fmt)
    invokeRestart("muffleMessage")
  }
  warnings <- character()
  wHandler <- function(w) {
    w = trimws(w$message)
    if (!is.null(.logger)) .logger(w$message)
    # if (.console) futile.logger::flog.info(w$message)

    fmt = futile.logger::flog.layout()(futile.logger::WARN, w)
    if (.console) message(fmt)
    messages <<- c(messages, fmt)
    # warnings <<- c(warnings, w$message)
    invokeRestart("muffleWarning")
  }
  .add_trace <- function(e) {
    e$trace=capture.output(rlang::trace_back())
    fmt = futile.logger::flog.layout()(futile.logger::FATAL, e$message)
    messages <<- c(messages, fmt)
    fmt = futile.logger::flog.layout()(futile.logger::FATAL, paste(c("", e$trace), collapse = "\n"))
    messages <<- c(messages, fmt)
    signalCondition(e)
  }

  temp <- file()
  sink(temp)
  on.exit({
    sink()
    close(temp)
  })

  tryCatch({
    result <- withCallingHandlers(
      eval(code),
      warning = wHandler, message = mHandler, error = .add_trace,
      packageStartupMessage = function(c) tryInvokeRestart("muffleMessage")
    )
    output <- paste0(readLines(temp, warn = FALSE), collapse = "\n")
    list(result = result, output = output, warnings = warnings, messages = messages, error=NULL)
  },
    error = function(e){
      output <- paste0(readLines(temp, warn = FALSE), collapse = "\n")
      list(
        result = NULL, output = output, warnings = warnings, messages = messages,
        error_trace=e$trace, error=e$message
      )
    }
  )
}

safely2 <- function(f) {
  function(..., .logger=.logger, .console=TRUE) with_handlers(f(...), .logger=.logger, .console=.console)
}

funwrapper <- function(f) {
  return (function(..., .logger=NULL, .console=TRUE, .nowrap=FALSE) {
    if (.nowrap) {
      return(f(...))
    }
    res <- safely2(f)(..., .logger=.logger, .console=.console)
    jsonlite::toJSON(res, digits=2)
  })
}
