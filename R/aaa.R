suppressively <- function(f) {
  function(...) suppressPackageStartupMessages(f(...))
}


with_handlers <- function(code, .logger=NULL, .console=TRUE) {
  messages <- character()
  mHandler <- function(m) {
    if (!is.null(.logger)) .logger(m$message)

    fmt = futile.logger::flog.layout()(futile.logger::INFO, m$message)
    if (.console) message(fmt)
    messages <<- c(messages, fmt)
    invokeRestart("muffleMessage")
  }
  warnings <- character()
  wHandler <- function(w) {
    if (!grepl("^package .* was built under R", w$message)) {
      if (!is.null(.logger)) .logger(w$message)
      if (.console) futile.logger::flog.info(w$message)

      fmt = futile.logger::flog.layout()(futile.logger::WARN, w$message)
      if (.console) message(fmt)
      messages <<- c(messages, fmt)
      # warnings <<- c(warnings, w$message)
    }
    invokeRestart("muffleWarning")
  }
  .add_trace <- function(e) {
    e$trace=capture.output(rlang::trace_back())
    fmt = futile.logger::flog.layout()(futile.logger::FATAL, e$message)
    messages <<- c(messages, fmt)
    fmt = futile.logger::flog.layout()(futile.logger::FATAL, paste(e$trace, collapse = "\n"))
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
  return (function(..., .logger=NULL, .console=TRUE) {
    res <- safely2(f)(..., .logger=.logger, .console=.console)
    jsonlite::toJSON(res, digits=2)
  })
}
