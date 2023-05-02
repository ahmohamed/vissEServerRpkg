shortlist_genesets = function(nodes) {
  nodes |> dplyr::ungroup() |>
    dplyr::slice_head(n=50, by=Cluster) |>
    dplyr::mutate(rown=1:dplyr::n()) |> dplyr::group_by(Cluster) |>
    dplyr::filter(!any(rown > 330))
}

summarise_geneset_gpt = function(geneset_names) {
  prompt = paste("Your task is to generate a short summary from a list of genesets enriched in a biological condition.",
    "The list, delimited by triple backticks, contains one geneset per line and each line.",
    "Summarize the geneset list to produce the following format in JSON:",
    "<highlight>: list of 3 to 5 highlight sentences using miminal words, as bullet points",
    "<abstract>: academic description of the biological meanning of the results (max 200 words)",
    "",
    "```",
    geneset_names,
    "```")

  res = openai::create_chat_completion('gpt-3.5-turbo', messages = list(
    list(
      "role" = "system",
      "content" = "You are an academic biomedical researcher interpreting experiment results."
    ),
    list(
      "role" = "user",
      "content" = prompt
    )
  ), temperature = 0)


  jsonlite::fromJSON(res$choices$message.content)
}

get_gpt_summary = function(nodes) {
  if (nrow(nodes) < 20) {
    return (list(highlight=list(), abstract="Not enough genesets to generate an abstract"))
  }

  if (nrow(nodes) > 300) {
    nodes = shortlist_genesets(nodes)
  }
  geneset_names = paste(unlist(sub("\\w+ ", "", gsub("_", " ", tolower(nodes$Geneset)))), collapse = "\n")
  tryCatch(
    summarise_geneset_gpt(geneset_names),
    error = function(e) {
      list(highlight=list(), abstract="Failed to generate an abstract")
    }
  )
}
