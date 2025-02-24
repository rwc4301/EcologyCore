generate_ai <- function(prompt) {
  # Check for API key in environment
  api_key <- Sys.getenv("OPENAI_API_KEY")
  if (api_key == "") {
    return("Error: OpenAI API key not found. Please set the OPENAI_API_KEY environment variable.")
  }

  response <- tryCatch({
    POST(
      url = "https://api.openai.com/v1/chat/completions",
      add_headers(
        "Authorization" = paste("Bearer", api_key),
        "Content-Type" = "application/json"
      ),
      body = list(
        model = "gpt-4",
        messages = list(
          list(
            role = "user",
            content = prompt
          )
        ),
        temperature = 0.7  # Add some randomness to responses
      ),
      encode = "json"
    )
  }, error = function(e) {
    return(paste("Error generating summary:", e$message))
  })

  if (is.character(response)) {
    return(response)  # Return error message if request failed
  }

  content <- jsonlite::fromJSON(rawToChar(response$content))
  return(content$choices$message$content)
}
