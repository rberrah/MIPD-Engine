# Public read-only search tables. No patient values, credentials or arbitrary URLs leave the session.
eucast_get <- function(antibiotic = NULL, page = 1L) {
  query <- list()
  if (!is.null(antibiotic)) {
    if (length(antibiotic) != 1 || !grepl("^[0-9]{1,6}$", antibiotic) || !page %in% 1:10) stop("Invalid EUCAST search.")
    query <- list("search[method]" = "mic", "search[antibiotic]" = antibiotic, "search[species]" = "-1",
      "search[disk_content]" = "-1", "search[limit]" = "50", page = page)
  }
  response <- httr::GET("https://mic.eucast.org/search/", query = query,
    httr::timeout(12), httr::config(followlocation = FALSE, maxfilesize = 5000000),
    httr::user_agent("Pharmacometrie-Pratique/1.0 (interactive EUCAST MIC lookup)"))
  httr::stop_for_status(response)
  if (httr::status_code(response) != 200) stop("EUCAST search moved; open the official site.")
  content <- httr::content(response, as = "raw")
  if (length(content) > 5000000) stop("EUCAST response exceeds 5 MB.")
  xml2::read_html(content, options = c("RECOVER", "NONET", "NOERROR", "NOWARNING"))
}

eucast_choices <- function(document = eucast_get()) {
  options <- xml2::xml_find_all(document, "//select[@id='search_antibiotic']/option")
  ids <- xml2::xml_attr(options, "value"); labels <- trimws(xml2::xml_text(options))
  keep <- grepl("^[0-9]+$", ids) & nzchar(labels)
  if (!any(keep)) stop("EUCAST antimicrobial list unavailable or format changed.")
  stats::setNames(ids[keep], labels[keep])
}

eucast_parse <- function(document) {
  table <- xml2::xml_find_first(document, "//table[@id='search-results-table']")
  if (inherits(table, "xml_missing")) stop("No EUCAST MIC table for this search, or format changed.")
  headings <- trimws(xml2::xml_text(xml2::xml_find_all(table, "./thead/tr[1]/th")))
  positions <- which(grepl("^[0-9]+(\\.[0-9]+)?$", headings))
  mic <- as.numeric(headings[positions])
  if (length(mic) < 3 || any(diff(mic) <= 0) || !all(c("Distributions", "Observations", "(T)ECOFF") %in% headings))
    stop("EUCAST MIC table columns changed; no data imported.")
  rows <- xml2::xml_find_all(table, "./tbody/tr[td/a]")
  result <- lapply(rows, function(row) {
    cells <- trimws(xml2::xml_text(xml2::xml_find_all(row, "./td")))
    if (length(cells) != length(headings)) stop("EUCAST row structure changed.")
    counts <- suppressWarnings(as.numeric(cells[positions]))
    total <- suppressWarnings(as.numeric(cells[match("Observations", headings)]))
    if (any(!is.finite(counts) | counts < 0 | counts != floor(counts)) || !is.finite(total) ||
        sum(counts) != total || total <= 0) stop("EUCAST observation counts do not reconcile.")
    href <- xml2::xml_attr(xml2::xml_find_first(row, "./td[1]/a"), "href")
    path <- sub("\\?.*$", "", href)
    if (!grepl("^/search/show-registration/[0-9]+$", path)) stop("Unexpected EUCAST source link.")
    list(species = cells[1], data = data.frame(mic = mic, count = counts), observations = total,
      distributions = cells[match("Distributions", headings)], ecoff = cells[match("(T)ECOFF", headings)],
      confidence = cells[match("Confidence interval", headings)], url = paste0("https://mic.eucast.org", path))
  })
  if (!length(result)) stop("No MIC distribution available.")
  list(rows = result, more = length(xml2::xml_find_all(document, "//a[@rel='next']")) > 0)
}

eucast_load <- function(antibiotic, label, get = eucast_get) {
  rows <- list()
  for (page in 1:10) {
    parsed <- eucast_parse(get(antibiotic, page))
    rows <- c(rows, parsed$rows)
    if (!parsed$more) break
    if (page == 10) stop("EUCAST pagination limit reached; use the official website.")
  }
  names(rows) <- vapply(rows, `[[`, character(1), "species")
  if (anyDuplicated(names(rows))) stop("Ambiguous EUCAST species rows.")
  list(antibiotic = antibiotic, label = label, accessed = format(Sys.Date()), rows = rows)
}
