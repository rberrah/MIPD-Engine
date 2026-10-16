pd_observations_ui <- function(id, t) {
  ns <- NS(id)
  tagList(DTOutput(ns("table")), div(class = "workbench-downloads",
    actionButton(ns("add"), t("Ajouter", "Add"), icon = icon("plus")),
    actionButton(ns("remove"), t("Supprimer la selection", "Remove selected"), icon = icon("trash"))))
}

pd_observations_server <- function(id, columns, initial, oncology = FALSE, dose_history = FALSE) {
  moduleServer(id, function(input, output, session) {
    t <- function(fr, en) app_t(app_language_from_query(session$clientData$url_search %||% ""), fr, en)
    stores <- reactiveVal(setNames(list(initial), paste(names(initial), collapse = ",")))
    revision <- reactiveVal(0L)
    key <- reactive(paste(columns(), collapse = ","))
    rows <- reactive(stores()[[key()]] %||% initial[FALSE, , drop = FALSE][, FALSE, drop = FALSE])
    validate <- function(data, expected = columns()) {
      if (nrow(data) > (if (dose_history) 500 else 2000) || !identical(names(data), expected)) stop("Invalid table columns or row count.")
      for (name in setdiff(names(data), "endpoint")) {
        supplied <- data[[name]]
        data[[name]] <- suppressWarnings(as.numeric(supplied))
        if (any(is.na(data[[name]]) & !is.na(supplied) & nzchar(trimws(as.character(supplied)))) || any(is.infinite(data[[name]]))) stop(t("Les valeurs doivent etre numeriques et finies.", "Values must be numeric and finite."))
      }
      if (any(data[[1]] < 0, na.rm = TRUE) || oncology && (any(data$value <= 0, na.rm = TRUE) || any(!data$endpoint %in% c("tumor", "anc")))) stop("Invalid observation: time >= 0; oncology value > 0; endpoint tumor or anc.")
      if (dose_history && (any(data$amount <= 0 | data$amount > 1e5, na.rm = TRUE) || any(data$infusion < 0 | data$infusion > 168, na.rm = TRUE))) {
        stop(t("Dose strictement positive (maximum 100000) ; perfusion entre 0 et 168 h.", "Dose must be positive (maximum 100000); infusion between 0 and 168 h."))
      }
      data
    }
    replace <- function(data) {
      data <- validate(data)
      value <- isolate(stores()); value[[isolate(key())]] <- data; stores(value)
    }
    replace_for <- function(data, target_columns) {
      data <- validate(data, target_columns)
      value <- isolate(stores()); value[[paste(target_columns, collapse = ",")]] <- data; stores(value)
    }
    output$table <- renderDT({
      revision()
      data <- rows()
      if (!ncol(data)) data <- as.data.frame(setNames(lapply(columns(), function(name) if (name == "endpoint") character() else numeric()), columns()))
      labels <- c(time = if (oncology || dose_history) t("Jour", "Day") else t("Temps (h)", "Time (h)"), concentration = "Concentration", effect = t("Valeur observee", "Observed value"), endpoint = t("Mesure", "Endpoint"), value = t("Valeur", "Value"),
        amount = t("Dose (unite PK)", "Dose (PK unit)"), infusion = t("Perfusion (h)", "Infusion (h)"))
      endpoint <- match("endpoint", names(data)) - 1L
      callback <- "table.on('keydown', 'tbody input', function(e) { if (e.key === 'Enter') { e.preventDefault(); this.blur(); } });"
      definitions <- NULL
      if (oncology) {
        choices <- c(tumor = t("Tumeur", "Tumor"), anc = t("Neutrophiles (ANC)", "Neutrophils (ANC)"))
        definitions <- list(list(targets = endpoint, render = DT::JS(sprintf("function(data, type) {
          if (type !== 'display') return data;
          var select = document.createElement('select');
          select.className = 'pd-endpoint'; select.setAttribute('aria-label', %s);
          var labels = %s;
          Object.keys(labels).forEach(function(key) { select.add(new Option(labels[key], key, key === data, key === data)); });
          return select.outerHTML;
        }", jsonlite::toJSON(unname(labels['endpoint']), auto_unbox = TRUE), jsonlite::toJSON(as.list(choices), auto_unbox = TRUE)))))
        callback <- paste(callback, sprintf("table.on('click dblclick', 'select.pd-endpoint', function(e) { e.stopPropagation(); });
          table.on('change', 'select.pd-endpoint', function() {
            var index = table.cell($(this).closest('td')).index();
            Shiny.setInputValue(%s, {row: index.row + 1, value: this.value}, {priority: 'event'});
          });", jsonlite::toJSON(session$ns('endpoint_edit'), auto_unbox = TRUE)))
      }
      datatable(data, colnames = unname(labels[names(data)]), rownames = FALSE, selection = "multiple",
        editable = if (oncology) list(target = "cell", disable = list(columns = endpoint)) else "cell",
        callback = DT::JS(callback),
        options = list(pageLength = 8, lengthChange = FALSE, searching = FALSE, scrollX = TRUE, dom = "tip",
          columnDefs = definitions,
          language = list(emptyTable = if (dose_history) t("Aucune dose passee", "No past doses") else t("Aucune observation", "No observations"), info = t("_START_ a _END_ sur _TOTAL_", "_START_ to _END_ of _TOTAL_"), infoEmpty = "0",
            paginate = list(previous = t("Precedent", "Previous"), "next" = t("Suivant", "Next")))))
    }, server = FALSE)
    handle <- function(fn) tryCatch(fn(), error = function(e) { showNotification(conditionMessage(e), type = "error", duration = 7); revision(isolate(revision()) + 1L) })
    observeEvent(input$table_cell_edit, handle(function() {
      edit <- input$table_cell_edit; data <- rows()
      if (!edit$row %in% seq_len(nrow(data)) || !edit$col %in% (seq_len(ncol(data)) - 1L)) stop("Invalid table cell.")
      data[edit$row, edit$col + 1L] <- edit$value
      replace(data)
    }))
    observeEvent(input$endpoint_edit, handle(function() {
      edit <- input$endpoint_edit; data <- rows()
      if (!oncology || length(edit$row) != 1L || !edit$row %in% seq_len(nrow(data)) ||
          length(edit$value) != 1L || !edit$value %in% c("tumor", "anc")) stop("Invalid oncology endpoint.")
      data$endpoint[edit$row] <- edit$value
      replace(data)
    }))
    observeEvent(input$add, handle(function() {
      row <- as.data.frame(setNames(lapply(columns(), function(name) if (name == "endpoint") "tumor" else if (name %in% c("value", "effect", "amount")) NA_real_ else 0), columns()))
      replace(if (ncol(rows())) rbind(rows(), row) else row)
    }))
    observeEvent(input$remove, handle(function() {
      selected <- input$table_rows_selected
      if (length(selected)) replace(rows()[-selected, , drop = FALSE])
    }))
    session$onSessionEnded(function() stores(list()))
    list(data = rows, replace = replace, replace_for = replace_for)
  })
}
