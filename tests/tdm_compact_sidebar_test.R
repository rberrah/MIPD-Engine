# Regression guard: the TDM sidebar exposes decisions, not duplicated guides.
# This tests presentation contracts only. Clinical eligibility remains exercised
# by tests/tdm_drug_profiles_test.R and tests/tdm_guide_presets_test.R.
app <- readLines("app.R", encoding = "UTF-8", warn = FALSE)
between <- function(first, last) {
  a <- grep(first, app, fixed = TRUE)
  b <- grep(last, app, fixed = TRUE)
  stopifnot(length(a) == 1L, length(b) == 1L, a < b)
  paste(app[a:b], collapse = "\n")
}
profile <- between('  output$tdm_drug_profile_ui <- renderUI({',
  '  shiny::outputOptions(output, "tdm_drug_profile_ui", suspendWhenHidden = FALSE)')
preset <- between('  output$target_preset_ui <- renderUI({',
  '  shiny::outputOptions(output, "target_preset_ui", suspendWhenHidden = FALSE)')
guide <- between('  output$tdm_guide_preset_ui <- renderUI({',
  '  shiny::outputOptions(output, "tdm_guide_preset_ui", suspendWhenHidden = FALSE)')
stopifnot(
  sum(grepl('uiOutput("tdm_drug_profile_ui")', app, fixed = TRUE)) == 1L,
  sum(grepl('uiOutput("tdm_guide_preset_ui")', app, fixed = TRUE)) == 1L,
  !any(grepl('radioButtons("tdm_guide_target_choice"', app, fixed = TRUE)),
  grepl('selectInput("tdm_guide_target_choice"', guide, fixed = TRUE),
  grepl('stats::setNames("", tx("Choisir une référence"', guide, fixed = TRUE),
  grepl('selected = selected_id', guide, fixed = TRUE),
  grepl('uiOutput("tdm_guide_drug_link_ui")', profile, fixed = TRUE),
  grepl('tags$details(', profile, fixed = TRUE),
  grepl('tags$details(', preset, fixed = TRUE),
  !grepl('doseGrid', guide, fixed = TRUE),
  !grepl('Grille de doses déclarée par MIPD', guide, fixed = TRUE),
  grepl('tdm_drug_profile_grid_eligible(record, mode, confirmed, st)', profile, fixed = TRUE),
  grepl('tdm_drug_profile_target_eligible(record, selected, st)', guide, fixed = TRUE),
  grepl('tdm_guide_preset_grid_applicable(record, mode)', paste(app, collapse = "\n"), fixed = TRUE),
  grepl('Simulation uniquement', profile, fixed = TRUE),
  grepl('Concentrations en ng/mL', profile, fixed = TRUE)
)
cat("Compact TDM sidebar, explicit target choice and safety gates: PASS\n")
