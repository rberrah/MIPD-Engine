normalize_app_language <- function(value) {
  value <- tolower(as.character(value[[1]] %||% "fr"))
  if (identical(value, "en")) "en" else "fr"
}

app_language_from_query <- function(query = "") {
  query <- sub("^\\?", "", as.character(query %||% ""))
  fields <- strsplit(query, "&", fixed = TRUE)[[1]]
  pair <- fields[startsWith(fields, "lang=")][1]
  if (!length(pair) || is.na(pair)) return("fr")
  normalize_app_language(utils::URLdecode(sub("^lang=", "", pair)))
}

app_language_from_request <- function(request) {
  query <- tryCatch(request$QUERY_STRING, error = function(error) "")
  app_language_from_query(query)
}

APP_TRANSLATIONS <- c(
  "Non calculable" = "Not available",
  "Analyse" = "Analysis",
  "Bibliothèque" = "Library",
  "Réglages" = "Settings",
  "Méthode" = "Method",
  "Outil de recherche et d'enseignement" = "Research and teaching tool",
  "Les résultats ne remplacent ni la validation locale du modèle ni le jugement clinique." = "Results do not replace local model validation or clinical judgment.",
  "Configurer et lancer" = "Configure and run",
  "Configuration" = "Configuration",
  "Modèle" = "Model",
  "Atelier Lego / C++" = "Lego workshop / C++",
  "Modèle principal" = "Primary model",
  "Mode local : le C++ libre est autorisé et s'exécute avec les droits du processus R." = "Local mode: unrestricted C++ is enabled and runs with the permissions of the R process.",
  "Serveur public : seuls les modèles portant la spécification contrôlée de l'Atelier Lego sont acceptés. Le C++ libre reste bloqué." = "Public server: only models carrying a validated Lego workshop specification are accepted. Unrestricted C++ remains disabled.",
  "Voie d'administration" = "Administration route",
  "Activer le model averaging" = "Enable model averaging",
  "Seuls les modèles de la même molécule compatibles avec la voie sélectionnée sont proposés." = "Only models for the same drug and selected route are offered.",
  "Seuls les modèles de la même molécule, de la même voie et du même mode d'administration sont proposés." = "Only models for the same drug, route, and administration mode are offered.",
  "Pondération" = "Weighting",
  "Critère d'Akaike" = "Akaike criterion",
  "Log-vraisemblance" = "Log likelihood",
  "Valider le modèle" = "Validate model",
  "Cible et grille de doses" = "Target and dose grid",
  "Métrique" = "Metric",
  "Concentration minimale" = "Minimum concentration",
  "Concentration maximale" = "Maximum concentration",
  "% du temps au-dessus d'une concentration" = "% of time above a concentration",
  "Cmin et Cmax simultanément" = "Cmin and Cmax simultaneously",
  "Seuil de concentration" = "Concentration threshold",
  "Les bornes ci-dessous indiquent le pourcentage cible de l'intervalle posologique passé au-dessus du seuil." = "The bounds below define the target percentage of the dosing interval spent above the threshold.",
  "Les premières bornes s'appliquent à Cmin. Renseignez ensuite l'intervalle cible de Cmax." = "The first bounds apply to Cmin. Then enter the Cmax target range.",
  "Cmax basse" = "Lower Cmax",
  "Cmax haute" = "Upper Cmax",
  "Cmin basse" = "Lower Cmin",
  "Cmin haute" = "Upper Cmin",
  "Borne basse" = "Lower bound",
  "Borne haute" = "Upper bound",
  "Dose min" = "Minimum dose",
  "Dose max" = "Maximum dose",
  "Pas" = "Step",
  "Intervalles testés" = "Candidate intervals",
  "Durée de perfusion (h, 0 = bolus IV)" = "Infusion duration (h, 0 = IV bolus)",
  "Voie orale : durée de perfusion fixée à 0 h." = "Oral route: infusion duration is fixed at 0 h.",
  "Voie intramusculaire : perfusion externe fixée à 0 h; absorption définie par le modèle." = "Intramuscular route: external infusion is fixed at 0 h; absorption is defined by the model.",
  "J'ai lu le statut et j'accepte d'utiliser ce prototype sous ma responsabilité." = "I have read the status statement and accept responsibility for using this prototype.",
  "Lancer l'analyse" = "Run analysis",
  "Données" = "Data",
  "Historique thérapeutique" = "Treatment history",
  "Saisissez des heures relatives, des jours avant une référence ou des dates civiles. Le moteur normalise ensuite l'origine sur la première administration." = "Enter relative hours, days before a reference date, or calendar dates. The engine then sets the first administration as time zero.",
  "Importer JSON" = "Import JSON",
  "Aucun fichier" = "No file selected",
  "Exporter JSON" = "Export JSON",
  "Session uniquement" = "Session only",
  "Le fichier est lu, validé puis supprimé du stockage temporaire. Aucun identifiant, dossier patient ou code C++ n'est enregistré." = "The file is read, validated, then removed from temporary storage. No identifier, patient record, or C++ code is stored.",
  "Format temporel" = "Time format",
  "Heures relatives" = "Relative hours",
  "Jours avant une date de référence" = "Days before a reference date",
  "Dates et heures" = "Dates and times",
  "Date de référence" = "Reference date",
  "Heure de référence" = "Reference time",
  "Administrations" = "Administrations",
  "Ajouter une administration" = "Add administration",
  "Retirer la dernière" = "Remove last",
  "Concentrations observées" = "Observed concentrations",
  "Ajouter un prélèvement" = "Add sample",
  "Retirer le dernier" = "Remove last",
  "Ajustement" = "Fit",
  "Ajustement bayésien" = "Bayesian fit",
  "Créer le rapport" = "Create report",
  "Doses supplémentaires : poursuivre ou modifier" = "Additional doses: maintain or modify",
  "Sensibilité du model averaging" = "Model averaging sensitivity",
  "Modèles et pondérations" = "Models and weights",
  "Posologies" = "Dosage regimens",
  "Distribution prédictive du meilleur scénario" = "Predictive distribution for the best scenario",
  "Scénarios classés" = "Ranked scenarios",
  "Exporter CSV" = "Export CSV",
  "Modèle mrgsolve complet" = "Complete mrgsolve model",
  "Bibliothèque de modèles" = "Model library",
  "Chaque modèle est associé à l'article qui décrit ses paramètres et sa population. Toute nouvelle entrée est relue avant publication." = "Each model is linked to the article describing its parameters and population. Every new entry is reviewed before publication.",
  "Réglages de simulation" = "Simulation settings",
  "Ces paramètres s'appliquent à la prochaine analyse et sont inclus dans le rapport." = "These settings apply to the next analysis and are included in the report.",
  "Résolution et projection" = "Resolution and projection",
  "Pas de simulation" = "Simulation step",
  "Haute · 0,05 h" = "High · 0.05 h",
  "Standard · 0,1 h" = "Standard · 0.1 h",
  "Rapide · 0,25 h" = "Fast · 0.25 h",
  "Doses supplémentaires comparées" = "Additional doses compared",
  "Afficher chaque modèle sur l'ajustement" = "Show each model on the fit",
  "Distribution prédictive" = "Predictive distribution",
  "Réplications Monte Carlo" = "Monte Carlo replicates",
  "Intervalle prédictif" = "Prediction interval",
  "Composantes simulées" = "Simulated components",
  "Incertitude postérieure MAP" = "MAP posterior uncertainty",
  "Erreur résiduelle" = "Residual error",
  "Incertitude des horaires" = "Timing uncertainty",
  "Réestimations avec horaires incertains" = "Refits with uncertain times",
  "Concentrations BLQ" = "BLQ concentrations",
  "Exclure de l'ajustement" = "Exclude from fit",
  "LLOQ / 2 (exploratoire)" = "LLOQ / 2 (exploratory)",
  "Erreur résiduelle MAP" = "MAP residual error",
  "Erreur utilisée pour l'ajustement" = "Error used for fitting",
  "Conserver l'erreur du modèle" = "Keep model error",
  "Fixer un CV proportionnel" = "Set a proportional CV",
  "CV proportionnel fixé (%)" = "Fixed proportional CV (%)",
  "Le mode fixé remplace les SIGMA de l'ajustement par le CV proportionnel choisi, ou son équivalent logarithmique pour un modèle exponentiel. Il modifie le poids des concentrations dans MAP-BE." = "Fixed mode replaces fitting SIGMA values with the selected proportional CV, or its log-scale equivalent for an exponential model. It changes the weight of concentrations in MAP-BE.",
  "Pas de lien partageable" = "No shareable link",
  "Les données patient ne sont pas encodées dans l'URL afin d'éviter leur présence dans l'historique du navigateur, les journaux réseau ou les outils d'analytique." = "Patient data are not encoded in the URL, preventing their inclusion in browser history, network logs, or analytics tools.",
  "Méthode et limites" = "Methods and limitations",
  "Estimation individuelle" = "Individual estimation",
  "Par défaut, l'ajustement conserve les variances résiduelles du modèle publié. Le réglage CV proportionnel fixé remplace chaque paire SIGMA par une variance proportionnelle commune, ou sa variance logarithmique équivalente si le modèle utilise une erreur exponentielle. Ce choix change le poids des observations dans l'estimation MAP et doit être documenté." = "By default, fitting retains the published model residual variances. The fixed proportional CV setting replaces each SIGMA pair with a common proportional variance, or its equivalent log variance when the model uses exponential error. This choice changes observation weighting in MAP estimation and must be documented.",
  "Model averaging" = "Model averaging",
  "Exposition et scénarios" = "Exposure and scenarios",
  "Cibles pharmacocinétiques" = "Pharmacokinetic targets",
  "La cible temps au-dessus d'un seuil est calculée sur un intervalle posologique complet à l'état stationnaire et exprimée en pourcentage. Les franchissements du seuil sont interpolés entre deux points simulés. Pour la cible combinée, Cmin et Cmax suivent chacun une contrainte indépendante: intervalle, valeur minimale ou valeur maximale. Les deux contraintes doivent être respectées; aucune moyenne entre les deux critères n'est utilisée." = "The time-above-threshold target is calculated over one complete steady-state dosing interval and expressed as a percentage. Threshold crossings are interpolated between simulated points. For the combined target, Cmin and Cmax each follow an independent constraint: range, minimum value, or maximum value. Both constraints must be met; the two criteria are not averaged.",
  "Qualité des données" = "Data quality",
  "Distribution" = "Distribution",
  "Apprentissage automatique" = "Machine learning",
  "Sécurité" = "Security",
  "Statut" = "Status",
  "Avertissement obligatoire" = "Mandatory warning",
  "Les effets aléatoires individuels sont estimés par maximum a posteriori avec mapbayr. Les administrations, covariables datées et concentrations sont converties en événements NM-TRAN puis simulées avec mrgsolve." = "Individual random effects are estimated by maximum a posteriori with mapbayr. Administrations, dated covariates, and concentrations are converted to NM-TRAN events and simulated with mrgsolve.",
  "Chaque modèle analyse les mêmes données. Les prédictions sont moyennées avec des poids issus de la vraisemblance ou du critère d'Akaike. Le serveur refuse l'agrégation de modèles ne partageant pas la même molécule, la même voie et le même mode d'administration." = "Each model analyzes the same data. Predictions are averaged using likelihood- or Akaike-based weights. The server rejects combinations that do not share the same drug, route, and administration mode.",
  "La robustesse est explorée avec des poids égaux, AIC, log-vraisemblance et des analyses laissant successivement de côté chaque modèle. Une divergence des doses proposées doit conduire à revoir l'applicabilité du model averaging." = "Robustness is explored with equal, AIC, and log-likelihood weights and leave-one-model-out analyses. Divergent proposed doses should prompt a review of model averaging applicability.",
  "L'exposition historique sur les dernières 24 heures est distinguée de l'exposition à l'état stationnaire du dernier schéma. La projection compare le maintien de la dernière posologie à l'application du scénario classé en tête, sur le nombre de doses supplémentaires choisi dans Réglages." = "Historical exposure over the latest 24 hours is separated from steady-state exposure for the last regimen. The projection compares maintaining that regimen with applying the top-ranked scenario over the number of additional doses selected in Settings.",
  "Une administration déclarée à l'état stationnaire est initialisée à t = 0 avec ss = 1. Les doses de la projection sont ensuite ajoutées explicitement avec ss = 0; la comparaison des doses supplémentaires continue donc au-delà de cet état initial." = "An administration declared at steady state is initialized at t = 0 with ss = 1. Projection doses are then added explicitly with ss = 0, so the additional-dose comparison continues beyond that initial state.",
  "Pour un modèle Lego qui répartit une dose entre une ou plusieurs voies d'entrée, l'état stationnaire est approché par 50 administrations explicites avant t = 0. Cette phase de préchauffage évite d'appliquer séparément ss = 1 à chaque voie, ce qui serait incorrect pour des cinétiques parallèles ou non linéaires." = "For a Lego model that splits one dose across one or more input pathways, steady state is approximated with 50 explicit administrations before t = 0. This warm-up avoids applying ss = 1 separately to each pathway, which would be incorrect for parallel or nonlinear kinetics.",
  "Les scénarios sont d'abord classés sur leur prédiction moyenne. La probabilité d'atteindre la cible, de sous-exposition et de surexposition est ensuite simulée pour les douze scénarios moyens les plus proches, puis utilisée pour leur classement final." = "Scenarios are first ranked by their mean prediction. Target attainment, underexposure, and overexposure probabilities are then simulated for the twelve closest mean scenarios and used in the final ranking.",
  "Les doses oubliées sont exclues. Les horaires incertains peuvent faire l'objet de réestimations de sensibilité. Les valeurs BLQ sont exclues par défaut; l'imputation LLOQ/2 est uniquement exploratoire. La matrice, les unités et le domaine exact des covariables ne sont jamais convertis ou déduits automatiquement." = "Missed doses are excluded. Uncertain times can be assessed through sensitivity refits. BLQ values are excluded by default; LLOQ/2 imputation is exploratory only. Matrix, units, and exact covariate domains are never converted or inferred automatically.",
  "La distribution prédictive est exploratoire. Après un ajustement, elle repose sur un bootstrap paramétrique de l'estimation MAP, auquel peuvent s'ajouter l'erreur résiduelle, l'incertitude des horaires et l'incertitude entre modèles. Sans concentration exploitable, elle revient à une simulation populationnelle." = "The predictive distribution is exploratory. After fitting, it uses a parametric bootstrap of the MAP estimate, optionally adding residual error, timing uncertainty, and between-model uncertainty. Without a usable concentration, it becomes a population simulation.",
  "Le module expérimental adapte à chaque modèle de la bibliothèque la méthodologie de simulation publiée pour le tacrolimus (doi:10.1016/j.phrs.2021.105578). XGBoost apprend le logarithme du rapport entre l'AUC24 individuelle simulée et l'AUC24 populationnelle du même schéma. Il utilise les concentrations et horaires, leurs prédictions populationnelles, la dose, l'intervalle, la durée de perfusion et les covariables du modèle. Une décomposition locale DALEX explique la prédiction par rapport à un échantillon de référence entièrement synthétique." = "The experimental module adapts the simulation methodology published for tacrolimus (doi:10.1016/j.phrs.2021.105578) to each library model. XGBoost learns the log ratio between simulated individual AUC24 and population AUC24 under the same regimen. It uses concentrations and times, their population predictions, dose, interval, infusion duration, and model covariates. A local DALEX decomposition explains the prediction against a fully synthetic reference sample.",
  "Une administration déclarée à l'état stationnaire (ss = 1) et au moins deux concentrations dans un même intervalle posologique sont requises. Cet intervalle ne doit pas nécessairement être le dernier. Les deux concentrations les plus récentes de l'intervalle admissible le plus récent alimentent le prédicteur." = "An administration declared at steady state (ss = 1) and at least two concentrations in the same dosing interval are required. This does not have to be the latest interval. The two latest concentrations from the latest eligible interval feed the predictor.",
  "L'estimation ML est une AUC24 expérimentale. Elle est comparée à l'AUC24 MAP du schéma répété à l'état stationnaire, et non à l'intégrale MAP historique partielle lorsque moins de 24 h sont disponibles. En model averaging, chaque modèle produit sa propre AUC24 ML, puis les estimations sont agrégées avec les poids de l'analyse. Tous les modèles doivent disposer d'un artefact compatible et partager la même voie et le même mode d'administration." = "The ML estimate is an experimental AUC24. It is compared with steady-state MAP AUC24 for the repeated regimen, not with the partial historical MAP integral when fewer than 24 hours are available. With model averaging, each model produces its own ML AUC24 and the estimates are then aggregated using the analysis weights. Every model must have a compatible artifact and share the same route and administration mode.",
  "Chaque artefact reste lié au même modèle PK, à la même voie, au même mode d'administration, au même schéma de variables et aux empreintes exactes du fichier mrgsolve et du booster. Le statut recherche exige le respect des seuils préspecifiés en validation croisée répétée et sur un test interne non touché; un artefact évalué mais sous ces seuils reste signalé comme expérimental. La transportabilité vers un autre modèle PopPK est rapportée séparément." = "Each artifact remains tied to the same PK model, route, administration mode, feature schema, and exact fingerprints of the mrgsolve file and booster. Research status requires meeting prespecified thresholds in repeated cross-validation and on an untouched internal test set; an evaluated artifact below those thresholds remains labelled experimental. Transportability to another PopPK model is reported separately.",
  "Une variable hors des bornes empiriques d'entraînement déclenche un avertissement sans masquer l'estimation. Une donnée manquante, moins de deux concentrations dans un même intervalle, l'absence d'état stationnaire, un protocole incompatible ou une AUC invalide restent bloquants. Ces entraînements synthétiques ne constituent pas une validation clinique pour les autres molécules. L'estimation ML ne remplace pas les projections de dose MAP tant qu'une validation favorable sur des patients réels indépendants n'est pas documentée." = "A feature outside empirical training bounds triggers a warning without hiding the estimate. Missing data, fewer than two concentrations in one interval, no steady state, an incompatible protocol, or an invalid AUC remain blocking. These synthetic training runs do not constitute clinical validation for the other drugs. ML does not replace MAP dose projections until favorable validation in independent real patients is documented.",
  "Le serveur public ne compile jamais directement le C++ reçu. Pour un modèle Atelier Lego, il extrait une spécification JSON, la valide, régénère lui-même le code mrgsolve puis compile uniquement ce code contrôlé. Tout autre C++ reste refusé tant qu'il n'est pas exécuté dans un conteneur éphémère isolé." = "The public server never compiles submitted C++ directly. For a Lego workshop model, it extracts and validates a JSON specification, regenerates mrgsolve code itself, and compiles only that controlled code. Other C++ remains rejected unless it runs in an isolated ephemeral container.",
  "Les imports JSON sont traités dans la session Shiny et leur fichier temporaire est supprimé immédiatement après lecture. Les exports sont produits à la demande sans base de données." = "JSON imports are processed within the Shiny session and their temporary file is deleted immediately after reading. Exports are generated on demand without a database.",
  "Ce prototype est destiné à la recherche et à l'enseignement. Il n'est pas enregistré comme dispositif médical, ne garantit ni l'exactitude d'un résultat ni son applicabilité à un patient particulier et ne remplace pas le jugement clinique." = "This prototype is intended for research and teaching. It is not registered as a medical device, does not guarantee result accuracy or applicability to an individual patient, and does not replace clinical judgment.",
  "Toute décision de dose reste sous la responsabilité du professionnel de santé et exige la vérification de la voie, des horaires, des unités, de la population source, des covariables, des concentrations, de la fonction d'organe et des recommandations locales. Une validation indépendante et une gouvernance documentée sont nécessaires avant toute utilisation clinique." = "Every dose decision remains the healthcare professional's responsibility and requires verification of route, times, units, source population, covariates, concentrations, organ function, and local guidance. Independent validation and documented governance are required before clinical use.",
  "Si vous recherchez un véritable dispositif médical avec expertise, consultez " = "For an actual medical device supported by clinical expertise, see ",
  "ABIS du CHU de Limoges." = "ABIS at Limoges University Hospital.",
  "PMx Explain · moteur R mrgsolve/mapbayr · aucun dossier patient n'est persisté" = "PMx Explain · R mrgsolve/mapbayr engine · no patient record is persisted",
  "Temps relatif (h)" = "Relative time (h)",
  "Jours avant référence" = "Days before reference",
  "Heure (HH:MM)" = "Time (HH:MM)",
  "Date/heure" = "Date/time",
  "JJ/MM HH:MM ou JJ/MM/AAAA HH:MM" = "DD/MM HH:MM or DD/MM/YYYY HH:MM",
  "Dose (mg)" = "Dose (mg)",
  "Administrée" = "Administered",
  "Horaire incertain" = "Uncertain time",
  "Oubliée / non prise" = "Missed / not taken",
  "Incertitude horaire (± h)" = "Time uncertainty (± h)",
  "Intervalle (h)" = "Interval (h)",
  "Nombre" = "Count",
  "Voie orale" = "Oral route",
  "La perfusion est imposée à 0 h." = "Infusion is fixed at 0 h.",
  "Voie intramusculaire" = "Intramuscular route",
  "La perfusion externe est imposée à 0 h; l'absorption est définie par le modèle." = "External infusion is fixed at 0 h; absorption is defined by the model.",
  "Temps imposé à t = 0 h (ss = 1)." = "Time fixed at t = 0 h (ss = 1).",
  "Concentration" = "Concentration",
  "Sous la limite de quantification (BLQ)" = "Below the limit of quantification (BLQ)",
  "Matrice biologique" = "Biological matrix",
  "Non précisée" = "Unspecified",
  "Sérum" = "Serum",
  "Sang total" = "Whole blood",
  "Temps (h)" = "Time (h)",
  "Modèle principal" = "Primary model",
  "Aucune analyse disponible." = "No analysis available.",
  "Molécule" = "Drug",
  "Voie" = "Route",
  "Administration" = "Administration",
  "Implémentation" = "Implementation",
  "Article source" = "Source article",
  "Population de l'article" = "Article population",
  "Type" = "Type",
  "Usage TDM" = "TDM use",
  "Statut bibliographique" = "Reference status",
  "Détail" = "Details",
  "Dans cible" = "In target",
  "Poids" = "Weight",
  "Perfusion" = "Infusion",
  "Distance" = "Distance",
  "Sous-cible (%)" = "Below target (%)",
  "Dans cible (%)" = "In target (%)",
  "Sur-cible (%)" = "Above target (%)",
  "Réplications" = "Replicates"
  ,"Temps relatif manquant ou invalide." = "Relative time is missing or invalid."
  ,"Le nombre de jours doit être positif." = "The number of days must be non-negative."
  ,"L'heure doit respecter HH:MM." = "Time must use HH:MM."
  ,"Utilisez JJ/MM HH:MM, JJ/MM/AAAA HH:MM ou AAAA-MM-JJ HH:MM." = "Use DD/MM HH:MM, DD/MM/YYYY HH:MM, or YYYY-MM-DD HH:MM."
  ,"La dose doit être strictement positive." = "Dose must be strictly positive."
  ,"L'intervalle doit être positif ou nul." = "Interval must be non-negative."
  ,"Le steady state exige un intervalle strictement positif." = "Steady state requires a strictly positive interval."
  ,"Le nombre d'administrations doit être au moins 1." = "Administration count must be at least 1."
  ,"La perfusion doit être positive ou nulle." = "Infusion duration must be non-negative."
  ,"La perfusion dépasse l'intervalle entre deux doses." = "Infusion duration exceeds the dosing interval."
  ,"L'incertitude horaire doit être positive ou nulle." = "Time uncertainty must be non-negative."
  ,"Précisez une incertitude horaire pour cette dose." = "Enter a time uncertainty for this dose."
  ,"Cette dose sera exclue des administrations reçues." = "This dose will be excluded from received administrations."
  ,"La concentration doit être positive ou nulle." = "Concentration must be non-negative."
  ,"Une LLOQ strictement positive est requise pour une valeur BLQ." = "A strictly positive LLOQ is required for a BLQ value."
  ,"Précisez une incertitude horaire pour ce prélèvement." = "Enter a time uncertainty for this sample."
  ,"Matrice biologique non précisée : vérifiez sa compatibilité avec l'article." = "Biological matrix unspecified: check compatibility with the article."
  ,"Ligne valide" = "Valid row"
  ,"Saisie cohérente" = "Consistent input"
  ,"Aucune anomalie de forme détectée avant l'analyse." = "No structural issue detected before analysis."
  ,"Rapport TDM" = "TDM report"
  ,"Rapport d'analyse TDM" = "TDM analysis report"
  ,"AUC actuelle glissante" = "Current rolling AUC"
  ,"AUC24 ML expérimentale" = "Experimental ML AUC24"
  ,"C0 actuelle (MAP)" = "Current C0 (MAP)"
  ,"AUC0-24 à l'état stationnaire" = "Steady-state AUC0-24"
  ,"C0 à l'état stationnaire" = "Steady-state C0"
  ,"Recommandation" = "Recommendation"
  ,"Applicabilité et qualité" = "Applicability and quality"
  ,"Ajustement" = "Fit"
  ,"Ajustement pharmacocinétique" = "Pharmacokinetic fit"
  ,"Comparaison des doses supplémentaires" = "Additional-dose comparison"
  ,"Comparaison de posologies" = "Dosage-regimen comparison"
  ,"Concordance MAP-BE et ML" = "MAP-BE and ML agreement"
  ,"Comparaison des AUC24 MAP-BE et ML" = "MAP-BE and ML AUC24 comparison"
  ,"La recommandation de dose reste calculée par MAP-BE." = "The dose recommendation remains MAP-BE based."
  ,"Explication locale DALEX de l'AUC24 ML" = "Local DALEX explanation of ML AUC24"
  ,"Les contributions décrivent le calcul de XGBoost par rapport à une référence synthétique; elles ne sont pas causales." = "Contributions describe the XGBoost calculation relative to a synthetic baseline; they are not causal."
  ,"Décomposition locale DALEX de l'AUC24 ML" = "Local DALEX decomposition of ML AUC24"
  ,"Concentrations et covariables" = "Concentrations and covariates"
  ,"Références" = "References"
  ,"Traçabilité" = "Traceability"
  ,"Empreintes SHA-256 des modèles" = "Model SHA-256 fingerprints"
  ,"Hypothèses" = "Assumptions"
  ,"Prototype de recherche et d'enseignement" = "Research and teaching prototype"
  ,"Ce rapport ne constitue pas une prescription ni une recommandation clinique validée. Vérifiez les données, unités, horaires, voie d'administration, population source et limites du modèle. Toute décision reste sous la responsabilité du professionnel de santé." = "This report is not a prescription or a validated clinical recommendation. Verify data, units, times, administration route, source population, and model limitations. Every decision remains the healthcare professional's responsibility."
)

app_t <- function(lang, fr, en = NULL) {
  if (!identical(normalize_app_language(lang), "en")) return(fr)
  if (!is.null(en)) return(en)
  translated <- unname(APP_TRANSLATIONS[as.character(fr)])
  missing <- is.na(translated)
  translated[missing] <- as.character(fr)[missing]
  attributes(translated) <- attributes(fr)
  translated
}

localize_ui <- function(value, lang) {
  if (!identical(normalize_app_language(lang), "en") || is.null(value)) return(value)
  if (is.character(value)) return(app_t(lang, value))
  if (inherits(value, "html_dependency")) return(value)
  if (inherits(value, "shiny.tag")) {
    for (attribute in intersect(c("title", "aria-label", "placeholder"), names(value$attribs))) {
      value$attribs[[attribute]] <- app_t(lang, value$attribs[[attribute]])
    }
    value$children <- lapply(value$children, localize_ui, lang = lang)
    return(value)
  }
  if (is.list(value)) {
    for (index in seq_along(value)) value[index] <- list(localize_ui(value[[index]], lang))
  }
  value
}

catalog_choices_i18n <- function(drug = NULL, route = NULL, mode = NULL, lang = "fr", analysis_only = FALSE) {
  rows <- MODEL_CATALOG
  if (isTRUE(analysis_only)) {
    keep <- vapply(seq_len(nrow(rows)), function(index) model_analysis_eligible(rows[index, , drop = FALSE]), logical(1))
    rows <- rows[keep, , drop = FALSE]
  }
  if (!is.null(drug)) rows <- rows[rows$drug == drug, , drop = FALSE]
  if (!is.null(route)) {
    keep <- vapply(seq_len(nrow(rows)), function(index) model_supports_route(rows[index, , drop = FALSE], route), logical(1))
    rows <- rows[keep, , drop = FALSE]
  }
  if (!is.null(mode)) {
    keep <- vapply(seq_len(nrow(rows)), function(index) model_supports_administration_mode(rows[index, , drop = FALSE], route, mode), logical(1))
    rows <- rows[keep, , drop = FALSE]
  }
  drug_labels <- if (identical(normalize_app_language(lang), "en") && "drugEn" %in% names(rows)) rows$drugEn else rows$drug
  stats::setNames(rows$id, paste(drug_labels, rows$model, sep = " - "))
}
