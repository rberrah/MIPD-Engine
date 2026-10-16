# MIPD Engine

Application Shiny séparée du site statique. Elle exécute les modèles avec `mrgsolve`, individualise les paramètres avec `mapbayr` ou `Pmetrics::PM_model$map()` et compare des scénarios de dose.

Le déploiement public Shinyapps utilise `mrgsolve` et `mapbayr`. Pmetrics y est
désactivé, car la chaîne Rust/Cargo actuelle de Shinyapps.io ne construit pas
Pmetrics 3.2.6. Le backend Pmetrics reste disponible localement et en CI, sans
empêcher la publication des corrections du moteur principal.

### Artefacts Pmetrics

Un artefact JSON versionné associe la structure Pmetrics, son erreur polynomiale et les points de support pondérés. Le schéma est dans `pmetrics/schema/` et les artefacts livrés avec l'application dans `pmetrics/artifacts/`. Le serveur n'exécute pas de code R provenant d'un téléversement : seuls les gabarits structurels explicitement pris en charge et les artefacts du dépôt sont acceptés.

Le backend convertit les administrations et observations en `PM_data`, transmet les points de support à `PM_model$map()` par un fichier temporaire, puis extrait `postPoints`. La prédiction centrale utilise la moyenne postérieure des paramètres et les simulations d'incertitude échantillonnent les points de support avec leurs probabilités postérieures. Le dossier temporaire Pmetrics est supprimé à la fin de chaque appel; le pont mrgsolve compilé reste limité à la session et est supprimé à sa fermeture.

L'artefact `demo_one_comp_iv` est synthétique et sert uniquement au test technique. Il ne représente pas un modèle publié ou validé cliniquement.

## Démarrage local

Depuis la racine du dépôt:

```powershell
npm run tdm:index
npm run tdm:run
```

Ouvrir ensuite `http://127.0.0.1:3838/`.

Le script local active `ALLOW_CUSTOM_MODELS=true`, ce qui autorise le collage et la compilation de C++. Ne pas activer cette option sur un serveur public non isolé.

## Déploiement Shiny

Depuis la racine du dépôt, vérifier d'abord le manifeste public sans déployer :

```powershell
Rscript scripts/deploy_shiny.R
```

Le manifeste doit contenir `mrgsolve` et `mapbayr`, mais pas `Pmetrics`. Après
revue, `Rscript scripts/deploy_shiny.R --deploy` publie ce profil compatible
Shinyapps.io avec `ENABLE_PMETRICS=false`. Le mode
`--with-pmetrics` est réservé à un hébergement dont la chaîne R/Rust est
compatible ; il ne doit pas être utilisé sur Shinyapps.io actuellement.

## Fonctionnalités

### Ateliers Locaux Interactions et Pharmacodynamie

Ces ateliers sont des prototypes de recherche intégrés au moteur. Accès direct : `?view=ddi` et `?view=pd`. Sur le serveur public, les modèles de bibliothèque et Lego contrôlés sont disponibles ; le C++ arbitraire reste désactivé.
Le port du lanceur local est configurable avec `TDM_PORT` (3838 par défaut).

- Interactions : bibliothèque ou collage mrgsolve pour chaque molécule, instantanés TDM indépendants, paramètre cible. Sept relations : facteur constant, inhibition/stimulation Emax instantanées, inhibition réversible Ki, inhibition Hill, inhibition dépendante du temps (kinact/KI), induction avec renouvellement enzymatique (Emax/EC50/kdeg).
- Export DDI : C++ du modèle affecté avec entrée externe `DDI_CP`, et script R autonome simulant les deux modèles. Le C++ seul ne contient pas la PK de la molécule 2. Les paramètres individuels sont conservés dans le script R, pas dans les valeurs par défaut du C++.
- Pharmacodynamie : effet linéaire, Emax/Hill, compartiment d'effet et quatre réponses indirectes. PK mrgsolve libre, profil exponentiel ou concentrations connues; reprise possible de l'exposition TDM de la même session.
- Ajustement PD : moindres carrés individuels via `minpack.lm`; ODE via `deSolve`. Sélection des paramètres estimés, diagnostics, données synthétiques, CSV et rapport HTML; exports mrgsolve et R.

Installer au besoin `install.packages(c("deSolve", "minpack.lm"))`. Les concentrations et paramètres doivent utiliser des unités cohérentes. Les mécanismes DDI restent des hypothèses utilisateur; les ajustements PD ne sont ni une analyse de population ni une validation clinique.

Les pages du site sont des ateliers autonomes avec blocs, paramètres numériques et explications FR/EN, sans iframe. « Ouvrir dans le moteur » transmet une spécification versionnée par `postMessage`, avec contrôle de l'origine, de la fenêtre et validation R. L'accusé de réception vient du serveur après validation. Le C++ collé n'est jamais compilé automatiquement à la réception. Aucun code, paramètre individuel ou historique n'est placé dans l'URL ou le stockage navigateur.

La reprise TDM nécessite d'effectuer le TDM dans cette même session, via les onglets du moteur. Aucun dossier patient n'est transféré automatiquement entre fenêtres. Un téléchargement explicite peut contenir le profil d'exposition ou les paramètres individuels sélectionnés : son stockage relève de l'utilisateur.

Les mécanismes dynamiques DDI utilisent A(0)=1, puis la solution exacte du renouvellement enzymatique par pas de concentration constante sur la grille de 0,1 h. La récupération continue après l'arrêt du traitement. Le C++ exporté attend une entrée externe `DDI_ACTIVITY` pour ces deux mécanismes : il ne contient pas lui-même l'ODE enzymatique. Le script R autonome calcule l'activité puis simule les deux PK. Le plancher inhibiteur est de 1 %. Pas de fraction métabolisée, de conversion libre/totale, de PBPK ni de couplage réciproque. Référence méthodologique : [ICH M12](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/m12-drug-interaction-studies).

### Parcours Oncologie

- PK IV d'exemple ou PK mrgsolve libre (bibliothèque, code personnel, Lego), avec voie et compartiment d'administration explicites. Les cycles restent en jours; le code PK peut utiliser heures ou jours. Les doses utilisent l'unité du modèle et un facteur explicite convertit les concentrations PK vers les unités PD (mg/L en oncologie). Les paramètres et covariables PK restent fixes pendant l'ajustement PD. Ce n'est pas encore une bibliothèque anticancéreuse validée.
- Croissance exponentielle, logistique ou Gompertz; destruction Emax avec perte empirique de sensibilité `exp(-RES*t)`. Taille exprimée en SLD (mm), pas en volume.
- Myélosuppression optionnelle : prolifération, trois transits, ANC et rétrocontrôle, `ktr=4/MTT`. Structure adaptée de [Friberg 2002](https://doi.org/10.1200/JCO.2002.02.140), avec inhibition linéaire bornée à 1; tous les états initiaux valent ANC0. Les paramètres par défaut sont illustratifs, non tirés d'un protocole de médicament.
- Observations CSV `time,endpoint,value`, endpoints `tumor`/`anc`, ajustement individuel de 1 à 5 paramètres par moindres carrés pondérés sur log(valeur). Écarts-types fixes et explicites par endpoint. Pas de MAP-BE ni d'incertitude prédictive.
- Toutes les observations utilisées précèdent ou coïncident avec la décision. L'historique réel des doses reste identique entre scénarios. Seules la prochaine dose et les suivantes changent, sans réinitialisation des états : fraction de dose, report et intervalle. Un historique vide permet une simulation depuis le début du traitement.
- Comparaison déterministe du ratio tumoral final, nadir ANC futur et temps sous un plancher exploratoire; aucune recommandation ni classement clinique. Résolution des métriques : 0,05 jour, sans classification RECIST ni survie.
- Exports CSV, C++ (PK/PD pour l'exemple IV; PD à entrée CP externe pour une PK libre), R autonome et rapport HTML. Les données restent en session sauf téléchargement explicite.

Tests supplémentaires : `Rscript tdm-engine/tests/ddi_mechanisms_test.R`, `Rscript tdm-engine/tests/oncology_test.R`, `Rscript tdm-engine/tests/onco_server_test.R`. Ils couvrent notamment les solutions analytiques DDI, les six exports C++ oncologiques, la récupération de paramètres sur données synthétiques et la conservation du passé entre scénarios.

Les pages statiques comprennent aussi des figures interactives : relation DDI à l'équilibre et réponse à une exposition constante avec récupération, comparaison tumorale avec/sans traitement, relation PD à l'équilibre et délai Ce. Les aperçus locaux sont pédagogiques, pas des prédictions patient. Les équations utilisent les paramètres de l'atelier; les hypothèses simplificatrices sont indiquées à côté des figures. Test : `node scripts/test_workshop_curves.mjs`.

Le graphique tumoral est maintenant unique : contrôle sans traitement en pointillés et traitement en continu. Le site calcule un aperçu pour la PK IV d'exemple; pour une PK libre, il attend le résultat R transmis par le moteur après comparaison, sans substituer une PK simplifiée. Le contrôle sans traitement exclut toutes les doses depuis t=0; les deux scénarios traités partagent leur historique jusqu'à la décision. Les conditions initiales, paramètres PK et paramètres PD sont séparés dans les interfaces.

Les ateliers JSON exportés par le site sont importables dans les onglets DDI et PD du moteur (450 kB maximum). L'import emploie la même validation que le pont `postMessage`, supprime le fichier temporaire et ne compile jamais automatiquement du C++ personnel. Les observations PD et oncologiques utilisent des tableaux éditables, avec ajout/suppression et import CSV. Une ligne ajoutée a une valeur vide, à compléter avant ajustement. Aucun stockage automatique hors session.

Tests PK libre : `Rscript tdm-engine/tests/pd_pk_test.R` (PK orale à deux compartiments, conversions heures/jours et concentrations, contrôle sans traitement, conservation du passé). Le code C++ oncologique exporté en mode PK libre est un module PD à entrée CP externe; le script R exporté compile et couple la PK sélectionnée. En PD générale, le script R contient le profil d'exposition calculé, figé pour les conditions exportées.

Tests : `Rscript tdm-engine/tests/pd_engine_test.R`, `Rscript tdm-engine/tests/ddi_builder_test.R` et `Rscript tdm-engine/tests/ddi_server_test.R`.

### TDM

- sélection d'un modèle de la bibliothèque ou import/collage d'un modèle produit par l'Atelier Lego;
- compilation C++ libre uniquement en exécution locale ou sur une infrastructure isolée;
- historique de plusieurs administrations et concentrations;
- saisie en heures relatives, en jours avant une date de référence ou en dates/heures calendaires;
- administrations répétées par nombre de doses ou directement à l'état stationnaire avec `ss = 1`;
- sélection explicite de la voie IV ou orale selon les compartiments exécutables du modèle;
- covariables générées depuis `$PARAM @covariates` et saisies à l'heure de chaque prélèvement;
- estimation MAP bayésienne avec `mapbayr`;
- model averaging AIC ou log-vraisemblance, limité à une même molécule, une même voie et un même mode d'administration;
- affichage de l'AUC0-24 et de la C0 actuelles à partir des paramètres postérieurs;
- simulation de profils individualisés;
- comparaison de la poursuite de la dernière posologie avec l'application de la recommandation sur plusieurs doses futures;
- classement de scénarios par AUC24, Cmin, Cmax, pourcentage d'un intervalle posologique au-dessus d'une concentration ou contraintes simultanées et indépendantes sur Cmin et Cmax (intervalle, minimum ou maximum);
- distribution prédictive Monte Carlo configurable avec variabilité interindividuelle et erreur résiduelle;
- export CSV de la grille posologique;
- génération d'un rapport HTML autonome sans identifiant patient;
- import et export JSON des administrations, observations, covariables et cibles.

## Confidentialité et absence de persistance

- aucun champ d'identité patient n'est demandé ni exporté;
- le JSON importé est limité à 1 Mo, validé puis supprimé immédiatement du dossier temporaire;
- les données restent dans la mémoire de la session Shiny et sont libérées à sa fermeture;
- le code C++ personnalisé n'est jamais inclus dans l'export patient;
- aucune donnée patient n'est placée dans une URL partageable, un stockage navigateur ou une base de données;
- les modèles Lego sont régénérés côté serveur puis compilés dans un dossier propre à la session, déchargés et supprimés à sa fermeture;
- les modèles publiés sont livrés avec l'application et aucun modèle soumis par un visiteur n'est ajouté automatiquement au serveur.

## Flux Atelier Lego vers TDM

L'atelier Lego génère un modèle mrgsolve compatible avec le contrat `mapbayr`: tags `[ADM]` et `[OBS]`, effets aléatoires, `OMEGA`, `SIGMA` et sortie `DV`. Il prend en charge les entrées instantanées ou d'ordre zéro, `Tlag`, les voies parallèles avec fractions de dose, les transferts de Michaelis-Menten et l'élimination paramétrée par `k` ou `CL`. Il permet aussi d'ajouter jusqu'à dix covariables continues ou catégorielles simples. Une covariable continue suit `paramètre × (covariable / référence)^β`; une covariable catégorielle suit `paramètre × exp(β)` pour la modalité comparée et conserve la valeur typique pour la référence. Aucune expression C++ libre n'est acceptée par ce constructeur.

Le code généré contient une spécification JSON versionnée. L'action **Ouvrir dans TDM** l'envoie à la fenêtre Shiny avec `postMessage`; un copier-coller du code complet conserve aussi cette spécification. Le serveur valide les types, identifiants, bornes, compartiments et transferts, puis régénère lui-même un code mrgsolve équivalent. Il ne compile jamais directement le texte C++ reçu.

Le moteur ouvre ensuite le mode Atelier Lego / C++ et demande toujours une validation explicite avant l'analyse. Les administrations, concentrations observées et autres données du patient sont saisies uniquement dans le TDM. Pour une dose répartie par le modèle Lego, le moteur approche l'état stationnaire par 50 administrations de préchauffage plutôt que par des événements `ss = 1` séparés sur chaque voie.

## Déploiement

Depuis la racine du dépôt, `Rscript scripts/deploy_shiny.R` vérifie la liste des
fichiers publics et les dépendances. Ajouter `--deploy` publie cette liste sur
`https://tdmhub.shinyapps.io/MIPD_Engine/`. Le paquet exclut les PDF, les exemples
locaux, les tests et les fichiers non suivis par Git.

GitHub Pages ne peut pas exécuter R. Déployer ce dossier séparément sur une infrastructure Shiny ou dans un conteneur R, puis fournir l'URL au build Svelte:

```powershell
$env:PUBLIC_TDM_ENGINE_URL = "https://tdm.votre-domaine.fr"
npm run build
```

Pour un serveur public:

```text
ALLOW_CUSTOM_MODELS=false
```

Avec `ALLOW_CUSTOM_MODELS=false`, les modèles Lego contrôlés restent utilisables publiquement, mais tout C++ arbitraire est refusé. La compilation libre exige un service isolé: conteneur éphémère, système de fichiers en lecture seule, absence de réseau sortant, limites CPU/mémoire/temps et aucune donnée patient persistante.

## Test d'intégration

```powershell
Rscript tdm-engine/tests/smoke_test.R
Rscript tdm-engine/tests/safe_lego_test.R
npm run tdm:validate-library
```

Le smoke test estime trois modèles vancomycine avec des covariables variant dans le temps, calcule l'AUC0-24/C0 actuelles, leurs poids AIC, compare les doses futures et vérifie la distribution prédictive. Le test Lego compile un modèle contrôlé en mode public et vérifie le refus du C++ libre et des spécifications invalides. La validation de bibliothèque compile les 46 fichiers et vérifie le contrat `mapbayr` de chacun.
