# MIPD Shiny — dossiers TDM JSON par molécule

## Source de vérité et couverture

`integration/tdm-drugs/<drugKey>.json` contient **un dossier pour chacune des 23 molécules** de la bibliothèque, indexant les **46 modèles**. Il expose notamment :

- `modelIds`, `concentrationUnit` et `doseUnit` (reconcilés avec `models/catalog.json`) ;
- `monitoring`, `sourceReferences` et `conditionalPdObjectives` extraits du TDM Guide PMx fusionné ;
- `guideTargetPresets` par mode, avec `target_metric`, `target_low`, `target_high`, `target_unit`, contexte, référence et lien ;
- `simulationPresets` : les cinq anciens préréglages sourcés (vancomycine, linézolide, voriconazole) ne sont plus codés en dur dans R ;
- `guideTargetPresets[mode].doseGrid` : bornes `dose_min`, `dose_max`, `dose_step`, `intervals`, `infusion`, spécialité et RCP quand une grille a été recherchée ;
- `guideHref` : navigation vers la molécule du TDM Guide, sans donnée patient.

Des champs restent volontairement vides pour les molécules dépourvues de référence numérique transposable ou de grille documentée. **Ce ne sont pas des doses « universelles » manquantes à deviner.**

## Utilisation dans Shiny

Le modèle choisi détermine le dossier JSON. L'app affiche l'indice PK/PD, les unités, les sources et les cibles contextualisées. Les cinq préréglages préexistants restent chargés automatiquement depuis JSON. Pour les autres molécules, les anciens champs **ne sont plus réutilisés** lors d'un changement de modèle ; il faut choisir explicitement une référence TDM ou saisir des bornes de simulation.

Une grille issue d'un RCP est **documentaire tant que sa spécialité n'est pas confirmée**. Elle n'est injectée dans les doses proposées **que sur un clic explicite**, après contrôle de concordance modèle–formulation. Le statut scientifique reste `PENDING` et la simulation n'est pas une prescription. Les modèles dont la spécialité est incompatible sont bloqués.

### Tacrolimus

Le modèle Woillard 2011 utilise `ST = 1` pour **Prograf immédiat 12 h**, `ST = 0` pour **Advagraf LP 24 h** (article original). Kamp 2024 correspond à **Envarsus LP 24 h**, qui **n'est pas interchangeable** avec Advagraf ou Prograf. La grille RCP Prograf `0,5–14 mg/prise, pas 0,5 mg, intervalle 12 h` est seulement un espace exploratoire, et elle n'est accessible que si `ST = 1` et Prograf a été confirmé. Aucune grille chiffrée n'est inventée pour Advagraf/Envarsus.

La concentration calculée est **déjà en ng/mL** : `CENT_TAC [mg] / V1 [L] × 1000`. Une division supplémentaire par 1 000 produirait des valeurs fausses. Les cibles C₀ du guide sont liées à la greffe, au délai post-transplantation, au traitement associé, à la formulation et à la méthode de dosage ; l'AUC peut être préférable pour apprécier l'exposition, mais n'est pas équivalente à une AUC24 unique.

## Requalification et déploiement

`node tests/tdm_drug_profiles.test.mjs` et `Rscript tests/tdm_drug_profiles_test.R` vérifient couverture, références, unités, sécurité et séparation des formulations. La suite CI sur GitHub doit rester entièrement verte.

**Une fusion du dépôt ne met pas Shinyapps à jour.** Le site public exige un déploiement rsconnect autorisé et la validation du SHA publié. La bibliothèque est destinée à la recherche et à l'enseignement, sans validation pour une décision de dose clinique.
