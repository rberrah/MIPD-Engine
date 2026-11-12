# DRAFT — Recommandations STP automatiques par modèle MIPD (46/46)

**Intégration du code uniquement (opt-in, recherche/enseignement).** La fusion dans `main` n'autorise ni déploiement Shinyapps ni activation clinique. Les cibles et grilles restent non applicables sans validation explicite des conditions de population, méthode analytique et formulation.

## Accès direct au TDM Guide depuis Shiny

Dans **Cible et grille de doses**, ainsi que sur chaque fiche de modèle, un lien
ouvre la page PMx Explain directement filtrée sur la molécule sélectionnée :
`https://rberrah.github.io/pharmacometrie/cibles-tdm/?drug=tac&lang=fr`
pour le tacrolimus. L'URL contient uniquement l'identifiant canonique de la
molécule et la langue, **jamais** de dose, de cible appliquée ni de donnée patient.
Le lien est proposé même sans préréglage chiffré, mais pas pour les cinq molécules
encore absentes du guide. La source PMx est désormais issue de la PR #237 fusionnée sur `main`
(commit `fc3fe8729f852033eeb5858b05e7a089de5071db`) ; ne considérer la navigation filtrée comme
publique qu'après déploiement vérifié des Pages. Aucune application de cible
ou de grille n'est déclenchée
par l'ouverture du lien.

## Deux sources autonomes

1. PMx Explain / TDM Guide produit un snapshot `pmx.tdm.model-references.preview.v1` avec source SHA256, règles STP, plages scientifiques, objectifs PK/PD et provenance bibliographique. **Une molécule, plusieurs références possibles**.
2. MIPD Engine possède `models/catalog.json` et les 46 modèles PK ; ce fichier reste sa source de vérité. `integration/tdm-guide-drug-crosswalk.preview.v1.json` traduit seulement les **23 clés de molécules** en identifiants canoniques du guide. Aucun modèle n'est copié dans PMx.
3. `node scripts/preview_tdm_guide_recommendations.mjs --snapshot /tmp/pmx-tdm-preview.json --output /tmp/mipd-tdm-preview.json` génère **exactement une fiche par modèle**, automatiquement. Pas d'accès réseau, aucun enregistrement patient ni calcul de dose.

## Résultat attendu pour chaque modèle

- Identification du modèle, molécule, population, voie/mode, statut simulation vs MAP.
- Critère privilégié du guide (AUC, résiduelle, pic, fT>CMI, etc.), prélèvement, source et date implicite par SHA.
- Références numériques originales possibles **sans choisir la « meilleure » arbitrairement**.
- Objectifs PK/PD dépendant de la CMI, distinguant fraction libre et totale. Pour perfusion continue, une règle IV intermittente est marquée `different_administration_mode`.
- Statut `published_target_candidates`, `contextual_pkpd_only`, `society_tdm_context_only`, `no_qualified_stp_reference` ou `not_indexed_in_guide`.
- Aucun modèle n'est interprétable cliniquement : `interpretationAllowed:false`, `automaticDoseAdjustmentAllowed:false`, `patientContextMatched:false`.

## Réconciliation actuelle

Le catalogue MIPD a **46 modèles PK de 23 principes actifs**. La couverture existante du guide est de **18 principes actifs** ; **5 clés MIPD non encore documentées dans le guide** : amiodarone, cobicistat, diltiazem, érythromycine et ritonavir (principaux modèles associés aux DDI). Leur présence déclenche une fiche informative « non référencé » et **jamais une cible thérapeutique inventée**.

Pour les molécules présentes dans le guide, plusieurs modèles requièrent des exclusions **propres à la population** : ciclosporine et formulations, tacrolimus IR vs Envarsus, amoxicilline adultes brûlés vs obèses vs endocardite, gentamicine dialyse vs autres, céfépime ICU/épuration, vancomycine continue vs intermittente. La cartographie par principe actif **n'implique pas** que chacune de ses cibles soit applicable à chacun de ses modèles.

## Garde-fous structurels d'applicabilité

`integration/tdm-guide-applicability-rules.preview.v1.json` déclare, de façon versionnée et relisible, les seules comparaisons que les métadonnées publiées permettent aujourd'hui. Une règle ne peut qu'**ajouter un avertissement** : elle ne déclare jamais une référence applicable, ne choisit aucune valeur et ne touche à aucune dose.

| Dimension | Règle | Effet |
| --- | --- | --- |
| Organe greffé | Pour les immunosuppresseurs déclarés (`tac`, `cyclosporine`, `everolimus`, `sirolimus`, `mpa`), l'indication du guide est comparée aux `populationTags` du modèle | `organ_indication_mismatch` affiché sur la référence et compté dans la bannière |
| Mode d'administration | Une référence de type résiduelle sur un modèle **exclusivement** en perfusion continue | `trough_reference_under_continuous_infusion_only_model` |
| Identité analytique | `analyteId`, `analyteBasis`, `matrixId`, `formulationId`, `assayMethod`, `samplingTime` absents du guide | listés comme **non vérifiables**, jamais comme concordants |

Le périmètre est volontairement restreint : d'autres molécules réutilisent les mêmes codes d'indication avec un autre sens (pour la digoxine, `heart` désigne l'insuffisance cardiaque et non une greffe). Un organe non déclaré par le modèle est traité comme *non vérifiable*, pas comme une incompatibilité.

Sur le snapshot réel, **8 références** sont signalées incompatibles (cible FDA greffe hépatique présentée sur des modèles de greffe rénale pour le tacrolimus et l'évérolimus, cibles greffe cardiaque de l'évérolimus) et **10 modèles** portent un objectif PK/PD d'un autre mode d'administration. Ces nombres sont verrouillés en CI : toute variation impose une nouvelle relecture.

État des métadonnées du guide à ce jour : `assayMethod` 0/141, `formulationId` 1/141, `analyteId`/`matrixId`/`samplingTime` 10/141, et **131/141 références sans aucune des six dimensions**. C'est la limite principale : tant que ces champs ne sont pas publiés par le guide, la distinction tacrolimus IR vs Envarsus reste **non vérifiable par machine** et demeure une relecture humaine.

## Préréglages cliniques dans « Cible et grille de doses »

Deux artefacts versionnés sont committés et lus directement par l'application, sans fichier local ni variable d'environnement :

| Fichier | Contenu | Généré par |
| --- | --- | --- |
| `integration/tdm-guide-recommendations.v1.json` | les 46 fiches STP affichées sous chaque carte de modèle | `scripts/preview_tdm_guide_recommendations.mjs` |
| `integration/tdm-guide-clinical-presets.v1.json` | les cibles publiées utilisables comme préréglage (9 clés `drugKey::mode`, 16 cibles) | `scripts/build_tdm_guide_presets.mjs` |
| `integration/tdm-guide-dose-grids.v1.json` | les grilles de doses **déclarées par MIPD**, sourcées RCP, `signoff: PENDING` | écrit à la main, relu en CI |

La CI régénère les deux artefacts depuis le commit PMx épinglé et **refuse toute dérive**.

### Ce qui peut devenir un préréglage

Une référence du guide n'alimente le panneau que si **toutes** ces conditions sont vraies, sinon elle reste de la documentation dans la fiche STP :

1. sa métrique désigne **la même grandeur** qu'une métrique de l'application (`appMetricEquivalence`). `AUC₀–₁₂` et `C₂h` sont **exclus explicitement** : une AUC sur 12 h n'est pas l'AUC sur 24 h optimisée par l'application, et une C2 n'est ni une résiduelle ni une Cmax ;
2. c'est un intervalle à deux bornes exploitables ;
3. aucune incompatibilité structurelle de contexte avec ce modèle (la cible FDA greffe hépatique du tacrolimus est donc **absente** des préréglages des modèles rénaux).

### Règle d'application

- **Plusieurs cibles qualifiées** (tacrolimus 3, ciclosporine 3, évérolimus 3, sirolimus 2) : le panneau affiche la liste, `autoApply` est faux et **rien n'est écrit** avant que le lecteur ne choisisse et ne clique.
- **Une seule cible qualifiée** (vancomycine ×2 voies, linézolide ×2 voies, voriconazole) : `autoApply` est vrai. Ces cinq clés possèdent **déjà** un préréglage curaté dans `R/clinical_presets.R` ; l'artefact n'introduit donc **aucune écriture automatique nouvelle**, et un test le vérifie.
- Le préréglage curaté MIPD **garde la priorité**. Un désaccord est affiché, jamais résolu en silence : le voriconazole en est le cas réel, guide 2–6 mg/L contre préréglage MIPD 1–5,5 mg/L (Ashbee 2014).
- Les grilles de doses ne viennent **jamais** du guide, qui ne publie aucune posologie. Elles sont déclarées par MIPD depuis les RCP (Prograf, Neoral, Certican, Rapamune) et bornent seulement l'espace de recherche.

### Deux garde-fous d'application

**1. Cible dépendante de la méthode de dosage.** `assayMethodRequiredForApplication` dans les règles liste les molécules dont l'intervalle publié n'a de sens qu'avec le dosage qui l'a produit. Le **sirolimus** y figure : le RCP Rapamune indique que les concentrations obtenues par chromatographie et par immunoessai **ne sont pas interchangeables** et que la cible doit être ajustée à la méthode employée. Le guide ne publiant pas `assayMethod`, les deux références sirolimus (4–12 et 12–20 ng/mL) restent **affichées et sourcées** mais **non applicables** : le bouton refuse et explique pourquoi. Le jour où le guide publie la méthode, la même référence devient applicable sans changer une ligne de code. Tacrolimus, ciclosporine et évérolimus rejoindront cette liste quand leur `assayMethod` sera publié.

**2. Grille liée à une formulation.** `integration/tdm-guide-model-formulations.v1.json` déclare la formulation de chaque modèle concerné, strictement d'après sa population publiée, et chaque grille déclare le périmètre qu'elle couvre (`requiresFormulationIds`). Une grille n'est appliquée à un modèle que si la formulation déclarée de ce modèle figure dans ce périmètre. **La signature seule ne suffit pas.**

Conséquence vérifiée par test : même si la grille tacrolimus passait à `APPROVED`, elle resterait **non applicable aux deux modèles tacrolimus existants** —

| Modèle | Formulation déclarée | Verdict |
| --- | --- | --- |
| `tacrolimus_kamp_envarsus_ddi` | Envarsus, libération prolongée 1×/j | hors périmètre d'une grille en deux prises |
| `tacrolimus_woillard_ddi` | « Prograf **ou** Advagraf » | formulation indéterminée, donc bloquée |
| `ciclosporine_press_ddi`, `everolimus_*`, `sirolimus_djebli_ddi` | aucune spécialité nommée par la source | bloquées |

`everolimus_terheine_ddi` mérite une mention : sa population mélange oncologie et transplantation, et l'évérolimus y est dosé à un ordre de grandeur supérieur (environ 10 mg/j contre 0,75 mg deux fois par jour). La grille Certican ne couvre pas cette population.

Une formulation non déclarée ou ambiguë **bloque** la grille au lieu de retomber sur le cas le plus courant. Un test de contrôle positif vérifie qu'en déclarant un modèle comme libération immédiate pure, la grille approuvée lui devient bien applicable : le garde-fou discrimine, il ne refuse pas tout par principe.

### Une grille non validée n'est jamais écrite

`signoff: PENDING` sur les quatre grilles, et **aucune n'est appliquée** : le bouton n'écrit que la métrique et les bornes de la cible publiée. `dose_min`, `dose_max`, `dose_step`, les intervalles testés et la durée de perfusion restent inchangés. La grille est affichée avec sa source et la mention « non appliquée ». Seul le passage à `signoff: APPROVED` autorise son écriture ; c'est le seul interrupteur, il est déclaré dans `applicationPolicy` et vérifié des deux côtés (Node et R).

### Contraintes d'étiquetage représentées

Deux bornes ont été corrigées après relecture du RCP, parce qu'une grille plausible n'est pas une grille licite :

- **Sirolimus : pas de 1 mg, pas 0,5 mg.** Le RCP Rapamune interdit d'utiliser des multiples du comprimé de 0,5 mg en remplacement des autres dosages, ce comprimé n'étant pas pleinement bioéquivalent aux 1, 2 et 5 mg sur la Cmax. Un pas de 0,5 mg aurait supposé exactement cette substitution. Le RCP ajoute que les concentrations obtenues par chromatographie et par immunoessai **ne sont pas interchangeables**, ce qui rend l'absence de `assayMethod` dans le guide d'autant plus décisive.
- **Tacrolimus : intervalle 12 h uniquement.** La posologie de 0,20–0,30 mg/kg/j citée est celle de Prograf à libération immédiate **en deux prises**. L'intervalle de 24 h a été retiré : Advagraf et Envarsus sont des formes à libération prolongée en une prise, dosées différemment, et le catalogue de modèles ne déclare pas la formulation. Le champ `excludesFormulations` le dit, et l'interface l'affiche.

Limite assumée et documentée dans `knownLimitations` : les grilles sont indexées par molécule et voie, **pas** par formulation ni par protocole de transplantation, et le poids n'est pas un paramètre. Une grille par formulation demanderait que le catalogue MIPD déclare la formulation de chaque modèle.

### Désactivation

Les fiches STP sont actives par défaut. `ENABLE_TDM_GUIDE_PREVIEW=false` les désactive entièrement ; `TDM_GUIDE_PREVIEW_JSON` permet de pointer un autre fichier local pour expérimenter.

## Passage ultérieur à une interface Shiny

Un composant **lecture seule** placé sur la fiche du modèle affichera « Indice STP conseillé / prélèvement / sources / limites du modèle ». L'utilisateur pourra ensuite préciser indication, formulation, schéma, âge, statut greffe, méthode, analyte, fraction libre et CMI. Le filtre de références restera **fail-closed** tant que toutes les dimensions nécessaires ne sont pas validées. Ne pas transformer une C0 en AUC sans méthode qualifiée.

Quand des modèles supplémentaires de Limoges seront ajoutés, CI échouera si une nouvelle `drugKey` n'est pas recensée dans le crosswalk. Ceci assure une ligne de recommandation pour **chaque nouveau modèle** sans bloquer l'ajout de plusieurs modèles d'une même molécule.

## Contrôles

`node tests/tdm_guide_adapter_preview.test.mjs` vérifie 46/46 modèles, 23 clés, 13 MAP, cinq fiches non indexées, cohérence statut/voie, absence de posologie, checksum de provenance et fail-closed à l'ajout d'une nouvelle molécule.

**Conception exclusivement**. Aucun code d'app Shiny modifié, aucun couplage runtime et aucun lien vers un résultat de patient.

## Prototype visuel R/Shiny en local (opt-in, non déployé)

L'interface Shiny dispose maintenant d'un panneau déroulant **lecture seule** par modèle, y compris les
modèles de simulation du catalogue, si et seulement si les variables suivantes sont définies
en session locale avant `shiny::runApp()` :

```bash
# Dans PMx Explain (main après fusion de la PR #237) :
node scripts/export_tdm_model_references.mjs --output /tmp/pmx-tdm-guide-preview.json
# Dans MIPD Engine (branche de la PR #39) :
node scripts/preview_tdm_guide_recommendations.mjs --snapshot /tmp/pmx-tdm-guide-preview.json --output /tmp/mipd-tdm-preview.json
ENABLE_TDM_GUIDE_PREVIEW=true TDM_GUIDE_PREVIEW_JSON=/tmp/mipd-tdm-preview.json Rscript run_local.R
```

Par défaut `ENABLE_TDM_GUIDE_PREVIEW=false` : aucun chargement du snapshot, aucun changement
visible pour les utilisateurs. Le panneau montre indice PK/PD, prélèvement, références candidates
(sources et contexte), objectifs conditionnels ; aucune plage n'est déclarée applicable au patient.
Le module refuse les fichiers incomplets, patients, références non verrouillées ou avec drapeaux
d'interprétation clinique actifs. Le panneau **ne modifie jamais les sliders de cible ou de dose**.

Tests R de présentation et de rejet de snapshots non sécurisés : `Rscript tests/tdm_guide_preview_test.R`.
Validation de cohérence, de SHA et de 46 modèles : tests Node de la PR #39.

## Contrôle de provenance historique renforcé

La CI source standalone vérifiait la présence locale d'objets Git historiques pour 51 fichiers
adaptés. Après nettoyage d'anciennes branches, sept objets sont devenus indisponibles dans le clone
MIPD **bien que leur contenu corresponde toujours, SHA Git identique, au commit PMx original
`2f97f109c3f4af257776a595f6546c9cd5f2321d`**.

La vérification conserve les empreintes immuables `migration/source-manifest.json`.
Si un objet existe dans le clone MIPD, son test initial s'applique. Sinon le test recherche
le même fichier dans un checkout de ce commit source exact et vérifie **l'égalité stricte
de son SHA Git**. Le test échoue si la source manque, si un octet diffère, ou si une entrée
n'est pas vérifiable. Aucun original n'est remplacé, ni statut scientifique modifié.


## Procédure de test local sous Windows (PowerShell)

Conserver **les deux dépôts dans deux dossiers voisins**. Les PR sont `draft` et ne sont jamais fusionnées dans `main`.

**PMx Explain, terminal 1 :**

```powershell
git fetch origin pull/235/head:tdm-preview-local
git switch tdm-preview-local
npm ci
$env:BASE_PATH="/pharmacometrie"
npm run dev
```

Ouvrir `http://localhost:5173/pharmacometrie/cibles-tdm/?lang=fr`.
Rechercher « tacrolimus » puis cliquer sur **Voir les modèles MIPD (2)**.
La page `/pharmacometrie/tdm/?drug=tacrolimus&lang=fr#modeles` doit montrer
deux modèles. Vancomycine : trois modèles. Lithium : aucun lien.
Les modèles non MAP sont toujours visibles dans la bibliothèque, mais non exécutables en MAP.

**MIPD Engine, terminal 2 :**

```powershell
git fetch origin pull/39/head:mipd-tdm-preview-local
git switch mipd-tdm-preview-local
node tests/tdm_guide_adapter_preview.test.mjs
cd ../rberrah.github.io
node scripts/export_tdm_model_references.mjs --output ../pmx-tdm-refs.preview.json
cd ../MIPD-Engine
node scripts/preview_tdm_guide_recommendations.mjs --snapshot ../pmx-tdm-refs.preview.json --output ../mipd-tdm.preview.json
$env:ENABLE_TDM_GUIDE_PREVIEW="true"
$env:TDM_GUIDE_PREVIEW_JSON=(Resolve-Path "../mipd-tdm.preview.json").Path
Rscript run_local.R
```

Ouvrir `http://127.0.0.1:3838/` : le bloc STP bibliographique apparaît sous
les cartes du catalogue et le modèle actif, **sans changer aucune cible ni dose**.
Préinstaller les dépendances R du dépôt. Adapter les chemins si les dépôts ne sont
pas voisins. Avec la variable `ENABLE_TDM_GUIDE_PREVIEW` absente, Shiny reste inchangé.

Tests complémentaires : PMx `npm run check` puis
`npx playwright test tests/e2e/tdm-targets.spec.js --config playwright.pmx.config.js` ;
MIPD `Rscript tests/tdm_guide_preview_test.R`.
