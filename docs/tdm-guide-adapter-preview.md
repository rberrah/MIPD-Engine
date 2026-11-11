# DRAFT — Recommandations STP automatiques par modèle MIPD (46/46)

**NE PAS FUSIONNER, DÉPLOYER OU ACTIVER dans Shiny.** Cette PR se limite à un adaptateur purement local et à ses tests.

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
# Dans PMx Explain (branche de la PR #235) :
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
