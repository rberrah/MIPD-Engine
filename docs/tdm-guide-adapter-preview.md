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
