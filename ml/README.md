# Entrainement PopPK + XGBoost de l'AUC24

Ce dossier ne contient aucune donnee patient. Le pipeline entraine un predicteur XGBoost distinct pour chaque couple modele mrgsolve / mode d'administration a partir de profils entierement simules. Il adapte a l'ensemble de la bibliotheque le principe valide pour le tacrolimus par Woillard et al. (doi:10.1016/j.phrs.2021.105578), sans extrapoler cette validation clinique aux autres molecules.

## Benchmark apparie MAP-BE / XGBoost

Le comparatif public utilise `--paired-benchmark`, implemente dans
`paired_auc_benchmark.R`. Pour chaque patient du test, `mapbayr` ajuste les ETA sur les **memes
concentrations observees** que celles donnees a XGBoost : C0 seul, puis C0+C1.

- C1 = Tinf + 1 h pour une perfusion intermittente.
- Pour l'oral et le bolus, Tinf de reference = 0 : C1 = 1 h apres la dose.
- Pour la perfusion continue, la pompe reste active a l'etat stationnaire et les
  prelevements sont effectues a t=0 et t=1 h. Aucun bolus fictif n'est ajoute.
- Hors perfusion continue, C0 est simule a tau moins 1 microseconde : il s'agit
  du creux pre-dose du cycle equivalent a l'etat stationnaire, pas du pic IV a t=0.
- L'AUC24 vraie est integree par trapezes sur le profil individuel sans bruit
  (pas de 0.1 h). Pour les dosages limites, la structure SIGMA publiee est
  conservee mais chaque ecart-type residuel est divise par 10. Ce scenario
  rapproche une erreur proportionnelle publiee de 10 % du scenario a 1 % etudie
  par Berrah et al., sans imposer la meme erreur additive a des unites differentes.
- Les ETA suivent la matrice OMEGA sans troncature a deux ecarts-types. Les
  covariables suivent `training-populations.json`; elles restent constantes dans
  chaque profil. Les 5 % d'AUC24 vraies les plus basses et les 5 % les plus
  hautes sont retirees separement pour chaque modele et mode d'administration.
- Pour Woillard, la formulation est liee a l'intervalle : `ST=1` pour Prograf
  q12h et `ST=0` pour Advagraf q24h.
- XGBoost apprend directement l'AUC24, sans cible logarithmique ni AUC
  populationnelle parmi ses entrees. Dose, intervalle, concentrations, horaires
  et covariables sont disponibles. C1 et ses differences sont absents du modele C0.
- Le partage 75/25 des patients est commun aux deux plans et intervient apres le
  filtre. La validation croisee a dix groupes est limitee aux 75 % d'apprentissage.
  ESM4 de Woillard publie le meilleur reglage a deux dosages (`tree_depth=1`,
  `min_n=40`, `mtry=7`, `learning rate=0.0261`), mais pas la grille complete.
  La grille predefinie pour un ou deux dosages inclut ce point et teste
  `max_depth={1,2,4}` avec `eta={0.0261,0.05}`, `min_child_weight=40`,
  `mtry=min(7,p)`, jusqu'a 1000 tours avec arret anticipe apres 30 tours. La
  profondeur 9 publiee pour le modele a trois dosages n'est pas transposee a ces
  plans plus parcimonieux. La selection minimise la RMSE absolue.

Avec `e = (AUC estimee - AUC vraie) / AUC vraie` : biais = `100*mean(e)`,
RMSE relative = `100*sqrt(mean(e^2))`, AUC bien predites = proportion avec
`abs(e) <= 0.20`. La page utilise les memes patients estimables par les deux
methodes pour les biais/RMSE; le taux a +/-20 % porte sur tout le jeu de test,
les echecs comptant comme mal predits. Les effectifs et echecs sont affiches.

Depuis `tdm-engine` :

```powershell
Rscript tests/paired_benchmark_test.R
Rscript ml/train_models_xgboost.R --paired-benchmark --smoke --workers=4 --report=../.svelte-kit/paired-smoke.json
Rscript ml/train_models_xgboost.R --paired-benchmark --n=9000 --workers=4 --report=ml/validation/paired-auc-benchmark.json
```

Avec `--n=9000`, Woillard conserve 9 000 profils simules avant filtrage pour
rester comparable au protocole source. Les autres couples modele/mode utilisent
3 000 profils avant filtrage. Apres le filtre et le partage, le pool de test
contient environ 2 025 profils Woillard et 675 profils par autre couple. Un
sous-echantillon fixe et non utilise pour le reglage fournit le comparatif
apparie : 1 098 profils pour Woillard, comme le test simule de l'article, et 300
pour chaque autre couple. Cette precision est affichee dans le rapport et limite
le cout des deux ajustements MAP individuels par profil.

`--base`, `--drug`, `--mode` et `--seed` restent disponibles. Les checkpoints
synthetiques sont locaux et ignores par Git; leurs empreintes dependent des
scripts, modeles et configurations. Seul le JSON agrege est publie sur le site,
jamais les observations simulees individuelles.

Cette extension reprend l'apprentissage direct, le partage et le principe de
selection des hyperparametres de Woillard 2021, pas son protocole integral.
L'article utilisait uniquement Prograf q12h, neuf doses precises, C0+C3 ou
C0+C1+C3, l'AUC0-12, des horaires legerement perturbes, une erreur additive de
0.1 microgramme/L avec une erreur proportionnelle de 0.01 %, ainsi que des
filtres propres au tacrolimus. Le benchmark utilise C0 ou C0+C1, l'AUC24, un
filtre central generique, Prograf q12h et Advagraf q24h, et des distributions de
covariables propres a chaque modele. Les scores ne sont donc pas censes etre
identiques. Les predicteurs du benchmark sont distincts des artefacts Shiny
historiques ci-dessous; leurs scores ne sont pas une validation de ces artefacts.

## Artefacts Shiny historiques : principe

Pour chaque patient virtuel, le script :

1. tire les covariables et effets aleatoires dans le domaine defini pour le modele;
2. simule un schema a l'etat stationnaire avec mrgsolve et l'erreur residuelle du modele;
3. simule un C0 pre-dose et une concentration une heure apres la fin de perfusion; pour une administration orale ou un bolus, la fin de perfusion est fixee a 0;
4. calcule l'AUC24 individuelle vraie et l'AUC24 populationnelle du meme schema;
5. entraine XGBoost sur `log(AUC24 vraie / AUC24 populationnelle)`;
6. reconstruit l'estimation par `AUC24 populationnelle * exp(prediction)`;
7. evalue le modele par validation croisee repetee, jeu de test interne non touche et, lorsqu'un autre modele compatible existe, transportabilite PopPK simulee;
8. produit un fond synthetique pour l'explication locale DALEX.

Cet ancien controle interne rapportait l'ecart entre l'AUC individuelle vraie et l'AUC obtenue avec les effets aleatoires fixes a zero (`POP_AUC24`), avant et apres correction par XGBoost. Il est conserve uniquement pour qualifier les artefacts historiques de l'application. Ce n'est ni le comparateur du benchmark public ni une estimation intitulee "AUC sans ML" : la comparaison scientifique publiee sur cette page est MAP-BE versus XGBoost, avec les memes patients et les memes dosages.

Les variables comprennent les covariables du modele, la dose, l'intervalle, la duree de perfusion, les horaires et concentrations, les predictions populationnelles correspondantes et les rapports observe/predit. Les domaines de dose et d'intervalle sont explicites dans `training-regimens.json`.

La validation historique des artefacts compare, sur les memes patients virtuels et le meme jeu de test, deux plans : `C0 seul` et `C0 + C(fin de perfusion + 1 h)`. Pour l'oral et le bolus, le second temps vaut donc 1 h apres la dose. Pour une perfusion intermittente de duree `Tinf`, il vaut `Tinf + 1 h`. Le C0 est simule juste avant la dose suivante. `zero_re(omega)` annule uniquement l'OMEGA interne de mrgsolve; les ETA sont fournis par `idata` et la SIGMA du modele reste appliquee aux concentrations.

La perfusion continue est exclue de cette ancienne comparaison d'artefacts. Elle est incluse dans le benchmark apparie ci-dessus selon la convention demandee : perfusion maintenue, prelevements a t=0 et t=1 h.

Cette ancienne validation ne reproduit pas le protocole de Woillard et al. L'article tacrolimus simulait 9 000 profils riches a l'etat stationnaire, avec neuf doses administrees toutes les 12 h, une concentration toutes les 30 minutes et une AUC0-12 cible. Le modele a deux prelevements utilisait C0 et C3; le modele a trois prelevements utilisait C0, C1 et C3. Dans le jeu de test simule, les RMSE relatives publiees etaient respectivement de 4.60 % et 2.61 %. Les artefacts historiques utilisent une AUC24, 1 000 profils par artefact et C0+C1 pour l'oral; leurs resultats ne sont donc pas une reproduction de l'article.

Les covariables ne sont jamais tirees dans une plage generique commune a toute la bibliotheque. `training-populations.json` contient les informations propres aux populations sources : bornes, moyenne/ecart-type ou proportions publiees. Une variable positive documentee uniquement par sa moyenne et son ecart-type peut utiliser une loi log-normale ajustee sur ces deux moments; cette hypothese est inscrite dans le fichier. Quand l'article ne documente pas suffisamment la population, la valeur de reference du modele est conservee; le pipeline n'invente pas une plage. Les covariables liees peuvent etre reconstruites conjointement (par exemple la creatinine a partir d'une clairance Cockcroft-Gault cible).

Les covariables dependantes du temps sont maintenues constantes pendant chaque profil simule. Leur trajectoire longitudinale n'est donc pas encore validee par ce benchmark.

Dans l'application, l'AUC24 ML reste separee de l'estimation MAP-BE. Elle ne modifie ni les trajectoires, ni les simulations de doses, ni la recommandation MAP-BE.

## Model averaging

Il n'existe pas d'artefact supplementaire propre au model averaging. Chaque modele selectionne produit son AUC24 ML, puis l'application applique les memes poids d'averaging que pour l'analyse pharmacometrique. L'agregation est disponible uniquement si tous les modeles ajustes disposent d'un artefact compatible et partagent la meme molecule, la meme voie et le meme mode d'administration.

## Utilisation

Depuis `tdm-engine` :

```powershell
# Verification rapide sans publier
Rscript ml/train_models_xgboost.R --smoke --base=all

# Evaluation d'un modele ou d'une molecule
Rscript ml/train_models_xgboost.R --n=1000 --base=vanco_pkjust --mode=IV_INTERMITTENT
Rscript ml/train_models_xgboost.R --n=1000 --drug=Vancomycine

# Entrainement et publication de toute la bibliotheque
Rscript ml/train_models_xgboost.R --n=1000 --base=all --publish --report=ml/validation/all-models.csv
```

Options disponibles :

- `--base=all` ou une liste d'identifiants separes par des virgules;
- `--drug=all`, une cle de molecule ou son nom;
- `--mode=all`, `ORAL`, `IV_INTERMITTENT` ou `IV_CONTINUOUS`;
- `--n=1000` pour l'effectif par couple modele/mode;
- `--seed=20260906` pour reproduire un entrainement;
- `--publish` pour ecrire les RDS et mettre a jour `registry.json`;
- `--report=...csv` pour conserver les metriques synthetiques.

La publication est refusee en mode `--smoke` ou avec moins de 1 000 profils par couple modele/mode. Pour ajouter un modele, il faut d'abord l'ajouter au catalogue, definir chaque schema d'administration pris en charge dans `training-regimens.json` et documenter sa population dans `training-populations.json`, puis relancer le script sur son identifiant. Aucun fichier patient n'est lu ou ecrit.

## Contrat et niveaux de preuve

Chaque artefact est lie a l'identifiant du modele, la molecule, la voie, le mode d'administration, le schema de variables et l'empreinte SHA-256 exacte du fichier mrgsolve. Le manifeste enregistre egalement la graine, l'effectif, les hyperparametres, les versions logicielles, le domaine d'entrainement et les metriques.

- `experimental` : artefact evalue en interne mais au moins un seuil de performance prespecifie n'est pas atteint;
- `research` : validation croisee repetee et jeu de test interne conformes aux seuils;
- validation clinique : toujours absente tant qu'une validation favorable sur des patients reels independants de la molecule concernee n'est pas documentee.

Une valeur hors du domaine empirique declenche un avertissement sans bloquer l'affichage. Pour l'artefact principal, un C0 pre-dose et une concentration post-dose dans l'intervalle suivant sont requis. L'absence d'etat stationnaire, un mode d'administration incompatible, une empreinte differente ou une AUC invalide restent bloquants. L'explication DALEX utilise exclusivement un echantillon synthetique stocke separement avec sa propre empreinte.

## Fichiers publies

Pour chaque couple modele/mode, `--publish` produit :

```text
artifacts/<modele>-<mode>-c0-post-auc24-xgb-v4.rds
artifacts/<modele>-<mode>-c0-post-auc24-xgb-v4-dalex-background.rds
```

Le premier RDS contient le booster XGBoost. Le second contient au maximum 200 profils synthetiques servant de reference DALEX. `registry.json` est le seul index charge par l'application.
