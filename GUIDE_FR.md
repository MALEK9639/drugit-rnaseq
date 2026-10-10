# Guide pas à pas (semaines 1 et 2)

Tout ce qui est marqué **[À FAIRE PAR TOI]** ne peut pas être fait à ta place (installation sur ton PC, compte GitHub, envoi du mail).

## Semaine 1

**Jour 1: environnement [À FAIRE PAR TOI]**
1. Installer R puis RStudio, puis Git.
2. Créer un compte GitHub, puis un dépôt `drugit-rnaseq`.
3. Dans un terminal, depuis ce dossier :
   ```
   git init
   git add .
   git commit -m "Initial project structure and week 1 documents"
   git branch -M main
   git remote add origin https://github.com/<ton-compte>/drugit-rnaseq.git
   git push -u origin main
   ```
4. Dans le README, remplir les dates de stage.

**Jour 2: théorie.** Lis `docs/week1_rnaseq_basics_note.md`, puis **réécris-le avec tes mots** (consigne du stage).

**Jours 3-4: choix du dataset.** `metadata/dataset_comparison.md` contient 3 candidats + 1 rejeté. Les cases **[TO VERIFY]** sont à vérifier sur la page GEO (accession, nombre d'échantillons, présence de la matrice de comptes bruts dans la section Download). Je n'ai pas pu ouvrir les pages GEO automatiquement.

**Jour 5: approbation.** Complète `docs/dataset_selection_form.md` (DOI de la publication, date), puis envoie `docs/email_dataset_approval.md` **[À FAIRE PAR TOI]**. Ne lance l'analyse principale (semaine 3) qu'après réponse.

## Semaine 2

1. `source("scripts/00_install_packages.R")` (une seule fois).
2. `source("scripts/01_download_or_import_data.R")`: télécharge comptes et métadonnées. Si le téléchargement échoue, prends le fichier à la main sur GEO et mets-le dans `data/raw/`.
3. `source("scripts/02_clean_metadata.R")`: **lis les tables affichées** (classes de tissu, patients). Si l'étiquetage automatique est faux, corrige les règles de la section 1.
4. `source("scripts/03_quality_control.R")`: produit 3 figures et `results/week2_qc_summary.md`.
5. Regarde chaque figure et écris une phrase d'interprétation. Ne supprime aucun échantillon sans l'accord écrit du superviseur.
6. Remplis `docs/week2_update_email_template.md` avec tes vrais résultats et envoie-le.

## Si une erreur apparaît
Utilise le format du PDF: Goal / What I tried / Code / Expected / Observed / Full error / Hypothesis. Colle-moi l'erreur complète et je t'aide à corriger.
