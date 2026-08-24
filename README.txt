Points d'amélioration :
---------------------

- La moyenne simple utilisée pour calculer le score de classement peut pénaliser fortement certaines colonnes. Il serait préférable de définir un score pondéré reflétant mieux la pertinence catalytique.

- Intégrer des descripteurs plus robustes, par exemple un score de conservation phylogénétique des résidus du site catalytique, afin de prioriser les candidats conservant ce site.

- Réaliser le BLAST en mode remote ou permettre le choix entre base locale et base distante.

- Remplacer ou compléter Swiss-Prot par RefSeq/NR pour augmenter la couverture des hits BLAST.

- Exécuter Boltz-2 sur GPU plutôt que sur CPU lorsque le matériel le permet.

- Séparer les descripteurs génériques applicables à tout ligand des descripteurs spécifiques à la nature chimique du ligand.

- Externaliser les commandes shell longues dans des scripts .sh

- Ajouter une documentation et des commentaires détaillés dans l'ensemble du workflow.

- Compléter la création d'environnements conda dédiés à chaque étape du pipeline afin d'isoler les dépendances, éviter les conflits de versions et améliorer la reproductibilité du workflow.


Commandes utiles :
----------------
snakemake --use-conda --cores all -p







