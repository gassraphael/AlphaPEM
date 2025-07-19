#!/bin/bash

# **Configuration:**
# Select correspond au nombre de serveurs, et ncpus au nombre de CPUS
#PBS -l select=1:ncpus=32     

# Durée maximale d'exécution du programme sur le cluster
#PBS -l walltime=90:00:00            
#PBS -o out
#PBS -e err

# Nom du job
#PBS -N AlphaPEM_calibration         
#PBS -V

# Adresse email pour être contacté à la fin du job.
#PBS -M raphael.gass@univ-reunion.fr 
# (b)eginning, (e)nd and (a)bortion. Un email est envoyé pour chacun de ces états.
#PBS -m bea                          


# **Infos utiles affichées à l'écran:**
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------


# Se placer dans le répertoire de soumission
cd "$PBS_O_WORKDIR"

# Identifier le dossier courant (où est lancé qsub) et son parent
FOLDER_ROOT=$(realpath .)                     # dossier courant (ex: .../project/subfolder)
FOLDER_NAME=$(basename "$FOLDER_ROOT")        # nom du dossier courant (ex: subfolder)
PROJECT_ROOT=$(realpath ..)      	      # dossier parent (ex: .../project)
PROJECT_NAME=$(basename "$PROJECT_ROOT")      # nom du dossier parent (ex: project)
SCRIPT_RELATIVE_PATH="${FOLDER_NAME}/${file}" # Chemin relatif du script à exécuter depuis la racine du projet

# Création du répertoire de travail temporaire (car le répertoire local est limité en taille)
export PBS_TMPDIR=/gpfs/scratch/$USER/$PBS_JOBID
mkdir -p "$PBS_TMPDIR"

# Copier tout le dossier parent (contenant l'ensemble du projet) dans le répertoire temporaire
cp -r "$PROJECT_ROOT" "$PBS_TMPDIR"

# Se déplacer dans la copie du projet
cd "$PBS_TMPDIR/$PROJECT_NAME"

# Chargement des modules
module purge
module load mpi/intel/2019_update2 # conserver ces versions, car ce sont celles installées sur le cluster.
module load tools/python/3.7.2     # conserver ces versions, car ce sont celles installées sur le cluster.

# Création de l’environnement virtuel
unset PYTHONPATH # pour s'assurer que les paquets installés dans l'environnement seront ceux utilisés.
python3 -m venv env
source env/bin/activate

# Installation des paquets permettant de faire fonctionner le programme
python3 -m pip install --upgrade --force-reinstall --no-cache-dir pip setuptools wheel
python3 -m pip install --upgrade --force-reinstall --no-cache-dir numpy scipy matplotlib colorama pygad

# Ajouter la racine du projet au PYTHONPATH pour permettre les imports relatifs
export PYTHONPATH=.

# Exécuter le script Python chargé avec qsub
mpiexec -np 1 python3 "$SCRIPT_RELATIVE_PATH"

# Copie des résultats de l'exécution vers le dossier d'origine et supprime le répertoire temporaire 
cd "$PBS_TMPDIR/$PROJECT_NAME"
cp -ru . "$PROJECT_ROOT"
cd "$PBS_O_WORKDIR"
rm -rf "$PBS_TMPDIR"
