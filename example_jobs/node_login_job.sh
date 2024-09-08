# Interactive Job Script 0.2
#
# This gives you a command prompt running in the kerenh queue.
# Resources:
# - Bash shell
# - 8 GB RAM
# - 1 CPU core
#####################################
#!/bin/bash
#PBS -S /bin/bash
#PBS -q itaym
#PBS -N clustInteractive
#PBS -l select=ncpus=16:mem=16gb
#PBS -I
