# sync github to nandi
#!/bin/bash

# What to backup. 
backup_files="/home/ben/projects/dehancer/*"

# Where to backup to.
dest="bsiranos@nandi:/users/bsiranos/projects/dehancer"

scp -r $backup_files $dest

