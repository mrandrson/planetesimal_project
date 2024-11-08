#!/bin/bash

PASSWORD="thisismynewpassword"
rm -rf particle_energies
mkdir particle_energies

for i in {1..100}
do
  sshpass -p "$PASSWORD" scp -r rbanderson@jakar.utep.edu:/home/rbanderson/projectfiles/shuscripts/outputdata/particle_dat$i.csv ./particle_energies/
done

