#!/bin/bash

PASSWORD="AnnaEvin1214!"

sshpass -p "$PASSWORD" scp -r rbanderson@jakar.utep.edu:/home/rbanderson/projectfiles/data1/simulation_data.h5 .
sshpass -p "$PASSWORD" scp -r rbanderson@jakar.utep.edu:/home/rbanderson/projectfiles/data2/simulation_data.h5 .
sshpass -p "$PASSWORD" scp -r rbanderson@jakar.utep.edu:/home/rbanderson/projectfiles/data3/simulation_data.h5 .
