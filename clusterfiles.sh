#!/bin/bash

PASSWORD="AnnaEvin1214!"

sshpass -p "$PASSWORD" scp -r rbanderson@jakar.utep.edu:/home/rbanderson/projectfiles/data1 .
sshpass -p "$PASSWORD" scp -r rbanderson@jakar.utep.edu:/home/rbanderson/projectfiles/data2 .
sshpass -p "$PASSWORD" scp -r rbanderson@jakar.utep.edu:/home/rbanderson/projectfiles/data3 .
