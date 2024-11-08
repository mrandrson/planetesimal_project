#!/bin/bash

PASSWORD="thisismynewpassword"

sshpass -p "$PASSWORD" scp -r rbanderson@jakar.utep.edu:/home/rbanderson/projectfiles/data3/simulation_data.h5 .
