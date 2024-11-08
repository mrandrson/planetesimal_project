#!/bin/bash

PASSWORD="thisismynewpassword"

sshpass -p "$PASSWORD" scp -r rbanderson@jakar.utep.edu:/home/rbanderson/projectfiles/simpledata/simpledata.h5 .
