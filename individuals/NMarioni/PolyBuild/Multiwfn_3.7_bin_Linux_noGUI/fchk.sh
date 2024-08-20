#!/bin/bash

file_name=$1

formchk -3 ${file_name}_syn.chk ${file_name}_syn.fchk
formchk -3 ${file_name}_iso.chk ${file_name}_iso.fchk