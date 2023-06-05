#!/bin/bash
python gen_barcodes.py 2> "errorlog.log"

ls | grep "\.yaml$" | while read line; do
	if [[ -e $line.svg ]]; then
		echo "DAG exists already"
	else
		snakemake all --configfile $line --dag | dot -Tsvg > $line.svg
	fi
	snakemake all --configfile $line --rerun-incomplete --cores 24 
done