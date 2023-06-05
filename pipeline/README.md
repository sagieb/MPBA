# Snakemake based pipline:<br>
## 1. put all the raw fastqfiles from ####geo_accession in the raw_data folder:<br> 
## 2. create environment:<br> 
`conda create -n mp_py python=3.8.3`<br> 
`conda activate mp_py`<br>
`pip install -r requirements.txt`<br>
`python -m ipykernel install --user --name mp_py --display-name "Python (mp_py)"`<br>
`conda install -c bioconda adapterremoval=2.3.1`<br>
## 3. run the pipeline
`chmod +x run.sh`<br>
`./run.sh`<br>
## 4. run the script build_csv_files to generate the final csv files
`python build_csv_files `<br>

