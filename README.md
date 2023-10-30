# MPBA Pipeline (Snakemake based):<br>
## 1. put all the raw fastqfiles from GSE234455 in the pipeline/raw_data folder:<br> 
## 2. create environment:<br> 
`conda create -n MPBA_pl python=3.8.3`<br> 
`conda activate MPBA_pl`<br>
`pip install -r ./pipeline/requirements.txt`<br>
`conda install -c bioconda adapterremoval=2.3.1`<br>
## 3. run the pipeline
`chmod +x ./pipeline/run.sh`<br>
`./pipeline/run.sh`<br>
## 4. run the script build_csv_files to generate the final csv files
`python ./pipeline/build_csv_files.py `<br>

# Data analysis
1. To run the code first install the python requirements:<br>
`conda create -n mpba_ana python=3.8.3`<br>
`conda activate mpba_ana`<br>
`pip install -r ./requirements.txt`<br>
`python -m ipykernel install --user --name mpba_ana --display-name "mpba_ana"`<br>
2. When you run figure scripts use jupyter notebook and make sure you're within the 'figure_code' directory 
