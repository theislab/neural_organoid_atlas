Running the shell script `./integration_script.sh` executes all integration approaches. 

To obtain the merged adata object which contains all integrated embeddings, you can either download the full object from the link in the README of this repo or alternatively use the code in `../notebooks/03_load_annotations.ipynb` to add it to the adata object.

Running the python script `./metrics_script.py` computes the scIB metrics for all integrations using the merged object.

Note that the filepaths need to be modified accordingly.
