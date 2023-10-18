#!/usr/bin/env python
# coding: utf-8

import argparse
import scanpy as sc
import anndata as ad

import scvi
import torch
torch.set_float32_matmul_precision("medium")
import warnings
warnings.filterwarnings('ignore')


def main(args):
    adata = sc.read(args.adata_path)
    adata.X = adata.layers["counts_lengthnorm"].A.copy()

    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")

    scvi_model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    scvi_model.train()
    scvi_model.save(f"{args.output_dir}/scvi_model.pt")

    scvi_latent = scvi_model.get_latent_representation()

    out_ad = ad.AnnData(
        X=scvi_latent
    )

    out_ad.obs.index = adata.obs.index

    out_ad.write_h5ad(f"{args.output_dir}/{args.output_filename}.h5ad")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run scVI script with command-line arguments.')
    parser.add_argument('--adata_path', type=str, default="/storage/data/data/organoid_atlas/230510_09_organoids_cleaned_hvg.h5ad",
                        help='Path to the input AnnData file (.h5ad).')
    parser.add_argument('--output_dir', type=str, default="/storage/data/data/organoid_atlas/output/integrations_21_05",
                        help='Path to the output directory.')
    parser.add_argument('--output_filename', type=str, default='scvi_results',
                        help='Name of the output file (default: "scvi_results").')

    args = parser.parse_args()

    main(args)
