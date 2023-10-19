#!/usr/bin/env python
# coding: utf-8

import argparse
import os

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
    vae = scvi.model.SCVI.load(os.path.join(args.output_dir, "scvi_model.pt"), adata=adata)
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        labels_key=args.cell_type_key,
        unlabeled_category="unknown",
    )
    lvae.train()

    scanvi_latent = lvae.get_latent_representation(adata)

    out_ad = ad.AnnData(
        X=scanvi_latent
    )
    out_ad.obs.index = adata.obs.index

    out_ad.write_h5ad(f"{args.output_dir}/{args.output_filename}.h5ad")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run scANVI script with command-line arguments.')
    parser.add_argument('--adata_path', type=str,
                        default="/storage/data/data/organoid_atlas/230510_09_organoids_cleaned_hvg.h5ad",
                        help='Path to the input AnnData file (.h5ad).')
    parser.add_argument('--output_dir', type=str, default="/storage/data/data/organoid_atlas/output/integrations_21_05",
                        help='Path to the output directory.')
    parser.add_argument('--cell_type_key', type=str, default='snapseed_pca_rss_level_1',
                        help='cell type key (default: "snapseed_pca_rss_level_1").')
    parser.add_argument('--output_filename', type=str, default='scanvi_results',
                        help='Name of the output file (default: "scanvi_results").')

    args = parser.parse_args()

    main(args)
