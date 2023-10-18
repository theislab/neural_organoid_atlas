#!/usr/bin/env python
# coding: utf-8

import argparse
import numpy as np
import scanpy as sc
import anndata as ad

from scarches.models.scpoli import scPoli
import warnings
warnings.filterwarnings('ignore')


def main(args):
    adata = sc.read(args.adata_path)
    adata.X = adata.layers["counts_lengthnorm"].A.copy()

    scpoli_model = scPoli(
        adata=adata,
        unknown_ct_names=["unknown"],
        condition_key="batch",
        cell_type_keys=args.cell_type_keys,
        embedding_dim=5,
        hidden_layer_sizes=[1024]
    )

    early_stopping_kwargs = {
        "early_stopping_metric": "val_prototype_loss",
        "mode": "min",
        "threshold": 0,
        "patience": 20,
        "reduce_lr": True,
        "lr_patience": 13,
        "lr_factor": 0.1,
    }

    scpoli_model.train(
        unlabeled_prototype_training=False,
        n_epochs=7,
        pretraining_epochs=5,
        early_stopping_kwargs=early_stopping_kwargs,
        eta=10,
        alpha_epoch_anneal=100
    )

    if args.store_model:
        scpoli_model.save(f"{args.output_dir}/{args.output_filename}.pkl")
    if args.store_batch_embeddings:
        scpoli_model.get_conditional_embeddings().write_h5ad(
            f"{args.output_dir}/{args.output_filename}_conditional_embeddings.h5ad")

    scpoli_latent = scpoli_model.get_latent(
        adata.X,
        adata.obs["batch"].values,
        mean=True
    )

    out_ad = ad.AnnData(
        X=scpoli_latent
    )

    out_ad.obs.index = adata.obs.index

    out_ad.write_h5ad(f"{args.output_dir}/{args.output_filename}.h5ad")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run scPoli script with command-line arguments.')
    parser.add_argument('--adata_path', type=str,
                        default="/storage/data/data/organoid_atlas/230510_09_organoids_cleaned_hvg.h5ad",
                        help='Path to the input AnnData file (.h5ad).')
    parser.add_argument('--output_dir', type=str, default="/storage/data/data/organoid_atlas/output/integrations_21_05",
                        help='Path to the output directory.')
    parser.add_argument('--cell_type_keys', type=str, nargs='+', default=['snapseed_pca_rss_level_1'],
                        help='List of cell type keys')
    parser.add_argument('--output_filename', type=str, default='scpoli_123_hierarchical',
                        help='Name of the output file')
    parser.add_argument('--store_model', action='store_true')
    parser.add_argument('--store_batch_embeddings', action='store_true')

    args = parser.parse_args()

    main(args)
