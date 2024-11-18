import warnings

warnings.filterwarnings("ignore")

import os
import sys

import numpy as np
import pandas as pd

import scvi
from scvi.model.utils import mde

import scanpy as sc

from interfaces import vae_interface


if __name__ == "__main__":
    args = vae_interface()

    # PATHS
    PROJECT_DIR = os.path.join(
        "/home/fleckj/projects/atlas/data/results/atlasv3_scanvi_scarches/v7"
    )
    OUT_DIR = os.path.join(PROJECT_DIR, args.annot)
    REF_PATH = "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_common_hv2k.h5ad"
    QUERY_PATH = "/home/fleckj/projects/atlas/scratch/v3/230510_phase3_hv2k_braun.h5ad"

    # PARAMS
    BATCH_SIZE = 1024
    EPOCHS_PRETRAIN = 500
    EPOCHS_FINETUNE = 100

    print("Loading reference & query data")
    adata_ref = sc.read(REF_PATH)
    print(adata_ref)

    adata_q = sc.read(QUERY_PATH)
    print(adata_q)

    # Subset to common features
    adata_q.varm = dict()
    adata_q.obs["batch"] = adata_q.obs["bio_sample"].astype(str).copy()

    adata_ref.varm = dict()
    adata_ref.obs["batch"] = adata_ref.obs["Donor"].astype(str).copy()

    ##########################
    #### SCVI PRETRAINING ####
    ##########################
    print("Setting up anndata")
    scvi.model.SCVI.setup_anndata(
        adata_ref,
        batch_key="batch",
        layer="counts",
    )

    try:
        print("Loading pretrained reference model")
        vae = scvi.model.SCVI.load(
            os.path.join(PROJECT_DIR, "ref_scvi/model.pt"),
            adata=adata_ref,
        )

    except ValueError as err:
        print(err)
        print("No pretrained reference model found, starting from scratch")

        # Set up with params recommended for scARCHES mapping
        arches_params = dict(
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
            dropout_rate=0.2,
        )

        vae = scvi.model.SCVI(
            adata_ref,
            n_latent=20,
            n_layers=2,
            n_hidden=256,
            **arches_params,
        )

        print("Fitting reference model")
        vae.train(
            batch_size=BATCH_SIZE,
            max_epochs=EPOCHS_PRETRAIN,
            early_stopping=True,
        )

        out_dir = os.path.join(PROJECT_DIR, "ref_scvi")
        os.makedirs(out_dir, exist_ok=True)
        vae.save(
            os.path.join(out_dir, "model.pt"),
            overwrite=True,
            save_anndata=False,
        )

    ###########################
    #### SCANVI FINETUNING ####
    ###########################
    try:
        print("Loading finetuned reference model")
        lvae = scvi.model.SCANVI.load(
            os.path.join(OUT_DIR, "ref_scanvi/model.pt"),
            adata=adata_ref,
        )

    except ValueError as err:
        print(err)
        print("No finetuned reference model found, starting from scratch")
        # Use snapseed labels as scANVI integration labels
        adata_ref.obs[args.annot] = adata_ref.obs[args.annot].astype(str)
        lvae = scvi.model.SCANVI.from_scvi_model(
            vae,
            adata=adata_ref,
            labels_key=args.annot,
            unlabeled_category="Unknown",
        )

        print("Finetuning model")
        lvae.train(
            batch_size=BATCH_SIZE,
            max_epochs=EPOCHS_FINETUNE,
            n_samples_per_label=100,
            early_stopping=True,
        )

        out_dir = os.path.join(OUT_DIR, "ref_scanvi")
        os.makedirs(out_dir, exist_ok=True)
        lvae.save(
            os.path.join(out_dir, "model.pt"),
            overwrite=True,
            save_anndata=False,
        )

    print("Computing and plotting reference MDE")
    adata_ref.obsm["X_scANVI"] = lvae.get_latent_representation(give_mean=True)
    adata_ref.obsm["X_umap"] = mde(adata_ref.obsm["X_scANVI"])

    ref_mde = pd.DataFrame(adata_ref.obsm["X_umap"])
    ref_mde.index = adata_ref.obs.index
    ref_mde.to_csv(os.path.join(OUT_DIR, "ref_mde.tsv"), sep="\t")

    PLOT_DIR = os.path.join(OUT_DIR, "plots")
    os.makedirs(PLOT_DIR, exist_ok=True)

    # Plot MDE with snapseed and batch labels
    p = sc.pl.umap(
        adata_ref,
        color=[args.annot, "batch"],
        show=False,
    )
    p[0].figure.savefig(os.path.join(PLOT_DIR, "scanvi_ref_mde.png"))

    ###########################
    #### REFERENCE MAPPING ####
    ###########################

    if args.remove_query_labels or args.annot not in adata_q.obs.columns:
        adata_q.obs[args.annot] = "Unknown"
        id_string = "_nolabs"
    else:
        adata_q.obs[args.annot] = adata_q.obs[args.annot].astype(str)
        # Subset query adata to common annotation labels
        common_labels = np.intersect1d(
            adata_ref.obs[args.annot].unique(),
            adata_q.obs[args.annot].unique(),
        )
        adata_q = adata_q[adata_q.obs[args.annot].isin(common_labels), :]
        id_string = ""

    scvi.model.SCANVI.prepare_query_anndata(adata_q, lvae)
    vae_q = scvi.model.SCANVI.load_query_data(adata_q, lvae)

    try:
        print("Loading query model")
        out_dir = os.path.join(OUT_DIR, "query_scarches" + id_string)
        vae_q = scvi.model.SCANVI.load(
            os.path.join(out_dir, "model.pt"),
            adata=adata_q,
        )

    except ValueError as err:
        print(err)
        print("No query model found, starting from scratch")
        print("Training query model")
        vae_q.train(
            batch_size=BATCH_SIZE,
            max_epochs=EPOCHS_FINETUNE,
            plan_kwargs=dict(weight_decay=0.0),
            check_val_every_n_epoch=10,
        )

        out_dir = os.path.join(OUT_DIR, "query_scarches" + id_string)
        os.makedirs(out_dir, exist_ok=True)
        vae_q.save(
            os.path.join(out_dir, "model.pt"),
            overwrite=True,
            save_anndata=False,
        )

    print("Computing and plotting query + reference projection MDE")
    adata_q.obs["predictions"] = vae_q.predict()

    adata_q.obs["ref"] = "query"
    adata_ref.obs["ref"] = "reference"

    adata_full = adata_q.concatenate(adata_ref)
    # For some reason batch col becomes messed up when concatenating,
    # so we need to fix it
    adata_full.obs["batch"] = (
        adata_q.obs["batch"].tolist() + adata_ref.obs["batch"].tolist()
    )
    print(adata_full)
    adata_full.obsm["X_scARCHES"] = vae_q.get_latent_representation(adata_full)
    adata_full.obsm["X_umap"] = mde(adata_full.obsm["X_scARCHES"])

    full_mde = pd.DataFrame(adata_full.obsm["X_umap"])
    full_mde.index = adata_full.obs.index
    full_mde.to_csv(os.path.join(OUT_DIR, f"full{id_string}_mde.tsv"), sep="\t")

    full_latent = pd.DataFrame(adata_full.obsm["X_scARCHES"])
    full_latent.index = adata_full.obs.index
    full_latent.to_csv(os.path.join(OUT_DIR, f"full{id_string}_latent.tsv"), sep="\t")

    adata_full.obs.to_csv(os.path.join(OUT_DIR, "full_meta.tsv"), sep="\t")
    adata_ref.obs.to_csv(os.path.join(OUT_DIR, "ref_meta.tsv"), sep="\t")
    adata_q.obs.to_csv(os.path.join(OUT_DIR, "query_meta.tsv"), sep="\t")

    PLOT_DIR = os.path.join(OUT_DIR, "plots")
    os.makedirs(PLOT_DIR, exist_ok=True)

    # Plot MDE with snapseed and batch labels
    p = sc.pl.umap(
        adata_full,
        color=[args.annot, "batch", "predictions", "ref"],
        show=False,
    )
    p[0].figure.savefig(os.path.join(PLOT_DIR, f"scanvi_full{id_string}_mde.png"))

    if args.compute_full_umap:
        print("Computing full UMAP")
        sc.pp.neighbors(adata_full, use_rep="X_scARCHES")
        sc.tl.umap(adata_full, min_dist=0.1)

        full_umap = pd.DataFrame(adata_full.obsm["X_umap"])
        full_umap.index = adata_full.obs.index
        full_umap.to_csv(os.path.join(OUT_DIR, f"full{id_string}_umap.tsv"), sep="\t")

        # Plot UMAP with snapseed and batch labels
        p = sc.pl.umap(
            adata_full,
            color=[args.annot, "batch", "predictions", "ref"],
            show=False,
        )
        p[0].figure.savefig(os.path.join(PLOT_DIR, f"scanvi_full{id_string}_umap.png"))

    if args.compute_query_umap:
        print("Computing query UMAP")
        adata_q.obsm["X_scARCHES"] = vae_q.get_latent_representation(adata_q)
        sc.pp.neighbors(adata_q, use_rep="X_scARCHES")
        sc.tl.umap(adata_q, min_dist=0.1)

        query_umap = pd.DataFrame(adata_q.obsm["X_umap"])
        query_umap.index = adata_q.obs.index
        query_umap.to_csv(os.path.join(OUT_DIR, f"query{id_string}_umap.tsv"), sep="\t")

        # Plot UMAP with snapseed and batch labels
        p = sc.pl.umap(
            adata_q,
            color=[args.annot, "predictions"],
            show=False,
        )
        p[0].figure.savefig(os.path.join(PLOT_DIR, f"scanvi_query{id_string}_umap.png"))

    if args.compute_ref_umap:
        print("Computing reference UMAP")
        adata_ref.obsm["X_scARCHES"] = vae_q.get_latent_representation(adata_ref)
        sc.pp.neighbors(adata_ref, use_rep="X_scARCHES")
        sc.tl.umap(adata_ref, min_dist=0.1)

        ref_umap = pd.DataFrame(adata_ref.obsm["X_umap"])
        ref_umap.index = adata_ref.obs.index
        ref_umap.to_csv(os.path.join(OUT_DIR, f"ref{id_string}_umap.tsv"), sep="\t")

        # Plot UMAP with snapseed and batch labels
        p = sc.pl.umap(
            adata_ref,
            color=[args.annot, "batch"],
            show=False,
        )
        p[0].figure.savefig(
            os.path.join(PLOT_DIR, f"scanvi_reference{id_string}_umap.png")
        )
