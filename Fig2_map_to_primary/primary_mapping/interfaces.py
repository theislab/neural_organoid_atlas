import argparse


def vae_interface():
    parser = argparse.ArgumentParser(description="Fits autoencoder")

    parser.add_argument(
        "--h5ad",
        type=str,
        dest="H5AD",
    )

    parser.add_argument(
        "--ref_h5ad",
        type=str,
        dest="REF_H5AD",
    )

    parser.add_argument("-o", "--output", type=str, dest="output", default="output")

    parser.add_argument("-m", "--model", type=str, dest="model")

    parser.add_argument(
        "-b",
        "--batch",
        type=str,
        dest="batch",
    )

    parser.add_argument(
        "--annotation",
        type=str,
        dest="annot",
    )

    parser.add_argument(
        "-e",
        "--epochs",
        type=int,
        dest="epochs",
        default=500,
    )

    parser.add_argument(
        "--n_latent",
        type=int,
        dest="n_latent",
        default=20,
    )

    parser.add_argument(
        "--batch_size",
        type=int,
        dest="batch_size",
        default=1024,
    )

    parser.add_argument(
        "--cat_covars",
        nargs="+",
        type=str,
        dest="cat_covars",
    )

    parser.add_argument(
        "--cont_covars",
        nargs="+",
        type=str,
        dest="cont_covars",
    )

    parser.add_argument(
        "--disentangling_method", type=str, dest="dis_method", default="tc"
    )

    parser.add_argument(
        "--disentangling_weight", type=float, dest="dis_weight", default="1e-5"
    )

    parser.add_argument("--kld_weight", type=float, dest="kld_weight", default="1e-5")

    parser.add_argument(
        "--subset_hvg",
        action="store_true",
        dest="subset_hvg",
    )

    parser.add_argument(
        "--scarches",
        action="store_true",
        dest="scarches",
    )

    parser.add_argument(
        "--scanvi",
        action="store_true",
        dest="scanvi",
    )

    parser.add_argument(
        "--new_batch_key",
        type=str,
        dest="new_batch_key",
    )

    parser.add_argument("--color", nargs="+", type=str, dest="color")

    parser.add_argument("--reduction", type=str, dest="reduction", default="umap")

    parser.add_argument("--save_anndata", action="store_true", dest="save_anndata")

    parser.add_argument("--pretrain", action="store_true", dest="pretrain")

    parser.add_argument(
        "--compute_full_umap", action="store_true", dest="compute_full_umap"
    )

    parser.add_argument(
        "--compute_query_umap", action="store_true", dest="compute_query_umap"
    )

    parser.add_argument(
        "--compute_ref_umap", action="store_true", dest="compute_ref_umap"
    )

    parser.add_argument(
        "--remove_query_labels", action="store_true", dest="remove_query_labels"
    )

    args = parser.parse_args()

    args.color = args.color if args.color else []
    args.color = args.color + [args.batch] if args.batch else args.color
    args.color = args.color + args.cat_covars if args.cat_covars else args.color
    args.color = args.color + args.cont_covars if args.cont_covars else args.color
    args.color = args.color + [args.annot] if args.annot else args.color

    return args


import argparse
