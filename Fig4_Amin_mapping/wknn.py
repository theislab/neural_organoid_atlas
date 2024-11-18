from scipy import sparse
from typing import Optional, Union, Mapping, Literal
import warnings
import pandas as pd
import numpy as np

import cudf
from cuml.neighbors import NearestNeighbors


def nn2adj(nn, n1=None, n2=None):
    if n1 is None:
        n1 = nn[1].shape[0]
    if n2 is None:
        n2 = np.max(nn[1].flatten())

    df = pd.DataFrame(
        {
            "i": np.repeat(range(nn[1].shape[0]), nn[1].shape[1]),
            "j": nn[1].flatten(),
            "x": nn[0].flatten(),
        }
    )
    adj = sparse.csr_matrix(
        (np.repeat(1, df.shape[0]), (df["i"], df["j"])), shape=(n1, n2)
    )
    return adj


def build_nn(ref, query=None, k=100):
    if query is None:
        query = ref

    model = NearestNeighbors(n_neighbors=k)
    model.fit(ref)
    knn = model.kneighbors(query)

    adj = nn2adj(knn, n1=query.shape[0], n2=ref.shape[0])
    return adj


def build_mutual_nn(dat1, dat2=None, k1=100, k2=None):
    if dat2 is None:
        dat2 = dat1
    if k2 is None:
        k2 = k1

    adj_12 = build_nn(dat1, dat2, k=k2)
    adj_21 = build_nn(dat2, dat1, k=k1)

    adj_mnn = adj_12.multiply(adj_21.T)
    return adj_mnn


def get_transition_prob_mat(dat, k=50, symm=True):
    adj = build_nn(dat, k=k)
    if symm:
        adj = adj + adj.transpose()
    prob = sparse.diags(1 / np.array(adj.sum(1)).flatten()) @ adj.transpose()
    return prob


def random_walk_with_restart(init, transition_prob, alpha=0.5, num_rounds=100):
    init = np.array(init).flatten()
    heat = init[:, None]
    for i in range(num_rounds):
        heat = init[:, None] * alpha + (1 - alpha) * (
            transition_prob.transpose() @ heat
        )
    return heat


def get_wknn(
    ref,  # the ref representation to build ref-query neighbor graph
    query,  # the query representation to build ref-query neighbor graph
    ref2=None,  # the ref representation to build ref-ref neighbor graph
    k: int = 100,  # number of neighbors per cell
    query2ref: bool = True,  # consider query-to-ref neighbors
    ref2query: bool = True,  # consider ref-to-query neighbors
    weighting_scheme: Literal[
        "n", "top_n", "jaccard", "jaccard_square"
    ] = "jaccard_square",  # how to weight edges in the ref-query neighbor graph
    top_n: Optional[int] = None,
    return_adjs: bool = False,
):
    adj_q2r = build_nn(ref=ref, query=query, k=k)

    adj_r2q = None
    if ref2query:
        adj_r2q = build_nn(ref=query, query=ref, k=k)

    if query2ref and not ref2query:
        adj_knn = adj_q2r.T
    elif ref2query and not query2ref:
        adj_knn = adj_r2q
    elif ref2query and query2ref:
        adj_knn = ((adj_r2q + adj_q2r.T) > 0) + 0
    else:
        warnings.warn(
            "At least one of query2ref and ref2query should be True. Reset to default with both being True."
        )
        adj_knn = ((adj_r2q + adj_q2r.T) > 0) + 0

    if ref2 is None:
        ref2 = ref
    adj_ref = build_nn(ref=ref2, k=k)
    num_shared_neighbors = adj_q2r @ adj_ref.T
    num_shared_neighbors_nn = num_shared_neighbors.multiply(adj_knn.T)

    wknn = num_shared_neighbors_nn.copy()
    if weighting_scheme == "top_n":
        if top_n is None:
            top_n = k // 4 if k > 4 else 1
        wknn = (wknn > top_n) * 1
    elif weighting_scheme == "jaccard":
        wknn.data = wknn.data / (k + k - wknn.data)
    elif weighting_scheme == "jaccard_square":
        wknn.data = (wknn.data / (k + k - wknn.data)) ** 2

    if return_adjs:
        adjs = {"q2r": adj_q2r, "r2q": adj_r2q, "knn": adj_knn, "r2r": adj_ref}
        return (wknn, adjs)
    else:
        return wknn


def transfer_labels(ref_adata, query_adata, wknn, label_key="celltype"):
    scores = pd.DataFrame(
        wknn @ pd.get_dummies(ref_adata.obs[label_key]),
        columns=pd.get_dummies(ref_adata.obs[label_key]).columns,
        index=query_adata.obs_names,
    )
    scores["best_label"] = scores.idxmax(1)
    scores["best_score"] = scores.max(1)
    return scores
