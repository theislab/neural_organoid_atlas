import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import NMF
from sklearn.decomposition import non_negative_factorization
from scipy.stats import rankdata

import scanpy as sc
import igraph
import anndata as ad

# Run NMF on the input matrix X
def run_NMF(X, model = None, W = None, H = None, n_components = 10, max_iter = 200, return_model = False, **kwargs):
  if model is None or (not isinstance(model, NMF)):
    if H is None and W is None:
      model = NMF(n_components = int(n_components), max_iter = int(max_iter), **kwargs)
      W = model.fit_transform(X)
      H = model.components_
    else:
      W, H, n_iter = non_negative_factorization(X, W = W, H = H, n_components = H.shape[0], max_iter = int(max_iter), **kwargs)
  else:
    W = model.transform(X)
    H = model.components_
  
  if return_model and model is not None:
    return(model)
  else:
    return({"W" : W, "H" : H})

# Run PAGA given the path to the h5ad file
def run_PAGA_only(path, rep, k, grp, overwrite = True):
  h5ad = sc.read_h5ad(path)
  sc.pp.neighbors(h5ad, n_neighbors=k, use_rep=rep)
  sc.tl.paga(h5ad, groups=grp)
  if overwrite:
    h5ad.write_h5ad(path)
  return h5ad

# Calculate correlation between columns in matrix X and those in matrix Y, with X and Y being dense or sparse
def corSparse(X, Y = None):
  if Y is None: Y = X
  
  n = X.shape[0]
  muX = np.ravel(X.mean(0))
  muY = np.ravel(Y.mean(0))
  covmat = ( X.T.dot(Y) - (n * muX[:,np.newaxis].dot(muY[:,np.newaxis].T)) ) / (n-1)
  sdvecX = np.ravel(np.sqrt(((X.power(2)).sum(0) - n*(muX**2)) / (n-1)) if sparse.issparse(X) else np.sqrt(((X**2).sum(0) - n*(muX**2)) / (n-1)))
  sdvecY = np.ravel(np.sqrt(((Y.power(2)).sum(0) - n*(muY**2)) / (n-1)) if sparse.issparse(Y) else np.sqrt(((Y**2).sum(0) - n*(muY**2)) / (n-1)))
  cormat = covmat / sdvecX[:,np.newaxis].dot(sdvecY[:,np.newaxis].T)
  
  return cormat

# For each column of the input matrix X, convert the values into the ranks (for calculation of Spearman correlation, for example)
def rankMatrix(X):
  if sparse.issparse(X):
    idx_row, idx_col, dat = sparse.find(X)
    df = pd.DataFrame({'i' : idx_row, 'j' : idx_col, 'x' : dat})
    df['r'] = np.concatenate([ rankdata(x) + (df.shape[0]-len(x)) - (1+(df.shape[0]-len(x)))/2 for x in np.split(df.x, np.unique(df.j, return_index=True)[1][1:]) ])
    ranked = sparse.csr_matrix((df['r'], (df['i'], df['j'])), shape = X.shape)
  else:
    ranked = rankdata(X, method = "average", axis=0)
  
  return ranked

# Similar to rankMatrix but to use "dense" method for ranking non-zero values
def rankMatrix_dense(X):
  if sparse.issparse(X):
    idx_row, idx_col, dat = sparse.find(X)
    df = pd.DataFrame({'i' : idx_row, 'j' : idx_col, 'x' : dat})
    df['r'] = np.concatenate([ rankdata(x, method = "dense") for x in np.split(df.x, np.unique(df.j, return_index=True)[1][1:]) ])
    ranked = sparse.csr_matrix((df['r'], (df['i'], df['j'])), shape = X.shape)
  else:
    ranked = rankdata(X, method = "average", axis=0)
  
  return ranked

# Similar to rankMatrix, but only ranking non-zero values with zero values remaining zero
def rankMatrix_nonzero(X):
  if sparse.issparse(X):
    idx_row, idx_col, dat = sparse.find(X)
    df = pd.DataFrame({'i' : idx_row, 'j' : idx_col, 'x' : dat})
    df['r'] = np.concatenate([ rankdata(x) for x in np.split(df.x, np.unique(df.j, return_index=True)[1][1:]) ])
    ranked = sparse.csr_matrix((df['r'], (df['i'], df['j'])), shape = X.shape)
  else:
    ranked = rankdata(X, method = "average", axis=0)
  
  return ranked


# AggreCell algorithm to obtain the pseudocell index, given the anndata as the input
def generate_aggrecell_idx(adata, n_neighbors = 20, use_rep = None, n_pcs = 20, ratio = 0.1, initdist = 1, maxdist = 2, mergedist = 2, minsize = 3, weighted = False, seed = 123, constraint_on_grp = None, thres_grp_connectivities = 0.5):
  ## re-calculate the knn network
  sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep = use_rep)
  
  ## extract the knn network
  adj_knn = adata.obsp['connectivities']
  adj_knn.setdiag(adj_knn.max())
  dist_knn = adata.obsp['distances']
  
  ## if required, constrain the knn network on cells in the same cluster or connected cluster
  if constraint_on_grp is not None and constraint_on_grp in adata.obs.keys():
    mat_grp_ident = sparse.csr_matrix((np.repeat(1, adata.shape[0]), (np.arange(adata.shape[0]), adata.obs[constraint_on_grp].cat.codes)), shape = [adata.shape[0], len(adata.obs[constraint_on_grp].cat.categories)])
    sc.tl.paga(adata, groups=constraint_on_grp)
    adj_knn_grp = ((adata.uns['paga']['connectivities'] >= thres_grp_connectivities) + sparse.diags(np.repeat(1, len(adata.obs[constraint_on_grp].cat.categories)))).astype('bool')
    constr = mat_grp_ident.dot(adj_knn_grp).dot(mat_grp_ident.T)
    adj_knn = adj_knn.multiply(constr)
  
  ## build knn graph
  tup_knn = adj_knn.nonzero()
  tup_knn = [(adata.obs_names[tup_knn[0][i]], adata.obs_names[tup_knn[1][i]], dist_knn[tup_knn[0][i], tup_knn[1][i]]+1e-5) for i in range(len(tup_knn[0]))]
  graph = igraph.Graph.TupleList(tup_knn, weights=weighted)
  avg_weight = dist_knn[dist_knn.nonzero()].mean() if graph.is_weighted() else 1
  
  ## initial pseudocell territories and randomly assign cells being close to multiple pseudocell capitals
  N = adata.shape[0]
  node_index = np.ones(N) * -1
  if seed is not None: np.random.seed(seed)
  base_filter = np.random.uniform(0, 1, len(node_index)) < ratio
  shortestDistance = pd.DataFrame(graph.distances(source=np.where(base_filter)[0]))
  initdist = initdist if graph.is_weighted() or initdist >= 1 else 1
  selected = shortestDistance <= initdist * avg_weight
  for i in range(sum(base_filter)):
    selected.iloc[:, np.where(base_filter)[0][i]] = False
    selected.iloc[i, np.where(base_filter)[0][i]] = True
  belongToMulti = selected.apply(lambda x: np.sum(x) > 1, axis=0)
  for i in np.where(belongToMulti)[0]:
    x = selected.iloc[:,i]
    if seed is not None: np.random.seed(seed + i)
    selected.iloc[:,i] = np.array(range(len(x))) == np.where(x)[0][np.random.rand(np.sum(x)).argmax()]
  
  ## eliminate too small pseudocells and assign cells without pseudocell assignment but close to at least one
  selected = selected.loc[selected.sum(axis=1) >= minsize, :]
  if maxdist > initdist:
    belongToNowhereButNotFar = np.where(selected.mean(0) == 0 & shortestDistance.apply(lambda x: np.min(x[x != 0]) <= maxdist * avg_weight, axis=0))[0]
    candidates = [np.where(shortestDistance.iloc[:, i] == np.min(shortestDistance.iloc[:, i]))[0] for i in belongToNowhereButNotFar]
    for i in range(len(candidates)):
      x = candidates[i]
      if seed is not None: np.random.seed(seed * i + 11)
      candidates[i] = x[np.random.rand(len(x)).argmax()]
    newSel = pd.DataFrame([np.repeat([False, True, False], [i, 1, np.sum(base_filter) - i - 1]) for i in candidates]).transpose()
    newSel.columns = belongToNowhereButNotFar
    selected.update(newSel)
  
  ## generate output data frame
  node_index = list(selected.apply(lambda x: np.where(x)[0][0] if sum(x) > 0 else -1, axis=0))
  node_id = list(graph.get_vertex_dataframe().iloc[:, 0])
  aggrecell_idx = dict(zip(node_id, node_index))
  aggrecell_idx = pd.DataFrame({'cell': adata.obs_names, 'idx_aggrecell': [aggrecell_idx[x] for x in adata.obs_names]})
  aggrecell_idx.iloc[aggrecell_idx.iloc[:,1] == -1,1] = None
  
  return aggrecell_idx

## summarize single-cell data to pseudocell-level
def group_vec_to_ident_mat(group, norm = True):
  i = np.where(group.notnull())[0]
  j = group[i].astype(int)
  mat_ident = sparse.csr_matrix((np.repeat(1, len(i)), (i, j)), shape=(len(group), len(np.unique(j))))
  if norm:
      num_cells_per_group = np.ravel(mat_ident.sum(axis=0))
      mat_ident = mat_ident @ sparse.diags(1 / num_cells_per_group)
  
  return mat_ident

def summarize_numeric_matrix(mat, group, use_mean = True): # mat.shape = cell*feature; len(group) = cell
  mat_ident = group_vec_to_ident_mat(group, norm = use_mean)
  mat_summ = mat_ident.transpose() @ mat
  return mat_summ

def summarize_categorical_vector(fac, group): # len(fac) = cell; len(group) = cell
  mat_ident = group_vec_to_ident_mat(group, norm=False)
  mat_fac = pd.get_dummies(fac)
  mat_fac_summ = mat_ident.transpose() @ mat_fac
  if fac.isna().sum() == len(fac):
    return np.repeat(np.nan, mat_fac_summ.shape[0])
  fac_summ = mat_fac.columns[mat_fac_summ.argmax(axis = 1)]
  return fac_summ

def summarize_dataframe(df, group): # df.shape = cell*meta; len(group) = cell
  mats = list()
  # numeric columns
  if (df.select_dtypes("number").shape[1] > 0):
    mat_num = pd.DataFrame(summarize_numeric_matrix(df.select_dtypes("number"), group))
    mat_num.columns = df.select_dtypes("number").columns
    mats.append(mat_num)
  # boolean columns
  if (df.select_dtypes("bool").shape[1] > 0):
    mat_bool_raw = df.select_dtypes("bool")
    mat_bool = pd.concat([pd.DataFrame(summarize_categorical_vector(mat_bool_raw[x], group) == True) for x in mat_bool_raw.columns], axis=1)
    mat_bool.columns = mat_bool_raw.columns
    mats.append(mat_bool)
  # categorical columns
  if (df.select_dtypes(exclude=[ "number", "bool"]).shape[1] > 0):
    mat_categorical_raw = df.select_dtypes(exclude=[ "number", "bool"])
    mat_categorical = pd.concat([pd.DataFrame(summarize_categorical_vector(mat_categorical_raw[x], group)) for x in mat_categorical_raw.columns], axis=1)
    mat_categorical.columns = mat_categorical_raw.columns
    mats.append(mat_categorical)
  
  df_summ = pd.concat(mats, axis = 1)
  df_summ = df_summ[df.columns]
  return df_summ

def generate_aggrecells(adata, aggrecell_idx, prefix = "aggrcell_", X_is_count = False, counts_layer = ['counts'], **kwargs):
  expr_aggrecell = summarize_numeric_matrix(adata.X, group=aggrecell_idx.iloc[:, 1], use_mean = not X_is_count)
  
  meta_aggrecell = summarize_dataframe(adata.obs, group=aggrecell_idx.iloc[:, 1])
  meta_aggrecell.index = [prefix + str(x) for x in meta_aggrecell.index]
  
  ad_aggrecells = ad.AnnData(expr_aggrecell)
  ad_aggrecells.obs_names = meta_aggrecell.index
  ad_aggrecells.var_names = adata.var_names.copy()
  ad_aggrecells.obs = meta_aggrecell
  ad_aggrecells.var = adata.var.copy()
  
  for layer in adata.layers.keys():
    ad_aggrecells.layers[layer] = summarize_numeric_matrix(adata.X, group=aggrecell_idx.iloc[:, 1], use_mean = layer in counts_layer)
  
  return ad_aggrecells
