# ClusterEnsembles.py
#   Author: Takehiro Sano
#   Contact: tsano430@gmail.com
#   License: MIT License
import os
import warnings
from typing import Optional
import numpy as np
import pymetis
import kahypar
from scipy import sparse
from sklearn.metrics import pairwise_distances, normalized_mutual_info_score
from sklearn.utils.extmath import safe_sparse_dot


def create_hypergraph(base_clusters):
    """Create the incidence matrix of base clusters' hypergraph

    Parameter
    ----------
    base_clusters: labels generated by base clustering algorithms

    Return
    -------
    H: incidence matrix of base clusters' hypergraph
    """
    H = []
    len_bcs = base_clusters.shape[1]

    for bc in base_clusters:
        bc = np.nan_to_num(bc, nan=float('inf'))
        unique_bc = np.unique(bc)
        len_unique_bc = len(unique_bc)
        bc2id = dict(zip(unique_bc, np.arange(len_unique_bc)))
        tmp = [bc2id[bc_elem] for bc_elem in bc]
        h = np.identity(len_unique_bc, dtype=int)[tmp]
        if float('inf') in bc2id.keys():
            h = np.delete(h, obj=bc2id[float('inf')], axis=1)
        H.append(sparse.csc_matrix(h))

    return sparse.hstack(H)


def to_pymetis_format(adj_mat):
    """Transform an adjacency matrix into the pymetis format

    Parameter
    ---------
    adj_mat: adjacency matrix 

    Returns
    -------
    xadj, adjncy, eweights: parameters for pymetis
    """
    xadj = [0]
    adjncy = []
    eweights = []
    n_rows = adj_mat.shape[0]
    adj_mat = adj_mat.tolil()

    for i in range(n_rows):
        row = adj_mat.getrow(i)
        idx_row, idx_col = row.nonzero()
        val = row[idx_row, idx_col]
        adjncy += list(idx_col)
        eweights += list(val.toarray()[0])
        xadj.append(len(adjncy))

    return xadj, adjncy, eweights


def cspa(base_clusters, nclass):
    """Cluster-based Similarity Partitioning Algorithm (CSPA)

    Parameters
    ----------
    base_clusters: labels generated by base clustering algorithms
    nclass: number of classes 

    Return
    -------
    celabel: consensus clustering label obtained from CSPA
    """
    H = create_hypergraph(base_clusters)
    S = H * H.T

    xadj, adjncy, eweights = to_pymetis_format(S)
    membership = pymetis.part_graph(
        nparts=nclass, xadj=xadj, adjncy=adjncy, eweights=eweights)[1]
    celabel = np.array(membership)

    return celabel


# def hgpa(base_clusters, nclass, random_state):
    """HyperGraph Partitioning Algorithm (HGPA)

    Parameters
    ----------
    base_clusters: labels generated by base clustering algorithms
    nclass: number of classes 
    random_state: used for reproducible results

    Return
    -------
    celabel: consensus clustering label obtained from HGPA
    """
    # Create hypergraph for kahypar
    H = create_hypergraph(base_clusters)
    n_nodes, n_nets = H.shape

    node_weights = [1] * n_nodes
    edge_weights = [1] * n_nets

    hyperedge_indices = [0]
    hyperedges = []
    HT = H.T
    for i in range(n_nets):
        h = HT.getrow(i)
        idx_row, idx_col = h.nonzero()
        hyperedges += list(idx_col)
        hyperedge_indices.append(len(hyperedges))

    hypergraph = kahypar.Hypergraph(
        n_nodes, n_nets, hyperedge_indices, hyperedges, nclass, edge_weights, node_weights)

    # Settings for kahypar
    context = kahypar.Context()
    config_path = './kahypar_config/km1_kKaHyPar_sea20.ini'
    context.loadINIconfiguration(config_path)
    if random_state is not None:
        context.setSeed(random_state)
    context.setK(nclass)
    context.setEpsilon(0.03)
    context.suppressOutput(True)

    # Hypergraph partitioning
    kahypar.partition(hypergraph, context)

    celabel = np.empty(n_nodes, dtype=int)
    for i in range(n_nodes):
        celabel[i] = hypergraph.blockID(i)

    return celabel


def mcla(base_clusters, nclass, random_state):
    """Meta-CLustering Algorithm (MCLA)

    Parameters
    ----------
    base_clusters: labels generated by base clustering algorithms
    nclass: number of classes 
    random_state: used for reproducible results

    Return
    -------
    celabel: consensus clustering label obtained from MCLA
    """
    np.random.seed(random_state)

    # Construct Meta-graph
    H = create_hypergraph(base_clusters)
    n_cols = H.shape[1]

    W = sparse.identity(n_cols, dtype=float, format='lil')
    for i in range(n_cols):
        hi = H.getcol(i)
        norm_hi = (hi.T * hi)[0, 0]
        for j in range(n_cols):
            if i < j:
                hj = H.getcol(j)
                norm_hj = (hj.T * hj)[0, 0]
                inner_prod = (hi.T * hj)[0, 0]
                W[i, j] = inner_prod / (norm_hi + norm_hj - inner_prod)
                W[j, i] = W[i, j]
    W *= 1e3
    W = W.astype(int)

    # Cluster Hyperedges
    xadj, adjncy, eweights = to_pymetis_format(W)
    membership = pymetis.part_graph(
        nparts=nclass, xadj=xadj, adjncy=adjncy, eweights=eweights)[1]

    # Collapse Meta-clusters
    meta_clusters = sparse.dok_matrix(
        (base_clusters.shape[1], nclass), dtype=float).tolil()
    for i, v in enumerate(membership):
        meta_clusters[:, v] += H.getcol(i)

    # Compete for Objects
    celabel = np.empty(base_clusters.shape[1], dtype=int)
    for i, v in enumerate(meta_clusters):
        v = v.toarray()[0]
        celabel[i] = np.random.choice(np.nonzero(v == np.max(v))[0])

    return celabel


def hbgf(base_clusters, nclass):
    """Hybrid Bipartite Graph Formulation (HBGF) 

    Parameters
    ----------
    base_clusters: labels generated by base clustering algorithms
    nclass: number of classes 

    Return
    -------
    celabel: consensus clustering label obtained from HBGF
    """
    A = create_hypergraph(base_clusters)
    rowA, colA = A.shape
    W = sparse.bmat([[sparse.dok_matrix((colA, colA)), A.T],
                    [A, sparse.dok_matrix((rowA, rowA))]])
    xadj, adjncy, _ = to_pymetis_format(W)
    membership = pymetis.part_graph(
        nparts=nclass, xadj=xadj, adjncy=adjncy, eweights=None)[1]
    celabel = np.array(membership[colA:])
    return celabel


def create_connectivity_matrix(base_clusters):
    """Create the connectivity matrix

    Parameter
    ---------
    base_clusters: labels generated by base clustering algorithms

    Return
    ------
    M: connectivity matrix
    """
    n_bcs, len_bcs = base_clusters.shape
    M = np.zeros((len_bcs, len_bcs))
    m = np.zeros_like(M)

    for bc in base_clusters:
        for i, elem_bc in enumerate(bc):
            m[i] = np.where(elem_bc == bc, 1, 0)
        M += m

    M /= n_bcs
    return sparse.csr_matrix(M)


def orthogonal_nmf_algorithm(W, nclass, random_state, maxiter):
    """Algorithm for bi-orthogonal three-factor NMF problem

    Parameters
    ----------
    W: given matrix 
    random_state: used for reproducible results
    maxiter: maximum number of iterations

    Return
    -------
    Q, S: factor matrices
    """
    np.random.seed(random_state)

    n = W.shape[0]
    Q = np.random.rand(n, nclass).reshape(n, nclass)
    S = np.diag(np.random.rand(nclass))

    for _ in range(maxiter):
        # Update Q
        WQS = safe_sparse_dot(W, np.dot(Q, S), dense_output=True)
        Q = Q * np.sqrt(WQS / (np.dot(Q, np.dot(Q.T, WQS)) + 1e-8))
        # Update S
        QTQ = np.dot(Q.T, Q)
        WQ = safe_sparse_dot(W, Q, dense_output=False)
        QTWQ = safe_sparse_dot(Q.T, WQ, dense_output=True)
        S = S * np.sqrt(QTWQ / (np.dot(QTQ, np.dot(S, QTQ)) + 1e-8))

    return Q, S


def nmf(base_clusters, nclass, random_state, maxiter=200):
    """NMF-based consensus clustering

    Parameters
    ----------
    base_clusters: labels generated by base clustering algorithms
    nclass: number of classes 
    random_state: used for reproducible results
    maxiter: maximum number of iterations 

    Return
    -------
    celabel: consensus clustering label obtained from NMF
    """
    M = create_connectivity_matrix(base_clusters)
    Q, S = orthogonal_nmf_algorithm(M, nclass, random_state, maxiter)
    celabel = np.argmax(np.dot(Q, np.sqrt(S)), axis=1)
    return celabel


def calc_objective(base_clusters, consensus_cluster):
    """Calculate the objective function value for cluster ensembles

    Parameters
    ----------
    base_clusters: labels generated by base clustering algorithms
    consensus_cluster: consensus clustering label

    Return
    -------
    objv: objective function value
    """
    objv = 0.0
    for bc in base_clusters:
        idx = np.isfinite(bc)
        objv += normalized_mutual_info_score(
            consensus_cluster[idx], bc[idx], average_method='geometric')
    objv /= base_clusters.shape[0]
    return objv


def cluster_ensembles(
        base_clusters: np.ndarray,
        nclass: Optional[int] = None,
        solver: str = 'hbgf',
        random_state: Optional[int] = None,
        verbose: bool = False) -> np.ndarray:
    """Generate a single consensus cluster using base clusters obtained from multiple clustering algorithms

    Parameters
    ----------
    base_clusters: labels generated by base clustering algorithms
    nclass: number of classes
    solver: cluster ensembles solver to use
    random_state: used for 'hgpa', 'mcla', and 'nmf'. Pass a nonnegative integer for reproducible results.
    verbose: whether to be verbose

    Return
    -------
    celabel: consensus clustering label 
    """
    if nclass is None:
        nclass = -1
        for bc in base_clusters:
            len_unique_bc = len(np.unique(bc[~np.isnan(bc)]))
            nclass = max(nclass, len_unique_bc)

    if verbose:
        print('Cluster Ensembles')
        print('    - number of classes:', nclass)
        print('    - solver:', solver)
        print('    - length of base clustering labels:',
              base_clusters.shape[1])
        print('    - number of base clusters:', base_clusters.shape[0])

    if not (isinstance(nclass, int) and nclass > 0):
        raise ValueError(
            'Number of class must be a positive integer; got (nclass={})'.format(nclass))

    if not ((random_state is None) or isinstance(random_state, int)):
        raise ValueError(
            'Number of random_state must be a nonnegative integer; got (random_state={})'.format(random_state))

    if isinstance(random_state, int):
        random_state = abs(random_state)

    if solver == 'cspa':
        if base_clusters.shape[1] > 5000:
            warnings.warn(
                '`base_clusters.shape[1]` is too large, so the use of another solvers is recommended.')
        celabel = cspa(base_clusters, nclass)
    # elif solver == 'hgpa':
    #     celabel = hgpa(base_clusters, nclass, random_state)
    elif solver == 'mcla':
        celabel = mcla(base_clusters, nclass, random_state)
    elif solver == 'hbgf':
        celabel = hbgf(base_clusters, nclass)
    elif solver == 'nmf':
        celabel = nmf(base_clusters, nclass, random_state)
    elif solver == 'all':
        if verbose:
            print('    - ANMI:')
        # ce_solvers = {'hgpa': hgpa, 'mcla': mcla, 'hbgf': hbgf}
        ce_solvers = { 'mcla': mcla, 'hbgf': hbgf}
        if base_clusters.shape[1] <= 5000:
            ce_solvers['cspa'] = cspa
            ce_solvers['nmf'] = nmf
        best_objv = None
        for name, ce_solver in ce_solvers.items():
            if ce_solver == cspa or ce_solver == hbgf:
                label = ce_solver(base_clusters, nclass)
            else:
                label = ce_solver(base_clusters, nclass, random_state)
            objv = calc_objective(base_clusters, label)
            if verbose:
                print('        -', name, ':', objv)
            if best_objv is None:
                best_objv = objv
                best_solver = name
                celabel = label
            if best_objv < objv:
                best_objv = objv
                best_solver = name
                celabel = label
        if verbose:
            print('    - Best solver:', best_solver)
    else:
        raise ValueError(
            "Invalid solver parameter: got '{}' instead of one of ('cspa', 'mcla', 'hbgf', 'nmf', 'all')".format(solver))

    return celabel