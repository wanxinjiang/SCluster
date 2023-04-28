import numpy as np
from scipy.sparse import coo_matrix

__all__ = ['bcubed', 'conditional_entropy', 'contingency_matrix', 'der',
           'goodman_kruskal_tau', 'mutual_information']


EPS = np.finfo(float).eps


def contingency_matrix(ref_labels, sys_labels):
    """Return contingency matrix between ``ref_labels`` and ``sys_labels``."""
    ref_classes, ref_class_inds = np.unique(ref_labels, return_inverse=True)
    sys_classes, sys_class_inds = np.unique(sys_labels, return_inverse=True)
    n_frames = ref_labels.size
    # Following works because coo_matrix sums duplicate entries. Is roughly
    # twice as fast as np.histogram2d.
    cmatrix = coo_matrix(
        (np.ones(n_frames), (ref_class_inds, sys_class_inds)),
        shape=(ref_classes.size, sys_classes.size),
        dtype=np.int)
    cmatrix = cmatrix.toarray()
    return cmatrix, ref_classes, sys_classes


def bcubed(ref_labels, sys_labels, cm=None):
    """Return B-cubed precision, recall, and F1.

    The B-cubed precision of an item is the proportion of items with its
    system label that share its reference label (Bagga and Baldwin, 1998).
    Similarly, the B-cubed recall of an item is the proportion pf items
    with its reference label that share its system label. The overall B-cubed
    precision and recall, then, are the means of the precision and recall for
    each item.

    Parameters
    ----------
    ref_labels : ndarray, (n_frames,)
        Reference labels.

    sys_labels : ndarray, (n_frames,)
        System labels.

    cm : ndarray, (n_ref_classes, n_sys_classes)
        Contingency matrix between reference and system labelings. If None,
        will be computed automatically from ``ref_labels`` and ``sys_labels``.
        Otherwise, the given value will be used and ``ref_labels`` and
        ``sys_labels`` ignored.
        (Default: None)

    Returns
    -------
    precision : float
        B-cubed precision.

    recall : float
        B-cubed recall.

    f1 : float
        B-cubed F1.

    References
    ----------
    Bagga, A. andBaldwin, B. (1998). "Algorithmsfor scoring coreference
    chains." Proceedings of LREC 1998.
    """
    if cm is None:
        cm, _, _ = contingency_matrix(ref_labels, sys_labels)
    cm = cm.astype('float64')
    cm_norm = cm / cm.sum()
    precision = np.sum(cm_norm * (cm / cm.sum(axis=0)))
    recall = np.sum(cm_norm * (cm / np.expand_dims(cm.sum(axis=1), 1)))
    f1 = 2*(precision*recall)/(precision + recall)
    return  f1