from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from subscript import *

def SVM_Classifier(X, y, indep=None, fold=5, kernel='rbf', degree=3, gamma=0.0125,
                   coef0=0, C=32):
    default_params = {'degree': degree, 'gamma': gamma, 'coef0': coef0, 'C': C}

    classes = sorted(list(set(y)))
    svms = []
    cvs = np.zeros((X.shape[0], len(classes) + 1))
    folds = StratifiedKFold(fold).split(X, y)

    inds = np.array([])
    if indep.shape[0] != 0:
        inds = np.zeros((len(indep), len(classes) + 1))
        inds[:, 0] = indep[:, 0]

    prediction_result_cv = []

    for trained, valided in folds:
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]

        model = svm.SVC(C=default_params['C'], kernel=kernel, degree=default_params['degree'],
                        gamma=default_params['gamma'], coef0=default_params['coef0'], probability=True,
                        random_state=1)

        svc = model.fit(train_X, train_y)
        svms.append(svc)
        proba_ = svc.predict_proba(valid_X)
        cvs[valided, 0], cvs[valided, 1:] = valid_y, proba_

        print('save the sample label and prediction result to ndarray')
        # save the sample label and prediction result to ndarray
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, proba_
        prediction_result_cv.append(tmp_result)

        print('begin independent')
        # independent
        if indep.shape[0] != 0:
            inds[:, 1:] += svc.predict_proba(indep[:, 1:])

    header = 'C=%f\tgamma=%s' % (default_params['C'], default_params['gamma'])
    if indep.shape[0] != 0:
        inds[:, 1:] /= fold
    return header, prediction_result_cv, inds


def RF_Classifier(X, y, indep=None, fold=5, n_trees=100, out='RF_output'):
    """
    Parameters:
    ----------
    :param X: 2-D ndarray
    :param y: 1-D ndarray
    :param indep: 2-D ndarray, the first column is labels and the rest are feature values
    :param fold: int, default 5
    :param n_trees: int, number of trees, default: 5
    :param out:
    :return:
        info: str, the model parameters
        cross-validation result: list with element is ndarray
        independent result: ndarray, the first column is labels and the rest are prediction scores.
    """
    classes = sorted(list(set(y)))
    if indep.shape[0] != 0:
        indep_out = np.zeros((indep.shape[0], len(classes) + 1))
        indep_out[:, 0] = indep[:, 0]

    prediction_result_cv = []
    prediction_result_ind = np.array([])
    if indep.shape[0] != 0:
        prediction_result_ind = np.zeros((len(indep), len(classes) + 1))
        prediction_result_ind[:, 0] = indep[:, 0]

    folds = StratifiedKFold(fold).split(X, y)
    for i, (trained, valided) in enumerate(folds):
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = RandomForestClassifier(n_estimators=n_trees, bootstrap=False)
        rfc = model.fit(train_X, train_y)
        scores = rfc.predict_proba(valid_X)
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, scores
        prediction_result_cv.append(tmp_result)
        # independent
        if indep.shape[0] != 0:
            prediction_result_ind[:, 1:] += rfc.predict_proba(indep[:, 1:])
    if indep.shape[0] != 0:
        prediction_result_ind[:, 1:] /= fold
    header = 'n_trees: %d' % n_trees
    return header, prediction_result_cv, prediction_result_ind


def LR_Classifier(X, y, indep=None, fold=5, out='LR_output'):
    classes = sorted(list(set(y)))
    prediction_result_cv = []
    indep_out = np.zeros([])
    if indep.shape[0] != 0:
        indep_out = np.zeros((indep.shape[0], len(classes) + 1))
        indep_out[:, 0] = indep[:, 0]

    folds = StratifiedKFold(fold).split(X, y)
    for trained, valided in folds:
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = LogisticRegression(C=1.0, random_state=0).fit(train_X, train_y)
        scores = model.predict_proba(valid_X)
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, scores
        prediction_result_cv.append(tmp_result)
        # independent
        if indep.shape[0] != 0:
            indep_out[:, 1:] += model.predict_proba(indep[:, 1:])
    if indep.shape[0] != 0:
        indep_out[:, 1:] /= fold
    header = ''
    return header, prediction_result_cv, indep_out
