from ace.cace import cace
from sklearn.metrics import r2_score
from sklearn.cross_validation import KFold,ShuffleSplit
import copy
import numpy as np
import pandas
import sys

def _encode_categories(s):
    vals = s.unique()
    rd = dict(zip(np.concatenate(([0], vals)), range(len(vals)+1)))
    s = s.fillna(0).apply(lambda x: rd[x])
    return (s.values.astype(float), len(vals))


def _impute_na(df,col):
    is_null = pandas.isnull(df[col])
    if np.any(is_null):
        m = df[col].median()
        if np.isnan(m):
            m=0
        df['__NA_IND__'] = is_null.astype(float)
        df[col] = df[col].fillna(m)


def r2_scorew(a, b, w):
    ''' default metric doesn't support weights do
        drop them for r2_score '''
    return r2_score(a, b)


def ace(x, target, cat_cols, cv=True, K=0, metric=r2_scorew, metric_dir=-1,
        weight=None, balance_weights=False):
    """ Calculate ACE for each column of a DataFrame
    against the target variable

    Parameters
    ----------
    x : DataFrame
        data
    target : str
        response variable name
    cat_cols : list
        contains names of categorical columns
    cv : bool
        If True ACE score is cacluated on validation set(s).
    K : float
        Parameter for Buhlmann credibility. -1 means automatic.
    metric : func
        Function for scoring predictions
    metric_dir : int
        -1 = higher is better, 1 = lower is better
    weight : string
        name of weight col in x
    balance_weights : bool
        If True, make the sum of weights for target classes equal

    Returns
    -------
    list of lists
        List where each item is a list containing the correlations between each column
        of `x` and `target` and the correlation of a mean-only prediction and
        the target
    """

    tol = 0.01
    np.random.seed(1234)

    # drop rows with N/A targets
    x = x[x[target].notnull()]
    x.reset_index(drop=True,inplace=True)
    rows = x.shape[0]
    pos_scale = 1
    if target in cat_cols and weight is not None and balance_weights:
        pos_scale = x[weight][x[target]==0].sum() / x[weight][x[target]==1].sum()

    # Choose # of CV folds
    if not cv or rows > 16000:
        folds = 1
    elif rows > 8000:
        folds = 2
    else:
        folds = 3
    #limit size of dataset to 25000 rows
    if rows > 25000:
        x = x.iloc[np.random.choice(np.arange(rows),25000,False),:]
    unique_Y = len(np.unique(x[target]))

    # inpute NA and encode categries as intergers
    for col in x:
        if col != target:
            is_cat = col in cat_cols
            if x[col].dtype == 'object' or is_cat:
                x[col], cc = _encode_categories(x[col])
            else:
                _impute_na(x,col)

    # Create data sets
    train_data = []
    test_data = []
    wtrain = []
    if not cv:
        train_data.append(x)
    elif folds == 1:
        train,test = next(iter(ShuffleSplit(x.shape[0],1,0.5,random_state=1234)))
        train_data.append(x.iloc[train])
        test_data.append(x.iloc[test])
    else:
        kf = iter(KFold(x.shape[0],folds,shuffle=True,random_state=1234))
        for f in range(folds):
            train, test = next(kf)
            train_data.append(x.iloc[train])
            test_data.append(x.iloc[test])

    rsqs = []
    # loop over columns in data
    for col in x:
        if col == weight and x.shape[1]==3 or col == '__NA_IND__':
            continue
        # no score for target var or if target is constant
        if col == target or unique_Y < 2: # or (col != 'A_1' and col != "A.1"):
            rsqs.append([0,1])
            continue
        yhat = []
        act = []
        ym = []
        vweight = []
        for f in range(folds):
            if '__NA_IND__' in x:
                # NOTE: The column order matters.  Having the N/A indicator first seems to work better.
                xtrain = train_data[f][['__NA_IND__',col]]
            else:
                xtrain = train_data[f][col]
            is_cat = col in cat_cols
            if is_cat:
                if '__NA_IND__' in x:
                    l = np.array([5, 5, 1])
                else:
                    l = np.array([5, 1])
                adj = cc * 1.0 / xtrain.shape[0]
            else:
                if '__NA_IND__' in x:
                    l = np.array([5, 1, 1])
                else:
                    l = np.array([1, 1])
                adj = 0
            if target in cat_cols:
                if '__NA_IND__' in x:
                    l[2] = 5
                else:
                    l[1] = 5
            else:
                if '__NA_IND__' in x:
                    l[2] = 1
                else:
                    l[1] = 1
            ytrain = train_data[f][target]
            if weight is not None:
                wtrain = train_data[f][weight].values
                if balance_weights and target in cat_cols:
                    wtrain[ytrain==1] = wtrain[ytrain==1] * pos_scale
                wtrain = np.reshape(wtrain,(-1,1)).astype(float)
            else:
                wtrain = None
            ca = cace()
            if '__NA_IND__' in x:
                xt = xtrain.values.astype(float)
            else:
                xt = np.reshape(xtrain.values,(-1,1)).astype(float)
            rsq = ca.fit(xt, ytrain.values.astype(float),
                         l.tolist(), weights=wtrain)[0]
            if cv:
                if '__NA_IND__' in x:
                    xtest = test_data[f][['__NA_IND__', col]].values
                else:
                    xtest = np.reshape(test_data[f][col].values, (-1,1)).astype(float)
                ytest = test_data[f][target]
                if weight is not None:
                    wtest = test_data[f][weight]
                    ym_ = np.average(ytrain.values.flatten(), weights=wtrain.flatten())
                else:
                    ym_ = ytrain.mean()
                ym.extend(np.repeat([ym_],ytest.shape[0]))
                adj = 0
                pred = ca.predict(xtest, K)
                if target in cat_cols:
                    pred = np.clip(pred, 0.000001, 0.999999)
                yhat.extend(pred)
                act.extend(ytest.values.tolist())
                if weight is not None:
                    vweight.extend(wtest.tolist())
                #np.savetxt('/home/glen/test'+col+'_'+str(f)+'_'+str(K)+'.out',np.column_stack((xtest,ca.predict(np.reshape(xtest.values,(-1,1)).astype(float), K),ytest)),fmt='%f')
        if cv:
            act = np.array(act)
            yhat = np.array(yhat)
            ym = np.array(ym)
            if vweight == []:
                vweight = None
            null_score = metric(act, ym, vweight)
            rsq = metric(act, yhat, vweight)
        else:
            if weight is None:
                weights = None
            else:
                weights = train_data[0][weight].values.flatten()
            null_score = metric(train_data[0][target],
                                np.repeat(np.average(train_data[0][target].values.flatten(),
                                                     weights=weights),
                                          train_data[0].shape[0]), weights)

        rsqs.append([rsq - adj, null_score])
    return rsqs
