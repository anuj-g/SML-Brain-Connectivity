#Jagat 11/04/13
#Runs 5 fold cross validation using linear SVM 
#68.5% accuracy
#Args: X file csv, y file
import sklearn
import numpy as np

from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn import grid_search

from sklearn.cross_validation import StratifiedKFold
from sklearn.cross_validation import ShuffleSplit
from sklearn.cross_validation import LeaveOneOut

from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB

from numpy import *
from pylab import *
from matplotlib import *
from sklearn.preprocessing import scale
from sklearn.metrics import *
import pylab as pl
import sys
import uuid


    
def getXyFromCsv(X_csv, y_csv=None):
    X = matrix(np.genfromtxt(X_csv, delimiter=","))
    if y_csv == None:
        y = array(X[:, 0].T)[0]
        X = X[:, 1:]
    else:
        y = np.genfromtxt(y_csv, delimiter=",")
    return array(X), y

def gridCVParams(numFeats, numSamples):
    paramMap = {}
    cvals = []
    cvalstart = 1/1024.0
    while cvalstart < 1024:
        cvals.append(cvalstart)
        cvalstart = cvalstart * 2;

    C_range = 2.0 ** np.arange(-5, 5)
    gamma_range = 2.0 ** np.arange(-7, 3)

    paramMap['LinSVC'] = {'C':cvals}
    paramMap['SVC'] = {'C':cvals}
    paramMap['DT'] = {'compute_importances':[True], 'min_samples_leaf':arange(2, numSamples-1, 1) }
    paramMap['RBF'] = {'C':C_range, 'gamma': gamma_range }
    return paramMap


def cv(X, y, clf, gridSearch=False, gridSearchParams=None, score_type=None, numFolds=10, leave_one_out=False):
    
    if(gridSearch):
        print "Performing a Grid Search"
        clf = grid_search.GridSearchCV(clf, gridSearchParams, refit=True)
        clf.fit(X, y)
        print "Best Estimator"
        print clf.best_estimator_

    if (leave_one_out == True):
        #print "Using Leave One Out Cross Validation"
        skf = LeaveOneOut(len(y))
    else:
        #print "Using Stratified K Fold"
        skf = StratifiedKFold(y, numFolds)

    accuracies = sklearn.cross_validation.cross_val_score(clf, X, y, cv=skf, n_jobs=-1, score_func=score_type)
    if score_type == None:

       if leave_one_out == False:
          print "Mean", np.mean(accuracies), "Standard Deviation", np.std(accuracies)
       else:
          print "Mean", np.mean(accuracies)

       return np.mean(accuracies)
    else:
      return accuracies

def performClassification(X, y, clf, classifierKey='DT', gridSearch=False, params = None, score_type=None, loo=True):
    numSamps = shape(y)
    numFeats = shape(X)[1]
    if(params == None):
        params = gridCVParams(numFeats, numSamps[0])[classifierKey]
    acc = cv(scale(X), y, clf, gridSearch=gridSearch, \
            gridSearchParams=params, score_type=score_type, leave_one_out=loo)
    return acc

def luck_classifier(y):
    num_ones = numpy.sum(y)
    num_zeros = len(y) - sum(y)
    return float(np.maximum(num_ones, num_zeros))/len(y)

def permutation_test_helper(X, y, clf, gridSearch=False, gridSearchParams=None, score_type=None, numFolds=10, leave_one_out=False, plotfilename=""):
    
    if(gridSearch):
        print "Performing a Grid Search"
        clf = grid_search.GridSearchCV(clf, gridSearchParams, refit=True)
        clf.fit(X, y)
        print "Best params"
        print clf.best_params_

    if (leave_one_out == True):
        print "Using Leave One Out Cross Validation"
        skf = LeaveOneOut(len(y))
    else:
        print "Using Stratified K Fold"
        skf = StratifiedKFold(y, numFolds)

    score, permutation_scores, pvalue = sklearn.cross_validation.permutation_test_score(clf, X, y, zero_one_score, cv=skf, n_permutations=1000, n_jobs = 1)

    print "Classification Score %s (pvalue: %s)" % (score, pvalue)

    # View histogram of permutation scores
    pl.hist(permutation_scores, 20, label='Permutation scores')
    ylim = pl.ylim()
    pl.plot(2 * [score], ylim, '--g', linewidth=3,
            label='Classification Score'
            ' (pvalue %s)' % round(pvalue, 3))

    luck_acc = luck_classifier(y)
    pl.plot(2 * [luck_acc], ylim, '--k', linewidth=3, label='Luck')

    pl.ylim(ylim)
    pl.legend()
    pl.xlabel('Score')
    if (plotfilename == ""):
        plotfilename = str(uuid.uuid4())
    pl.savefig(plotfilename)
    pl.close()

def perform_permutation_test(X, y, clf, classifierKey='DT', gridSearch=False, params = None, score_type=None, loo=True, filename=""):
    numSamps = shape(y)
    numFeats = shape(X)[1]
    if(params == None):
        params = gridCVParams(numFeats, numSamps[0])[classifierKey]
    permutation_test_helper(scale(X), y, clf, gridSearch=gridSearch, \
            gridSearchParams=params, score_type=score_type, leave_one_out=loo, plotfilename = filename)
