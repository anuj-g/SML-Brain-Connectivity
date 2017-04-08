from sklearn import *
from classify_utils import *
from plot_decision_tree import *
import numpy as np
import pylab as pl

#Get the data
CCFX, CCFy = getXyFromCsv('./csvs/set2/CCFCombined_SVM.csv')
EBCX, EBCy = getXyFromCsv('./csvs/set2/EBCCombined_SVM.csv')
PPFX, PPFy = getXyFromCsv('./csvs/set2/PPFCombined_SVM.csv')
ConnX, ConnY = getXyFromCsv('./csvs/set2/SexCombined.csv')

CombX = CCFX[:, [24, 54, 67]]
CombX = column_stack((CombX, EBCX[:,[840]]))
CombX = column_stack((CombX, PPFX[:, [17, 60, 67]]))
CombX = column_stack((CombX, ConnX[:, [2307, 4722]]))
CombY = ConnY

print "Number of features", shape(CombX)

CRANGE = 2.0 ** np.arange(-5, 5)
GAMMARANGE = 2.0 ** np.arange(-7, 3)

NUMROWS = len(CRANGE)
NUMCOLS = len(GAMMARANGE)

R_CCF  = np.zeros((NUMROWS, NUMCOLS))
R_EBC  = np.zeros((NUMROWS, NUMCOLS))
R_PPF  = np.zeros((NUMROWS, NUMCOLS))
R_CONN = np.zeros((NUMROWS, NUMCOLS))
R_COMB = np.zeros((NUMROWS, NUMCOLS))

for c in np.arange(0, len(CRANGE)):
    for g in np.arange(0, len(GAMMARANGE)):
        print c, g
        cval = CRANGE[c]
        gval = GAMMARANGE[g]
        acc = performClassification(CCFX[:,[24, 54, 67]], CCFy, SVC(kernel='rbf', C = cval, gamma = gval), 'RBF', False)
        R_CCF[c][g]  = acc

        acc = performClassification(EBCX[:,[840]], EBCy, SVC(kernel='rbf', C = cval, gamma = gval), 'RBF', False)
        R_EBC[c][g]  = acc

        acc = performClassification(PPFX[:,[17, 60, 67]], PPFy, SVC(kernel='rbf', C = cval, gamma = gval), 'RBF', False)
        R_PPF[c][g]  = acc

        acc = performClassification(ConnX[:,[2307, 4722]], ConnY, SVC(kernel='rbf', C = cval, gamma = gval), 'RBF',  False)
        R_CONN[c][g] = acc

        acc = performClassification(CombX, CombY, SVC(kernel='rbf', C = cval, gamma = gval), 'RBF', False)
        R_COMB[c][g] = acc

w, h = pl.figaspect(1.)
pl.figure(figsize = (w,h))
pl.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
pl.title('C and Gamma for CCF')
pl.imshow(R_CCF, interpolation = 'nearest', cmap = pl.cm.spectral)
pl.xlabel('gamma')
pl.ylabel('C')
pl.colorbar()
pl.xticks(np.arange(NUMCOLS), GAMMARANGE, rotation = 45)
pl.yticks(np.arange(NUMROWS), CRANGE, rotation = 45)
pl.savefig('hyperparameters_estimate_plots/CCF_RBF_SET2.png')

pl.figure(figsize = (8,8))
pl.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
pl.title('C and Gamma for EBC')
pl.imshow(R_EBC, interpolation = 'nearest', cmap = pl.cm.spectral)
pl.xlabel('gamma')
pl.ylabel('C')
pl.colorbar()
pl.xticks(np.arange(NUMCOLS), GAMMARANGE, rotation = 45)
pl.yticks(np.arange(NUMROWS), CRANGE, rotation = 45)
pl.savefig('hyperparameters_estimate_plots/EBC_RBF_SET2.png')

pl.figure(figsize = (8,8))
pl.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
pl.title('C and Gamma for PPF')
pl.imshow(R_PPF, interpolation = 'nearest', cmap = pl.cm.spectral)
pl.xlabel('gamma')
pl.ylabel('C')
pl.colorbar()
pl.xticks(np.arange(NUMCOLS), GAMMARANGE, rotation = 45)
pl.yticks(np.arange(NUMROWS), CRANGE, rotation = 45)
pl.savefig('hyperparameters_estimate_plots/PPF_RBF_SET2.png')


pl.figure(figsize = (8,8))
pl.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
pl.title('C and Gamma for CONN')
pl.imshow(R_CONN, interpolation = 'nearest', cmap = pl.cm.spectral)
pl.xlabel('gamma')
pl.ylabel('C')
pl.colorbar()
pl.xticks(np.arange(NUMCOLS), GAMMARANGE, rotation = 45)
pl.yticks(np.arange(NUMROWS), CRANGE, rotation = 45)
pl.savefig('hyperparameters_estimate_plots/CONN_RBF_SET2.png')

pl.figure(figsize = (8,8))
pl.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
pl.title('C and Gamma for COMB')
pl.imshow(R_COMB, interpolation = 'nearest', cmap = pl.cm.spectral)
pl.xlabel('gamma')
pl.ylabel('C')
pl.colorbar()
pl.xticks(np.arange(NUMCOLS), GAMMARANGE, rotation = 45)
pl.yticks(np.arange(NUMROWS), CRANGE, rotation = 45)
pl.savefig('hyperparameters_estimate_plots/COMB_RBF_SET2.png')

pl.show()

print "All Done!! Have a great day"
