from sklearn import *
from classify_utils import *
from plot_decision_tree import *
import numpy as np

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

print "Using decision trees"

X_DEC = np.arange(2, 20)
CCF_ACC_DEC = []
EBC_ACC_DEC = []
PPF_ACC_DEC = []
CONN_ACC_DEC = []
COMB_ACC_DEC = []

#Use 0 based index while indexing
print "Using CCF 25, 55 68"
for min_leaves in X_DEC:
    acc = performClassification(CCFX[:,[24, 54, 67]], CCFy, DecisionTreeClassifier(min_samples_leaf=min_leaves), 'DT', False)
    CCF_ACC_DEC.append(acc)
  
    print "Using EBC 841"
    acc = performClassification(EBCX[:,[840]], EBCy, DecisionTreeClassifier(min_samples_leaf=min_leaves), 'DT', False)
    EBC_ACC_DEC.append(acc)

    print "Using PPF 18 and 61 and 68"
    acc = performClassification(PPFX[:,[17, 60, 67]], PPFy, DecisionTreeClassifier(min_samples_leaf=min_leaves), 'DT', False)
    PPF_ACC_DEC.append(acc)
      

    print "Using only edge connectivity between node 33 and 68"
    acc = performClassification(ConnX[:,[2307, 4722]], ConnY, DecisionTreeClassifier(min_samples_leaf=min_leaves), 'DT', False)
    CONN_ACC_DEC.append(acc)

    print "Combining all features"
    acc = performClassification(CombX, CombY, DecisionTreeClassifier(min_samples_leaf=min_leaves), 'DT', False)
    COMB_ACC_DEC.append(acc)

print "----------------------------------------------------------------------------------------------------------------------------------"


print "Using SVM with linear kernel"
X_SVM = [1.0/1024, 1.0/512, 1.0/256, 1.0/128, 1.0/64, 1.0/32.0, 1.0/16, 1.0/8, 1.0/4, 1.0/2.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0]

CCF_ACC_SVM  = []
EBC_ACC_SVM  = []
PPF_ACC_SVM  = []
CONN_ACC_SVM = []
COMB_ACC_SVM = []

for c in X_SVM:
    #Use 0 based index while indexing
    print "Using CCF 25, 55 68"
    acc = performClassification(CCFX[:,[24, 54, 67]], CCFy, SVC(kernel='linear', C = c), 'LinSVC', False)
    CCF_ACC_SVM.append(acc)

    print "Using EBC 841"
    acc = performClassification(EBCX[:,[840]], EBCy, SVC(kernel='linear', C = c), 'LinSVC', False )
    EBC_ACC_SVM.append(acc)

    print "Using PPF 18 and 61 and 68"
    acc = performClassification(PPFX[:,[17, 60, 67]], PPFy, SVC(kernel='linear', C = c), 'LinSVC', False)
    PPF_ACC_SVM.append(acc)

    print "Using only edge connectivity between node 33 and 68"
    acc = performClassification(ConnX[:,[2307, 4722]], ConnY, SVC(kernel='linear', C = c), 'LinSVC', False)
    CONN_ACC_SVM.append(acc)

    print "Combining all features"
    acc = performClassification(CombX, CombY, SVC(kernel='linear', C = c), 'LinSVC', False)
    COMB_ACC_SVM.append(acc)

#Plot min leaves for decision tree
label_str =['CCF', 'EBC', 'PPF' ,'CONN', 'COMB']
pl.figure()
pl.title('Acc versus Min leaves for decision trees')
pl.plot(X_DEC, CCF_ACC_DEC, X_DEC, EBC_ACC_DEC, X_DEC, PPF_ACC_DEC, X_DEC, CONN_ACC_DEC, X_DEC, COMB_ACC_DEC)
pl.legend(label_str)

#Plot C values for Linear SVM
label_str =['CCF', 'EBC', 'PPF' ,'CONN', 'COMB']
pl.figure()
pl.title('Acc versus C for SVM(Linear)')
pl.xscale('log',basex=2)
pl.plot(X_SVM, CCF_ACC_SVM, X_SVM, EBC_ACC_SVM, X_SVM, PPF_ACC_SVM, X_SVM, CONN_ACC_SVM, X_SVM, COMB_ACC_SVM)
pl.legend(label_str)
pl.grid(True)
pl.show()

print "All Done!! Have a great day"
