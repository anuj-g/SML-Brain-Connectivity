from sklearn import *
from classify_utils import *
from plot_decision_tree import *

#Get the data
CCFX, CCFy = getXyFromCsv('./csvs/set1/CCFCombined_SVM.csv')
EBCX, EBCy = getXyFromCsv('./csvs/set1/EBCCombined_SVM.csv')
PPFX, PPFy = getXyFromCsv('./csvs/set1/PPFCombined_SVM.csv')
ConnX, ConnY = getXyFromCsv('./csvs/set1/SexCombined.csv')

CombX = CCFX[:, [24, 54, 67]]
CombX = column_stack((CombX, EBCX[:,[840]]))
CombX = column_stack((CombX, PPFX[:, [17, 60, 67]]))
CombX = column_stack((CombX, ConnX[:, [2307, 4722]]))
CombY = ConnY

print "Number of features", shape(CombX)

print "Using decision trees"

#Use 0 based index while indexing
print "Using CCF 25, 55 68"
print performClassification(CCFX[:,[24, 54, 67]], CCFy, DecisionTreeClassifier(min_samples_leaf=17), 'DT', False)

print "Using EBC 841"
performClassification(EBCX[:,[840]], EBCy, DecisionTreeClassifier(min_samples_leaf=14), 'DT', False)

print "Using PPF 18 and 61 and 68"
performClassification(PPFX[:,[17, 60, 67]], PPFy, DecisionTreeClassifier(min_samples_leaf=14), 'DT', False)

print "Using only edge connectivity between node 33 and 68"
performClassification(ConnX[:,[2307, 4722]], ConnY, DecisionTreeClassifier(min_samples_leaf=14), 'DT', False)

print "Combining all features"
performClassification(CombX, CombY, DecisionTreeClassifier(min_samples_leaf=15), 'DT', False)

print "----------------------------------------------------------------------------------------------------------------------------------"
print "Using SVM with Linear kernel"

#Use 0 based index while indexing
print "Using CCF 25, 55 68"
performClassification(CCFX[:,[24, 54, 67]], CCFy, SVC(kernel='linear', C = 1.0/4.0), 'SVC', False)

print "Using EBC 841"
performClassification(EBCX[:,[840]], EBCy, SVC(kernel='linear', C = 4.0), 'SVC', False)

print "Using PPF 18 and 61 and 68"
performClassification(PPFX[:,[17, 60, 67]], PPFy, SVC(kernel='linear', C = 1.0/32.0), 'SVC', False)

print "Using only edge connectivity between node 33 and 68"
performClassification(ConnX[:,[2307, 4722]], ConnY, SVC(kernel='linear', C = 1.0/8.0), 'SVC',  False)

print "Combining all features"

performClassification(CombX, CombY, SVC(kernel='linear', C = 1.0/16.0), 'SVC', False)

print "----------------------------------------------------------------------------------------------------------------------------------"
print "Using SVM with RBF kernel"

#Use 0 based index while indexing
print "Using CCF 25, 55 68"
performClassification(CCFX[:,[24, 54, 67]], CCFy, SVC(kernel='rbf',C = 0.5, gamma = 1.0), 'RBF', False)

print "Using EBC 841"
performClassification(EBCX[:,[840]], EBCy, SVC(kernel='rbf', C = 0.5, gamma = 0.25), 'RBF', False)

print "Using PPF 18 and 61 and 68"
performClassification(PPFX[:,[17, 60, 67]], PPFy, SVC(kernel='rbf', C = 4.0, gamma = 1.0/128.0), 'RBF', False)

print "Using only edge connectivity between node 33 and 68"
performClassification(ConnX[:,[2307, 4722]], ConnY, SVC(kernel='rbf', C = 4.0, gamma = 4.0), 'RBF',  False)

print "Combining all features"
performClassification(CombX, CombY, SVC(kernel='rbf', C = 0.25, gamma = 1.0/64.0), 'RBF', False)

print "Generating the decision tree plot"
drawDecTree(DecisionTreeClassifier(criterion='entropy', min_samples_leaf=15), CombX, CombY, 'dectree', featNames=['CCF25', 'CCF55', 'CCF68', 'EBC841', 'PPF18', 'PPF61',
'PPF68', 'CONN33_68', 'CONN68_33'], label='All_Features')

print "All Done!! Have a great day"
