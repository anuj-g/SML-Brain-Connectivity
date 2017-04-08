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
perform_permutation_test(CCFX[:,[24, 54, 67]], CCFy, DecisionTreeClassifier(min_samples_leaf=17), 'DT', False, filename = 'perm_test_plots/CCF_25_55_68_dec.png')

print "Using EBC 841"
perform_permutation_test(EBCX[:,[840]], EBCy, DecisionTreeClassifier(min_samples_leaf=14), 'DT', False, filename = 'perm_test_plots/EBC_841_dec.png')

print "Using PPF 18 and 61 and 68"
perform_permutation_test(PPFX[:,[17, 60, 67]], PPFy, DecisionTreeClassifier(min_samples_leaf=14), 'DT', False, filename = 'perm_test_plots/PPF_18_61_68_dec.png')

print "Using only edge connectivity between node 33 and 68"
perform_permutation_test(ConnX[:,[2307, 4722]], ConnY, DecisionTreeClassifier(min_samples_leaf=14), 'DT', False, filename = 'perm_test_plots/CONN_33_68_dec.png')

print "Combining all features"
perform_permutation_test(CombX, CombY, DecisionTreeClassifier(min_samples_leaf = 15), 'DT', False, filename = 'perm_test_plots/All_dec.png')

print "----------------------------------------------------------------------------------------------------------------------------------"
print "Using SVM with linear kernel"

#Use 0 based index while indexing
print "Using CCF 25, 55 68"
perform_permutation_test(CCFX[:,[24, 54, 67]], CCFy, SVC(kernel='linear', C = 1.0/4.0), 'LinSVC', False, filename = 'perm_test_plots/CCF_25_55_68_svm_linear.png')

print "Using EBC 841"
perform_permutation_test(EBCX[:,[840]], EBCy, SVC(kernel='linear', C = 4.0), 'LinSVC', False , filename = 'perm_test_plots/EBC_841_svm_linear.png')

print "Using PPF 18 and 61 and 68"
perform_permutation_test(PPFX[:,[17, 60, 67]], PPFy, SVC(kernel='linear', C = 1.0/32.0), 'LinSVC', False, filename = 'perm_test_plots/PPF_18_61_68_svm_linear.png')

print "Using only edge connectivity between node 33 and 68"
perform_permutation_test(ConnX[:,[2307, 4722]], ConnY, SVC(kernel='linear', C = 1.0/8.0), 'LinSVC', False, filename = 'perm_test_plots/CONN_33_68_svm_linear.png')

print "Combining all features"
perform_permutation_test(CombX, CombY, SVC(kernel='linear', C = 1.0/16.0), 'LinSVC', False, filename = 'perm_test_plots/All_svm_linear.png')

print "----------------------------------------------------------------------------------------------------------------------------------"
print "Using SVM with RBF kernel"

#Use 0 based index while indexing
print "Using CCF 25, 55 68"
perform_permutation_test(CCFX[:,[24, 54, 67]], CCFy, SVC(kernel='rbf', C = 1.0, gamma = 0.5), 'RBF', False, filename = 'perm_test_plots/CCF_25_55_68_svm_rbf.png')

print "Using EBC 841"
perform_permutation_test(EBCX[:,[840]], EBCy, SVC(kernel='rbf', C = 0.5, gamma = 0.25), 'RBF', False , filename = 'perm_test_plots/EBC_841_svm_rbf.png')

print "Using PPF 18 and 61 and 68"
perform_permutation_test(PPFX[:,[17, 60, 67]], PPFy, SVC(kernel='rbf', C = 4.0, gamma = 1.0/128.0), 'RBF', False, filename = 'perm_test_plots/PPF_18_61_68_svm_rbf.png')

print "Using only edge connectivity between node 33 and 68"
perform_permutation_test(ConnX[:,[2307, 4722]], ConnY, SVC(kernel='rbf', C = 4.0, gamma = 4.0), 'RBF', False, filename = 'perm_test_plots/CONN_33_68_svm_rbf.png')

print "Combining all features"
perform_permutation_test(CombX, CombY, SVC(kernel='rbf', C = 0.25, gamma = 1.0/64.0), 'RBF', False, filename = 'perm_test_plots/All_svm_rbf.png')

print "All Done!! Have a great day"
