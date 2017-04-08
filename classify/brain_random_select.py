from sklearn import *
from classify_utils import *
from plot_decision_tree import *
import numpy as np
import numpy.random as rng
import pylab as pl

#Code for moving average:http://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
def movingaverage(interval, window_size):
    window = numpy.ones(int(window_size))/float(window_size)
    return numpy.convolve(interval, window, 'same')


#Get the data
CCFX, CCFy = getXyFromCsv('./csvs/set2/CCFCombined_SVM.csv')
EBCX, EBCy = getXyFromCsv('./csvs/set2/EBCCombined_SVM.csv')
PPFX, PPFy = getXyFromCsv('./csvs/set2/PPFCombined_SVM.csv')
ConnX, ConnY = getXyFromCsv('./csvs/set2/SexCombined.csv')

CCF_ACC_DEC  = []
EBC_ACC_DEC  = []
PPF_ACC_DEC  = []
CONN_ACC_DEC = []

CCF_ACC_SVM_LINEAR  = []
EBC_ACC_SVM_LINEAR  = []
PPF_ACC_SVM_LINEAR  = []
CONN_ACC_SVM_LINEAR = []

CCF_ACC_SVM_RBF  = []
EBC_ACC_SVM_RBF  = []
PPF_ACC_SVM_RBF  = []
CONN_ACC_SVM_RBF = []

numIterations = 100
for i in np.arange(numIterations):
	CCF_SELECTED = rng.choice(70, 3) #Select 3 random nodes for CCF features
	EBC_SELECTED = rng.choice(4900, 1) #Select a random edge for EBC
	PPF_SELECTED = rng.choice(70, 3) # Select 3 random nodes for PPF

	#Select a bidirectional edge
	v1 = rng.randint(0,70)
	v2 = rng.randint(0,70) 
	e1 = v1 * 70 + v2;
	e2 = v2 * 70 + v1
	CONN_SELECTED = [e1, e2]

	CombX = CCFX[:, CCF_SELECTED]
	CombX = column_stack((CombX, EBCX[:, EBC_SELECTED]))
	CombX = column_stack((CombX, PPFX[:, PPF_SELECTED]))
	CombX = column_stack((CombX, ConnX[:, CONN_SELECTED]))
	CombY = ConnY

	print "Number of features", shape(CombX)

	print "Using decision trees"

	#Use 0 based index while indexing
	print "Using CCF", CCF_SELECTED
	acc = performClassification(CCFX[:,CCF_SELECTED], CCFy, DecisionTreeClassifier(), 'DT', gridSearch=True)
        CCF_ACC_DEC.append(acc)

	print "Using EBC", EBC_SELECTED
	acc = performClassification(EBCX[:,EBC_SELECTED], EBCy, DecisionTreeClassifier(), 'DT', gridSearch=True)
        EBC_ACC_DEC.append(acc)

	print "Using PPF", PPF_SELECTED
	acc = performClassification(PPFX[:, PPF_SELECTED], PPFy, DecisionTreeClassifier(), 'DT', gridSearch=True)
        PPF_ACC_DEC.append(acc)

	print "Using only edge connectivity", CONN_SELECTED
	acc = performClassification(ConnX[:, CONN_SELECTED], ConnY, DecisionTreeClassifier(), 'DT', gridSearch=True)
        CONN_ACC_DEC.append(acc)

	print "----------------------------------------------------------------------------------------------------------------------------------"
	print "Using SVM with linear kernel"

	#Use 0 based index while indexing
	print "Using CCF", CCF_SELECTED
	acc = performClassification(CCFX[:,CCF_SELECTED], CCFy, SVC(kernel='linear'), 'LinSVC', gridSearch=True)
        CCF_ACC_SVM_LINEAR.append(acc)

	print "Using EBC", EBC_SELECTED
	acc = performClassification(EBCX[:,EBC_SELECTED], EBCy, SVC(kernel='linear'), 'LinSVC', gridSearch=True)
        EBC_ACC_SVM_LINEAR.append(acc)

	print "Using PPF", PPF_SELECTED
	acc = performClassification(PPFX[:,PPF_SELECTED], PPFy, SVC(kernel='linear'), 'LinSVC', gridSearch=True)
        PPF_ACC_SVM_LINEAR.append(acc)

	print "Using only edge connectivity", CONN_SELECTED
	acc = performClassification(ConnX[:,CONN_SELECTED], ConnY, SVC(kernel='linear'), 'LinSVC', gridSearch=True)
        CONN_ACC_SVM_LINEAR.append(acc)
	print "----------------------------------------------------------------------------------------------------------------------------------"

	print "Using SVM with RBF kernel"

	#Use 0 based index while indexing
	print "Using CCF", CCF_SELECTED
	acc = performClassification(CCFX[:,CCF_SELECTED], CCFy, SVC(kernel='rbf'), 'RBF', gridSearch=True)
        CCF_ACC_SVM_RBF.append(acc)

	print "Using EBC", EBC_SELECTED
	acc = performClassification(EBCX[:,EBC_SELECTED], EBCy, SVC(kernel='rbf'), 'RBF', gridSearch=True)
        EBC_ACC_SVM_RBF.append(acc)

	print "Using PPF", PPF_SELECTED
	acc = performClassification(PPFX[:,PPF_SELECTED], PPFy, SVC(kernel='rbf'), 'RBF', gridSearch=True)
        PPF_ACC_SVM_RBF.append(acc)

	print "Using only edge connectivity", CONN_SELECTED
	acc = performClassification(ConnX[:,CONN_SELECTED], ConnY, SVC(kernel='rbf'), 'RBF', gridSearch=True)
        CONN_ACC_SVM_RBF.append(acc)
	print "----------------------------------------------------------------------------------------------------------------------------------"

luck = luck_classifier(CCFy)
luck_y = luck * np.ones(numIterations)

X_AXIS = np.arange(numIterations)
Y_DEC  = np.zeros((numIterations, 5))
Y_SVM_LINEAR  = np.zeros((numIterations, 5))
Y_SVM_RBF  = np.zeros((numIterations, 5))

windowSize = 20
Y_DEC[:,0] = movingaverage(CCF_ACC_DEC,  windowSize)
Y_DEC[:,1] = movingaverage(EBC_ACC_DEC,  windowSize)
Y_DEC[:,2] = movingaverage(PPF_ACC_DEC,  windowSize)
Y_DEC[:,3] = movingaverage(CONN_ACC_DEC, windowSize)
Y_DEC[:,4] = luck_y

Y_SVM_LINEAR[:,0] = movingaverage(CCF_ACC_SVM_LINEAR,  windowSize)
Y_SVM_LINEAR[:,1] = movingaverage(EBC_ACC_SVM_LINEAR,  windowSize)
Y_SVM_LINEAR[:,2] = movingaverage(PPF_ACC_SVM_LINEAR,  windowSize)
Y_SVM_LINEAR[:,3] = movingaverage(CONN_ACC_SVM_LINEAR, windowSize)
Y_SVM_LINEAR[:,4] = luck_y

Y_SVM_RBF[:,0] = movingaverage(CCF_ACC_SVM_RBF,  windowSize)
Y_SVM_RBF[:,1] = movingaverage(EBC_ACC_SVM_RBF,  windowSize)
Y_SVM_RBF[:,2] = movingaverage(PPF_ACC_SVM_RBF,  windowSize)
Y_SVM_RBF[:,3] = movingaverage(CONN_ACC_SVM_RBF, windowSize)
Y_SVM_RBF[:,4] = luck_y

labelstr = ['CCF', 'EBC', 'PPF', 'CONN', 'Luck']

pl.figure()
pl.title('Performance of random features on classification tasks using Decision Trees')
pl.plot(X_AXIS, Y_DEC)
pl.legend(labelstr)

pl.figure()
pl.title('Performance of random features on classification tasks using SVM(Linear kernel)')
pl.plot(X_AXIS, Y_SVM_LINEAR)
pl.legend(labelstr)

pl.figure()
pl.title('Performance of random features on classification tasks using SVM(RBF kernel)')
pl.plot(X_AXIS, Y_SVM_RBF)
pl.legend(labelstr)
pl.show()

print "All Done!! Have a great day"
