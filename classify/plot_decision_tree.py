from numpy import *
from sklearn import tree
from random import randint
from sklearn.tree import DecisionTreeClassifier
from classify_utils import *
import StringIO, pydot

def drawDecTree(decTree, X, y, outdir, label=randint(100), featNames=None):
    decTree.fit(X, y)
    print decTree.feature_importances_
    dot_data = StringIO.StringIO()
    tree.export_graphviz(decTree, out_file=dot_data, feature_names=featNames)
    graph = pydot.graph_from_dot_data(dot_data.getvalue())
    graph.write_png(outdir +  "/" + str(label) + "_graph" + ".png")

#drawDecTree(DecisionTreeClassifier(), X, y, 'dectreeimg'))
#drawDecTree(DecisionTreeClassifier(compute_importances=True, criterion='entropy', min_samples_leaf=105, max_depth=5), X, y, 'dectreeimg', featNames=featNames)
