# Powered by Python 3.6

# To cancel the modifications performed by the script
# on the current graph, click on the undo button.

# Some useful keyboards shortcuts : 
#   * Ctrl + D : comment selected lines.
#   * Ctrl + Shift + D  : uncomment selected lines.
#   * Ctrl + I : indent selected lines.
#   * Ctrl + Shift + I  : unindent selected lines.
#   * Ctrl + Return  : run script.
#   * Ctrl + F  : find selected text.
#   * Ctrl + R  : replace selected text.
#   * Ctrl + Space  : show auto-completion dialog.

from tulip import tlp
from math import *

# The updateVisualization(centerViews = True) function can be called
# during script execution to update the opened views

# The pauseScript() function can be called to pause the script execution.
# To resume the script execution, you will have to click on the "Run script " button.

# The runGraphScript(scriptFile, graph) function can be called to launch
# another edited script on a tlp.Graph object.
# The scriptFile parameter defines the script name to call (in the form [a-zA-Z0-9_]+.py)

# The main(graph) function must be defined 
# to run the script on the current graph
BASESIZE=1.0
EDGE_THRESHOLD=1.0

#PARTIE 1 : fonctions
def setNodeLabelAndSize(lab,locus,size,viewSize,node):
  lab[node] = locus[node]
  viewSize[node] = tlp.Size(size,size*0.2,0.) #* (graph.deg(n)+1)

def setEdgesNodesColors(graph, n, pos, neg, color, shape, tgtShape):
  color[n] = tlp.Color.Azure
  for e in graph.getInEdges(n):
    if pos[e]:
      if neg[e]: #regulation positive et négative
        color[e] = tlp.Color.Black
        tgtShape[e] = tlp.EdgeExtremityShape.Diamond
      else: #regulation positive sur le node n
        color[e] = tlp.Color.Green
        tgtShape[e] = tlp.EdgeExtremityShape.Arrow
    else:
      if neg[e]: #regulation negative sur le node n
        color[e] = tlp.Color.Red
        tgtShape[e] = tlp.EdgeExtremityShape.Cross
      else: #aucune regulation
        color[e] = tlp.Color.Black

def applyForce(graph,layout,force="FM^3 (OGDF)"):
  param={"Unit edge length":1}
  graph.applyLayoutAlgorithm(force,layout,param)
  
#PARTIE 2 : fonctions
def getExpressionXY(graph,nX,nY,mean_X,mean_Y):  
  """
  Returns the value of " sum( val_X * val_Y) " (@see getPearsonValue formula )
  """
  temp = 0.0
  for i in range(1,18):
    temp += (graph.getDoubleProperty("tp{} s".format(i))[nX] - mean_X) * (graph.getDoubleProperty("tp{} s".format(i))[nY] - mean_Y)
  return temp

def getSquaredValue(graph,n,mean):
  """
  Returns " sum( val_N²) ", N being either the X or Y node values (@see getPearsonValue)
  """
  temp = 0.0
  for i in range(1,18):
    temp += pow(graph.getDoubleProperty("tp{} s".format(i))[n] - mean,2)
  return temp

def getMeanExpression(graph,n):  
  """
  Returns the mean of the expression levels from a given gene
  """
  tp_mean = 0.0
  for i in range(1,18):
    tp_mean += graph.getDoubleProperty("tp{} s".format(i))[n]
  tp_mean /= 17.0
  return tp_mean
  
def getPearsonValue(graph,nX,nY):
  """
  Pearson correlation inspired from https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
  The formula used for the pearson correlation is the following : 
  P = sum( val_X * val_Y ) / sqrt( sum(val_X²) )* sqrt( sum(val_Y²) )
  where : 
    i is the considered time point ( from 1 to 17 )
    val_X = x[i] - mean(x)
    x[i] being the Expression value of X for the time point i
    mean(x) being the mean of all expression levels for x
    analogously for val_Y with Y instead of X
    sum() is the sum of given value for each i betwen 1 and 17
  """
  mean_X = getMeanExpression(graph,nX)
  mean_Y = getMeanExpression(graph,nY)
  val_XY = getExpressionXY(graph,nX,nY,mean_X,mean_Y)
  val_X_squared = getSquaredValue(graph,nX,mean_X)
  val_Y_squared = getSquaredValue(graph,nY,mean_Y)
  if val_X_squared == 0.0 or val_Y_squared == 0.0 :
    return 0.0
  return val_XY / ( sqrt(val_X_squared) * sqrt(val_Y_squared) )

def setEdgesWeight(graph,n,other,poids):
  poids[graph.addEdge(n,other)]= getPearsonValue(graph,n,other)
  
#PARTIE 3 : fonction 
def placeHeatMapLine(gr,pos,nodeList):
  size = gr.getSizeProperty('viewSize')
  layout = gr.getLayoutProperty('viewLayout')
  metric = gr.getDoubleProperty('viewMetric')
  color = gr.getColorProperty('viewColor')
  """size[nodeList[pos]] = tlp.Size(1.,1.,0.)
  layout[nodeList[pos]] = tlp.Coord(0.,pos,0.)
  color[nodeList[pos]] = tlp.Color(255,255,255)"""
  for i in range(1,18):
    node = gr.addNode()
    size[node] = tlp.Size(1.,1.,0.)
    layout[node] = tlp.Coord(i,pos,0.)
    metric[node] = gr.getDoubleProperty("tp{} s".format(i))[nodeList[pos]]
  #Delete the nodes coming from the original graph
  gr.delNode(nodeList[pos])

#MAIN
def main(graph): 
  Locus = graph.getStringProperty("Locus")
  Negative = graph.getBooleanProperty("Negative")
  Positive = graph.getBooleanProperty("Positive")
  tp1_s = graph.getDoubleProperty("tp1 s")
  tp10_s = graph.getDoubleProperty("tp10 s")
  tp11_s = graph.getDoubleProperty("tp11 s")
  tp12_s = graph.getDoubleProperty("tp12 s")
  tp13_s = graph.getDoubleProperty("tp13 s")
  tp14_s = graph.getDoubleProperty("tp14 s")
  tp15_s = graph.getDoubleProperty("tp15 s")
  tp16_s = graph.getDoubleProperty("tp16 s")
  tp17_s = graph.getDoubleProperty("tp17 s")
  tp2_s = graph.getDoubleProperty("tp2 s")
  tp3_s = graph.getDoubleProperty("tp3 s")
  tp4_s = graph.getDoubleProperty("tp4 s")
  tp5_s = graph.getDoubleProperty("tp5 s")
  tp6_s = graph.getDoubleProperty("tp6 s")
  tp7_s = graph.getDoubleProperty("tp7 s")
  tp8_s = graph.getDoubleProperty("tp8 s")
  tp9_s = graph.getDoubleProperty("tp9 s")
  viewBorderColor = graph.getColorProperty("viewBorderColor")
  viewBorderWidth = graph.getDoubleProperty("viewBorderWidth")
  viewColor = graph.getColorProperty("viewColor")
  viewFont = graph.getStringProperty("viewFont")
  viewFontSize = graph.getIntegerProperty("viewFontSize")
  viewIcon = graph.getStringProperty("viewIcon")
  viewLabel = graph.getStringProperty("viewLabel")
  viewLabelBorderColor = graph.getColorProperty("viewLabelBorderColor")
  viewLabelBorderWidth = graph.getDoubleProperty("viewLabelBorderWidth")
  viewLabelColor = graph.getColorProperty("viewLabelColor")
  viewLabelPosition = graph.getIntegerProperty("viewLabelPosition")
  viewLayout = graph.getLayoutProperty("viewLayout")
  viewMetric = graph.getDoubleProperty("viewMetric")
  viewRotation = graph.getDoubleProperty("viewRotation")
  viewSelection = graph.getBooleanProperty("viewSelection")
  viewShape = graph.getIntegerProperty("viewShape")
  viewSize = graph.getSizeProperty("viewSize")
  viewSrcAnchorShape = graph.getIntegerProperty("viewSrcAnchorShape")
  viewSrcAnchorSize = graph.getSizeProperty("viewSrcAnchorSize")
  viewTexture = graph.getStringProperty("viewTexture")
  viewTgtAnchorShape = graph.getIntegerProperty("viewTgtAnchorShape")
  viewTgtAnchorSize = graph.getSizeProperty("viewTgtAnchorSize")
  updateVisualization(centerViews = True)
  #
  clusterized = False
  
  #PARTIE 1
  applyForce(graph, viewLayout, "FM^3 (OGDF)")
  for n in graph.getNodes():
    setNodeLabelAndSize(viewLabel,Locus,BASESIZE,viewSize,n)
    setEdgesNodesColors(graph, n, Positive, Negative, viewColor, 
    viewShape, viewTgtAnchorShape)
  
  graphTmp = graph.addCloneSubGraph("Parties 2 et 3")
  
  """
  #PARTIE 2
  #creation du nouveau graphe  
  graphCopy = graphTmp.addCloneSubGraph("partitionnement")
  graphCopy.delEdges(graphCopy.getEdges())
  nodeList = graphCopy.nodes()
  #ajout des poids des arretes  
  poids = graphCopy.getDoubleProperty("poids");
  for i in range(len(nodeList)):
    for j in range(len(nodeList[i+1::])):
      setEdgesWeight(graphCopy,nodeList[i],nodeList[j],poids)
  #suppression des arretes "superflues"
  for e in graphCopy.getEdges():
    poids[e] = 1-poids[e] #placing values between 0 and 2
    if poids[e] > 0.5 and poids[e] < 1.5:
      graphCopy.delEdge(e)
  #clustering  
  params = tlp.getDefaultPluginParameters('MCL Clustering', graphCopy)
  params["weights"]=poids
  clusterValue = graphCopy.getDoubleProperty('clusterValue')
  #clusterized becomes True if the algorithm processes properly
  clusterized = graphCopy.applyDoubleAlgorithm('MCL Clustering', clusterValue, params)
  """
  
  #PARTIE 3
  #creation d'un nouveau graphe pour la heatmap
  graphHeat = graphTmp.addCloneSubGraph("heat map")
  graphHeat.delEdges(graphHeat.getEdges())  
  nodeListHeat = graphHeat.nodes()
  #Tri des nodes en fonction des clusters
  if clusterized :
    nodeListHeat.sort(key=lambda x: clusterValue[x], reverse=False)
  #placement des nodes pour la Heat Map
  for i in range(len(nodeListHeat)):
    placeHeatMapLine(graphHeat,i,nodeListHeat)
  colorMappingParams = tlp.getDefaultPluginParameters('Color Mapping', graphHeat)
  colorMappingParams['color scale'] = tlpgui.ColorScalesManager.getColorScale('BiologicalHeatMap')
  colorMappingParams['input property'] = graphHeat.getDoubleProperty("viewMetric")
  success = graphHeat.applyColorAlgorithm('Color Mapping', colorMappingParams)
