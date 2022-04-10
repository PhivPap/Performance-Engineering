from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker

import numpy as np
import sys
import pylab
import re
import json
import math

#Define the colors of the benchmark icons
#The returns are codes for different colors
def get_color(point):
  if point["label"] == "Sequential":
    return "C3"
  elif point["label"] == "Parallel":
    return "C4"
    print("need new color for \"" +  point["label"] + "\"")
    sys.exit()

#Define the icons for the benchmark
#The returns are codes for different colors
def get_marker(point):
  if point["prop"][1] == "Sequential":
    return "X"
  elif point["prop"][1] == "Parallel":
    return "*"
  else:
    print("need new marker for \"" +  point["prop"][1] + "\"")
    sys.exit()

#Define the size of the benchmark icons
def get_markersize(point):
  base=30
  if point["prop"][1] == "Sequential":
    return base*1.2
  elif point["prop"][1] == "Parallel":
    return base*1.2
  else:
    print("need new markersize for \"" +  point["prop"][1] + "\"")
    sys.exit()

#Defines statical vertical lines for thresholds (This is not necessarly necessary at least for now)
def str_benchmark(benchmark):
  if benchmark == "taylor-green3D" or benchmark == "taylor-green2D":
    return "Taylor Green"
  elif benchmark == "couette-flow":
    return "Couette Flow"

#Canvas 
fig = plt.figure(figsize=(10, 40))
ax = plt.subplot(1,1,1)

def set_size(w,h, ax=None):
  """ w, h: width, height in inches """
  if not ax: ax=plt.gca()
  l = ax.figure.subplotpars.left
  r = ax.figure.subplotpars.right
  t = ax.figure.subplotpars.top
  b = ax.figure.subplotpars.bottom
  figw = float(w)/(r-l)
  figh = float(h)/(t-b)
  ax.figure.set_size_inches(figw, figh)

ax.grid(color="#dddddd", zorder=-1)
#Variable names
ax.set_xlabel("Arithmetic Intensity (FLOP/Byte)", fontsize=15)
ax.set_ylabel("Performance (GFLOP/s)", fontsize=15)

#Architecture Specifics
cpu_roofs = [
  {"name" : "Max single-core CPU performance ",    "val" : 38.4}
]

mem_bottlenecks = [
  {"name" : "Max Theoretical Bandwidth",     "val" : 59}
]

#Plot Settings
# Roofline window: Rect (x0, x1, y0, y1)
window_rect = [0, 10, 0, 45]

#Threshold Location
AI_v = {"taylor-green3D" : 6, "couette-flow" : 3}

#Data Points
datapoints = [
    {"type" : "taylor-green3D", "label" : "Sequential", "prop" : ["D3Q27", "Sequential", 64,  192,  0], "GFLOPs" : 10},
    {"type" : "couette-flow", "label" : "Parallel", "prop" : ["D3Q27", "Parallel", 64,  192,  0], "GFLOPs" : 5}
    ]

max_roof = cpu_roofs[0]["val"]
max_slip = mem_bottlenecks[0]["val"]

for roof in cpu_roofs:
  if roof["val"] > max_roof:
    max_roof = roof["val"]

for slip in mem_bottlenecks:
  print(slip)


#LOGARITHMIC SCALE
  y = [0, max_roof]
  x = [yy/slip["val"] for yy in y]
  ax.loglog(x, y, linewidth=1.0,
    linestyle='-.',
    marker="2",
    color="grey",
    zorder=10)

#LINEAR SCALE
  # y = [0, max_roof]
  # x = [yy/slip["val"] for yy in y]
  # ax.plot(x, y, linewidth=1.0,
  #   linestyle='-.',
  #   marker="2",
  #   color="grey",
  #   zorder=10)

  ax.annotate(slip["name"] + ": " + str(slip["val"]) + " GB/s", (window_rect[0]*1.1, window_rect[0]*slip["val"]*1.35),
    rotation="23.3",
    fontsize=11,
    ha="left", va='bottom',
    color="grey")
    
  if slip["val"] > max_slip:
    max_slip = slip["val"]

for roof in cpu_roofs:
  print(roof)

#LOGARITHMIC SCALE
  x = [roof["val"]/max_slip, window_rect[1]*10]
  ax.loglog(x, [roof["val"] for xx in x], linewidth=1.0,
    linestyle='-.',
    marker="2",
    color="grey",
    zorder=10)

#LINEAR SCALE
  # x = [roof["val"]/max_slip, window_rect[1]*10]
  # ax.loglog(x, [roof["val"] for xx in x], linewidth=1.0,
  #   linestyle='-.',
  #   marker="2",
  #   color="grey",
  #   zorder=10)

  ax.text(
    window_rect[1]/1.1, roof["val"]*1.1,
    roof["name"] + ": " + str(roof["val"]) + " GFLOPs",
    ha="right",
    fontsize=11,
    color="grey")

for benchmark in AI_v:
  AI = AI_v[benchmark]
  print(benchmark)
  print(AI)

  plt.axvline(x=AI, dashes=[10, 10, 3, 10], linewidth=0.4, color="#aaaaaa")

  ax.text(
    AI/1.15, window_rect[2]*1.24,
    str_benchmark(benchmark),
    fontsize=12,
    rotation="90",
    va="bottom",
    color="k")

i = 0
for point in datapoints:
  AI = AI_v[point["type"]]

  ax.scatter(AI, point["GFLOPs"], marker=get_marker(point), s=get_markersize(point), color=get_color(point),
    zorder=100)

  if i < 4:
    ax.plot([], [], linestyle="", marker=get_marker(point), ms=get_markersize(point)/6, color="black", label=point["prop"][1])
  if i % 4 == 0:
    ax.plot([], [], linestyle="", marker="s", color=get_color(point), label=point["label"])

  i += 1
set_size(6,3)

ax.set_xlim(window_rect[0], window_rect[1])
ax.set_ylim(window_rect[2], window_rect[3])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(0,0)))).format(y)))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(0,0)))).format(y)))


handles,labels = ax.get_legend_handles_labels()
handles = [handles[0], handles[1]]
labels = [labels[0], labels[1]]
plt.figlegend(handles=handles,  labels=labels, loc="upper center", bbox_to_anchor=(0.7, 0.97), ncol=4, fontsize=10)

handles,labels = ax.get_legend_handles_labels()
handles = [handles[2], handles[1]]
labels = [labels[2], labels[1]]
plt.figlegend(handles=handles,  labels=labels, loc="lower right", bbox_to_anchor=(0.96, 0.38), fontsize=11)

plt.tight_layout()
plt.show()

pp = PdfPages("finalversion")
pp.savefig(fig)
pp.close()




