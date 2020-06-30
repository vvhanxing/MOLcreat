import numpy as np
import matplotlib.pyplot as plt
import csv

with open("gdb9.sdf_gdb_129135removed_.csv","r") as csvfile:
        reader=csv.reader(csvfile)
       
        data = list(reader)[1:133884]
        gap_data = [ float( x[8]) for x in data]
        lumo_data = [ float( x[4]) for x in data]
#print(gap_data[0])  
# Fixing random state for reproducibility
#np.random.seed(19680801)

# some random data
x =  gap_data[:133882]
y = lumo_data[:133882]

Dx = 200

def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, alpha=0.1)

    # now determine nice limits by hand:

    
    xmax = np.max(np.abs(x))
    ymax = np.max(np.abs(y))
    xmin = np.min( np.array(x) )
    ymin = np.min( np.array(y) )

    binwidth_x = (xmax-xmin)/Dx
    binwidth_y = (ymax-ymin)/Dx
    
    lim_x = (int(xmax/binwidth_x) + 1) * binwidth_x
    lim_y = (int(ymax/binwidth_y) + 1) * binwidth_y
    



    

    bins_x = np.arange(xmin, lim_x + binwidth_x, binwidth_x)
    bins_y = np.arange(ymin, lim_y + binwidth_y, binwidth_y)
    ax_histx.hist(x, bins=bins_x)
    ax_histy.hist(y, bins=bins_y, orientation='horizontal')







left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005


rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a square Figure
fig = plt.figure(figsize=(8, 8))

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# use the previously defined function
scatter_hist(x, y, ax, ax_histx, ax_histy)

plt.show()
