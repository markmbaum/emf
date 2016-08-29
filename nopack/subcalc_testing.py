import cProfile
import pstats
import sys

sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/")
sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/emf/")

from emf import subcalc
import pandas as pd

def get_coords(df, close):
    x = list(df['x'].values)
    y = list(df['y'].values)
    if(close):
        x.append(x[0])
        y.append(y[0])
    return({'x': x, 'y': y})

center = lambda v: float(max(v) + min(v))/2.

data,grid = subcalc.read_REF(r"G:\Projects\215045-Canal_Station_EMF\Models\SUBCALC\08-22-16_update\REF_GRID.REF")

mod = subcalc.Model(data, grid)

fig, ax, CS = subcalc.contour_plot(mod)

fn = r"G:\Projects\215045-Canal_Station_EMF\Models\SUBCALC\08-22-16_update\coordinates.csv"
df = pd.read_csv(fn, index_col = 0)
line = get_coords(df.ix[:7], False)
rail = get_coords(df.ix[8:12], False)
switch = get_coords(df.ix[13:19], True)
site = get_coords(df.ix[20:23], True)
h1 = get_coords(df.ix[24:27], True)
h2 = get_coords(df.ix[28:33], True)
h3 = get_coords(df.ix[34:37], True)
h4 = get_coords(df.ix[38:], True)

alpha = 0.8

ax.plot(line['x'], line['y'], 'k--', alpha = alpha, zorder = -1)

ax.plot(rail['x'], rail['y'], 'k', alpha = alpha, zorder = -1)
ax.text(min(rail['x']), max(rail['y']), 'Railway', ha = 'left', va = 'bottom',
        alpha = alpha, zorder = -1)

ax.plot(switch['x'], switch['y'], 'k', alpha = alpha, zorder = -1)
ax.text(center(switch['x']), center(switch['y']), 'Switchyard Boundary',
        ha = 'center', va = 'center', alpha = alpha, zorder = -1)

ax.plot(site['x'], site['y'], 'k', alpha = alpha, zorder = -1)
ax.text(center(site['x']), center(site['y']), 'Plant Boundary',
        ha = 'center', va = 'center', alpha = alpha, zorder = -1)

ax.plot(h1['x'], h1['y'], 'k', alpha = alpha, zorder = -1)
ax.plot(h2['x'], h2['y'], 'k', alpha = alpha, zorder = -1)
ax.plot(h3['x'], h3['y'], 'k', alpha = alpha, zorder = -1)
ax.plot(h4['x'], h4['y'], 'k', alpha = alpha, zorder = -1)
ax.text(1600, 300, 'Nearby Buildings',
        ha = 'center', va = 'center', alpha = alpha, zorder = -1)


mod.show()
