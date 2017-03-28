import matplotlib.pyplot as plt

from emf import fields as fld

sb = fld.load_template('practice_xcs')

n = []
keys = ['Babsmax', 'Bpermax', 'Eabsmax', 'Epermax']
eps = dict(zip(keys, [[] for i in range(len(keys))]))
for xs in sb:
    fn = 'compare_DAT/' + xs.sheet.upper() + '.DAT'
    try:
        pan = xs.compare_DAT(fn, False)
    except(fld.EMFError):
        pass
    else:
        n.append(len(xs))
        p = pan['percent-difference']
        a = pan['absolute-difference']
        eps['Babsmax'].append(a['Bmax'].abs().max())
        eps['Bpermax'].append(p['Bmax'].abs().max())
        eps['Eabsmax'].append(a['Emax'].abs().max())
        eps['Epermax'].append(p['Emax'].abs().max())

fig, axs = plt.subplots(2,2)
for i in range(2):
    for j in range(2):
        ax = axs[i][j]
        k = keys[2*i + j]
        ax.plot(n, eps[k], 'k.')
        ax.set_title(k)

plt.show()
