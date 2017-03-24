from emf import fields as fld
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 1)
for ax in axs:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.75, box.height])
ax1, ax2 = axs

for k in ['Bmax', 'Emax']:
    for i in range(1,4):

        bfn = 'Envsmpl' + str(i)
        df1 = fld.read_o01('enviro-files/' + bfn + '.o01')
        xs, df2 = fld.read_o02('enviro-files/' + bfn + '.o02')

        x = df1.index.values
        max1 = df1[k].values
        max2 = df2[k].values

        ax1.clear()
        ax1.plot(x, max2 - max1, 'k', label='Absolute Difference')
        ax1.plot(x, 100*(max2 - max1)/max1, 'r', label='Percent Difference')
        ax1.legend(loc='center left', bbox_to_anchor=(1.025, 0.5))
        ax1.grid()
        ax1.set_title(bfn)
        ax1.set_ylabel('Differences')

        ax2.clear()
        ax2.plot(x, max2, label='emf.fields results')
        ax2.plot(x, max1, label='ENVIRO results')
        yl = ax2.get_ylim()
        ax2.plot(xs.x, xs.y*0.25*(yl[1]/max(xs.y)), 'ko', label='Conductors')
        ax2.legend(loc='center left', bbox_to_anchor=(1.025, 0.5), numpoints=1)
        ax2.grid()
        ax2.set_ylabel(k)
        ax2.set_xlabel('Distance $(ft)$')

        plt.savefig('plots/' + bfn + '-' + k)
