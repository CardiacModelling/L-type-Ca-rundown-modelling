import matplotlib
import matplotlib.pyplot as plt

# Set font
matplotlib.rc('font', family='arial', size = 8)

fig = plt.figure(figsize=(3.4, 2.6))
ax = fig.add_subplot(111)

ax.hlines(-90, xmin = -10, xmax = 0, color = 'tab:blue')
ax.hlines(0, xmin = 0, xmax = 120, color = 'tab:blue')
ax.hlines(-90, xmin = 120, xmax = 130, color = 'tab:blue')
ax.vlines(0, ymin = -90, ymax = 0, color = 'tab:blue')
ax.vlines(120, ymin = -90, ymax = 0, color = 'tab:blue')
ax.text(128, -95, '/', size = 20, color = 'tab:blue')
ax.text(135, -95, '/', size = 20, color = 'tab:blue')
ax.hlines(-90, xmin = 137, xmax = 147, color = 'tab:blue')

# label the times
ax.text(125, -85, r'$\mathrm{t}=\mathrm{t}_{hold}$')

ax.set_ylim(-100, 10)

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
plt.tight_layout()
plt.savefig('figures/figure3.pdf')
plt.close()