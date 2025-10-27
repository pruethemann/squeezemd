
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

dist = pd.read_csv('S-12_distances_GLU92_ARG578.csv')

# Convert steps to time in ns
dist['time'] /= 50

print(dist.seed.unique())

print(dist['seed'].value_counts()
)

sns.lineplot(data=dist,
             x='time',
             y='distance',
             errorbar='ci',
            linewidth=0.5,   # or lw=0.5
             )

plt.ylabel('Distance in Angstrom')
plt.xlabel('Time in ns')

plt.savefig("C_terminal_binding_avg.svg")
plt.show()


sns.lineplot(data=dist,
             x='time',
             y='distance',
            palette=sns.color_palette('Blues', n_colors=8),
            linewidth=0.5,
             hue='seed')

plt.ylabel('Distance in Angstrom')
plt.xlabel('Time in ns')

plt.savefig("C_terminal_binding_all.svg", format='svg')
plt.show()

