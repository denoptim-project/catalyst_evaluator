import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.axisartist as axisartist
from matplotlib.font_manager import FontProperties

# Set common font
plt.rcParams['font.size'] = 8
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['mathtext.bf'] = 'Arial:bold'

# Data for the bar chart
objects = ['$\mathbf{H0}$', '$\mathbf{HI}$', '$\mathbf{HII\'}$', '$\mathbf{HII}$', '$\mathbf{HC1^{Ph}}$']
categories = ['Yield, 1h', 'Yield, 4h']

data = np.array([
    [1, 12],
    [5, 21],
    [24, 33],
    [29, 40],
    [80, 90]
])

# Set up figure and axis
fig = plt.figure(figsize=(3.18, 2.3))
ax = fig.add_subplot(axes_class=axisartist.Axes)

# Bar width
bar_width = 0.35

# X-axis positions for the groups
x_positions = np.arange(len(objects))

# Plot the bars for each category
for i, category in enumerate(categories):
    ax.bar(x_positions + i * bar_width, data[:, i], bar_width, label=category)

# Adjusting tick marks
ax.axis["top","right"].toggle(all=False)
ax.axis[:].major_ticks.set_tick_out(True)

# Set labels and title
ax.set_xlabel('')
ax.set_ylabel('%yield')
ax.set_title('')

# Set object labels on the x-axis and make them bold
ax.set_xticks(x_positions + bar_width / 2)
ax.set_xticklabels(objects, fontweight='bold')

# Add legend
ax.legend()
plt.yticks([0, 20, 40, 60, 80, 100])
plt.yticks(fontname='Arial')

plt.tight_layout()
plt.savefig("figures/yieldBarChart.svg", format="svg", bbox_inches="tight")

# Show the plot
plt.show()
