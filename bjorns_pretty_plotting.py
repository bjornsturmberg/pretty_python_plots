"""
An amalgamation of plotting routines to produce pretty plots with matplotlib.
Original author: Bjorn Sturmberg <https://github.com/bjornsturmberg>.

This file are released under the CC0 license / public domain dedication.
We would appreciate credit if you use or redistribute this file,
but do not impose any legal restrictions.

To the extent possible under law, the persons who associated CC0 with
bjorns_pretty_plotting have waived all copyright and related or neighboring
rights to bjorns_pretty_plotting.

You should have received a copy of the CC0 legalcode along with this
work. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

"""

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
# If you want to use the new default matplotlib colourmaps, which you do!
# See https://bids.github.io/colormap/ and place files in this directory.
import colormaps as cmaps
# import csv


# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

bluey = tableau20[0]
redy = tableau20[6]
cyany = tableau20[18]
purpley = tableau20[8]
orangey = tableau20[2]
greeny = tableau20[4]
grey = tableau20[14]
pinky = tableau20[12]


def dB_function(dB):
    """ Convert dB into linear units """
    linear = 10**(dB/10.)
    return linear


def eV_function(energies):
    """ Convert energy in eV into wavelengths in nm """
    Plancks_h = 6.62606957*1e-34   # Planck's constant
    speed_c = 299792458            # Speed of light in vacuum
    charge_e = 1.602176565*1e-19   # Charge of an electron
    wls = Plancks_h*speed_c*1e9/(energies*charge_e)
    return wls


# The first thing we're plotting
nu_points = 1000
x_array = np.linspace(300, 1500, nu_points)
y_array = np.sin(x_array*np.pi/180)**2
y_array_dB = 10*np.log10(y_array)

# Or we might want to read in some data...

# # How to read in a txt file
# data = np.loadtxt('%s.txt' % data_name)
# wls = data[:, 0]
# Re_n = data[:, 1]
# Im_n = data[:, 2]

# # How to read in a csv with robustness to missing entries.
# with open('ExtractedForBjorn_rev3.csv', 'rb') as csvfile:
#     spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
#     for row in spamreader:
#         sim_wl_p.append(row[1])
#         sim_r_p.append(row[2])
# sim_wl_p = np.array([float(x) for x in sim_wl_p if any(c.isdigit() for c in x)])
# sim_r_p = np.array([float(x) for x in sim_r_p if any(c.isdigit() for c in x)])


# Layout
fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0      # Aesthetic ratio
scale = 2
fig_width = scale*fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean            # height in inches
fig_size = [fig_width, fig_height]

font_label = 14
font_ticks_major = 12
font_ticks_minor = 8
font_legend = 8
font_text = 10

linesstrength = 1.5

# Make figure
plt.clf()
fig = plt.figure(figsize=fig_size)

h_ratio = [1, 1.2]
w_ratio = [.8, 1]
gs = gridspec.GridSpec(2, 2, width_ratios=w_ratio, height_ratios=h_ratio,
                       wspace=0.7, hspace=0.6)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0, 1])
ax3 = plt.subplot(gs[1, :])

# Simple line plot
ax1.plot(x_array, y_array, '-', color=bluey, linewidth=linesstrength,
         label=r'leg label')
# Add some markers
ax1.plot(x_array[1::20], y_array[1::20], 'o', markerfacecolor='none',
         markeredgecolor=purpley, markersize=3,
         markeredgewidth=linesstrength-0.1)
# Major manipulation
major_xticks = np.arange(min(x_array), max(x_array)+1, 400)
minor_xticks = np.arange(min(x_array), max(x_array)+1, 50)
ax1.set_xticks(major_xticks)
ax1.set_xticks(minor_xticks, minor=True)
major_yticks = np.arange(min(y_array), max(y_array)+1, .5)
minor_yticks = np.arange(min(y_array), max(y_array)+1, 0.05)
ax1.set_yticks(major_yticks)
ax1.set_yticks(minor_yticks, minor=True)
ax1.set_ylabel(r"fan$_{c}$$^{y}$", fontsize=font_label)
ax1.set_xlabel("x (nm)", fontsize=font_label)
ax1.tick_params(axis='both', which='major', labelsize=font_ticks_major)
ax1.tick_params(axis='both', which='minor', labelsize=font_ticks_minor)
# Let's have another x-axis in eV
ax1_2 = ax1.twiny()
new_tick_values = np.array([4, 2, 1])
new_tick_locations = eV_function(new_tick_values)
ax1_2.set_xticks(new_tick_locations)
new_tick_labels = ["%.1f" % z for z in new_tick_values]
ax1_2.set_xticklabels(new_tick_labels, fontsize=font_ticks_major)
ax1.set_xlim((min(x_array), max(x_array)))
ax1_2.set_xlim((min(x_array), max(x_array)))
ax1_2.set_xlabel("x (eV)", fontsize=font_label)
# Let's have another y-axis in dB
ax1_3 = ax1.twinx()
new_tick_values = np.array([-5, -3, -1, 0])
new_tick_locations = dB_function(new_tick_values)
ax1_3.set_yticks(new_tick_locations)
new_tick_labels = ["%.1f" % z for z in new_tick_values]
ax1_3.set_yticklabels(new_tick_labels, fontsize=font_ticks_major)
ax1.set_ylim((min(y_array), max(y_array)+0.2))
ax1_3.set_ylim((min(y_array), max(y_array)+0.2))
ax1_3.set_ylabel(r"dB", fontsize=font_label)

handles, labels = ax1.get_legend_handles_labels()
axbox = ax1.get_position()
x_value = +.3    # Offset by eye
y_value = +.23
lgd = ax1.legend(handles, labels, loc=(axbox.x0 + x_value, axbox.y0 + y_value),
                 prop={'size':font_legend}, handlelength=2.8, frameon=False)

ax1.text(500, 1.1, '(a)', color='k', ha='center', va='center',
         fontsize=font_text)


# This time actually plot some data in dB
line = ax2.plot(x_array, y_array, '-', color=bluey, linewidth=linesstrength)
seq = [4, 2, 4, 2]
line[0].set_dashes(seq)
# Major manipulation
major_xticks = np.arange(min(x_array), max(x_array)+1, 200)
minor_xticks = np.arange(min(x_array), max(x_array)+1, 50)
ax2.set_xticks(major_xticks)
ax2.set_xticks(minor_xticks, minor=True)
major_yticks = np.arange(min(y_array), max(y_array)+1, .2)
minor_yticks = np.arange(min(y_array), max(y_array)+1, 0.05)
ax2.set_yticks(major_yticks)
ax2.set_yticks(minor_yticks, minor=True)
# ax2.set_ylabel(r"fan$_{c}$$^{y}$", fontsize=font)
ax2.set_xlabel("x (nm)", fontsize=font_label)
ax2.tick_params(axis='both', which='major', labelsize=font_ticks_major)
ax2.tick_params(axis='both', which='minor', labelsize=font_ticks_minor)
# Let's have another x-axis in eV
ax2_2 = ax2.twiny()
new_tick_values = np.array([4, 2, 1])
new_tick_locations = eV_function(new_tick_values)
ax2_2.set_xticks(new_tick_locations)
new_tick_labels = ["%.1f" % z for z in new_tick_values]
ax2_2.set_xticklabels(new_tick_labels, fontsize=font_ticks_major)
ax2.set_xlim((min(x_array), max(x_array)))
ax2_2.set_xlim((min(x_array), max(x_array)))
ax2_2.set_xlabel("x (eV)", fontsize=font_label)
ax2_3 = ax2.twinx()
line = ax2_3.plot(x_array, y_array_dB, '--', color=redy,
                  linewidth=linesstrength)
seq = [8, 6, 16, 6]
line[0].set_dashes(seq)
new_tick_values = np.array([-50, -30, -10, 0])
new_tick_locations = np.array([-50, -30, -10, 0])
ax2_3.set_yticks(new_tick_locations)
new_tick_labels = ["%.1f" % z for z in new_tick_values]
ax2_3.set_yticklabels(new_tick_labels, fontsize=font_ticks_major)
ax2.set_ylim((min(y_array), max(y_array)))
ax2_3.set_ylim((min(y_array_dB), max(y_array_dB)))
ax2_3.set_ylabel(r"dB", fontsize=font_label)
ax2_3.yaxis.label.set_color(redy)
ax2_3.yaxis.set_label_position("right")
ax2_3.tick_params(axis='y', colors=redy)
ax2.spines['right'].set_color(redy)


x_min = -10
x_max = 10
y_min = 0.0
y_max = 10
num_x = 21
num_y = 11
x_list = np.linspace(x_min, x_max, num_x)
y_list = np.linspace(y_max, y_min, num_y)
X2, Y2 = np.meshgrid(range(num_x), range(num_y))

mat = np.zeros((num_y, num_x))
for x in range(num_x):
    for y in range(num_y):
        mat[y, x] = x_list[x]**2 + y_list[y]**2


ax3 = plt.subplot(gs[1, :])
levels = np.linspace(np.min(mat), np.max(mat), 11)

cont_show = ax3.contour(X2, Y2, mat, levels, origin='lower', linewidths=2)

mat_show = ax3.matshow(mat, cmap=cmaps.inferno, aspect='equal',
                       vmin=np.min(mat), vmax=np.max(mat))
# Alternative colour scheme = viridis
ax3.xaxis.set_ticks_position('bottom')
ax3.xaxis.set_ticks_position('both')
x_ints = 5
y_ints = 3
ax3.set_xticks(np.linspace(0, num_x-1, x_ints))
ax3.set_yticks(np.linspace(0, num_y-1, y_ints))
ax3.set_xticklabels(["%3.0f" % i for i in np.linspace(x_min, x_max, x_ints)])
ax3.set_yticklabels(["%3.0f" % i for i in np.linspace(y_max, y_min, y_ints)])

# # manual_locations = [(-20, 40), (10, 50), (30, 50), (40, 50), (50, 50), (50, 30)]
# plt.clabel(cont_show, inline=True, fontsize=12, fmt='%2.0f')#, manual=manual_locations)
fmt = {}
strs = ['first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh',
        '8th', '9th', '10th', '11th']
for l, s in zip(cont_show.levels, strs):
    fmt[l] = s
plt.clabel(cont_show, cont_show.levels[::2], inline=True, fontsize=12, fmt=fmt)

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "5%", pad="3%")
cbar1 = plt.colorbar(mat_show, cax=cax, extend='neither', alpha=1,
                     ticks=levels[::2])
cbar1.solids.set_edgecolor("face")  # get rid of white lines in cbar
cbar1.set_label(r'I dunno', fontsize=font_label)
# cbar1.ax.set_yticklabels(['0', '50', '100'])
cbar1.set_clim(np.min(mat), np.max(mat))
cbar1.add_lines(cont_show)
ax3.set_xlabel(r"bottom", fontsize=font_label)
ax3.set_ylabel(r"top", fontsize=font_label)
plt.xticks(fontsize=font_ticks_major)
plt.yticks(fontsize=font_ticks_major)

# Add an inset
inset_font_label = 8
fig_ratio_x = fig_size[0]
fig_ratio_y = fig_size[1]
inset_ratio = 0.025
inset_x = inset_ratio*fig_ratio_x
inset_y = inset_ratio*fig_ratio_y
ax_inset = fig.add_axes([0.5, 0.18, inset_x, inset_y])
line = ax_inset.plot(x_array, y_array, 'k', linewidth=linesstrength*0.5)
ax_inset.set_xlabel(r"$\lambda$ (nm)", fontsize=inset_font_label, color='w')
ax_inset.set_ylabel(r"Abs (%)", fontsize=inset_font_label, color='w')
plt.xticks(fontsize=font_ticks_minor, color='w')
plt.yticks(fontsize=font_ticks_minor, color='w')
plt.xlim(520, 640)
ax_inset.set_xticks([520, 560, 600, 640])
ax_inset.set_yticks([0.0, 0.4, 0.8])


# plt.subplots_adjust(hspace = -.1)
plt.savefig('pretty', bbox_inches='tight')
plt.close()
