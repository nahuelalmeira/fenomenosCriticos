import matplotlib as mpl

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
          '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
markers = ['s', 'o', 'v', '^', '*', '<', '>', '8', '.', '',
		   's', 'o', 'v', '^', '*', '<', '>', '8', '.', '']

rc_font_size = 30
rc_label_size = 28
rc_legend_size = 26
mpl.rcParams['figure.figsize'] = (10, 6)
#mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['legend.fontsize'] = rc_legend_size
mpl.rcParams['savefig.transparent'] = True
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['axes.labelsize'] = rc_font_size
mpl.rcParams['xtick.labelsize'] = rc_label_size
mpl.rcParams['ytick.labelsize'] = rc_label_size
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['patch.linewidth'] = 1.5