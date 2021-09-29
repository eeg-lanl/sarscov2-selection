## plotting functions for EstaVoir simulations
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import scipy.stats as sts
from matplotlib.transforms import blended_transform_factory

from evpytools import auxiliary as aux

def hybrid_trans(ax):
    """
    This function returns a blended/hybrid
    transformation w.r.t. subplot ax.
    x-coord: use the data for positioning
    y-coord: use relative position w.r.t y-axis
    """
    return blended_transform_factory(ax.transData, ax.transAxes)


def y_fmt(x, y):
    """
    manual formatting of semi-scientific notation
    """
    if x == 0:
        return '0'
    prefix = ''
    if x < 0:
        x = -x
        prefix = '-'
    e = int(np.log10(x))
    b = x/10**e
    return f"${prefix}{b:0.1f}\\times 10^{{{e}}}$"


def make_hybrid_axs(fig, threshold):
    gs = GridSpec(5,1)
    hx1 = fig.add_subplot(gs[0:4,:])
    lx1 = fig.add_subplot(gs[4:5,:], sharex=hx1)
    gs.update(wspace=0.025, hspace=0) # set the spacing between axes. 

    lx1.set_ylim(0, threshold)

    hx1.get_xaxis().set_visible(False)
    hx1.spines['bottom'].set_visible(False)
    hx1.set_yscale('log')
    return lx1, hx1


def make_hybrid_axs_grid(fig, threshold, shape, ratio=(4,1), wspace=0.025):
    """same as make_hybrid_axs, but then for multiple subplots"""
    top, bottom = ratio
    space = 1
    rows = top + bottom + space
    gs = GridSpec(rows * shape[0], shape[1])
    gs.update(wspace=wspace, hspace=0) # set the spacing between axes. 
    axs = [] ## to-be-returned
    for i in range(shape[0]):
        for j in range(shape[1]):
            hx = fig.add_subplot(gs[i*rows:i*rows+top,j])
            lx = fig.add_subplot(gs[i*rows+top:i*rows+top+bottom,j], sharex=hx)
            lx.set_ylim(0, threshold)
            hx.get_xaxis().set_visible(False)
            hx.spines['bottom'].set_visible(False)
            hx.set_yscale('log')
            axs.append((lx, hx))
    return axs


def plot_realization(times, states, threshold, varnames):
    fig = plt.figure(figsize=(14,5))
    lx1, hx1 = make_hybrid_axs(fig, threshold)

    for ax in (hx1, lx1):
        for varname in varnames:
            vals = np.array([state[varname] for state in states])
            if len(vals.shape) == 1:
                ## just plot the one time series and add a label
                ax.plot(times, vals, label=f"${varname}$")    
            else:
                ## only add a label to the first graph, not the others
                ax.plot(times, vals[:,0], label=f"${varname}$")
                if vals.shape[1] > 1:
                    ax.plot(times, vals[:,1:])
    lx1.legend(loc=4)    

    lx1.set_ylim(0, threshold)
    hx1.set_ylim(threshold, hx1.get_ylim()[1])
    hx1.set_yticks([10**k for k in range(int(np.log10(threshold))+1, int(np.log10(hx1.get_ylim()[1]))+1)])

    lx1.set_xlabel("time since infection (days)")
    varnamestring = ", ".join([f"${varname}$" for varname in varnames])
    fig.text(0.07, 0.5, varnamestring, rotation=90, va='center')
    return fig, (lx1, hx1)


def violin_plot(ax, xs, datas, color="violet", facecolor=None, alpha=1, midline=True, **kwargs):
    if facecolor is None: facecolor = color
    dist = np.max(xs)-np.min(xs)
    w = 0.15*max(dist,1.0) / len(xs)
    for x, d in zip(xs, datas):
        try:
            k = sts.gaussian_kde(d) #calculates the kernel density          
            m = k.dataset.min() #lower bound of violin
            M = k.dataset.max() #upper bound of violin
            y = np.arange(m, M, (M-m)/100) # support for violin
            v = k.evaluate(y) #violin profile (density curve)
            v = v/v.max()*w #scaling the violin to the available space
            if midline:
                ax.fill_betweenx(y, x, v+x, facecolor=facecolor, color=color, alpha=alpha, **kwargs)
                ax.fill_betweenx(y, x, -v+x, facecolor=facecolor, color=color, alpha=alpha, **kwargs)
            else:
                ax.fill_betweenx(y, -v+x, v+x, facecolor=facecolor, color=color, alpha=alpha, **kwargs)
        except:
            if len(d) > 1:
                ax.plot([-w+x, w+x], [np.mean(d), np.mean(d)], color=color, alpha=alpha, **kwargs)
            else:
                print("warning (violin_plot): no data in column " + str(x))

def range_plot(ax, t, l, u, color='k', dt=1, **kwargs):
    """
    a single range plot
    """
    ax.plot([t, t], [l, u], color=color, **kwargs)
    ax.plot([t-dt/2, t+dt/2], [u, u], color=color, **kwargs)
    ax.plot([t-dt/2, t+dt/2], [l, l], color=color, **kwargs)
    
    
def range_plots(ax, ts, ls, us, color='k', dt=1, **kwargs):
    """
    a series of range plots
    """
    for t, l, u in zip(ts, ls, us):
        range_plot(ax, t, l, u, color=color, dt=dt, **kwargs)
    

def pfilter_boxplot(ax, ts, ranges, pred, filt, color='k', dt=3, mask=None, **kwargs):
    """
    Range plots with predicted and filtered mean/median.
    The user can pass a masking array [bool] 'mask' 
    to prevent plotting some boxplots.
    """
    if mask is not None:
        ts = [x for x, m in zip(ts, mask) if not m]
        ranges = [x for x, m in zip(ranges, mask) if not m]
        pred = [x for x, m in zip(pred, mask) if not m]
        filt = [x for x, m in zip(filt, mask) if not m]
    range_plots(ax, ts, *aux.unzip(ranges), color=color, dt=dt, zorder=2, **kwargs)
    ax.scatter(ts, pred, color=color, marker=0, zorder=2)
    ax.scatter(ts, filt, color=color, marker=1, zorder=2)   

    
def simple_boxplot(ax, pos, data, color='k', color_med='red', p=50, horizontal=False, **kwargs):
    uppers = [np.nanpercentile(x, 50+p/2) for x in data]
    downers = [np.nanpercentile(x, 50-p/2) for x in data]
    medians = [np.nanmedian(x) for x in data]
    if horizontal:
        ax.scatter(uppers, pos, marker='|', color=color, **kwargs)
        ax.scatter(downers, pos, marker='|', color=color, **kwargs)
        ax.scatter(medians, pos, marker='|', color=color_med, **kwargs)
        for p, d, u in zip(pos, downers, uppers):
            ax.plot([d, u], [p, p], color=color, **kwargs)        
    else: ## vertical
        ax.scatter(pos, uppers, marker='_', color=color, **kwargs)
        ax.scatter(pos, downers, marker='_', color=color, **kwargs)
        ax.scatter(pos, medians, marker='_', color=color_med, **kwargs)
        for p, d, u in zip(pos, downers, uppers):
            ax.plot([p, p], [d, u], color=color, **kwargs)
