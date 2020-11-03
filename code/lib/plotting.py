import string, re, math, itertools, glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import palettable
import seaborn as sns
import scipy.stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import statsmodels.api as sm

from .config import *


# save and close plotted stuff
def plot_save(datadir, fidx, tag):
    plt.savefig("%sout%d-%s.pdf" % (datadir, fidx, tag))
    print("Printed %sout%d-%s.pdf" % (datadir, fidx, tag))
    plt.close()

# plot zipf
def plot_zipf(sizes, normed=True, ax=None,
              scalex=1.0, scaley=1.0, **kwargs):
    if ax is None:
        ax = plt.gca()

    size_unique, size_count = np.unique(sizes, return_counts=True)
    if normed:
        size_count = size_count.astype(float)
        size_count /= np.sum(size_count)
    ax.plot(scalex*size_unique, scaley*size_count, **kwargs)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel("Clone Size Frequency")
    ax.set_xlabel("Clone Size")

def label_axes(fig_or_axes, labels=string.ascii_uppercase,
               labelstyle=r'%s',
               xy=(-0.1, 0.95), xycoords='axes fraction', **kwargs):
    """
    Walks through axes and labels each.
    kwargs are collected and passed to `annotate`

    Parameters
    ----------
    fig : Figure or Axes to work on
    labels : iterable or None
        iterable of strings to use to label the axes.
        If None, lower case letters are used.

    loc : Where to put the label units (len=2 tuple of floats)
    xycoords : loc relative to axes, figure, etc.
    kwargs : to be passed to annotate
    """
    # re-use labels rather than stop labeling
    defkwargs = dict(fontweight='bold')
    defkwargs.update(kwargs)
    labels = itertools.cycle(labels)
    axes = fig_or_axes.axes if isinstance(fig_or_axes, plt.Figure) else fig_or_axes
    for ax, label in zip(axes, labels):
        ax.annotate(labelstyle % label, xy=xy, xycoords=xycoords,
                    **defkwargs)

# plot CDF
def plot_rankfrequency(data, ax=None,
                       normalize_x=True, normalize_y=False,
                       log_x=True, log_y=True,
                       scalex=1.0, scaley=1.0, **kwargs):
    """
    Plot rank frequency plots. 

    data: count data to be plotted
    ax: matplotlib Axes instance
    normalize_x: if True (default) plot relative frequency, if False plot raw counts
    normalize_y: if False (default) plot rank, if True plot cumulative probability
    """
    if ax is None:
        ax = plt.gca()
    data = np.asarray(data)
    data = data[~np.isnan(data)]
    if normalize_x:
        data = data/np.sum(data)
    sorted_data = np.sort(data)  # Or data.sort(), if data can be modified
    # Cumulative counts:
    if normalize_y:
        norm = sorted_data.size
    else:
        norm = 1
    ret = ax.step(sorted_data[::-1]*scalex, scaley*np.arange(sorted_data.size)/norm, **kwargs)
    if log_x:
        ax.set_xscale('log')
    if log_y:
        ax.set_yscale('log')
    if normalize_x:
        ax.set_xlabel('Normalized clone size')
    else:
        ax.set_xlabel('Clone size')
    if not normalize_y:
        ax.set_ylabel('Clone size rank')
    return ret

def plot_insetcolorbar(vmin, vmax, cmap, step=0.1, label=None, ax=None):
    """
    Plot an inset colorbar based on a dummy axes
    """
    if ax is None:
        ax = plt.gca()
    fig, dummyax = plt.subplots()
    # make dummy plot for colorbar
    levels = np.arange(vmin, vmax+step, step)
    CS = dummyax.contourf([[0,0],[1,0]], levels, cmap=cmap)
    plt.close(fig)
    cax = inset_axes(ax, width="30%", height="3%", loc='upper right')
    cbar = plt.colorbar(CS, orientation='horizontal', cax=cax, ticks=[vmin, vmax])
    if label:
        cbar.set_label(label)
    return cax, cbar

def plot_referencescaling(ax=None, x=[4e-5, 4e-2], factor=1.0, color='k', exponent=-1.0, label=True, **kwargs):
    """
    Plot a reference power law scaling with slope -1.

    kwargs are passed to ax.plot
    """
    if ax is None:
        ax = plt.gca()
    x = np.asarray(x)
    ax.plot(x, factor*x**exponent, color=color, **kwargs)
    if label:
        xt = scipy.stats.gmean(x)
        xt = xt*1.05
        yt = factor*xt**exponent *1.05
        ax.text(xt, yt, '%g'%exponent, va='bottom', ha='left', color=color)

def statsmodels_regression(x, y):
    x = sm.add_constant(x)
    model = sm.OLS(y,x)
    results = model.fit()
    return model, results

def plot_regression(x, y, ax=None,
                    logy=False, p_cutoff=0.05, fit_slope=True,
                    extend=0.0, ci=95, plot_ci=True,
                    fittype='bootstrap',
                    fittransform=None,
                    data_label='',
                    label=None,
                    **kwargs):
    """Plot a linear regression analysis.
   
    logy: log-transform data before fitting
    p_cutoff: significance cutoff for showing the fitted slope
    fit_slope: fit slope if True else rate
    fittype : in bootstrap, scipy, statsmodels
    extend: by how much to extend the fitting function beyond the fitted values
    """
    if fittype not in ['bootstrap', 'scipy', 'statsmodels']:
        raise Exception('Invalid argument')
    if label is None:
        if fittype == 'bootstrap':
            label = '{0:.0f} [{1:.0f}, {2:.0f}]'
        elif fittype == 'scipy':
            label = '${0:.0f} \pm {1:.0f}$'
        elif fittype == 'statsmodels':
            label = '${0:.0f} \pm {2:.0f} x + {1:.0f} \pm {3:.0f}$'
    if ax is None:
        ax = plt.gca()
    l, = ax.plot(x, y, 'o', label=data_label, **kwargs)
    if logy:
        y = np.log(y)

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    if fittype == 'bootstrap':
        if fit_slope:
            def robust_linregress(x, y):
                try:
                    res = scipy.stats.linregress(x, y)
                    return res
                except FloatingPointError:
                    return [np.nan]*5
            fit_parameter = 1/slope
            bootstrap = sns.algorithms.bootstrap(x, y, func=lambda x, y: 1/robust_linregress(x, y)[0], n_boot=10000)
        else:
            fit_parameter = slope
            bootstrap = sns.algorithms.bootstrap(x, y, func=lambda x, y: scipy.stats.linregress(x, y)[0], n_boot=10000)
        low, high = sns.utils.ci(bootstrap)
        print('fit', fit_parameter, 'std', np.std(bootstrap),
              'p', p_value, 'low', low, 'high', high)
        label = label.format(fit_parameter, low, high)
    elif fittype == 'scipy':
        if fit_slope:
            label = label.format(1/slope, std_err/slope**2, r_value**2)
        else:
            label = label.format(slope, std_err, r_value**2)
    elif fittype == 'statsmodels':
        x_fit = x.copy()
        if not fittransform is None:
            x_fit = fittransform(x_fit)
        model, results = statsmodels_regression(x_fit, y)
        label = label.format(results.params[1], results.params[0], results.bse[1], results.bse[0])

    print(label)

    x_fit = np.linspace(min(x)-extend, max(x)+extend, 400)
    y_fit = intercept+slope*x_fit
    if logy:
        y_fit = np.exp(y_fit)
        ax.set_yscale('log')
    ax.plot(x_fit, y_fit, c=l.get_color(),
            label=label if p_value<p_cutoff else 'NS', **kwargs)

    # error band for plot
    if plot_ci:
        def reg_func(_x, _y):
            return np.linalg.pinv(_x).dot(_y)
        X = np.c_[np.ones(len(x)), x]
        grid = np.c_[np.ones(len(x_fit)), x_fit]
        yhat = grid.dot(reg_func(X, y))
        beta_boots = sns.algorithms.bootstrap(X, y, func=reg_func,
                                    n_boot=10000).T
        yhat_boots = grid.dot(beta_boots).T
        err_bands = sns.utils.ci(yhat_boots, ci, axis=0)
        if logy:
            err_bands = np.exp(err_bands)
        ax.fill_between(x_fit, *err_bands, facecolor=l.get_color(), alpha=.3)

    return slope, intercept, r_value**2

def _split(number):
    """ Split a number in python scientific notation in its parts.
        
        @return value and exponent of number

    """
    return re.search(r'(-?[0-9].[0-9]*)(?:e\+?)(-?[0-9]*)', number).groups()

def str_quant(u, uerr, scientific=False):
    """ Make string representation in nice readable format
    
        >>> str_quant(0.0235, 0.0042, scientific = True)
        '2.4(5) \\\cdot 10^{-2}'
        >>> str_quant(1.3, 0.4)
        '1.3(4)'
        >>> str_quant(8.4, 2.3)
        '8(3)'
        >>> str_quant(-2, 0.03)
        '-2.00(3)'
	>>> str_quant(1432, 95, scientific = True)
	'1.43(10) \\\cdot 10^{3}'
	>>> str_quant(1402, 95, scientific = True)
	'1.40(10) \\\cdot 10^{3}'
        >>> str_quant(6.54, 0.14)
        '6.54(14)'
        >>> str_quant(0.8, 0.2, scientific=False)
        '0.8(2)'
        >>> str_quant(45.00, 0.05, scientific=False)
        '45.00(5)'

    """
    # preformatting
    number = format(float(u), "e")
    error = format(float(uerr), "e")
    numberValue, numberExponent = _split(number) 
    errorValue, errorExponent = _split(error)
    numberExponent, errorExponent = int(numberExponent), int(errorExponent)    

    # Precision = number of significant digits
    precision = numberExponent - errorExponent
    # make error
    if errorValue.startswith("1"):
        precision += 1
        errorValue = float(errorValue) * 10  # roundup second digit
    error = int(math.ceil(float(errorValue))) # roundup first digit

    # number digits after point (if not scientific)
    nDigitsAfterPoint = precision - numberExponent
    # make number string
    if scientific:
        number = round(float(numberValue), precision)
        if precision == 0:
            number = int(number)
    else:
        number = round(float(numberValue) * 10**numberExponent, nDigitsAfterPoint)
        if nDigitsAfterPoint == 0:
            number = int(number)
    numberString = str(number)

    # pad with 0s on right if not long enough
    if "." in numberString:
        if scientific:
            length = numberString.index(".") + precision + 1
            numberString = numberString.ljust(length, "0")
        else:
            length = numberString.index(".") + nDigitsAfterPoint + 1
            numberString = numberString.ljust(length, "0")
    
    if scientific and numberExponent != 0:
        outputString = "%s(%d) \cdot 10^{%d}" % (numberString, error, numberExponent)
    else:
        outputString = "%s(%d)" % (numberString, error)

    return outputString

from scipy.interpolate import interpn

def density_scatter(x, y, ax=None, sort=True, bins=20, trans=None, **kwargs):
    """
    Scatter plot colored by 2d histogram
    """
    x = np.asarray(x)
    y = np.asarray(y)
    if ax is None :
        ax = plt.gca()
    if trans is None:
        trans = lambda x: x
    data , x_e, y_e = np.histogram2d(trans(x), trans(y), bins=bins)
    z = interpn(( 0.5*(x_e[1:] + x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1]) ),
                data, np.vstack([trans(x),trans(y)]).T,
                method="splinef2d", bounds_error=False)

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter(x, y, c=z, **kwargs)
    return ax

def plot_zeroinsertion_aging(df_enrichments, name,
        minrank=1, maxrank_fitting=9, maxrank_plotting=9, agebinsize=10.0,
        misspecification_error=2e-3, alpha=1.17):
    """
    Plotting zeroinsertions as a function of aging.
    """
    agebins = np.arange(0.0, 81.0, agebinsize)
    bin_ts = agebins[:-1]+agebinsize/2
    bins = np.array([1, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000])
    bins = bins[:max(maxrank_fitting, maxrank_plotting)+1]
    binmids = 0.5*(bins[1:] + bins[:-1])
    #binmids[1] = 250
    print(binmids[minrank:maxrank_fitting])

    #df_enrichments['zeroInsertion500'] = (3*df_enrichments['zeroInsertion500'] + 2*df_enrichments['zeroInsertion200'])/5.0
    grouped = df_enrichments.groupby(pd.cut(df_enrichments['Age'], bins=agebins))

    meanfreq = grouped.agg('mean')
    meanfreq = np.array([list(meanfreq['zeroInsertion%s'%rank]) for rank in bins[1:]])
    semfreq = grouped.agg('sem')
    semfreq = np.array([list(semfreq['zeroInsertion%s'%rank]) for rank in bins[1:]])

    def prediction(rank, t, alpha, sigma, rank0, background_probability=0.0, max_probability=1.0):
        return background_probability+(max_probability-background_probability)*0.5*scipy.special.erfc((np.log(rank/rank0)/alpha+t*alpha*sigma**2)/(4*t*sigma**2)**.5)

    def func(params, alpha):
        taud, rank0, background_probability, max_probability = params
        sigma = (1/taud)**.5
        residuals = (meanfreq[minrank:maxrank_fitting]-prediction(binmids[minrank:maxrank_fitting][:, np.newaxis],
                                        bin_ts, alpha, sigma, rank0,
                                        background_probability, max_probability))/(semfreq[minrank:maxrank_fitting] + misspecification_error)
        return residuals.flatten()

    background_probability = meanfreq.min()
    max_probability = meanfreq.max()
    sigmasq = 0.05
    sigma = sigmasq**.5
    rank0 = 1.5e4
    #params = sigmasq, rank0
    #optparams, pcov, infodict = scipy.optimize.leastsq(lambda params, back_prob, max_prob, *arg: func([params[0],
    #                                                                        params[1], back_prob, max_prob], *arg),
    #                                                   params, args=(background_probability, max_probability, alpha,),
    #                                                   full_output=True)[:3]
    params = 1/sigmasq, rank0, background_probability, max_probability
    optparams, pcov, infodict = scipy.optimize.leastsq(func,
                                                       params, args=(alpha,),
                                                       full_output=True)[:3]
    s_sq = np.sum(infodict['fvec']**2) / (infodict['fvec'].size - optparams.size)
    print(s_sq)
    #pcov = pcov * s_sq
    optse = pcov.diagonal()**.5
    for p, pse in zip(optparams, optse):
        print(str_quant(p, pse, scientific=True))

    #sigmasq, rank0 = optparams
    tau_d, rank0, background_probability, max_probability = optparams
    sigma = (1/tau_d)**.5

    df_enrichments = df_enrichments[~df_enrichments['Age'].isna()]
    df_enrichments = df_enrichments.reset_index()

    fig, axes = plt.subplots(figsize=(4.42, 2.4), ncols=2, sharey=True)
    colors = np.asarray(palettable.matplotlib.Viridis_8.mpl_colors)
    marker = ["o", "v", "^", "<", ">", "1", "2", "3", "4", "x"]

    ax = axes[0]
    for i, t in enumerate(bin_ts):
        #x = binmids[minrank:maxrank_plotting]
        x = bins[minrank+1:maxrank_plotting+1]
        y = meanfreq[minrank:maxrank_plotting, i]
        yerr = semfreq[minrank:maxrank_plotting, i]
        ax.plot(x, y, '-', c=colors[i])
        ax.fill_between(x,
                    y-yerr, y+yerr, facecolor=colors[i], alpha=.5, edgecolor=None)
    for i in range(minrank, maxrank_plotting):
        #rank = binmids[i]
        rank = bins[i+1]
        for j, t in enumerate(bin_ts):
            x = rank
            y = meanfreq[i, j]
            ax.plot(x, y, marker[i-minrank],
                    c=colors[j], markersize=3,
                    label='%g-%g'%(t-5, t+5) if i == minrank else '')
     

    ax.set_xlabel('Clone size rank (binned)')
    ax.set_ylabel('Zero insertion clones')
    legend_kwargs = dict(ncol=2, fontsize='x-small', handlelength=1, title_fontsize='x-small',
              loc='upper right', bbox_to_anchor=(1.0, 1.1))
    ax.legend(title='Age in years (binned)', **legend_kwargs)
    ax.set_ylim(0.0, 0.09)
    ax.set_xscale('log')

    ax = axes[1]
    for i in range(minrank, maxrank_plotting):
        rank = binmids[i]
        for j in range(len(agebins)-1):
            x = (np.log(rank/rank0)/alpha+bin_ts[j]*alpha*sigma**2)/((4*sigma**2 * bin_ts[j])**.5)
            y = meanfreq[i, j]
            if j == minrank:
                ax.plot(x, y, marker[i-minrank],
                    c=colors[j],# label='%g'%bins[i+1],
                    zorder=-1)
            ax.errorbar(x, y, semfreq[i, j], fmt=marker[i-minrank],
                    c=colors[j], zorder=3)
    x = np.linspace(-3.2, 2.4)
    ax.plot(x, (max_probability-background_probability)*0.5*scipy.special.erfc(x)+background_probability,
            label='Theory',
            lw=2, c='k', zorder=2)
    ax.set_xlabel(#'Rescaled rank and age\n'+
                  r'$\left[\log(r/r^\star)+t/\tau_d\right] \, / \sqrt{4 t/\tau_d}$')
    #legend_kwargs.update(dict(loc='lower left', bbox_to_anchor=None))
    ax.locator_params(axis='x', nbins=10)
    #ax.legend(title='Clone size rank (binned)', **legend_kwargs)
    ax.legend(loc='upper right')
    ax.set_xlim(min(x), max(x))
    ax.set_ylim(0.0, 0.085)
    fig.tight_layout()

    fig.savefig(figure_directory+'%s.svg'%name)
    return fig

class OffsetHandlerTuple(mpl.legend_handler.HandlerTuple):
    """
    Legend Handler for tuple plotting markers on top of each other
    """
    def __init__(self, **kwargs):
        mpl.legend_handler.HandlerTuple.__init__(self, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):
        nhandles = len(orig_handle)
        perside = (nhandles - 1) / 2
        offset = height / nhandles
        handler_map = legend.get_legend_handler_map()
        a_list = []
        for i, handle1 in enumerate(orig_handle):
            handler = legend.get_legend_handler(handler_map, handle1)
            _a_list = handler.create_artists(legend, handle1,
                                             xdescent,
                                             offset*i+ydescent-offset*perside,
                                             width, height,
                                             fontsize,
                                             trans)
            a_list.extend(_a_list)
        return a_list


