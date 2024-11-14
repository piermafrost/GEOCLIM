import numpy as np
import netCDF4 as nc
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from itertools import cycle
import sys



##################
# GEOCLIM BASINS #
##################

# "historical" GEOCLIM order
HIST_BOX_INDEX = {'Pole': [[0,1], [7,8]], 'Mid-lat': [[2,3,4]], 'Epicont': [[5,6]]}
HIST_IND_SEDI = (1,1,0,1,1,1,1,1,1)


###################
# GEOCLIM OUTPUT: #
###################

# Get argument pass to script

output_files = sys.argv
# remove 1st argument: script name
del(output_files[0])



#############################
# BASIN GEOMETRY FUNCTIONS: #
#############################

def get_box_index(ncdat: nc.Dataset):
    '''
    Generate the dictionary of box indices from the box definition variables
    of the GEOCLIM netCDF dataset "ncdat"
    '''

    nbas = ncdat.dimensions['box'].size
    is_epic = ncdat.variables['indice_epicont'][:].astype(bool)
    is_pole = ncdat.variables['indice_polar'][:].astype(bool)
    is_surf = ncdat.variables['indice_surface'][:].astype(bool)

    box_index = { 'Pole': [], 'Mid-lat': [], 'Epicont': []}
    for j in range(nbas-1):
        if is_epic[j]:
            key = 'Epicont'
        elif is_pole[j]:
            key = 'Pole'
        else:
            key = 'Mid-lat'

        if is_surf[j]:
            box_index[key].append([j])
        else:
            box_index[key][-1].append(j)

    return box_index



######################
# DRAWING FUNCTIONS: #
######################


def plot_flux(ncdat, pdf, tindex=-1, figsize=(7,6), netflux=True):


    fig = plt.figure(figsize=figsize)
    ax = []
    # ax[0]: Carbon fluxes
    ax += [fig.add_axes([.03, .69, .65, .25])]
    # ax[1]: Oxygene fluxes
    ax += [fig.add_axes([.03, .49, .65, .12])]
    # ax[2]: Carbonate fluxes
    ax += [fig.add_axes([.03, .28, .65, .12])]
    # ax[3]: Sulfur fluxes
    ax += [fig.add_axes([.03, .07, .65, .12])]

    axtxt = fig.add_axes([.75, .02, .24, .90])
    axtxt.set_xlim([0, 1])
    axtxt.set_ylim([0, 1])
    axtxt.axis('off')

    fig.suptitle('Summary', fontweight='bold')

    def resizeheight(ax):
        h = ax.get_ylim()[1]
        ax.set_ylim([0, 1.03*h])

    def annotate_flux(ax, flxname, flxval, inout, ypos, fontsize=9, xpad=0.02):
        '''
        put name and flux value on horizontal barplot (all on a same axis).
        inout = +1 if name is wanted outside, -1 if wanted inside.
        xpad is relative to x width (i.e., =1 for entire x axis)
        both fontsize and xpad may be an array
        '''
        # prelim
        n = np.size(flxname)
        if np.size(fontsize)==1:
            fontsize = n*[fontsize]

        if np.size(xpad)==1:
            xpad = xpad*np.ones(n)

        xpad = xpad*np.diff(ax.get_xlim())

        # Annotation
        for name, val, io, y, fsz, xp in zip(flxname, flxval, inout, ypos, fontsize, xpad):
            absval = abs(val)
            if absval < 10:
                valtxt = '{:.2f}'.format(absval)
            elif absval < 100:
                valtxt = '{:.1f}'.format(absval)
            else:
                valtxt = '{:.0f}'.format(absval)

            orient = np.sign(val*io)
            if orient >= 0:
                ax.annotate(name,
                            (val+xp, y+0.17),
                            va='center', ha='left', fontsize=fsz)
                ax.annotate(valtxt,
                            (val+xp, y-0.17),
                            va='center', ha='left', fontsize=fsz)
            else:
                ax.annotate(name,
                            (val-xp, y+0.17),
                            va='center', ha='right', fontsize=fsz)
                ax.annotate(valtxt,
                            (val-xp, y-0.17),
                            va='center', ha='right', fontsize=fsz)



    # Carbon fluxes
    #--------------

    flxname = ['degass', '+sulf', 'sil wth', 'bas wth', 'ker wth', 'OC bur']
    flx = ['tot_CO2_degassing', 'sil_sulfwth_flux', 'sil_wth_flux', 'bas_wth_flux', 'ker_wth_flux', 'org_C_tot_bur_flux']
    col = ['crimson', 'goldenrod', 'steelblue', 'rebeccapurple', 'indianred', 'seagreen']
    y = [0, 0, 0, 0, 1, 1]

    flxval = [1e-12*ncdat[fname][tindex] for fname in flx]
    flxval[1] += flxval[2] # sum sulfuric + carbonic silicate weathering

    # sources: +; sinks: -
    for k in [1, 2, 3, 5]:
        flxval[k] = -flxval[k]

    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    ax[0].barh(y, flxval, color=col, edgecolor='black', linewidth=0.5)
    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    if netflux:
        net = [round(sum(flxval[:2]), 1), round(sum(flxval[4:]), 1)]
        ax[0].vlines(net[0], y[0]+0.4, y[0]-0.4, linestyles='--', colors='black', linewidths=1) 
        ax[0].vlines(net[1], y[4]+0.4, y[3]-0.4, linestyles='--', colors='black', linewidths=1) 
        annotate_flux(ax[0], flxname+['', ''], flxval+net, [-1, +1, -1, +1, -1, -1, +1, +1], y+[y[0], y[4]])
    else:
        annotate_flux(ax[0], flxname, flxval, [-1, +1, -1, +1, -1, -1], y)

    ax[0].set_xlabel('Tmol(C)/a', fontsize=10)


    # Oxygen fluxes - oxidative weathering and marine burial
    #-------------------------------------------------------

    flxname = ['+sw', 'ker wth', 'OC bur']
    flx = ['pyrite_wth_flux', 'ker_wth_flux', 'orgCbur_O2_flux']
    col = ['goldenrod', 'indianred', 'seagreen']
    y = len(col)*[0]

    flxval = [-(15/8)*1e-12*ncdat[flx[0]][tindex], -1e-12*ncdat[flx[1]][tindex], 1e-12*ncdat[flx[2]][tindex,:].sum()]
    flxval[0] += flxval[1] # sum sulfide and kerogen oxidation

    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    ax[1].barh(y, flxval, color=col, edgecolor='black', linewidth=0.5)
    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    if netflux:
        net = round(flxval[0]+flxval[2], 1)
        ax[1].vlines(net, y[0]+0.4, y[0]-0.4, linestyles='--', colors='black', linewidths=1) 
        annotate_flux(ax[1], flxname+[''], flxval+[net], [-1, -1, -1, +1], y+[y[0]])
    else:
        annotate_flux(ax[1], flxname, flxval, [-1, -1, -1], y)

    ax[1].set_xlabel('Tmol(O$_2$)/a', fontsize=10)


    # Carbonate weathering and precipitation
    #---------------------------------------

    flxname = ['+sulf', '+sil', 'carb wth', 'tot', 'carb prec (ner.)']
    flx = ['sil_sulfwth_flux', 'carb_sulfwth_flux', 'sil_wth_flux', 'carb_wth_flux', 'carb_pel_tot_bur_flux', 'carb_ner_tot_bur_flux']
    col = ['goldenrod', 'seashell', 'orchid', 'mediumblue', 'skyblue']
    y = len(col)*[0]

    flxval = [1e-12*ncdat[fname][tindex] for fname in flx]
    flxval[-2] += flxval[-1] # sum pelagic + neritic carbonate burial
    flxval[2] += flxval[3] # sum carbonate and silicate weathering
    flxval[0] += flxval[1] + flxval[2] # sum sulfuric and carbonic weathering
    del(flxval[1]) # do not show separately carbonate and silicate sulfuric weathering

    # sources: +; sinks: -
    for k in [-2, -1]:
        flxval[k] = -flxval[k]

    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    ax[2].barh(y, flxval, color=col, edgecolor='black', linewidth=0.5)
    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    if netflux:
        net = round(flxval[0]+flxval[3], 1)
        ax[2].vlines(net, y[0]+0.4, y[0]-0.4, linestyles='--', colors='black', linewidths=1) 
        annotate_flux(ax[2], flxname+[''], flxval+[net], [+1, -1, -1, -1, -1, +1], y+[y[0]],
                    xpad=[0.005, 0.015, 0.02, 0.02, 0.02, 0.02],
                    fontsize=[8, 8, 9, 9, 9, 9])
    else:
        annotate_flux(ax[2], flxname, flxval, [+1, -1, -1, -1, -1], y,
                    xpad=[0.005, 0.015, 0.02, 0.02, 0.02],
                    fontsize=[8, 8, 9, 9, 9])

    ax[2].set_xlabel('Tmol(Ca)/a', fontsize=10)


    # Sulfide weathering and sulfate-reduction
    #-----------------------------------------

    flxname = ['sulfide wth', 'sulfate-red']
    flx = ['pyrite_wth_flux', 'sulf_red_flux']
    col = ['goldenrod', 'sienna']
    y = len(col)*[0]

    flxval = [1e-12*ncdat[flx[0]][tindex], -1e-12*ncdat[flx[1]][tindex,:].sum()]

    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    ax[3].barh(y, flxval, color=col, edgecolor='black', linewidth=0.5)
    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    if netflux:
        net = round(sum(flxval), 1)
        ax[3].vlines(net, y[0]+0.4, y[0]-0.4, linestyles='--', colors='black', linewidths=1) 
        annotate_flux(ax[3], flxname+[''], flxval+[net], [-1, -1, +1], y+[y[0]])
    else:
        annotate_flux(ax[3], flxname, flxval, [-1, -1], y)

    ax[3].set_xlabel('Tmol(S)/a', fontsize=10)


    # Customization
    #--------------

    for axi in ax:
        ylim = axi.get_ylim()
        axi.vlines(0, *ylim, colors='grey', linewidths=0.5)
        axi.set_ylim(ylim)
        axi.tick_params(axis='x', labelsize=7)
        axi.get_yaxis().set_visible(False)
        for wch in ['top', 'left', 'right']:
            axi.spines[wch].set_visible(False)


    # Main fluxes
    #------------

    axtxt.annotate('Discharge:\n  {:.3f} Sv'.format((1e-6/(365.2422*24*60*60))*ncdat['discharge'][tindex]),
                   (.05, 1.0), ha='left', va='top', fontsize=12)
    axtxt.annotate('TSS flux:\n  {:.3f} Gt/yr'.format(1e-12*ncdat['TSS'][tindex]),
                   (.05, .90), ha='left', va='top', fontsize=12)
    axtxt.annotate('sedim. flux:\n  {:.3f} Gt/yr'.format(1e-12*ncdat['sedim_flux'][tindex,:].sum()),
                   (.05, .80), ha='left', va='top', fontsize=12)
    axtxt.annotate('P weath:\n  {:.1f} Gmol/yr'.format(1e-9*ncdat['P_wth_flux'][tindex]),
                   (.05, .70), ha='left', va='top', fontsize=12)
    axtxt.annotate('Net bioprod:\n  {:.1f} Tmol/yr'.format(1e-12*ncdat['org_C_bioprod_flux'][tindex,:].sum()),
                   (.05, .60), ha='left', va='top', fontsize=12)
    axtxt.annotate('Carb pel prod:\n  {:.1f} Tmol/yr'.format(1e-12*ncdat['carb_bioprod_flux'][tindex,:].sum()),
                   (.05, .50), ha='left', va='top', fontsize=12)


    # CO2, O2, and SO4
    #-----------------

    axtxt.annotate('pCO$_2$:\n  {:.1f} µatm'.format(ncdat['pCO2_atm'][tindex]),
                   (.05, .35), ha='left', va='top', fontsize=12)
    axtxt.annotate('pO$_2$:\n  {:.4f} atm'.format(ncdat['pO2_atm'][tindex]),
                   (.05, .25), ha='left', va='top', fontsize=12)
    axtxt.annotate('SO$_4$$^2$$^-$:\n  {:.3f} mol/m$^3$'.format(ncdat['SO4_glob'][tindex]),
                   (.05, .15), ha='left', va='top', fontsize=12)


    # Printing/Saving
    #----------------

    pdf.savefig()
    plt.close('all')




############################################################################
#--------------------------------------------------------------------------#
############################################################################



def plot_ocechem(ncdat, pdf, tindex=-1, figsize=(7,4)):

    colors = mpl.rcParams['axes.prop_cycle']
    markers = ['o', 's', 'd', '*']


    try:
        box_index = get_box_index(ncdat)
    except KeyError as err:
        print('could not retrieve box information:')
        print(err)
        print('assume historical geoclim box order')
        box_index = HIST_BOX_INDEX

    zmax = np.max([len(idx) for lst in box_index.values() for idx in lst]) - 1

    ########################################################
    # Oceanic profiles (polar, Mid-lat and epicontinental) #
    ########################################################

    for specie in ['alkalinity', 'DIC', ['H2CO3', 'HCO3', 'CO3'], 'pH', 'omega', 'temperature', 'Ca', 'PIC', 'POC', 'O2', 'PO4', 'POP', 'SO4', 'DIC_d13C', 'Sr', 'Sr_iso_ratio']:

        if type(specie) is not list:
            listspe = [specie]
        else:
            listspe = specie

        fig, ax = plt.subplots(ncols=3, figsize=figsize, sharex=True, sharey=True)
        fig.subplots_adjust(left=0.09, right=0.96, bottom=0.25, top=0.92, wspace=0.1)
        axtxt = fig.add_axes([.01, .01, .75, .12])
        axleg = fig.add_axes([.77, .02, .17, .17])
        for axk in [axtxt, axleg]:
            axk.set_xlim([0, 1])
            axk.set_ylim([0, 1])
            axk.axis('off')

        for k,spe in enumerate(listspe):
            units = ncdat[spe].units
            if 'd13C' in spe and units == '-':
                fact = 1e3
                units = '‰'
            elif spe in ['PO4', 'POP', 'Sr']:
                fact = 1e3
                units = 'm'+units
            else:
                fact = 1

            for axk, (key, lst) in zip(ax, box_index.items()):

                axk.set_title(key, fontsize=13)

                if len(lst) > 1 and k==0:
                    label = lambda idx: 'b. {0:}–{1:}'.format(idx[0]+1, idx[-1]+1)
                else:
                    label = lambda idx: None

                for idx, col in zip(lst, colors):
                    idpth = np.arange(0, -len(idx), -1)
                    if key == 'Pole':
                        idpth[1:] -= 1

                    axk.plot(fact*ncdat[spe][tindex,idx], idpth, '-', color=col['color'], marker=markers[k], markersize=5, label=label(idx))

                if len(lst) != 1:
                    axk.legend(fontsize=8)

        for axk in ax:
            axk.set_yticks(range(-zmax, 1, 1))
            axk.set_yticklabels([])
            axk.tick_params(axis='x', labelsize=7)

        ax[0].set_ylabel('$\longleftarrow$ depth', fontsize=11)

        if listspe[0] == 'H2CO3':
            ax[0].set_xscale('log')

        if len(listspe) > 1:
            for k,spe in enumerate(listspe):
                axleg.plot([], [], '-', color='grey', marker=markers[k], markersize=5, label=spe)

            axleg.legend(fontsize=9, bbox_to_anchor=(0, 0, 1, 1))

        name = listspe[0]
        for spe in listspe[1:]:
            name += ', '+spe

        axtxt.annotate(name+' ('+units+')', (0.1, 0.5), fontsize=12, fontweight='bold', ha='left', va='bottom')


        # Mean value:
        #------------

        if type(specie) is str:

            if specie == 'alkalinity':
                globname = 'alk_glob'
            elif specie == 'Sr_iso_ratio':
                globname = 'Sr_iso_rat_glob'
            else:
                globname = specie+'_glob'

            try:
                meanval = fact*ncdat[globname][tindex]
                axtxt.annotate('Mean value: {:.4G} '.format(meanval)+units, (0.99, 0.5), ha='right', va='bottom')
            except IndexError:
                pass


        # Save page
        #----------

        pdf.savefig()
        plt.close('all')



    ###################################
    # Lysocline depths (special case) #
    ###################################

    listspe = ['diag_lys_dpth_calc', 'diag_lys_dpth_arag']
    listnam = ['calcite', 'aragonite']

    fig, ax = plt.subplots(ncols=3, figsize=figsize, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.09, right=0.96, bottom=0.25, top=0.92, wspace=0.1)
    ax[0].invert_yaxis()
    axtxt = fig.add_axes([.01, .01, .75, .12])
    axleg = fig.add_axes([.76, .03, .15, .15])
    for axk in [axtxt, axleg]:
        axk.set_xlim([0, 1])
        axk.set_ylim([0, 1])
        axk.axis('off')

    n = len(listspe)
    for (k,spe),ls in zip(enumerate(listspe), ['-', ':']):
        units = ncdat[spe].units

        for axk, (key, lst) in zip(ax, box_index.items()):

            axk.set_title(key, fontsize=12)

            if len(lst) > 1 and k==0:
                label = lambda idx: 'b. {0:}–{1:}'.format(idx[0]+1, idx[-1]+1)
            else:
                label = lambda idx: None

            for idx, col in zip(lst, colors):
                axk.hlines(ncdat[spe][tindex,idx[0]], 0, 1, linestyles=ls, colors=col['color'], label=label(idx))
                axk.annotate('{:.2f}'.format(ncdat[spe][tindex,idx[0]]), (0.05+k/n,ncdat[spe][tindex,idx[0]]), va='bottom')

            if len(lst) != 1:
                axk.legend(fontsize=8)

    ax[0].set_xlim([0, 1])
    for axk in ax:
        axk.set_xticks([])
        axk.set_xticklabels([])

    for axk in ax[1:]:
        axk.tick_params(axis='y', labelright=False, labelleft=False)

    ax[0].set_ylabel('depth ({:})'.format(units))

    for spe, ls in zip(listnam, ['-', ':']):
        axleg.hlines([], 0, 1, linestyles=ls, colors='grey', label=spe)

    axleg.legend(fontsize=9, bbox_to_anchor=(0, 0, 1, 1))

    axtxt.annotate('Lysocline depth ('+units+')', (0.1, 0.5), fontsize=12, fontweight='bold', ha='left', va='bottom')


    # Save page
    #----------

    pdf.savefig()
    plt.close('all')



    # === #



    # --------------------------------------------------------------------- #
    # Ordering boxes, by depth, and following "epicont -> mid-lat -> poles" #
    # --------------------------------------------------------------------- #

    cmap = mpl.colormaps['Set2']
    get_cmap_val = lambda n: [np.array(cmap(x)[:-1]) for x in np.arange(min(n,8))/7]

    order = []
    xpos  = []
    label = []
    color = []

    x = 0

    for key, lab in zip(['Epicont', 'Mid-lat', 'Pole'], ['epic.', 'm-lat', 'pole']):
        ncol = len(box_index[key])
        categories_colors = get_cmap_val(ncol)
        zmax = np.max([len(idx) for idx in box_index[key]])
        k = -1
        go_on = True
        while go_on:
            colorcycle = cycle(categories_colors)
            x += 0.5
            k += 1
            go_on = False
            for idx in box_index[key]:
                if len(idx) > k:
                    go_on = True
                    order.append(idx[k])
                    label.append(lab+' b#{:}'.format(1+idx[k]))
                    color.append((1-k/zmax)*next(colorcycle))
                    xpos.append(x)
                    x += 1

    ###############################################
    # Prod. fluxes + Bur. flux & Sedim. variables #
    ###############################################

    def_color_list = ['maroon', 'goldenrod', 'slateblue']
    #
    #color = ['goldenrod', 'saddlebrown', 'skyblue', 'royalblue', 'navy', 'turquoise', 'teal', 'plum', 'darkorchid']

    for var in ['org_C_bioprod_flux', 'carb_bioprod_flux', 'sedim_rate', 'sedim_flux', 'carb_bur',
                'org_C_bur_flux', 'burial_efficiency', 'P_bur', 'sulf_red_flux']:

        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([.09, .28, .88, .65])
        axtxt = fig.add_axes([.01, 0, .98, .11])
        axtxt.set_xlim([0, 1])
        axtxt.set_ylim([0, 1])
        axtxt.axis('off')

        color_list = def_color_list

        if var == 'P_bur':
            value = [1e-9*f[v][tindex,order] for v in ['P_bur_flux_orgC', 'P_bur_flux_phosph', 'P_bur_flux_hydro']]
            value[1] += value[2]
            value[0] += value[1]
            unit = 'Gmol/a'
            name = 'Phosphorus burial fluxes'
            legend = ['orgC-bound', 'phosphorite', 'Hydro Fe-bound']
        elif var == 'sulf_red_flux':
            value = 1e-12 * f[var][tindex,order]
            unit = 'Tmol/a'
            name = f[var].long_name
        elif var == 'sedim_rate':
            value = 1e5 * f[var][tindex,order]
            unit = 'cm/ka'
            name = f[var].long_name
        elif var in ['sedim_flux', 'org_C_bur_flux']:
            value = 100 * f[var][tindex,order] / f[var][tindex,:].sum()
            unit = '% of total'
            name = f[var].long_name
        elif var == 'carb_bur':
            value = [f[v][tindex,order] for v in ['carb_pel_bur_flux', 'carb_ner_bur_flux']]
            value[0] += value[1]
            value = [100*v/value[0].sum() for v in value]
            unit = '% of total'
            name = 'Carbonate burial fluxes'
            legend = ['pelagic', 'neritic']
            color_list = ['mediumblue', 'skyblue']
        else:
            value = f[var][tindex,order]
            unit = f[var].units
            name = f[var].long_name


        # Bar plot
        #---------

        if type(value) is list:
            for val,col in zip(value, color_list):
                ax.bar(xpos, val, width=0.8, color=col, edgecolor='black', linewidth=0.5)

            ax.legend(legend)
        else:
            ax.bar(xpos, value, width=0.8, color=color, edgecolor='black', linewidth=0.5)

        ax.set_xlim([xpos[0]-0.75, xpos[-1]+0.75])
        ax.set_xticks(xpos)
        ax.set_xticklabels(label, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel(unit)
        ax.set_title('depth $\longrightarrow$')

        axtxt.annotate(name, (0.1, 0.5), fontsize=12, fontweight='bold', ha='left', va='bottom')


        # Save page
        #----------

        pdf.savefig()
        plt.close('all')


############################################################################


for name in output_files:

    f = nc.Dataset(name)

    # remove "geoclim_output."
    k = name.find('geoclim_output')
    name = name[k+15:]
    # remove ".nc"
    name = name[:-3]

    #it = int(np.argwhere(f['time'][:]==700000))
    it = -1

    # figures' file name and pdf object
    pdf = PdfPages('final_ocean_chemistry--'+name+'.pdf')

    # Fluxes histograms
    plot_flux(f, pdf, tindex=it)

    # Ocean chemistry
    plot_ocechem(f, pdf, tindex=it)

    # Close figure file and dataset file
    pdf.close()
    f.close()

