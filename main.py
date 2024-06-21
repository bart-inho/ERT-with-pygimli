from ert_pygimli_functions import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

df = pd.read_excel('/Users/bart/Documents/Indépendant/Geophy/24IRCH_Modeling.xlsx')

method = 'slm' # schlumberger

WennerSchlum = pd.read_table('WennerSchlum.txt')
exportFile(WennerSchlum, 2.5, 'Array.ohm')


i = 0
# make a loop that itarates through the columns of the dataframe

for column_name, item in df.items():

    if i == 0:
        i += 1
        continue

    if i == 3:
        break
    
    rhomap, geom = implement_geometry(df.iloc[: ,i])


    model_id = column_name

    rho_min = df.iloc[2:5, i].min()
    rho_max = df.iloc[2:5, i].max()

    scheme = pg.load('Array.ohm', verbose=True, testAll=False, realName=None)
    mesh = create_mesh(scheme, geom, rhomap, quality = 30)
    data = generate_data(mesh, scheme, rhomap)

    filename = 'data_synth/ert_' + str(model_id) + '.dat'
    save_data(data, filename)

    manager = initialize_ERTMANAGER(path=f'data_synth/ert_{model_id}.dat')

    # mesh_inv = mt.createMesh(quality=25)
    inv = invert_data(data, mesh, manager)
    meshPD = generate_meshPD(manager)

    # manager.showResultAndFit()

    import matplotlib.pyplot as plt
    import matplotlib as mpl

    # Create figure and subplots with adjusted size
    fig = plt.figure(figsize=(12, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 2], height_ratios=[1, 1])

    axs = {
        'Left': fig.add_subplot(gs[:, 0]),            # Left subplot spanning both rows
        'TopRight': fig.add_subplot(gs[0, 1]),        # Top-right subplot
        'BottomRight': fig.add_subplot(gs[1, 1])# , sharex=fig.add_subplot(gs[0, 1]))  # Bottom-right subplot sharing x-axis
    }

    # Plot data
    pg.show(data, ax=axs['Left'], cMap="jet", logScale=False, cMin=rho_min, cMax=rho_max, colorBar=False)
    axs['Left'].set_title('Synthetic Data (WS)')

    # Plot initial model
    pg.show(mesh, rhomap, ax=axs['TopRight'], hold=True, cMap="jet", logScale=False, colorBar = False, cMin=rho_min, cMax=rho_max)
    axs['TopRight'].set_xlim([0, 180])
    axs['TopRight'].set_ylim([-40, 0])
    axs['TopRight'].set_title('True model')

    # Plot inverse model
    manager.showResult(ax=axs['BottomRight'], cMap="jet", logScale=False, cMin=rho_min, cMax=rho_max, colorBar=False)
    axs['BottomRight'].set_xlim([0, 180])
    axs['BottomRight'].set_ylim([-40, 0])
    axs['BottomRight'].set_title('Inversion unstructured grid')

    # Create an axes for the colorbar
    cbar_ax = fig.add_axes([.92, 0.15, 0.02, 0.4])

    # Create the colorbar
    cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=rho_min, vmax=rho_max)
    cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')

    # Add a label to the colorbar
    cb1.set_label('Resistivity (Ohm.m)', rotation=90, labelpad=10)

    fig.suptitle('Model ID : ' + column_name + ' - WS, 72 electrodes, 2.5m spacing')

    fig.savefig('figure_results/' + column_name + '.png')
    plt.close(fig)


    # fig = plt.figure(constrained_layout=True, figsize=(11, 5.5))
    # axs = fig.subplot_mosaic([['Left', 'TopRight'],['Left', 'BottomRight']],
    #                         gridspec_kw={'height_ratios':[1, 2]}, width_ratios=[1, 2])
    
    # # Plot data
    # pg.show(data, ax=axs['Left'], cMap="jet", logScale=False, cMin=rho_min, cMax=rho_max, colorBar = False)
    # axs['Left'].set_title('Synthetic Data (WS)')

    # # Plot initial model
    # pg.show(mesh, rhomap, ax=axs['TopRight'], hold=True, cMap="jet", logScale=False, colorBar = False, cMin=rho_min, cMax=rho_max)
    # axs['TopRight'].set_xlim([0, 180])
    # axs['TopRight'].set_ylim([-40, 0])
    # axs['TopRight'].set_title('True model')

    # # Plot inverse model
    # manager.showResult(ax=axs['BottomRight'], cMap="jet", logScale=False, cMin=rho_min, cMax=rho_max, colorBar = False)
    # # pg.show(inv, ax=axs['BottomRight'], cMap="jet", logScale=False, cMin=rho_min, cMax=rho_max, colorBar = False)
    # axs['BottomRight'].set_xlim([0, 180])
    # axs['BottomRight'].set_ylim([-40, 0])
    # axs['BottomRight'].set_title('Inversion unstructured grid')

    # axs['TopRight'].sharex(axs['BottomRight'])

    # # Create an axes for the colorbar
    # cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.6])

    # # Create the colorbar
    # cmap = mpl.cm.jet
    # norm = mpl.colors.Normalize(vmin=rho_min, vmax=rho_max)
    # cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')

    # # Add a label to the colorbar
    # cb1.set_label('Resistivity (Ohm.m)', rotation=90, labelpad=10)

    # fig.suptitle('Model ID : '+column_name+ ' - WS, 72 electrodes, 2.5m spacing')

    # # plt.tight_layout()

    # fig.savefig('figure_results/'+column_name+'.png')
    # plt.close(fig)

    # fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True, figsize = ((12, 7)))

    # pg.show(mesh, rhomap, ax=ax1, hold=True, cMap="jet", logScale=False, colorBar = False, orientation=None, cMin=rho_min, cMax=rho_max)
    # ax1.set_xlim([-0, 180])
    # ax1.set_ylim([-40, 0])

    # manager.showResult(ax=ax2, cMap="jet", logScale=False, cMin=rho_min, cMax=rho_max, colorBar = False)
    # # pg.show(meshPD, inv, ax=ax2, hold=True, cMap="jet", logScale=False, colorBar = False, cMin=rho_min, cMax=rho_max)
    # ax2.set_xlim([-0, 180])
    # ax2.set_ylim([-40, 0])

    # ax1.set_title('True model')
    # ax2.set_title('Inversion unstructured grid')

    # # Create an axes for the colorbar
    # cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.4])

    # # Create the colorbar
    # cmap = mpl.cm.jet
    # norm = mpl.colors.Normalize(vmin=rho_min, vmax=rho_max)
    # cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')

    # # Add a label to the colorbar
    # cb1.set_label('Resistivity (Ohm.m)', rotation=90, labelpad=10)

    # fig.suptitle('Model ID : '+column_name+ ' - WS, 72 electrodes, 2.5m spacing')

    # # plt.tight_layout()
    # # plt.axis('equal')
    # fig.savefig('figure_results/'+column_name+'.png')
    # plt.close(fig)
    # # break

    i += 1