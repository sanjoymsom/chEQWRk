####################################################################
# plotcsv.py
#
# Plotting file. Performs serpentinization with brines
#
# Sanjoy Som, December 2022
#######################################################################
# Initialization
#######################################################################
import numpy as np
import pandas as pd
import os.path
import sys
import eq36python as eq
from itertools import islice
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
#######################################################################
# Main plotting function
#######################################################################
def setup_and_plot_brines(plot_style,cwdpath,csv_destination,
  plot_lit=False,ax1=None, lt='k',default_label=None): #lt = line type
  allfiles  = [f for f in os.listdir(csv_destination)]
  allbrines = [x for x in allfiles if 'rb-0' in x]
  evapfile  = [x for x in allfiles if 'curve' in x]
  fig,ax1 = plt.subplots()
  #plot evaporation curve
  #open file
  try:
    df = pd.read_csv(csv_destination+evapfile[0]).fillna(0)
    df = df.set_index('Variables').T
  except:
    print(evapfile)
    print('Error: Cannot open evap file. It is either open elsewhere or does not exist.');sys.exit()
  #check if var has been set
  if 'yvar' in locals(): 
    pass
  else:
    xvar,yvar = _ask_for_plotting_variables(df)
  #get carbonates
  if (xvar == 'carbonates_out') or (yvar == 'carbonates_out'):
    df = _get_carbonates(df)
  #and plot
  df_xvar,xunit,x_var,xscale = _get_plot_units(xvar,df)
  df_yvar,yunit,y_var,yscale = _get_plot_units(yvar,df)
  plt.xlabel(' '.join([x_var,xunit]))
  plt.ylabel(' '.join([y_var,yunit]))
  plt.plot(df_xvar,df_yvar,'k',lw=2)
  #Now add the brines on top of this evaporation curve
  #check if csvs exist
  if not allfiles:
    print('Error: No csvs to process. Create them first.')
    sys.exit()
  for filecsv in allbrines:
    #open file
    df = pd.read_csv(csv_destination+filecsv).fillna(0)
    df = df.set_index('Variables')
    df = df[df.columns[-1]] #select only last column of brines
    df = df.T
    #get carbonates
    if (xvar == 'carbonates_out') or (yvar == 'carbonates_out'):
      df = _get_carbonates(df)
    #get desired label
    label = 'Aw = '+filecsv.split('-')[1]
    #------ start plot styles --------------------------------------
    df_xvar,xunit,x_var,xscale = _get_plot_units(xvar,df)
    df_yvar,yunit,y_var,yscale = _get_plot_units(yvar,df)
    #and plot
    plt.xlabel(' '.join([x_var,xunit]))
    plt.ylabel(' '.join([y_var,yunit]))
    #yvar = __force_label(yvar)
    plt.plot(df_xvar,df_yvar,'x',label=label,markersize=10,markeredgewidth=3)
    ax1.set_yscale(yscale)
    ax1.set_xscale(xscale)
    if y_var == 'pH':
      ax1.set_ylim([2,12])

    ax1.set_xscale(xscale)
    #plt.plot(df[xvar],df[yvar],'x')
    #if plot_style == 0: # 1 variable vs another
    #  _plot_vars(xvar,yvar,df,ax1,label,lw=3)
  #legend
  # Shrink current axis by 20%
  box = ax1.get_position()
  ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
  # Put a legend to the right of the current axis
  ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5),frameon=False)
  plt.show()


def setup_and_plot_vars(plot_style,cwdpath,csv_destination, 
    plot_lit=False,plot_field=False,ax1=None, lt=None,
    default_label=None): #lt = line type
  allfiles = [f for f in os.listdir(csv_destination)]
  allfiles = [x for x in allfiles if 'r' in x]
  try:
    fig,ax1 = plt.subplots()
  except:
    print('Error: Cannot connect to display. X-Window is not on.')
    sys.exit()
  #check if csvs exist
  if not allfiles:
    print('Error: No csvs to process. Create them first.')
    sys.exit()
  for filecsv in allfiles:
    #open file
    df = pd.read_csv(csv_destination+filecsv).fillna(0)
    df = df.set_index('Variables').T
    #remove '_out' from mineral names
    #df.columns=[x.rstrip('_out') for x in df.columns]
    #identify which variables are of interest
    if 'yvar' in locals(): #check if var has been set
      pass
    else:
      xvar,yvar = _ask_for_plotting_variables(df)
    #get carbonates
    if (xvar == 'carbonates_out') or (yvar == 'carbonates_out'):
      df = _get_carbonates(df)
    #sort data by xvar
    df = df.sort_values(by=[xvar])
    #get desired label
    label = 'Aw = '+filecsv.split('-')[1]
    #------ start plot styles --------------------------------------
    if plot_style == 0: # 1 variable vs another
      _plot_vars(xvar,yvar,df,ax1,label,lw=3,lt=lt)
  #Add literature values if desired
  if (yvar == 'H2,aq') and (xvar == 'Temp'):
    _get_and_plot_literature_H2(cwdpath,plot_lit)
    #Get WR to set y-range for esthetics
    WR = float(filecsv.split('-')[4])
    if WR - 1. < 1e-10: #WR = 1
      pass
    elif WR - 0.1 < 1e-10: #WR = 0.1
      pass
    elif WR - 10. < 1e-10: #WR = 10
      ax1.set_ylim([1e-3,1e3])
      ax1.set_yscale('log')
    else:
      pass
  elif (yvar == 'pH') and (xvar == 'Temp'):
    _get_and_plot_literature_pH(cwdpath,plot_lit)
  elif (yvar == 'pH') and (xvar == 'Aw'):
    _get_and_plot_fielddata_pHAw(cwdpath,plot_field)
  #legend
  # Shrink current axis by 20%
  box = ax1.get_position()
  ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
  # Put a legend to the right of the current axis
  ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5),frameon=False)
  plt.show()

def mineral_superset(cwdpath):
  #Read all output files for plotting
  outpath = cwdpath + 'plotting_files/csv/'
  allfiles = [f for f in os.listdir(outpath)]
  #check if csvs exist
  if not allfiles:
    print('Error: No csvs to process. Create them first.')
    sys.exit()
  #create empty dataframe
  dfm = pd.DataFrame()
  #iterate through all the 6o files
  for csvfile in allfiles:
    df = pd.read_csv(outpath+csvfile)
    bool1 = df.iloc[:,0].str.endswith('out')
    df_minerals=df[bool1].iloc[:,0]
    #remove principal components of solid solutions
    df_minerals = df_minerals[~df_minerals.str.contains('_pcss_')]
    #and concatenate
    dfm = pd.concat([dfm,df_minerals]).drop_duplicates().reset_index(drop=True)
  #add a numerical value that will corespond to index of colormap
  dfm['colors'] = np.arange(int(len(dfm.index)))
  #remove '_out' from mineral name
  newcol = dfm[0].map(lambda x: x.rstrip("_out"))
  dfm = dfm.drop(labels=0,axis="columns")
  dfm.insert(0,"Variables", newcol)
  #organize alphabetically
  dfm = dfm.sort_values("Variables")
  #save to csv so that the next file can use the same colors
  dfm.to_csv('mineral.colors', index=True)
  os.system('mv mineral.colors plotting_files') 
  return dfm

def plot_minerals(dfm, cwdpath, visualize):
  #------- rename columns to make figure pretty -----------
  def _rename_columns(df):
    #translation dictionary
    map_names = {
      'IDEAL-OLIVINE': 'OLIVINE',
      'BRUCITE-SS'   : 'BRUCITE',
      'SERP-SS'      : 'SERPENTINE',
      'TALC-SS'      : 'TALC',
      'MOL.-MIX.-AMPH.': 'AMPHIBOLE'
    }
    #map dictionary to dataframe
    df['Phases'] = df['Variables']
    df = df.drop('Variables', axis=1)
    df = df.set_index('Phases')
    df = df.T.rename(columns=map_names).T
    return df
  #--------------------------------------------------------
  outpath = cwdpath + 'plotting_files/csv/'
  allfiles = [f for f in os.listdir(outpath)]
  for filecsv in allfiles:
    df = pd.read_csv(outpath+filecsv)
    df = df.fillna(0.)
    #filter df to extrac Temp and stack above minerals
    bool1 = df["Variables"].str.startswith('Temp')
    #filter df to keep only the minerals
    bool2 = df["Variables"].str.endswith('out')
    dfbool2 = df[bool2]
    newcol = df[bool2]["Variables"].map(lambda x: x.rstrip("_out"))
    dfbool2 = dfbool2.drop(labels="Variables",axis="columns")
    dfbool2.insert(0,"Variables", newcol)
    #rename, concat, and release
    df = pd.concat([df[bool1],dfbool2]); del bool2
    df = df.sort_values("Variables")
    #filter df using superset for minerals for this particular file
    # to assign colorcode from superset
    df_merge = pd.merge(dfm, df, on=['Variables'])
    #add back the temp row from df to the merged dataframe
    df = df_merge.append(df.iloc[-1])
    #request user input on minerals to plot
    df = _ask_mineral_index(df.reset_index(drop=True))
    #assign fixed colors
    colors = _assign_colors(df,dfm,visualize)
    #now that colors are assigned, drop the colors column
    df=df.drop(labels='colors',axis='columns')
    #filter colors and transform into colormap
    from matplotlib.colors import ListedColormap
    colors = ListedColormap(colors)
    #jiggle dataframe
    df = _rename_columns(df) #makes things pretty
    df = df.T.set_index('Temp')
    df = df.sort_values('Temp')
    #normalize
    df = df.divide(df.sum(axis=1), axis=0)
    #stack plot
    ax=df.plot.area(lw=0,colormap=colors,stacked=True,alpha=0.8)
    #niceties
    ax.set_xlabel(r'Temperature [$\circ$C]')
    ax.set_ylabel('Norm. minerals')
    ax.set_ylim([0, 1])
    ax.set_xlim([10, 400])
    #create labels
    label = 'Aw = '+filecsv.split('-')[1]
    ax.set_title(label)
    plt.show()

def plot_h2_by_mineral(cwdpath,csv_destination):
  allfiles = [f for f in os.listdir(csv_destination)]
  if not allfiles:
    print('Error: No csvs to process. Create them first.')
    sys.exit()
  #set filename to read
  filename6i = 'serp.6i'
  itexists = os.path.exists(cwdpath+'plotting_files/rock6i/'+filename6i)
  if itexists:
    print('Reading '+filename6i+' for rock olivine content.')
  else:
    print('rock file missing or misnamed. Check plot_h2_by_mineral in plotcsv.py.')
    sys.exit()
  #check if csvs exist
  if not allfiles:
    print('Error: No csvs to process. Create them first.')
    sys.exit()
  #Get Olivine moles
  minline = eq.findline_in36o('|Reactants',cwdpath+'plotting_files/rock6i/',filename6i)
  detect_olivine = False
  with open(cwdpath+'plotting_files/rock6i/'+filename6i) as f:
    lines = islice(f,minline+1,minline+8)
    for line in lines:
      if ('OLIVINE' or 'Olivine') in line:
        detect_olivine = True
      if 'remaining (moles)' in line:
        OLIVI = float(line.split(' ')[4][:-1]) # Moles
  #Ensure Olivine was in fact detected
  if not detect_olivine:
    print('Olivine was not detected. Check some things.')
    sys.exit()
  #---------------------------------------------------------------  
  #Analysis
  for filecsv in allfiles:
    #open file
    df = pd.read_csv(csv_destination+filecsv).fillna(0)
    df = df.set_index('Variables').T
    #Simplify naming
    CHRYS = df['CHRYSOTILE_pcss_out'] #Moles
    GREEN = df['GREENALITE_pcss_out'] #Moles
    CRONS = df['CRONSTEDTITE_pcss_out'] #Moles 
    MAGNE = df['MAGNETITE_out'] #Moles
    #Calculate Magnesium Number of Serpentine
    Mg_serp = 3.*CHRYS/((3.*CHRYS+GREEN)+4.*CRONS) #Mg/Mg+Fe
    #Calculate ratio of FeIII/Fetotal
    Fe3Fet = 2.*CRONS/(3.*GREEN + 4.*CRONS) #see data0.mbn for multipliers
    #Calculate H2_mineral
    H2_serpentine = df['SERP-SS_out'] * 3. * (1.-Mg_serp) * Fe3Fet * 0.5 / OLIVI #moles /Mole OLIVINE
    H2_serpentine = 1000. * H2_serpentine #Moles --> mMoles 
    H2_serpentine_mmolal = H2_serpentine * OLIVI / (df['Solv_Mass']/1000.) #mMoles --> mmolal 
    H2_magnetite  = 1000.*MAGNE / OLIVI #Moles --> mMoles / Mole OLIVINE
    H2_magnetite_mmolal  = 1000.*MAGNE / (df['Solv_Mass']/1000.) #mMoles --> mmolal
    #and plot
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    stack1 = ax1.stackplot(df.Temp,H2_serpentine,H2_magnetite)
    stack2 = ax2.plot(df.Temp,H2_serpentine_mmolal+H2_magnetite_mmolal,'k')
    #Make things pretty
    title = 'Aw = '+filecsv.split('-')[1]
    plt.title(title)
    ax1.set_xlabel(r'Temperature [$\circ$C]')
    ax1.set_ylabel(r'H$_2$ generated [mmol/mol olivine]')
    ax1.set_ylim([0, 54.5]) #limit here is chosen such that both stacks overlap
    ax2.set_ylabel(r'H$_2$ generated [mmolal]')
    ax2.set_ylim([0, 350])
    plt.show()

def plot_Fe_content_by_mineral(cwdpath,csv_destination):
  allfiles = [f for f in os.listdir(csv_destination)]
  #check if csvs exist
  if not allfiles:
    print('Error: No csvs to process. Create them first.')
    sys.exit()
  total_plots = 4
  #Analysis
  for p in range(total_plots):
    fig,ax1 = plt.subplots()
    for filecsv in allfiles:
      #open file
      df = pd.read_csv(csv_destination+filecsv).fillna(0)
      df = df.set_index('Variables').T
      #Simplify naming of Serpentine SS
      CHRYS = df['CHRYSOTILE_pcss_out']   #Moles
      GREEN = df['GREENALITE_pcss_out']   #Moles
      CRONS = df['CRONSTEDTITE_pcss_out'] #Moles
      MAGNE = df['MAGNETITE_out']         #Moles
      TSERP = df['SERP-SS_out']           #Moles
      #Simplify naming of Brucite SS
      MgBRU = df['BRUCITE_pcss_out']      #Moles
      FeBRU = df['Fe(OH)2_pcss_out']      #Moles
      TBRUC = df['BRUCITE-SS_out']        #Moles
      #Simplify naming of Olivine         #Moles
      MgOLI = df['FORSTERITE_pcss_out']   #Moles
      FeOLI = df['FAYALITE_pcss_out']     #Moles
      TOLIV = df['IDEAL-OLIVINE_out']     #Moles
      #Calculate X_Mg/Fe of Serpentine, Brucite and Olivine
      X_Fe_serp = 2.*CRONS/(2.*CRONS+3.*GREEN) #Fe3/Fe3+Fe2
      X_Mg_serp = 3.*CHRYS/((3.*CHRYS+GREEN)+4.*CRONS) #Mg/Mg+Fe
      X_Mg_bru = MgBRU/TBRUC
      X_Mg_oli = MgOLI/TOLIV
      #get desired label
      label = 'Aw = '+filecsv.split('-')[1]
      #and plot
      if p == 0:
        plt.plot(df.Temp,100.*X_Fe_serp,label=label,lw=2)
        plt.xlabel(r'Temperature [$\circ$C]')
        plt.ylabel(r'[Fe$^{3+}]$ / [Fe$^{3+}$+Fe$^{2+}$] $\times$ 100')
      elif p == 1:
        plt.plot(df.Temp,X_Mg_serp,label=label,lw=2)
        plt.xlabel(r'Temperature [$\circ$C]')
        plt.ylabel(r'X$_{Mg}$ Serpentine')
      elif p == 2:
        plt.plot(df.Temp,X_Mg_bru,label=label,lw=2)
        plt.xlabel(r'Temperature [$\circ$C]')
        plt.ylabel(r'X$_{Mg}$ Brucite')
      elif p == 3:
        plt.plot(df.Temp,X_Mg_oli,label=label,lw=2)
        plt.xlabel(r'Temperature [$\circ$C]')
        plt.ylabel(r'X$_{Mg}$ Olivine')
      else:
        print('Error: no plotting details listed.');sys.exit()
    #legend
    # Shrink current axis by 20%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5),frameon=False)
    plt.show()

#######################################################################
# Supporting function
#######################################################################
def _ask_for_plotting_variables(df,second_plot=False):
  #identify which variables are of interest
  yvar=input('Which Y-AXIS variable to plot? (or list, or carbonates_out): ')
  if (yvar == 'list') or (yvar == 'ls'):
    print(list(df.T.index))
    yvar=input('Which Y-AXIS variable to plot? ')
  #set defaults
  if not yvar: yvar = 'H2,aq'
  #check if second plot flag is passed. If not, ask for x-axis
  if not second_plot:
    xvar=input('Which X-AXIS variable to plot? ')
    if not xvar: xvar = 'Temp'
  else:
    xvar=[] #dummy
  return xvar,yvar

def _plot_vars(x_var,y_var,df,ax,label,lw=1,lt=None,default_label=None):
    import numpy as np
    import matplotlib.pyplot as plt

    #def __force_label(desired_label):
    #  if default_label is not None:
    #    desired_label = default_label
    #  else:
    #    pass
    #  return desired_label
    yvar_orig = y_var
    df_xvar,xunit,x_var,xscale = _get_plot_units(x_var,df)
    df_yvar,yunit,y_var,yscale = _get_plot_units(y_var,df)
    #and plot
    plt.xlabel(' '.join([x_var,xunit]))
    plt.ylabel(' '.join([y_var,yunit]))
    #yvar = __force_label(yvar)
    plt.plot(df_xvar,df_yvar,label=label,linewidth=lw,color=lt)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    #enter plot tayloring here
    #if yvar_orig == 'pH':
    #  plt.gca().set_ylim(3,13)

def _get_plot_units(varstr,df):
    if varstr == 'WR':
      var    = df[varstr]
      unit   = ''
      varstr = 'W:R'
      scale  = 'log'
    elif varstr == 'Temp':
      var    = df[varstr]
      unit   = r'[$\circ$C]'
      varstr = 'Temperature'
      scale  = 'linear'
    elif varstr == 'pH':
      var    = df[varstr]
      unit   = ''
      varstr = 'pH'
      scale  = 'linear'
    elif varstr == 'Aw':
      var    = df[varstr]
      unit   = ''
      varstr = 'Water Activity (Aw)'
      scale  = 'linear'
    elif varstr == 'Solv_Mass':
      var    = df[varstr]
      unit   = '[g]'
      varstr = 'Solvent mass'
      scale  = 'linear'
    elif varstr == 'TDS%':
      var    = df[varstr]
      unit   = '[%]'
      varstr = 'TDS'
      scale  = 'linear'
    elif varstr == 'TFor':
      var    = df[varstr]*1000.
      unit   = '[mmolal]'
      varstr = 'Total formate'
      scale  = 'linear'
    elif varstr == 'Xi':
      var    = df[varstr]
      unit   = ''
      varstr = r'Reaction progress $\xi$'
      scale  = 'linear'
    elif varstr.startswith('g'):
      var    = 10.**df[varstr]
      unit   = ''
      varstr = r'Activity coefficient: $\gamma$'+varstr[1:]
      scale  = 'linear'
    elif '_out' in varstr:
      try:
        var    = df[varstr]
      except:
        try:
          df[varstr] = df['Temp']
          df[varstr].iloc[:] = 0.00 
          var    = df[varstr]
        except:
          var = 0
          print('Caution: The Mineral '+varstr+' is likely suppressed or dropped. Setting to zero.')
      #mineral output is moles, so divide by solvent mass
      denom  = df['Solv_Mass']/1000. # g --> kg
      var = var/denom
      varstr = varstr[:-4].capitalize() #remove '_out'
      unit   = r'[Moles/kg. H$_2$O]'
      scale  = 'linear'
    else:
      try:
        var    = df[varstr]*1000.
      except:
        df[varstr] = df['Temp']
        df[varstr].iloc[:] = 0.00 
        var    = df[varstr]
      unit   = '[mmolal]'
      scale  = 'linear'
    return var,unit,varstr,scale

def _get_carbonates(df):
    carbonates_mbn=['ARAGONITE', 'ARTINITE','AZURITE','CALCITE','CERUSSITE',
                'DOLOMITE', 'DOLOMITE,DISORDERED','DOLOMITE,ORDERED',
                'HUNTITE','HYDROMAGNESITE','MAGNESITE','MALACHITE',
                'NESQUEHONITE','RHODOCHROSITE','SIDERITE','SMITHSONITE',
                'STRONTIANITE','WITHERITE','Carbonate-Calcite',
                'Calcite-SS']
    carbonates_ypf=['Siderite','Hydromagnesite','Artinite','Nesquehonite','Calcite',
                'Magnesite','Aragonite','Huntite','Dolomite']
    carbonates = carbonates_mbn + carbonates_ypf
    #add the suffix for consistency
    carbonates = [x+'_out' for x in carbonates]
    #create carbonates column
    df['carbonates_out'] = 0
    for solid in carbonates:
      try:
        df[solid]
        df['carbonates_out'] = df['carbonates_out'] + df[solid]
      except:
        continue
    return df

def _get_and_plot_literature_pH(cwdpath,Doit):
  if Doit:
    #open file
    mccollomfile = 'plotting_files/literature/McCollomBach2009_pH_Ol85_Opx10_Cpx5.csv'
    df2 = pd.read_csv(cwdpath+mccollomfile)
    plt.plot(df2.Temp, df2.pH,color='0', ls='--', lw=2, label='McCollom and Bach 2009')
    neutralphfile = 'plotting_files/reference/neutral_pH.csv'
    df3 = pd.read_csv(cwdpath+neutralphfile)
    plt.plot(df3.Temp, df3.pH,color='0', ls='-', lw=3, label='Neutral pH')

def _get_and_plot_literature_H2(cwdpath,Doit):
  if Doit:
    #open file
    kleinfile    = 'plotting_files/literature/klein2013_H2_Ol80_Opx10_Cpx10.csv'
    mccollomfile = 'plotting_files/literature/McCollomBach2009_H2_Ol85_Opx10_Cpx5.csv'
    df  = pd.read_csv(cwdpath+kleinfile)
    df1 = pd.read_csv(cwdpath+mccollomfile)
    plt.plot(df.Temp, df.H2,color='0', ls='-', lw=2, label='Klein et al. 2013')
    plt.plot(df1.Temp, df1.H2,color='0', ls='--', lw=2, label='McCollom and Bach 2009')

def _get_and_plot_fielddata_pHAw(cwdpath,Doit):
  if Doit:
    symbsize=10
    sbsw_s1, = plt.plot(0.697,6.97,color='#E05252',marker='>',markersize=symbsize,label='Site 1',linestyle='None')
    sbsw_s2t, = plt.plot(0.9415,8.16,color='#263F73',marker='o',markersize=symbsize,label='Site 2t',linestyle='None')
    sbsw_s2b, = plt.plot(0.9356,8.03,color='#263F73',marker='s',markersize=symbsize,label='Site 2b',linestyle='None')
    sbsw_s3, = plt.plot(0.6364,6.44,color='#E05252',marker='o',markersize=symbsize,label='Site 3',linestyle='None')
    sbsw_s4, = plt.plot(0.3923,5.3,color='#AACC66', marker='D',markersize=symbsize,label='Site 4',linestyle='None')
    sbsw_s5, = plt.plot(0.4863,5.70,color='#AACC66',marker='<',markersize=symbsize,label='Site 5',linestyle='None')
    sbsw_s6, = plt.plot(0.6819,6.83,color='#E05252',marker='s',markersize=symbsize,label='Site 6',linestyle='None')
    sbsw_s7, = plt.plot(0.7062,6.96,color='#E05252',marker='^',markersize=symbsize,label='Site 7',linestyle='None')
    sbsw_s8, = plt.plot(0.6760,6.68,color='#E05252',marker='v',markersize=symbsize,label='Site 8',linestyle='None')

def _ask_mineral_index(df):
    print(df["Variables"][:-1]) #part of the code, not a debug. leave it there.
    delvar = input('Enter indices of minerals to remove (or leave blank to include all): ')
    if not delvar: pass
    else:
      try:
        delvar = [int(x) for x in list(delvar.split(','))]
      except:
        print('Please use only indices (first column) separated by commas.')
        _ask_mineral_index(df)
      #remove selected minerals
      [df.drop(labels=i,axis=0,inplace=True) for i in delvar]
    return df

def _assign_colors(df,dfm,visualize):
  #automatic color selection (works okay, but less control on esthetics
  seed = 10 
  np.random.seed(seed)
  colorwheel = 'tab10'
  cmap = plt.get_cmap(colorwheel)
  #print('Colorwheel is '+colorwheel)
  num_colors = 20 
  values = np.linspace(0,1,num_colors)
  colors = cmap(np.random.permutation(values))
  #randomly re-arrange rows so that colors don't follow each other
  nrows = colors.shape[0]
  random_indices = np.random.permutation(nrows)
  colors = colors[random_indices]
  # manual override of color selection to control esthetics - selects 20 colors
  colors = __manual_colorwheel(visualize)
  #set the number of colors equal to the number of total minerals
  color_index = df['colors'].dropna().values.tolist()
  color_index = [int(x) for x in color_index]
  colors=colors[color_index]
  return colors

def __manual_colorwheel(visualize=False):
  colors = np.array([
    [0.83921569, 0.15294118, 0.15686275, 1.0],  #reddish
    [0.12156863, 0.46666667, 0.70588235, 1.0],  #blueish
    [0.17254902, 0.62745098, 0.17254902, 1.0],  #green
    [1.0, 0.49803922, 0.05490196, 1.0],         #orange
    [0.58039216, 0.40392157, 0.74117647, 1.0],  #purple
    [0.54901961, 0.3372549, 0.29411765, 1.0],   #brown
    [0.89019608, 0.46666667, 0.76078431, 1.0],  #pink
    [0.7372549, 0.74117647, 0.13333333, 1.0],   #yellow
    [0.09019608, 0.74509804, 0.81176471, 1.0],  #teal
    [0.0, 0.0, 0.0, 1.0],                       #Black
    [1.0, 0.0, 0.0, 1.0],                       #Red 
    [0.0, 1.0, 0.0, 1.0],                       #Green   
    [1.0, 1.0, 0.0, 1.0],                       #Yellow
    [0.0, 0.0, 1.0, 1.0],                       #Blue
    [1.0, 0.0, 1.0, 1.0],                       #Magenta
    [0.0, 1.0, 1.0, 1.0],                       #Cyan
    [0.5, 0.5, 0.5, 1.0],                       #Aqua
    [0.62745098, 0.32156863, 0.17647059, 1.0],  #Gray
    [0.32156863, 0.63921569, 0.23529412, 1.0],  #Sienna
    [0.80392157, 0.44313725, 0.3372549, 1.0]    #Tomato
  ])
  if visualize:
    # Create a color wheel plot
    angles = np.linspace(0, 2 * np.pi, len(colors), endpoint=False)
    angles = np.concatenate((angles, [angles[0]]))  # Close the loop
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(range(1, len(colors) + 1))
    ax.set_yticklabels([])

    # Plot the colors on the wheel
    for i, color in enumerate(colors):
      ax.fill_between([angles[i], angles[i + 1]], 0, 1, color=color)
  plt.show()  
  return colors
