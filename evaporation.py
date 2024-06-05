####################################################################
# evaporation.py
#
# Evaporation initializer. Performs serpentinization with brines
#
# Sanjoy Som, December 2022
#######################################################################
# Initialization
#######################################################################
import os
import sys
import itertools
import eq36python as eq
from itertools import islice
from scipy import interpolate
#######################################################################
# Common Functions
#######################################################################
def extract_path_files(file):
  cwdpath, eqpath, platform, brine6opath = [x for x in file]
  return cwdpath, eqpath, platform, brine6opath

def extract_wr_files(file):
  f3i, _, f6i, _, th_evap, db, _ = [x for x in file]
  return f3i, th_evap, f6i, db
#######################################################################
# Evaporation
#######################################################################
def evaporation(T,WR,xi,dxi,filename3i,tophalf,platform,db,cwdpath,\
      eqpath,drop_minerals=False):
  '''
  All the steps needed to evaporate seawater using EQ3/6
  '''
  #first, set the T and WR in the 3i file
  eq.edit_in3i('|Temperature',T,cwdpath,'Temp',filename3i)
  eq.edit_in3i('|Aq.',WR,cwdpath,'W:R',filename3i)
  #Now we are ready to run EQ3
  eq.run_eq(3,eqpath,cwdpath,filename3i,db,platform)
  #set Xi in 6i file (55.51 mol is max evaporation)
  eq.edit_in6i('|Maximum Xi',xi,cwdpath,'Ximax',tophalf)
  eq.edit_in6i('|--->|dXi(n)/dXi',dxi,cwdpath,'dxixi',tophalf)
  #Now update the temperature in the rock file in LINEAR TRACKING
  eq.edit_in6i('|'+13*' '+'Value (C',T,cwdpath,'Temp',tophalf)
  #now add the output of EQ3 under the rock file
  filename6i = 'brine.6i'
  eq.stitch(filename3i,tophalf,filename6i)
  if drop_minerals:
    #drop minerals at the end of the run
    eq.edit_iopt_in_6i(filename6i,cwdpath,7,'Do it')
  #And run EQ6
  conv=eq.run_eq(6,eqpath,cwdpath,filename6i,db,platform)
  return conv
#######################################################################
# Seawater evaporation with field data
#######################################################################
def evaporation_and_nothing_else(pH,T_start,WR,path_files,wr_files,drop_minerals = False):
  '''
  Evaporate surface (pH = 8.2) seawater with and without suppression of carbonates,
  and include pH & Aw data from Klempay et al. 2021.
  '''
  print(25*'-')#------------------------------------------------------
  print('EVAPORATION')
  cwdpath, eqpath, platform, _ = extract_path_files(path_files)
  filename3i, tophalf, filename6i, evap_db = extract_wr_files(wr_files)
  xi  = 1.0
  dxi = -55.25*WR
  #update pH to surface condition
  eq.edit_in3i('|H+',pH,cwdpath,'pH',filename3i)
  conv = evaporation(T_start,WR,xi,dxi,filename3i,tophalf,platform,evap_db,cwdpath,eqpath,drop_minerals)
  if not conv:
    print('ERROR: The evaporation function in evaporation.py failed to converge.')
    sys.exit()
  print(25*'-')#-------------------------------------------------------
  #Move files
  #rename and move to output folder
  suffix = input('Enter file suffix (e.g. carb_supp): ')
  newbrinefile = 'rb-'+suffix+'-evap'
  filename6o = 'brine.6o'
  os.system('cp '+filename6o+' '+newbrinefile+'.6o')
  os.system('mv '+newbrinefile+'.6o'+' plotting_files/raw6o')
#######################################################################
# Specific brine creation
#######################################################################
def make_brine(T_start,WR,Awset,path_files,wr_files,drop_minerals = False):
  '''
  Establishes the series of steps necessary to create a brine from evaporation
  of seawater. First, seawater is completely evaporated. Then, the output is 
  examined to see which xi corresponds to the Aw of interest, following which
  seawater is evaporated again up to that new xi to create the desired brine.
  '''
  print(25*'-')#------------------------------------------------------
  print('EVAPORATION')
  cwdpath, eqpath, platform, _ = extract_path_files(path_files)
  filename3i, tophalf, filename6i, evap_db = extract_wr_files(wr_files)
  xi  = 1.0
  dxi = -55.2*WR
  conv = evaporation(T_start,WR,xi,dxi,filename3i,tophalf,platform,evap_db,cwdpath,eqpath,drop_minerals)
  if not conv:
    print('ERROR: The evaporation function in evaporation.py failed to converge.')
    sys.exit()
  print(25*'-')#-------------------------------------------------------
  print('BRINE FORMATION')
  #Finding brine with desired Aw
  wantaw = float(Awset)
  #wantaw=0.5
  newxi = _find_xi_from_aw(cwdpath,filename6i,wantaw)
  conv=evaporation(T_start,WR,newxi,dxi,filename3i,tophalf,platform,evap_db,cwdpath,eqpath)
  if not conv:
    sys.exit()
  filename6o = 'brine.6o'
  Aw=eq.findvar_in36o(18*' '+'Activity of water',cwdpath,'Aw',filename6o)
  print(25*'-')#-------------------------------------------------------
  print('RESCALING')
  SolMass=eq.findvar_in36o(23*' '+'Solvent mass',cwdpath,'SolMass'\
    ,filename6o) #grams
  Scaling_factor = 1./(SolMass*1.e-3)*WR #kg-1
  #evaporate
  evaporation(T_start,WR*Scaling_factor,newxi,dxi*Scaling_factor,filename3i,tophalf\
    ,platform,evap_db,cwdpath,eqpath,drop_minerals)
  #Get brine data for on screen display
  Aw=eq.findvar_in36o(18*' '+'Activity of water',cwdpath,'Aw',filename6o)
  SolMass=eq.findvar_in36o(23*' '+'Solvent mass',cwdpath,'SolMass'\
    ,filename6o) #grams
  print(25*'-')
  print('Brine data:')
  print(25*'-')
  print('Water activity: '+str(Aw))
  msg = ('Solvent mass:',str(SolMass),'grams')
  print(' '.join(msg))
  print('filename: brine.6i')
  #Move files
  os.system('cp brine.6p brine.6i')
  os.system('mv brine.6p'+' output_files/just_the_brine.6p') 
  os.system('cp brine.6o'+' output_files/just_the_brine.6o') 
  #rename and move to output folder
  allinput = [Aw,T_start]
  newbrinefile = 'rb-'+'-'.join([str('%.2f'%x) for x in allinput])
  os.system('cp '+filename6o+' '+newbrinefile+'.6o')
  os.system('mv '+newbrinefile+'.6o'+' plotting_files/brine6o')
  return Aw

def _find_xi_from_aw(cwdpath,filename6i,wantaw):
  #first, extract the relevant section of the eq6 output file
  #to determine array size
  dumpfile='eq6out.log'
  import numpy as np
  from scipy import interpolate
  from itertools import islice
  text = '   step=    0'
  lstart = eq.findline_in36o(text,cwdpath,dumpfile)
  text = '   Have'
  lend = eq.findline_in36o(text,cwdpath,dumpfile)
  #need to search again the 6o file and not the csv file
  #due to not enough sig figs in the csv file. aw is very
  #sensitive to xi
  #seach for Xi corresponding to wanted aw
  with open(dumpfile) as f:
    lines = islice(f,lstart,lend-1)
    for line in lines:
      if line.startswith('   step'):
        line = line.split(',')
        iteraw = float(line[3][5:])
        iterxi = float(line[1][5:])
        if iteraw < wantaw:
          Awprev = float(prevline[3][5:])
          xiprev = float(prevline[1][5:])
          xp     = [iteraw,Awprev]
          yp     = [iterxi,xiprev]
          f      = interpolate.interp1d(xp,yp)
          newxi  = f(wantaw)
          break
        prevline = line
  return newxi
#######################################################################
# Test how well the brines turned out creation
#######################################################################
def test_brines(T_start,WR,cwdpath,csv_destination,path_files,wr_files,brine6opath):
  evap_curve_filename = 'rb-evaporation-curve.6o'
  print(25*'-')#------------------------------------------------------
  print('Building Evaporation Curve')
  #check if curve exists
  cwdpath, eqpath, platform, brine6opath = extract_path_files(path_files)
  test=os.path.isfile(cwdpath+'plotting_files/brine6o/'+evap_curve_filename)
  if test:
    print('Evaporation curve (6o) exists.')
  else:
    filename3i, tophalf, filename6i, evap_db = extract_wr_files(wr_files)
    xi  = 1.0
    dxi = -55.2
    conv = evaporation(T_start,WR,xi,dxi,filename3i,tophalf,platform,evap_db,
      cwdpath,eqpath)
    if conv:
      cmd = 'mv brine.6o '+cwdpath+brine6opath+evap_curve_filename
      os.system(cmd)
    else:
      print('ERROR: The evaporation function in test_brines in evaporation.py failed to converge.')
      sys.exit()
  #check if csvs with the rb prefix exist
  test=os.path.isfile(cwdpath+csv_destination+'rb-evaporation-curve.6o.csv')
  if test:
    print('Evaporation curve (csv) exists')
  else:
    eq.generate_tablefrom6o('evaporation',cwdpath,eqpath,brine6opath,\
      csv_destination)  
  test=[filename for filename in os.listdir(cwdpath+csv_destination) if filename.startswith('rb-0')]
  if test:
    print('Brine csvs exist.')
  else:
    eq.generate_tablefrom6o('evaporation',cwdpath,eqpath,brine6opath,\
      csv_destination) 
  print(25*'-')#------------------------------------------------------
