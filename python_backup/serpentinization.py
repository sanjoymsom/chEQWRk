#####################################################################
# serpentinization.py
#
# Main functions. Performs serpentinization with brines
#
# Sanjoy Som, December 2022
#######################################################################
# Initialization
#######################################################################
import os
import sys
import numpy as np
import itertools
import eq36python as eq
#######################################################################
# Common Functions
#######################################################################
def extract_path_files(file):
  cwdpath, eqpath, platform, _ = [x for x in file]
  return cwdpath, eqpath, platform

def extract_wr_files(file):
  f3i, th_rock, f6i, db, _, _, _ = [x for x in file]
  return f3i, th_rock, f6i, db

def _get_cpx_harzburgite(oli,OLIwt):
  '''
  Returns a fixed CPX composition if OLI is set to a single composition
  Otherwise, returns an array of CPX compositions to loop through
  '''
  if len(OLIwt) == 1:
    CPXwt = [2.5]
    print('CPX wt% set to: '+str(CPXwt[0]))
  else: #used for ternary
    opxcpx = 100 - oli
    if opxcpx <= 100:
      CPXwt = np.linspace(0,opxcpx,len(OLIwt)+10)
    elif opxcpx > 100:
      CPXwt = np.linspace(0,100,len(OLIwt))
  return CPXwt

def _generate_equidistant_rock_compositions(num_points):
    '''
    Generates equally spaced compositional points within
    the ultramafic ternary. Thank you chatGPT for the help. 
    '''
    compositions = []
    sqrt_n = int(np.sqrt(num_points))
    for i in range(sqrt_n + 1):
        for j in range(sqrt_n + 1 - i):
            k = sqrt_n - i - j
            compositions.append([round(i*100. / sqrt_n,6), round(j*100. / sqrt_n,6), round(k*100. / sqrt_n,6)])
    #append the 3 corners (end-members) of the ternary
    end_members = [[100.,0.,0.],[0.,100.,0.],[0.,0.,100.]]
    [compositions.append(em) for em in end_members]    
    Wtp = compositions
    return Wtp

def _check_olivine(tophalf):
  if 'Fo90 Olivine' in open(tophalf).read():
    OLIVINE = 'Fo90 Olivine'
    OLI_pure_mineral = True
  elif 'IDEAL OLIVINE' in open(tophalf).read():
    OLIVINE = 'IDEAL OLIVINE'
    OLI_pure_mineral = False
  else:
    print('ERROR: Olivine not detected in tophalf file. Add your olivine')
    print('to the _check_olivine function in serpentinization.py.');sys.exit()
  return OLIVINE, OLI_pure_mineral

#######################################################################
# Option 2: react water with rock
#######################################################################
def react_heating(T_start, T_array, WR_array, OLIwt, path_files, wr_files, ternary=False):
  '''
  Reacts water and rock over a loop of Temperature, W:R, and composition
  '''
  for T, WR, oli in itertools.product(T_array, WR_array, OLIwt):
    if ternary:
      density = OLIwt[0]
      rock_compositions = _generate_equidistant_rock_compositions(density) #Weight percentage (Wtp)
      for i, Wtp in enumerate(rock_compositions):
        oli = Wtp[0]
        opx = Wtp[1]
        cpx = Wtp[2]
        #print(f"Rock Composition {i+1}: {Wtp}")
        _react_heating(T_start,T,WR,oli,opx,cpx,OLIwt,Wtp,path_files,wr_files)
    else:  
      _iterate_over_CPX(T_start,T,WR,oli,OLIwt,path_files,wr_files)

def _iterate_over_CPX(T_start,T,WR,oli,OLIwt,path_files,wr_files):
    cwdpath, eqpath, platform = extract_path_files(path_files)
    filename3i, tophalf, filename6i, serp_db = extract_wr_files(wr_files)
    #Prep minerals in harzburgite
    oli=round(oli,2)
    CPXwt = _get_cpx_harzburgite(oli,OLIwt) #identify CPX value
    opxcpx = 100-oli
    for cpx in CPXwt:
      cpx = round(cpx,2)
      opx = round(opxcpx - cpx,2)
      Wtp = [oli,opx,cpx]
      if abs(round(sum(Wtp,2)) - 100.0) > 1e-5:
        pass
      else:
        print(sum(Wtp))
        print('Sum of not equal to 100: '+str(Wtp))
        sys.exit()
      _react_heating(T_start,T,WR,oli,opx,cpx,OLIwt,Wtp,path_files,wr_files)

def _react_heating(T_start,T,WR,oli,opx,cpx,OLIwt,Wtp,path_files,wr_files):
      cwdpath, eqpath, platform = extract_path_files(path_files)
      filename3i, tophalf, filename6i, serp_db = extract_wr_files(wr_files)    
      #set minerals in tophalf file
      OLIVINE, OLI_pure_mineral = _check_olivine(tophalf)
      eq.peridotite_in6i(Wtp,cwdpath,tophalf, OLIVINE, OLI_pure_mineral)
      #Prep 3i file
      eq.edit_in3i('|Temperature',T_start,cwdpath,'Temp',filename3i)
      eq.edit_in3i('|Aq.',WR,cwdpath,'W:R',filename3i)
      #conver pitzer file to extended debye-huckel file
      eq.convert_pitz2bdot(filename3i,cwdpath,serp_db)
      #run EQ3
      conv = eq.run_eq(3,eqpath,cwdpath,filename3i,serp_db,platform)
      #add rocks and run EQ6 at temperature lower-end
      eq.set_temp_option_in6i(tophalf,cwdpath,'( 0) Constant')
      eq.edit_in6i('|'+13*' '+'Value (C',T_start,cwdpath,'Temp',tophalf)
      eq.stitch(filename3i,tophalf,filename6i)
      #and run
      conv = eq.run_eq(6,eqpath,cwdpath,filename6i,serp_db,platform)
      if conv:
        #Rename the pickup file to a 6i, and enable heating
        heat6i = 'serp_heat.6i'
        cmd = ' '.join(['mv',filename6i[:-3]+'.6p',heat6i])
        os.system(cmd)
        eq.set_temp_option_in6i(heat6i,cwdpath,'( 1) Linear')
        eq.edit_in6i('|'+13*' '+'Base Value (C',T_start,cwdpath,'Temp',heat6i)
        eq.edit_in6i('|'+13*' '+'Derivative',T-T_start,cwdpath,'Temp',heat6i)
        eq.edit_in6i('|Starting Xi',0.0,cwdpath,'Ximin',heat6i)
        eq.edit_in6i('|Maximum Xi',1.0,cwdpath,'Ximax',heat6i)
        #and run
        conv = eq.run_eq(6,eqpath,cwdpath,heat6i,serp_db,platform)
        #Move files
        if conv:
          #rename and move to output folder
          allinput = [T_start,T,WR,oli,opx,cpx]
          newserpfile = 'r6'+'-'+'seawater-'+'-'.join([str('%.2f'%x) for x in allinput])
          os.system('cp '+heat6i[:-3]+'.6o'+' '+newserpfile+'.6o')
          os.system('mv '+newserpfile+'.6o'+' plotting_files/raw6o')
      else:
        pass #keep iterating even if there is a failure.
#######################################################################
# Option 3: react brine with rock
#######################################################################
def react_brine_heating(Aw, T_start, T_array, WR_array, OLIwt, path_files, \
  wr_files):
  '''
  Reacts water and rock over a loop of Temperature, W:R, and composition.
  The brine file comes from water evaporation to a specific Aw. See the
  evaporation.py file for implementation.
  '''
  cwdpath, eqpath, platform = extract_path_files(path_files)
  filename3i, tophalf, filename6i, serp_db = extract_wr_files(wr_files)
  filename6i = 'brine.6i'
  if Aw == 'NaCl(sat)':
      cmd = ' '.join(['mv',filename6i[:-3]+'.6p',filename6i])
      os.system(cmd)
  Wtp_prev = []
  for T, WR, oli in itertools.product(T_array, WR_array, OLIwt):
    #Prep minerals in harzburgite
    oli=round(oli,2)
    CPXwt = _get_cpx_harzburgite(oli,OLIwt)
    opxcpx = 100-oli
    for cpx in CPXwt:
      cpx = round(cpx,2)
      opx = round(opxcpx - cpx,2)
      Wtp = [oli,opx,cpx]
      if abs(round(sum(Wtp,2)) - 100.0) > 1e-5:
        pass
      else:
        print('Sum of not equal to 100: '+str(Wtp))
        sys.exit()
      #skip dublicates
      #if Wtp == Wtp_prev
      #  continue
      Wtp_prev = Wtp
      #set minerals in tophalf file
      OLIVINE, OLI_pure_mineral = _check_olivine(tophalf)
      eq.peridotite_in6i(Wtp,cwdpath,tophalf,OLIVINE, OLI_pure_mineral)
      #conver pitzer file to extended debye-huckel file
      eq.convert_pitz2bdot(filename6i,cwdpath,serp_db)
      #Swap the header of the 6i file to include the rocks
      eq.swapheader_in6i(cwdpath,tophalf,filename6i)
      #run EQ6 at temperature lower-end
      eq.set_temp_option_in6i(filename6i,cwdpath,'( 0) Constant')
      eq.edit_in6i('|'+13*' '+'Value (C',T_start,cwdpath,'Temp',filename6i)
      conv = eq.run_eq(6,eqpath,cwdpath,filename6i,serp_db,platform)
      if conv:
        pass
      else:
        continue
      #Rename the pickup file to a 6i, and enable heating
      heat6i = 'serp_heat.6i'
      cmd = ' '.join(['mv',filename6i[:-3]+'.6p',heat6i])
      os.system(cmd)
      eq.set_temp_option_in6i(heat6i,cwdpath,'( 1) Linear')
      eq.edit_in6i('|'+13*' '+'Base Value (C',T_start,cwdpath,'Temp',heat6i)
      eq.edit_in6i('|'+13*' '+'Derivative',T-T_start,cwdpath,'Temp',heat6i)
      eq.edit_in6i('|Starting Xi',0.0,cwdpath,'Ximin',heat6i)
      eq.edit_in6i('|Maximum Xi',1.0,cwdpath,'Ximax',heat6i)
      #and run
      conv = eq.run_eq(6,eqpath,cwdpath,heat6i,serp_db,platform)
      #Move files
      if conv:
        #rename and move to output folder
        if Aw == 'NaCl(sat)':
          allinput = [T_start,T,WR,oli,opx,cpx]
          prefix   = 'r6-'+'NaClsat-'
        else:
          allinput = [Aw,T_start,T,WR,oli,opx,cpx]
          prefix   = 'r6-'
        newserpfile = prefix+'-'.join([str('%.2f'%x) for x in allinput])
        os.system('cp '+heat6i[:-3]+'.6o'+' '+newserpfile+'.6o')
        os.system('mv '+newserpfile+'.6o'+' plotting_files/raw6o')
