#####################################################################
# eq36python.py
#
# EQ3/6 interface. Performs serpentinization with brines
#
# Sanjoy Som, December 2022
#######################################################################
# Initialization
#######################################################################
import os
import sys
import numpy as np
#######################################################################
# Main function (called from the outside)
#######################################################################
def peridotite_in6i(RC,path,filename,OLIVINE,OLI_pure_mineral):
  RC = [float(x) for x in RC] #fix type here. #RC=Rock Composition
  #Check if RC is a list
  if type(RC) is not list:
    print('Error: First element passed needs to be a list. Quitting')
    sys.exit()
  #check rock composition is 100%. % is weight percent
  if abs(sum([x for x in RC]) - 100.0) > 1e-5:
    print('Error: rock composition does not add to 100%. Quitting.')
    sys.exit()
  #Check state of Olivine and enter mineral data:
  if OLI_pure_mineral:
    mol_min = _calculate_mineral_moles_OLIpure(OLIVINE,RC,path,filename)
  else: #olivine is not a pure mineral, but a solid solution
    mol_min = _calculate_mineral_moles_OLIss(OLIVINE,RC,path,filename)
  #error out if a value is negative
  if any(x<0 for x in mol_min) == True:
    print('Neg. mol values; error in rock comp. passed to peridotite_in6i')
    sys.exit()
  if OLI_pure_mineral:
    _edit_inrockfile(OLIVINE,mol_min[0],path,filename,'pure mineral')
  else:
    _edit_inrockfile(OLIVINE,mol_min[0],path,filename,'solid solution')
  _edit_inrockfile('ORTHOPYROXENE',mol_min[1],path,filename,'solid solution')
  _edit_inrockfile('CLINOPYROXENE',mol_min[2],path,filename,'solid solution')

def edit_in3i(string,newval,path,element,filename):
    '''
    This function allows the editing of variables inside a 3i file.
    '''
    check=True
    print(25*'-')
    print("Editing "+filename+": "+element)
    inp = open(path+filename,'r')
    out = open(path+'text.temp','w')
    for line in inp:
        if string in line:
           line.strip()
           if element == 'Temp':
              newval = '%.5E'%newval
              newline = line[0:27]+newval+line[38:]
           elif element == 'Density':
              newval = '%.5E'%newval
              newline = line[0:27]+newval+line[38:]
           elif element == 'TDS':
              newval = '%.5E'%newval
              newline = line[0:32]+newval+line[43:]
           elif element == 'balance':
              restart = 32+len(newval)
              newline = line[0:32]+newval+line[restart:]
           elif element == 'Gas':
              import numpy as np
              newval = '%.5E'%np.log10(newval)
              if float(newval) < 0.:
                newline = line[0:50]+newval+line[62:]
              else:
                newline = line[0:51]+newval+line[62:]
           elif (element == 'Aqueous') and (line.startswith(string)):
              newval = '%.5E'%float(newval)
              newline = line[0:51]+newval+line[62:]
           elif element == 'pH':
              newval = '%.5E'%float(newval)
              newline = line[0:51]+newval+line[62:]
           elif element == 'Unsuppress':
              newval = 'Molality  ' #<-- leave those spaces
              newline = line[0:63]+newval+line[73:]
           elif element == 'Suppress':
              newval = 'Suppressed'
              newline = line[0:63]+newval+line[73:]
           elif element == 'W:R':
              newval = '%.5E'%newval
              newline = line[0:27]+newval+line[38:]
           elif element == 'pCO2':
              import numpy as np
              lognewval = np.log10(newval)
              newval = '%.5E'%lognewval
              check = _check_in3i_pco2setup(path,filename)
              if check:
                if lognewval < 0:
                  newline = line[0:50]+newval+line[62:]
                else:
                  newline = line[0:50]+' '+newval+line[62:]
              else:
                newline = line
           else:
               check = False
               newline = line
           out.write(newline)
        else:
           out.write(line)
    inp.close()
    out.close()
    cmd = ('cp text.temp '+filename)
    os.system(cmd)

def edit_in6i(string,newval,path,element,filename):
    '''
    This function allows the editing of variables inside a 6i file.
    '''
    print(25*'-')
    print("Editing "+filename+": "+element)
    inp = open(path+filename,'r')
    out = open(path+'text.temp','w')
    for line in inp:
      if line.startswith(string):
         line.strip()
         if element == 'Temp':
            newval = float(newval)
            newval = '%.5E'%newval
            newline = line[0:34]+newval+line[45:]
         elif element == 'Ximax':
            newval = '%.5E'%newval
            newline = line[0:28]+newval+line[39:]
         elif element == 'Ximin':
            newval = '%.5E'%newval
            newline = line[0:28]+newval+line[39:]
         elif element == 'dxixi':
            newval = '%.5E'%newval
            newline = line[0:33]+newval+line[45:]
         elif element == 'co2gas':
            newval = '%.6E'%newval
            newline = line[0:27]+newval+line[39:]

         else:
            print('Couldnt find requested 6i element during edit.')
            sys.exit()
         out.write(newline)
      else:
         out.write(line)
    inp.close()
    out.close()
    cmd = ('cp text.temp '+filename)
    os.system(cmd)

def findline_in36o(string,path,filename,EndFileSignal=True):
    num_lines = sum(1 for line in open(path+filename))
    #important that this is a REVERSED search!
    for line in reversed(open(path+filename).readlines()):
      if line.startswith(string):
         break
      else:
        num_lines=num_lines-1
    if num_lines == 0 and EndFileSignal:
      print('error in '+filename)
      print('findline_in36o failed. Check your string: '+string)
    return num_lines

def findvar_in36o(string,path,element,filename):
   for line in reversed(open(path+filename).readlines()):
     if line.startswith(string):
        line.strip()
        if element == 'Temp':
          val = line[15:22]
          break
        elif element == 'W:R':
          val = line[27:38]
        elif element == 'pH':
          val = -float(line.split()[4])
          break
        elif element == 'Aw':
          val = line[37:43]
          break
        elif element == 'Aqueous':
          val = line.split()[1]
          break
        elif element == 'Gamma':
          val = line.split()[3]
          break
        elif element == 'SS component':
          val = line[30:41]
        elif element == 'SS mineral':
          from itertools import islice
          minline = findline_in36o(string,path,filename)
          with open(path+filename) as f:
            lines = islice(f,minline+5,minline+6)
            for line in lines:
              if line.startswith('|->|Amount'):
               val = line[31:42]
        elif element == 'Ionic':
          val = line[37:48]
          break
        elif element == 'SolMass':
           val = line.split(' ')[-2]
           break
        else:
          print('Couldnt find requested 36o element.')
   #check if value is actually a float
   try:
     float(val)
   except:
     val = ''.join(c for c in val if (c.isdigit() or c == '.'))
   return float(val)

def run_eq(which_eq,eqpath,cwdpath,filename,datafilekey,platform):
    if (platform == 'linux2') or (platform == 'linux'):
       conv = _run_eq_linux(which_eq,eqpath,cwdpath,filename,datafilekey)
    else:
       print('platform not setup: '+platform)
    return conv

def set_temp_option_in6i(filename,path,model):
    '''
    Change Temperature model in 6i file
    '''
    inp = open(path+filename,'r')
    out = open(path+'text.temp','w')
    text ='|Temperature option (jtemp)'
    lstart = findline_in36o(text,path,filename)
    lcount=0 #line counter
    for line in inp:
      if (lcount > lstart-1 and lcount < lstart+8):
        if model in line:
          s=list(line);s[4]='x'
          line=''.join(s)
        else:
          s=list(line);s[4]=' '
          line=''.join(s)
      out.write(line)
      lcount=lcount+1
    inp.close()
    out.close()
    cmd = ('cp text.temp '+filename)
    os.system(cmd)

def stitch(f1,f2,f3):
    print(25*'-')
    print( "Stitching EQ6 file.")
    weq = f1[-2] #which EQ
    cmd = ('cat '+f2+' '+f1[:-2]+weq+'p'+' > '+f3)
    os.system(cmd)

def actcoefmodel_in36i(filename,path,model):
    '''
    Change activity coefficient model:
    model = "B-dot equation"
    model = "Pitzer\'s equations"
    '''
    inp = open(path+filename,'r')
    out = open(path+'text.temp','w')
    text ='|iopg(1)'
    lstart = findline_in36o(text,path,filename)
    lcount=0 #line counter
    for line in inp:
      if (lcount > lstart-1 and lcount < lstart+4):
        if model in line:
          s=list(line);s[4]='x'
          line=''.join(s)
        else:
          s=list(line);s[4]=' '
          line=''.join(s)
      out.write(line)
      lcount=lcount+1
    inp.close()
    out.close()
    cmd = ('cp text.temp '+filename)
    os.system(cmd)

def convert_pitz2bdot(filename,cwdpath):
  '''
  Converts a eq3 file from pitzer to bdot
  '''
  #Replace var names to enable ypf --> mbn
  efs = False #EndFileSignal flag
  _swapvarname_in6i(filename,cwdpath,'|Ca++','|Ca+2',efs)
  _swapvarname_in6i(filename,cwdpath,'|Mg++','|Mg+2',efs)
  _swapvarname_in6i(filename,cwdpath,'|Fe++','|Fe+2',efs)
  _swapvarname_in6i(filename,cwdpath,'|O2(aq)','|O2,aq ',efs)
  _swapvarname_in6i(filename,cwdpath,'|H2(aq)','|H2,aq ',efs)
  _swapvarname_in6i(filename,cwdpath,'|SiO2(aq)','|SiO2,aq ',efs)
  _swapvarname_in6i(filename,cwdpath,'|CO2(g)','|CO2,g ',efs)
  _swap_case_in6i(filename,cwdpath,'|Dolomite','upper',efs)
  _swap_case_in6i(filename,cwdpath,'|Halite','upper',efs)
  _swap_case_in6i(filename,cwdpath,'|Hematite','upper',efs)
  _swap_case_in6i(filename,cwdpath,'|Quartz','upper',efs)
  _swap_case_in6i(filename,cwdpath,'|Talc','upper',efs)
  _swap_case_in6i(filename,cwdpath,'|Sylvite','upper',efs)
  #Update Pitzer --> B-Dot
  model="B-dot equation"
  actcoefmodel_in36i(filename,cwdpath,model)  

def swapheader_in6i(path,tophalf,filename6i):
    #Do the swap
    string='* Start of the bottom'
    num_line = findline_in36o(string,path,filename6i)
    lines = open(path+filename6i).readlines()
    open('file6p.temp', 'w').writelines(lines[num_line-1:])
    stitch('file6p.temp',tophalf,filename6i)

def edit_iopt_in_6i(filename,path,iopt,model):
    '''
    Sets iopt(5): Clear the ES Solids Read from the INPUT File:
    :
    model = "Don't do it"
    model = "Do it"
    '''
    print(25*'-')
    inp = open(path+filename,'r')
    out = open(path+'temp.txt','w')
    iopttext ='|iopt('+str(iopt)+')'
    lstart = findline_in36o(iopttext,path,filename)
    lcount=0 #line counter
    for line in inp:
      if (lcount > lstart-1 and lcount < lstart+2):
        if model in line:
          s=list(line);s[4]='x'
          line=''.join(s)
        else:
          s=list(line);s[4]=' '
          line=''.join(s)
      out.write(line)
      lcount=lcount+1
    inp.close()
    out.close()
    cmd = ('cp temp.txt '+filename)
    os.system(cmd)
    print('NOTE: '+str(iopttext[1:])+' is set to: '+model)

def generate_tablefrom6o(process,cwdpath,eqpath,raw6opath,
  csv_destination):
  import numpy as np
  import pandas as pd
  from itertools import islice
  #process the table based on the type of files loaded
  if process == 'waterrock':
    key = 'r6'
  elif process == 'evaporation':
    key = 'rb'
  else:
    print('Process not defined. Error in generate_tablefrom6o.')
  #Read solunit file for units of solids
  try:
    f = open("solunit.temp", "r")
    solunit=f.read()
    f.close()
  except:
    solunit = 'Moles'
  #Read all output files for plotting
  outpath = cwdpath + raw6opath
  allfiles = np.array([f for f in os.listdir(outpath)])
  all6ofiles_bool = np.array([f.startswith(key) for f in os.listdir(outpath)])
  try:
    all6ofiles = allfiles[all6ofiles_bool]
  except:
    print('Error: No 6o files found.');sys.exit() 
  #iterate through all the 6o files
  for file6o in all6ofiles:
    print('translating '+file6o+' into csv...')
    #create empty dataframe
    dfm = pd.DataFrame({'Variables': []})
    iter=0
    #find number of datablocks
    for line in (open(outpath+file6o).readlines()):
      if line.startswith(20*' '+'Xi='):
        iter=iter+1
    blocklines=np.zeros(iter)
    #find total number of lines in file
    total = _find_totalines(outpath,file6o)
    #extract line numbers corresponding to datablocks
    iter=0
    for line in (open(outpath+file6o).readlines()):
      if line.startswith(20*' '+'Xi='):
        lstart = findline_in36o(line,outpath,file6o)
        blocklines[iter]=lstart
        iter = iter+1
    blocklines=np.append(blocklines, total)
    blocklines=np.unique(blocklines) #remove duplicates
    #extract datablocks themselves and save to temp file
    iter = 0
    for step in range(0,len(blocklines)-1):
      block = islice(open(outpath+file6o),int(blocklines[iter]-1),\
        int(blocklines[iter+1]))
      with open(cwdpath+'block.temp','w') as out:
        for lines in block:
          out.write(lines)
      #extract xi
      with open(cwdpath+'block.temp','r') as out:
        first_line = out.readline()
      xi=float(first_line.split()[1])
      #extract temperature, Ionic strength, water activity, and TDS
      with open(cwdpath+'block.temp','r') as out:
        lines = out.readlines()
        for line in lines:
          if line.startswith(' Temperature'):
            T=float(line.split()[1])
          if line.startswith('                  Activity of water'):
            Aw=float(line.split()[3])
          if line.startswith('                 Ionic strength'):
            I = float(line.split()[3])
          if line.startswith(' NBS'):
            pH = float(line.split()[3])
          if line.startswith('                    Solute fraction'):
            TDS_perc = float(line.split()[2])*100.
          if line.startswith('                       Solvent mass'):
            Solvent_Mass = float(line.split()[2])
      #place values in dataframe
      data = {'Variables':['Xi','Aw','Temp','I', 'pH', 'TDS%','Solv_Mass'],\
         'Step'+str(step):[xi,Aw,T,I,pH,TDS_perc,Solvent_Mass]}
      dfin = pd.DataFrame(data)
      #extract aqueous species
      aq_species,aq_values,aq_gammas= _findaqueous_in36o(cwdpath,'block.temp')
      aq_specgam = ['g'+x for x in aq_species]
      data = {'Variables':aq_species, 'Step'+str(step):aq_values}
      dfaq = pd.DataFrame(data)
      data = {'Variables':aq_specgam, 'Step'+str(step):aq_gammas}
      dfaqg = pd.DataFrame(data)
      #extract total aqueous HCO3-, Na+
      #tDIC = _find_totalaqueous('HCO3-',cwdpath,'block.temp')
      tDIC = 999
      tNa = _find_totalaqueous('Na+',cwdpath,'block.temp')
      data = {'Variables':['total_HCO3-','total_Na+'],
        'Step'+str(step):[str(tDIC),str(tNa)]}
      dftaq = pd.DataFrame(data)
      #extract total carbon
      tco2 = _findtotal_co2(cwdpath,'block.temp')
      data = {'Variables':['TCO2'], 'Step'+str(step):[tco2]}
      dftc = pd.DataFrame(data)
      #extract total formate
      tFor = _findtotal_formate(cwdpath,'block.temp')
      data = {'Variables':['TFor'], 'Step'+str(step):[tFor]}
      dftf = pd.DataFrame(data)
      #extract solid species
      sol_species,sol_values = _findsolids_in6o(eqpath,cwdpath,
        'block.temp',solunit,process)
      data = {'Variables':sol_species, 'Step'+str(step):sol_values}
      dfsol = pd.DataFrame(data)
      #concat results
      df = pd.concat([dfin,dfaq,dfaqg,dftaq,dftc,dftf,dfsol], ignore_index=True)
      #remove elements that have values of zero
      df = df.replace('0.0000E+00',np.nan)
      df = df.dropna()
      #merge with master dataframe
      dfm = pd.merge(dfm,df,on='Variables',how='outer')
      #iterate block counter
      iter = iter + 1
    #print dataframe to csv
    dfm.set_index('Variables').to_csv(file6o+'.csv',index=True)
    os.system('mv '+file6o+'.csv '+ csv_destination)
    print(file6o+'.csv created')
    del df;del dfm

#######################################################################
# Secondary function (called from inside this file)
#######################################################################

def _edit_inrockfile(string,newval,path,filename,mintype):
    newval = '%.5E'%newval
    print(25*'-')
    print("Editing mineral in rock file: "+string)
    text = '|Reactant        |'+string
    l_mineral = findline_in36o(text,path,filename)
    lines = open(path+filename,'r').readlines()
    #change "Amount remaining in moles"
    newline1= lines[l_mineral+5][0:31]+newval+\
       lines[l_mineral+5][42:]
    lines[l_mineral+5] = newline1
    #change dxi(n)/dxi"
    if mintype == 'pure mineral':
       line2edit = 21
    elif mintype == 'solid solution':
       line2edit = 28
    newline2 = lines[l_mineral+line2edit][0:34]+newval+\
       lines[l_mineral+line2edit][45:]
    lines[l_mineral+line2edit]= newline2
    out = open(path+filename,'w')
    out.writelines(lines)
    out.close()

def _run_eq_linux(which_eq,eqpath,cwdpath,filename,datafilekey):
    import subprocess
    weq=str(which_eq)
    runeq36 = 'csh '+eqpath+'scripts/runeq'
    print(25*'-')
    print("Running EQ"+weq+' on '+filename)
    #run EQ & move output back to the cwd
    print('File xfer: cwd --> eqwork.')
    cmd = ('cp '+filename+' '+eqpath+'work/')
    os.system(cmd)
    cmd='cd '+eqpath+';cd work/;'+runeq36+weq+' '+datafilekey+' '+filename+' > eq'+weq+'out.log'
    #cmd = 'cd "' + eqpath + '"; cd work/; ' + runeq36 + weq + ' ' + datafilekey + ' "' + filename + '" > eq' + weq + 'out.log'
    subprocess.run(cmd, shell=True)
    print('File xfer: eqwork --> cwd.')
    cmd = ('cd '+eqpath+'work/; cp '+filename[:-2]+weq+'o '+cwdpath)
    os.system(cmd)
    cmd = ('cd '+eqpath+'work/; cp '+filename[:-2]+weq+'p '+cwdpath)
    os.system(cmd)
    cmd = ('cd '+eqpath+'work/; cp eq'+weq+'out.log '+cwdpath)
    os.system(cmd)
    if weq == '6':
      cmd = ('cd '+eqpath+ 'work; cp '+filename[:-2]+'csv '+cwdpath)
      os.system(cmd)
    #check log file to ensure convergence
    DidEQConverge = _converge_eq('eq'+weq+'out.log')
    if DidEQConverge == False:
      print('EQ'+weq+' doesnt appear to have converged. Check log file.')
    else:
      print('EQ'+weq+' properly converged.')
    #check if file exists
    if not os.path.isfile(filename[:-2]+weq+'o'):
      print('EQ'+weq+' .'+weq+'o file did not make it back to cwd.')
      sys.exit()
    return DidEQConverge

def _converge_eq(filename):
    #check if file exists
    if not os.path.isfile(filename):
      print('log file not created. Something went wrong.')
      sys.exit()
    else:
      for line in reversed(open(filename).readlines()):
        if line.startswith('  The following input files') or\
           line.startswith(' The following input files'):
          print('(ignore possible setenv error above.)')
          statement = line.split()
          if statement[6] == 'without':
            DidEQConverge = True
          else:
            DidEQConverge = False
    return DidEQConverge

def _swapvarname_in6i(filename,path,oldvar,newvar,EndFileSignal=True):
    print("Swapping "+oldvar+' to '+newvar)
    if len(oldvar) != len(newvar):
      print('check that both old and new var strings are same length')
      sys.exit()
    l_swap = findline_in36o(oldvar,path,filename,EndFileSignal)
    while l_swap !=0:
      lines = open(path+filename,'r').readlines()
      newline=newvar+lines[l_swap-1][len(newvar):]
      lines[l_swap-1] = newline
      out = open(path+filename,'w')
      out.writelines(lines)
      out.close()
      l_swap = findline_in36o(oldvar,path,filename,EndFileSignal)

def _swap_case_in6i(filename,path,var,case,EndFileSignal=True):
    l_swap = findline_in36o(var,path,filename,EndFileSignal)
    while l_swap != 0:
      lines = open(path+filename,'r').readlines()
      if case == 'upper':
        newline = lines[l_swap-1].upper()
      else:
        newline = lines[l_swap-1].lower()
      lines[l_swap-1]=newline
      out = open(path+filename,'w')
      out.writelines(lines)
      out.close()
      l_swap = findline_in36o(var,path,filename,EndFileSignal)

def _find_totalines(path,filename):
    ''' courtesy stackoverflow'''
    f = open(path+filename)
    lines = 0
    buf_size = 1024*1024
    read_f = f.read

    buf = read_f(buf_size)
    while buf:
      lines += buf.count('\n')
      buf = read_f(buf_size)
    return lines

def _findtotal_co2(path,file36o):
  import numpy as np
  from itertools import islice
  from itertools import compress
  tco2_species = ['MgCO3,aq','CaCO3,aq','MgHCO3+','CaHCO3+',
    'HCO3-','CO3-2','NaHCO3,aq','CO2,aq']
  alt_tco2_species = ['MgCO3(aq)','CaCO3(aq)','MgHCO3+',
    'CaHCO3+','HCO3-','CO3--','NaHCO3(aq)','CO2(aq)']
  aq_species,aq_values,aq_gammas = \
    _findaqueous_in36o(path,file36o)
  #get boolean
  try:
    tco2_bool = [x in tco2_species for x in aq_species]
  except:
    tco2_bool = [x in alt_tco2_species for x in aq_species]
  #filter
  tco2_species_infile = list(compress(aq_species,tco2_bool))
  tco2_values_infile  = list(compress(aq_values,tco2_bool))
  tco2_values_infile  = [float(x) for x in tco2_values_infile]
  total_co2           = sum(tco2_values_infile)
  return str('%.5E'%total_co2)

def _findtotal_formate(path,file36o):
  import numpy as np
  from itertools import islice
  from itertools import compress
  formate_species = ['FORMATE,aq','Mg(For)+','Mg(For)2,aq','Na(For),aq',
    'Fe(For)+','Fe(For)2,aq']
  aq_species,aq_values,aq_gammas = \
    _findaqueous_in36o(path,file36o)
  #get boolean
  formate_bool = [x in formate_species for x in aq_species]
  #filter
  formate_species_infile = list(compress(aq_species,formate_bool))
  formate_values_infile  = list(compress(aq_values,formate_bool))
  formate_values_infile  = [float(x) for x in formate_values_infile]
  total_formate          = sum(formate_values_infile)
  return str('%.5E'%total_formate)

def _findaqueous_in36o(path,file36o):
  import numpy as np
  from itertools import islice
  #find aqueous species in 3o or 6o file
  text ='                --- Distribution of Aqueous Solute Species'
  lstart = findline_in36o(text,path,file36o)
  text='    Species with molalities less than'
  lend= findline_in36o(text,path,file36o)
  #create empty array where those will be populated
  nsize = lend-lstart-5
  aqueous_species = [None]*nsize
  try:
    aqueous_values = [None]*nsize
    aqueous_gammas = [None]*nsize
  except:
    print('Error: Something happened. Maybe close .o file')
    sys.exit()
  #populate array
  row = 0 #counter
  with open(path+file36o) as f:
    lines = islice(f,lstart+3,lend-2)
    for line in lines:
        line=line.strip().split()
        try:
          float(line[1])
          aqueous_species[row] = line[0]
          aqueous_values[row]  = line[1]
          aqueous_gammas[row]  = line[3] #log gamma
        except:
          aqueous_species[row] = line[0]+' '+line[1]
          aqueous_values[row]  = line[2]
          aqueous_gammas[row]  = line[4] #log gamma
        row = row+1
    #EQ3/6 has a bug in how it reports E-100 values, in that it omits
    #the 'E', which screws things up. Here, if -100 is detected, it is
    #removed
    if '-100' in aqueous_species[-1]:
      aqueous_species = aqueous_species[:-1]
      aqueous_values  = aqueous_values[:-1]
      aqueous_gammas  = aqueous_gammas[:-1]
  return aqueous_species, aqueous_values, aqueous_gammas

def _find_totalaqueous(var,path,file36o):
  import numpy as np
  from itertools import islice
  #find species in 36o file
  text = ' Species Accounting for 99% or More of Aqueous '+var
  lstart = findline_in36o(text,path,file36o)
  if file36o[-2] == '3':
    text = '                --- Aqueous Redox'
  else:
    text = '                     --- Summary'
  lend= findline_in36o(text,path,file36o)
  nsize = lend-lstart
  with open(path+file36o) as f:
    lines = islice(f,lstart+3,lend-3)
    for line in lines:
      if not line.isspace():
        line=line.split()
        if line[0] == 'Subtotal':
          total_aqval = float(line[1])
          percentage  = float(line[2])
          total_aqval = 100.*total_aqval/percentage
          break
  return str('%.5E'%total_aqval) #Molality

def _findsolids_in6o(eqpath,path,file36o, solunit,process='waterrock',get_solid_solution=False):
  if process == 'evaporation':
    text = '                          Mass, grams'
  else:
    text = '            --- Grand Summary of Solid'
  #now find the solids
  solid_species, solid_values = __findsolids_in6o(eqpath,path,file36o,\
    solunit, text)
  return solid_species, solid_values

def __findsolids_in6o(eqpath,path,file36o, solunit,end_text):
  import numpy as np
  from itertools import islice
  #find species in 36o file
  text = '                     --- Summary of Solid'
  lstart = findline_in36o(text,path,file36o)
  lend= findline_in36o(end_text,path,file36o)
  #create empty array where those will be populated
  nsize = lend-lstart
  solid_species = [None]*nsize
  try:
    solid_values = [None]*nsize
  except:
    print('Error: Something happened. Maybe close .o file')
    sys.exit()
  #populate array
  row = 0 #counter
  with open(path+file36o) as f:
    lines = islice(f,lstart+3,lend-3)
    for line in lines:
      if line.startswith(' None'):
        continue
      if not line.isspace():
        #identify solid solutions
        if line[2] == ' ':
          #Extract name and values. Append '_out' for later filter
          solid_species[row] = ('-').join(line.split()[:-4])+'_pcss_out'
          solid_values[row]  = line.split()[-3]
        else:
          #Extract name and values. Append '_out' for later filter
          solid_species[row] = ('-').join(line.split()[:-4])+'_out'
          solid_values[row]  = line.split()[-3]
        if solunit=='Grams':
          solid_values[row]  = line.split()[-2]
        row = row+1
  #remove [None] from vector
  solid_species = list(filter(None, solid_species))
  solid_values = list(filter(None, solid_values))
  return solid_species, solid_values

def _calculate_mineral_moles_OLIpure(OLIVINE,RC,path,filename):
  rockmass = 1000. #rock mass in grams
  pure_mins  = [OLIVINE]
  wgm_pm = [147.] #weight in g/mol
  solid_sols = ['ENSTATITE','FERROSILITE',
              'DIOPSIDE','HEDENBERGITE']
  wgm_ss = [100., 132., 217., 249]  #weight in g/mol
  #search for the mole fractions in 6i file
  element='SS component'
  mfr=[findvar_in36o('|--->|'+x,path,element,filename) for x in solid_sols]
  #calculate mineral density.
  minerals = ['ORTHOPYROXENE','CLINOPYROXENE']
  weig_min = np.zeros(len(minerals)) # weight  of minerals in g/mol
  c = 0 #array counter
  for i in range(len(solid_sols)):
    if i % 2 == 0: #if i is even
      weig_min[c] = mfr[i]*wgm_ss[i] + mfr[i+1]*wgm_ss[i+1] #g/mol
      c=c+1
  weig_min = np.hstack((wgm_pm,weig_min))
  #calculate mineral moles
  mol_min = np.divide(np.multiply(rockmass,[x/100. for x in RC]),weig_min)
  return mol_min

def _calculate_mineral_moles_OLIss(OLIVINE,RC,path,filename):
  rockmass = 1000. #rock mass in grams
  solid_sols = ['FAYALITE','FORSTERITE','ENSTATITE','FERROSILITE',
              'DIOPSIDE','HEDENBERGITE']
  wgm_ss = [204., 141., 100., 132., 217., 249]  #weight in g/mol
  #search for the mole fractions in 6i file
  element='SS component'
  mfr=[findvar_in36o('|--->|'+x,path,element,filename) for x in solid_sols]
  #calculate mineral density.
  minerals = [OLIVINE, 'ORTHOPYROXENE','CLINOPYROXENE']
  weig_min = np.zeros(len(minerals)) # weight  of minerals in g/mol
  c = 0 #array counter
  for i in range(len(solid_sols)):
    if i % 2 == 0: #if i is even
      weig_min[c] = mfr[i]*wgm_ss[i] + mfr[i+1]*wgm_ss[i+1] #g/mol
      c=c+1
  #calculate mineral moles
  mol_min = np.divide(np.multiply(rockmass,[x/100. for x in RC]),weig_min)
  return mol_min
