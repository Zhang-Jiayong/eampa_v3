######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
import os
from labels import labels
from pwscf_output import pwscf_output
from potential_functions import potential_functions
from f2py_lib.f_interp import interp

class configs:

  @staticmethod
  def load():
    try:
      configs_dir = g.inp['configs']['dir']
    except:  
      return False
    if(not os.path.isdir(configs_dir)):
      return False      
    
    # Make file list
    g.configs['config_files'] = []
    g.configs['config_files'] = configs.load_files(configs_dir, g.configs['config_files'])
    
    # Read files in
    configs.read()
    
    #print(g.configs)
    
    return True

  @staticmethod
  def load_files(path_in, files):
    for path in os.listdir(path_in):
      path_new = path_in + "/" + path
      if(os.path.isdir(path_new)):
        files = configs.load_files(path_new, files)
      else:
        type = configs.file_type(path_new)
        files.append([path_new,type])
    return files
    
  @staticmethod
  def read():
    # Loop through files and read in
    for file in g.configs['config_files']:      
      if(file[1] == 'std'):        
        configs.std(file[0])
      elif(file[1] == 'qe'):
        configs.qe(file[0])
        
        


  @staticmethod
  def std(file_path):   
  
    # Read content from file
    content = std.file_to_list(file_path)
    content = std.prep_data(content)
    content = std.remove_comments(content)
  
    # Split up into individual configs
    config_list = []    
    temp = []
    last_line = None
    for line in content:
      if(last_line != None and line != '' and last_line != '' and last_line[0] != "#" and line[0] == "#"):
        config_list.append(temp)
        temp = []
      if(line != ''):  
        temp.append(line)      
        last_line = line
    if(len(temp)>0):
      config_list.append(temp)
      
    # Read each config  
    for i in range(len(config_list)):
      configs.add_config(config_list[i], file_path, i)
  
  @staticmethod
  def make_config():
    return {
    'file_type': '',  
    'file_path': '',  
    'file_part': 0,  
    'alat': 0.0,  
    'uv_prim': numpy.zeros((3,3,),),
    'uv': numpy.zeros((3,3,),),
    'stress': numpy.zeros((3,3,),),
    'c': numpy.zeros((3,),dtype=numpy.int32,),
    'rcut': 0.0,
    'rverlet': 0.0,
    'mtemp': {},
    'm': {},
    'coord_count_prim': 0,
    'coords_label_prim': None,
    'coords_label_id_prim': None,
    'coords_prim': None,
    'forces_prim': None,
    'coord_count': 0,
    'coords_label': None,
    'coords_label_id': None,
    'coords': None,
    'forces': None,
    'energy_per_atom': None,
    'energy': None,
    'e': 0,
    'f': 0,
    's': 0,
    'n_atoms_prim': 0,
    'n_atoms': 0,
    'l_units': 'ANG',
    'e_units': 'EV',
    'f_units': 'EV/ANG',
    's_units': 'EV/ANG3',
    }    
  
  @staticmethod
  def add_config(content, file_path, i, config_type='STANDARD'):
    fd = configs.make_config()
    fd['file_type'] = config_type
    fd['file_path'] = file_path
    fd['file_part'] = i
    
    # FIRST READ
    for line in content:
      line = line.strip() 
      f = std.to_fields(line, ' ')
      epa_set = False

      if(len(f) > 1 and f[0].upper() == "#ALAT"):
        fd['alat'] = float(f[1])
      if(len(f) > 3 and f[0].upper() == "#X"):
        fd['uv_prim'][0,0] = float(f[1])
        fd['uv_prim'][0,1] = float(f[2])
        fd['uv_prim'][0,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#Y"):
        fd['uv_prim'][1,0] = float(f[1])
        fd['uv_prim'][1,1] = float(f[2])
        fd['uv_prim'][1,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#Z"):
        fd['uv_prim'][2,0] = float(f[1])
        fd['uv_prim'][2,1] = float(f[2])
        fd['uv_prim'][2,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#SX"):
        fd['s'] = 1
        fd['stress'][0,0] = float(f[1])
        fd['stress'][0,1] = float(f[2])
        fd['stress'][0,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#SY"):
        fd['s'] = 1
        fd['stress'][1,0] = float(f[1])
        fd['stress'][1,1] = float(f[2])
        fd['stress'][1,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#SZ"):
        fd['s'] = 1
        fd['stress'][2,0] = float(f[1])
        fd['stress'][2,1] = float(f[2])
        fd['stress'][2,2] = float(f[3])
      if(len(f) > 2 and f[0].upper() == "#M"):
        fd['mtemp'][f[1]] = f[2]      
      if(len(f) > 3 and f[0].upper() == "#C"):
        fd['c'][0] = int(f[1])
        fd['c'][1] = int(f[2])
        fd['c'][2] = int(f[3])
      if(len(f) > 1 and f[0].upper() == "#E"):
        fd['energy'] = float(f[1])
        fd['e'] = 1    
      if(len(f) > 1 and f[0].upper() == "#EPA"):
        fd['energy_per_atom'] = float(f[1])
        fd['e'] = 1   
      if(len(f) > 1 and f[0].upper() == "#RCUT"):
        fd['rcut'] = f[1]
        if(fd['rverlet'] == 0):
          fd['rverlet'] = f[1]        
      if(len(f) > 1 and f[0].upper() == "#RVERLET"):
        fd['rverlet'] = f[1]      
      if(len(f) > 1 and f[0].upper() == "#L_UNITS"):
        fd['l_units'] = f[1]  
      if(len(f) > 1 and f[0].upper() == "#E_UNITS"):
        fd['e_units'] = f[1]  
      if(len(f) > 1 and f[0].upper() == "#F_UNITS"):
        fd['f_units'] = f[1]
      if(len(f) > 1 and f[0].upper() == "#S_UNITS"):
        fd['s_units'] = f[1]
      if(len(f) >= 4 and f[0][0] != "#"):    
        count_coord = False
        if(len(f) >= 4):
          count_coord = True
        if(len(f) >= 7):
          fd['f'] = 1
        if(count_coord):
          fd['coord_count_prim'] = fd['coord_count_prim'] + 1
          
 
    # COORD SIZE
    c = numpy.identity(3)
    
    # EXPAND
    c[0,0] = float(fd['c'][0]) 
    c[1,1] = float(fd['c'][1]) 
    c[2,2] = float(fd['c'][2]) 
    fd['coord_count'] = (fd['c'][0] * fd['c'][1] * fd['c'][2]) * fd['coord_count_prim']
    fd['uv'] = numpy.matmul(c, fd['uv_prim'])
    
    # MAKE ARRAYS
    fd['coords_label_prim'] = [None] * fd['coord_count_prim']
    fd['coords_label_id_prim'] = numpy.zeros((fd['coord_count_prim'],), dtype=numpy.int32,)
    fd['coords_prim'] = numpy.zeros((fd['coord_count_prim'],3,),)
    if(fd['f'] == 1):
      fd['forces_prim'] = numpy.zeros((fd['coord_count_prim'],3,),)
    
    fd['coords_label'] = [None] * fd['coord_count']
    fd['coords_label_id'] = numpy.zeros((fd['coord_count'],), dtype=numpy.int32,)
    fd['coords'] = numpy.zeros((fd['coord_count'],3,),)
    if(fd['f'] == 1):
      fd['forces'] = numpy.zeros((fd['coord_count'],3,),)
        
    fd['n_atoms_prim'] = fd['coord_count_prim']
    fd['n_atoms'] = fd['coord_count']   
    

    # READ COORDS
    n = 0
    for line in content:
      line = line.strip() 
      f = std.to_fields(line, ' ')
      if(len(f) >= 4 and f[0][0] != "#"):    
        if(len(f) >= 4):
          label_str, label_id = labels.add(f[0])
          fd['coords_label_prim'][n] = label_str
          fd['coords_label_id_prim'][n] = label_id
          fd['coords_prim'][n,0] = float(f[1])
          fd['coords_prim'][n,1] = float(f[2])
          fd['coords_prim'][n,2] = float(f[3])
          
        if(len(f) >= 7):
          fd['forces_prim'][n,0] = float(f[4])
          fd['forces_prim'][n,1] = float(f[5])
          fd['forces_prim'][n,2] = float(f[6])
        n = n + 1
        
    # EXPAND COORDS
    m = 0
    for i in range(fd['c'][0]):
      for j in range(fd['c'][1]):
        for k in range(fd['c'][2]):
          for n in range(fd['coord_count_prim']):
            fd['coords_label'][m] = fd['coords_label_prim'][n]
            fd['coords_label_id'][m] = fd['coords_label_id_prim'][n]
            fd['coords'][m,0] = (i + fd['coords_prim'][n,0]) / fd['c'][0]
            fd['coords'][m,1] = (j + fd['coords_prim'][n,1]) / fd['c'][1]
            fd['coords'][m,2] = (k + fd['coords_prim'][n,2]) / fd['c'][2]
            if(fd['f'] == 1):
              fd['forces'][m,0] = fd['forces_prim'][n,0]
              fd['forces'][m,1] = fd['forces_prim'][n,1]
              fd['forces'][m,2] = fd['forces_prim'][n,2]
            m = m + 1
            
    # SORT ENERGY       
    if(fd['e'] == 1 and fd['energy_per_atom'] != None and fd['energy'] == None):
      try:
        fd['energy'] = fd['energy_per_atom'] * fd['n_atoms']
      except:
        pass            
    if(fd['e'] == 1 and fd['energy'] != None and fd['energy_per_atom'] == None):
      try:
        fd['energy_per_atom'] = fd['energy'] / fd['n_atoms_prim']
        fd['energy'] = fd['energy_per_atom'] * fd['n_atoms']
      except:
        pass

    
    # APPEND
    g.configs['configs'].append(fd)




  #
  #  FILE TYPE (standard, quantum espresso)
  #
        
  @staticmethod
  def file_type(file_path):
    content = std.file_to_list(file_path)
    
    # Check if standard file
    count = 0
    for line in content:
      fields = line.split(" ")
      if(fields[0].upper() == "#ALAT"):
        count = count + 1
      if(fields[0].upper() == "#X"):
        count = count + 1
      if(fields[0].upper() == "#Y"):
        count = count + 1
      if(fields[0].upper() == "#Z"):
        count = count + 1
    if(count >= 4):
      return 'std'
  
    # Check if pwscf/qe file
    count = 0
    for line in content:
      if(line.strip()[0:13] == "Program PWSCF"):
        count = count + 1
      if(line.strip()[0:27] == "bravais-lattice index     ="):
        count = count + 1
      if(line.strip()[0:27] == "kinetic-energy cutoff     ="):
        count = count + 1
      if(line.strip()[0:27] == "mixing beta               ="):
        count = count + 1
      if(line.strip()[0:27] == "Exchange-correlation      ="):
        count = count + 1
      if(line.strip()[0:9] == "JOB DONE."):
        count = count + 1
    if(count >= 4):
      return 'qe'   
      

  # Espresso files

  @staticmethod
  def qe(file_path):   
    
    qe = pwscf_output(file_path)
    xyz = qe.make_xyz()    
    n = 0
    for xyz_inner in xyz:
      n = n + 1
      configs.add_config(xyz_inner, file_path, n, config_type='QE')



    
  @staticmethod
  def output(): 
    if(g.outputs):     
      fh = open(g.dirs['output'] + '/' + 'configs.dat', 'w')
      
      n = 0
      for c in g.configs['configs']:
              
        n = n + 1
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('CONFIG ' + str(n) + '\n')
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('file_type  ' + str(c['file_type']) + ' \n')
        fh.write('file_path  ' + str(c['file_path']) + ' \n')
        fh.write('file_part  ' + str(c['file_part']) + ' \n')
        fh.write('alat  ' + str(c['alat']) + ' \n')
        fh.write('uv_prim   ' + str(c['uv_prim'][0,0]) + ' ' + str(c['uv_prim'][0,1]) + ' ' + str(c['uv_prim'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv_prim'][1,0]) + ' ' + str(c['uv_prim'][1,1]) + ' ' + str(c['uv_prim'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv_prim'][2,0]) + ' ' + str(c['uv_prim'][2,1]) + ' ' + str(c['uv_prim'][2,2]) + ' ' + ' \n')
        fh.write('uv        ' + str(c['uv'][0,0]) + ' ' + str(c['uv'][0,1]) + ' ' + str(c['uv'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv'][1,0]) + ' ' + str(c['uv'][1,1]) + ' ' + str(c['uv'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv'][2,0]) + ' ' + str(c['uv'][2,1]) + ' ' + str(c['uv'][2,2]) + ' ' + ' \n')
        fh.write('stress    ' + str(c['stress'][0,0]) + ' ' + str(c['stress'][0,1]) + ' ' + str(c['stress'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['stress'][1,0]) + ' ' + str(c['stress'][1,1]) + ' ' + str(c['stress'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['stress'][2,0]) + ' ' + str(c['stress'][2,1]) + ' ' + str(c['stress'][2,2]) + ' ' + ' \n')
        fh.write('c         ' + str(c['c'][0])  + ' ' + str(c['c'][1]) + ' ' + str(c['c'][2]) + ' ' + ' \n')
        fh.write('rcut      ' + str(c['rcut']) + ' \n')
        fh.write('rverlet   ' + str(c['rverlet']) + ' \n')
        fh.write('e         ' + str(c['e']) + ' \n')
        fh.write('f         ' + str(c['f']) + ' \n')
        fh.write('s         ' + str(c['s']) + ' \n')
        
        fh.write('\n')        
        fh.write('COORDS PRIM (' + str(c['coord_count_prim']) + ')\n')
        for k in range(c['coord_count_prim']):        
          fh.write(str(c['coords_label_prim'][k]) + '  ')  
          fh.write(str(c['coords_prim'][k,0]) + '  ')  
          fh.write(str(c['coords_prim'][k,1]) + '  ')  
          fh.write(str(c['coords_prim'][k,2]) + '  ')             
          if(c['f'] == 1):     
            fh.write(str(c['forces_prim'][k,0]) + '  ')  
            fh.write(str(c['forces_prim'][k,1]) + '  ')  
            fh.write(str(c['forces_prim'][k,2]) + '  ')           
          fh.write('\n')  
        fh.write('\n')    
        
        fh.write('\n')        
        fh.write('COORDS (' + str(c['coord_count']) + ')\n')
        for k in range(c['coord_count']):        
          fh.write(str(c['coords_label'][k]) + '  ')  
          fh.write(str(c['coords'][k,0]) + '  ')  
          fh.write(str(c['coords'][k,1]) + '  ')  
          fh.write(str(c['coords'][k,2]) + '  ')             
          if(c['f'] == 1):     
            fh.write(str(c['forces'][k,0]) + '  ')  
            fh.write(str(c['forces'][k,1]) + '  ')  
            fh.write(str(c['forces'][k,2]) + '  ')           
          fh.write('\n')  
        fh.write('\n')    
             
             
             
        fh.write('\n')               
        fh.write('\n')  
         

    fh.close()









  # Run after configs read in
  # After dft_energy_adjustments
    
  @staticmethod
  def complete(): 
    #print("Complete Configs")
    #print(g.configs['configs'])
    
    # CONVERT
    for i in range(len(g.configs['configs'])):
    
      g.configs['configs'][i]['alat'] = units.convert(g.configs['configs'][i]['l_units'], 'ang', g.configs['configs'][i]['alat'])
      g.configs['configs'][i]['l_units'] = 'ANG'
      
      #print(g.configs['configs'][i]['energy'], g.configs['configs'][i]['energy_per_atom'])
      
      g.configs['configs'][i]['energy_per_atom'] = units.convert(g.configs['configs'][i]['e_units'], 'ev', g.configs['configs'][i]['energy_per_atom'])
      g.configs['configs'][i]['energy'] = units.convert(g.configs['configs'][i]['e_units'], 'ev', g.configs['configs'][i]['energy'])
      g.configs['configs'][i]['e_units'] = 'EV'
      #print(g.configs['configs'][i]['energy'], g.configs['configs'][i]['energy_per_atom'])
      
      # print(g.configs['configs'][i]['f_units'])
      if(g.configs['configs'][i]['f'] == 1):
        for na in range(len(g.configs['configs'][i]['forces_prim'])):
          for nb in range(len(g.configs['configs'][i]['forces_prim'][na])):
            g.configs['configs'][i]['forces_prim'][na,nb] = units.convert(g.configs['configs'][i]['f_units'], 'EV/ANG', g.configs['configs'][i]['forces_prim'][na,nb])
      
      if(g.configs['configs'][i]['f'] == 1):
        for na in range(len(g.configs['configs'][i]['forces'])):
          for nb in range(len(g.configs['configs'][i]['forces'][na])):
            g.configs['configs'][i]['forces'][na,nb] = units.convert(g.configs['configs'][i]['f_units'], 'EV/ANG', g.configs['configs'][i]['forces'][na,nb])
     
     
     
      if(g.configs['configs'][i]['f'] == 1):
        for na in range(len(g.configs['configs'][i]['stress'])):
          for nb in range(len(g.configs['configs'][i]['stress'][na])):
            g.configs['configs'][i]['stress'][na,nb] = units.convert(g.configs['configs'][i]['s_units'], 'EV/ANG3', g.configs['configs'][i]['stress'][na,nb])
    
    
    # ADJUST ENERGY QE
    for i in range(len(g.configs['configs'])):
      if(g.configs['configs'][i]['file_type'] == 'QE'):
        #print(g.configs['configs'][i]['coord_count'], g.configs['configs'][i]['coord_count_prim'])
        for k in range(g.configs['configs'][i]['coord_count']):      
          e_adj = g.dft_energy_adjustments[g.configs['configs'][i]['coords_label_id'][k]]          
          g.configs['configs'][i]['energy'] = g.configs['configs'][i]['energy'] + e_adj['calc_apaev']
        g.configs['configs'][i]['energy_per_atom'] = g.configs['configs'][i]['energy'] / g.configs['configs'][i]['coord_count'] 
        #print(g.configs['configs'][i]['energy'], g.configs['configs'][i]['energy_per_atom'])  
          
          #print(e_adj['calc_apaev'])

        
        #c['coords_label'][k]
        #print('QE')
        
        
        

    #  g.dft_energy_adjustments[label_id]
  
  
  
  
  
  
  
        
      
###########################################################
# F2PY functions
###########################################################
  
  @staticmethod
  def efs_add_config():  
    for c in g.configs['configs']:
      efs.add_config( 
                     6.5,
                     c['alat'], 
                     c['uv'], 
                     c['coords_label_id'], 
                     c['coords'], 
                     c['energy'], 
                     c['forces'], 
                     c['stress'], 
                     c['f'], 
                     c['s']
                    )                 
    efs.make_nl()
  
  
  
  