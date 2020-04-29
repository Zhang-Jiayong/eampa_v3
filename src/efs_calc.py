######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from f2py_lib.f_efs import efs
from f2py_lib.f_bp import bp
import matplotlib.pyplot as plt
import time


class efs_calc:

  def run_energy():
  
    print("Calc Energy") 
  
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy()
    efs_calc.output_energy()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
  
  def run_energy_force():
  
    print("Calc Energy and Forces") 
    
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy_force() 
    efs_calc.output_energy()
    efs_calc.output_forces()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
  
  
  def run_energy_force_stress():
  
    print("Calc Energy, Forces and Stress") 
    
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy_force_stress() 
    efs_calc.output_energy()
    efs_calc.output_forces()
    efs_calc.output_stress()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
  


  def output_energy():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['results'] + '/' + 'config_energies.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('ENERGY RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    for n in range(efs.cc):    
      std.write_file_line(fh, 'Config ' + str(n+1) + ':', t_pad, efs.config_energy[n,:], f_pad)
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()
  
  
  def output_forces():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['results'] + '/' + 'config_forces.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('FORCE RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('\n')
    for n in range(efs.cc): 
      fh.write('##################\n')
      fh.write('Config ' + str(n) + '\n')
      fh.write('##################\n')
      a = efs.key[n, 0] - 1
      b = efs.key[n, 1]
      for l in range(a, b):
        std.write_file_line(fh, str(efs.labels[l]) + ':', t_pad, efs.config_forces[l,:], f_pad)
      fh.write('\n')
     
      
    fh.write('\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()



  def output_stress():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['results'] + '/' + 'config_stresses.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('STRESS RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    for n in range(efs.cc):    
      a = n * 3
      std.write_file_line(fh, 'Config ' + str(n+1) + ':', t_pad, efs.config_stresses[a,:], f_pad)
      std.write_file_line(fh, '', t_pad, efs.config_stresses[a+1,:], f_pad)
      std.write_file_line(fh, '', t_pad, efs.config_stresses[a+2,:], f_pad)
      fh.write('\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()
  

  def output_rss():
    fh = open(g.dirs['results'] + '/' + 'rss_configs.txt', 'w')
    fh.write('###############################################################\n')
    fh.write('ENERGY RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################\n')
    for n in range(efs.cc):
      fh.write(str(efs.config_energy[n,0]) + ' ')
      fh.write(str(efs.config_energy[n,1]) + ' ')
      fh.write(str(efs.config_energy[n,2]) + ' ')
      fh.write(str(efs.config_energy[n,3]) + ' ')
      fh.write(str(efs.config_energy[n,4]) + ' ')
      fh.write(str(efs.config_energy[n,5]) + ' ')
      fh.write('\n')
    fh.write('###############################################################\n')
    fh.close()
