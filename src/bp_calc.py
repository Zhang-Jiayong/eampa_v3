######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from f2py_lib.f_efs import efs
from f2py_lib.f_bp import bp
from f2py_lib.f_bp import polyfit

"""
1. init()
This prepares the arrays and clear any data

2. add configs for BP testing (e.g. BCC, FCC, and so on)
This makes the configs, ghost config and the neighbour list

3. add potentials
Update with the latest potential

4. calculate energy (no forces, no stresses)
Calculates and saves the energy for each


"""



class bp_calc:

  def run():  
  
    print("Calc Bulk Properties") 
    
    # Setup BP
    bp.init()
    potential.bp_add_potentials()
    b_props.bp_add()
    
    # Calculate
    bp.energy()
    bp.calculate_bp()    
    
    # Output to File
    b_props.bp_output()
    
    # Plots
    b_props.bp_eos_plot()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    