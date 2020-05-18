######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from f2py_lib.f_es import es



class es_calc:

  def run():  
  
    print("Energies: Surface, Vacany etc") 


    # Setup ES
    es.init()
    potential.es_add_potentials()
    es_calc.es_add()
    es.energy()
    for i in range(es.cc):
      print(es.config_energy[i,:])

  @staticmethod
  def es_add():  
    # Test values
    
    rcut = 6.5
    alat = 4.04
    label = 1
    type = 3
      
    # Add Config
    es_id = es.add_es_config(rcut, alat, label, type)





