################################################################
#    DFT Energy Adjustments
#
################################################################

import os
from std import std
from globals import globals as g
from units import units

"""
Inputs the actual cohesive energy of atoms and relaxed calculated energy so the total energy can be adjusted
"""

class e_adjust:

  @staticmethod
  def load():
    if('dft' not in g.inp.keys()):
      return None
    if('e_adjust' not in g.inp['dft'].keys()):
      return None
      
    dft_file = g.inp['dft']['e_adjust']
    if(os.path.isfile(dft_file)): 
      # Load csv to list
      fd = std.read_csv(dft_file)
      
      for row in fd:
        label_text = row[0].strip().upper()
        if(label_text != ''):
          label_str, label_id = labels.add(label_text)
          try:        
            relaxed_dft_ev = units.convert(row[3], "EV", float(row[2]) / int(row[1])) # relaxed per atom in eV
            coh_ev = units.convert(row[5], "EV", float(row[4]))
            apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
            g.dft_energy_adjustments[label_id] = {
                                                  'label_id': label_id,
                                                  'label_text': label_str,
                                                  'atom_count': int(row[1]),
                                                  'relaxed_energy': float(row[2]),
                                                  'relaxed_energy_unit': row[3],
                                                  'coh_energy': float(row[4]),
                                                  'coh_energy_unit': row[5],
                                                  'calc_relaxed_dft_ev': relaxed_dft_ev,
                                                  'calc_coh_ev': coh_ev,
                                                  'calc_apaev': apaev,
                                                 }
          except:
            pass  
          
      
      
      
      """
      g.rd['dft_energy_adjustments'] = {}
      
      
      dft_energy_adjustments
      for row in fd:
        label_text = row[0].strip().upper()
        if(label_text != '' and label_text not in g.rd['dft_energy_adjustments'].keys()):
          try:            
            relaxed_dft_ev = units.convert(row[3], "EV", float(row[2]) / int(row[1])) # relaxed per atom in eV
            coh_ev = units.convert(row[5], "EV", float(row[4]))
            apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
            g.rd['dft_energy_adjustments'][label_text] = {
                                                      'label': 0,
                                                      'label_text': label_text,
                                                      'atom_count': int(row[1]),
                                                      'relaxed_energy': float(row[2]),
                                                      'relaxed_energy_unit': row[3],
                                                      'coh_energy': float(row[4]),
                                                      'coh_energy_unit': row[5],
                                                      'calc_relaxed_dft_ev': relaxed_dft_ev,
                                                      'calc_coh_ev': coh_ev,
                                                      'calc_apaev': apaev,
                                                     }
          except:
            pass
      """