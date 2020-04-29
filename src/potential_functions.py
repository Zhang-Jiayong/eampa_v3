######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
import os

class potential_functions:


  # Lennard Jones Potential
  # p[0] = e
  # p[1] = rm
  @staticmethod
  def lennard_jones(r, p):
    return p[0] * ((p[1] / r)**12 - 2 * (p[1]/r)**6)
    
    
  # Morse Potential
  # p[0] = d
  # p[1] = a
  # p[2] = re
  @staticmethod
  def morse(r, p):
    return p[0] * (numpy.exp(-2 * p[1] * (r - p[2]) - 2 * numpy.exp (-p[1]*(r - p[2]))))

    
  # Buckingham Potential
  # p[0] = A
  # p[1] = B
  # p[2] = C
  @staticmethod
  def buckingham(r, p):
    return p[0] * numpy.exp(-1 * p[1] * r) - p[2] / r**6



  # Embedding Finnis-Sinclair
  @staticmethod
  def density_fs(r, p):
    return (r - p[0])**2




  
  # Embedding Finnis-Sinclair
  @staticmethod
  def embedding_fs(r, p):
    return p[0] * numpy.sqrt(r)



  
  # Embedding Mendelev, Han, Srolovitz, Ackland, Sun, Asta
  @staticmethod
  def embedding_1(r, p):
    return -1 * numpy.sqrt(r) + p[0] * r**2
    
    
    
    
    
    


    # MORSE
    # V(r) = d * (exp(-2 * a * (r - re) - 2 * exp (-a(r - re))))
    
    # POLY3
    # V(r) = a + br + cr^2 + dr^3