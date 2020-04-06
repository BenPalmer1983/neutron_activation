##########
# GLOBALS

import numpy

class g: 
  
  dirs = {
         'wd': 'wd',
         'log': 'wd/log',  
         'plots': 'wd/plots',  
         'plots_svg': 'wd/plots/svg', 
         'plots_eps': 'wd/plots/eps', 
         'xs': None,  
         'results': 'wd/results',  
         }
  
  times = {
          'start' : 0.0,
          'end' : 0.0,
          'duration' : 0.0,
          }
          
  inp = {}    
  
  
  mat = {}
  mat_tally = {}
  target = {
           'height': 0.01,
           'width': 0.01,
           'depth': 0.001,
           'atoms_per_m3': 0.0,
           'atoms': 0.0,
           }
           
  experiment = {
               'flux': 1.0e10,
               'i_time': 300.0,
               'c_time': 300.0,
               'i_points': 100,
               'c_points': 200,
               'time_line': None,
               'activity': None,
               'gammas': [],
               'gamma_energy': 0.0,
               }
               
               
  results_fh = None
  log_fh = None
         
  file_counter = 0 
         
  def file_name():
    globals.file_counter = globals.file_counter + 1
    name = "file_"
    file_counter_str = str(globals.file_counter)
    while(len(file_counter_str) < 6):
      file_counter_str = '0' + file_counter_str
    name = name + file_counter_str    
    return name
         
         