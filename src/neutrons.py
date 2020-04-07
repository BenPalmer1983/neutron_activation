######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from tendl import tendl
from isotopes import isotopes
from f2py.f_rng import rng
import matplotlib.pyplot as plt

class neutrons:

  def run():
    print("RUNNING")
    projectile_protons = 0
    projectile_neutrons = 1
    

    g.results_fh = open(g.dirs['results']+'/results.txt', 'w')
    
    #XS Dir
    g.dirs['xs'] = g.inp['xs_dir']['path']
    g.results_fh.write("XS Dir: " + str(g.dirs['xs']) + '\n')
    
    # LOAD MATERIAL
    g.mat = isotopes.make_material(g.inp['target']['isotopes'], g.inp['target']['rel_mass'], g.inp['target']['density'])
    
    if(g.inp['target']['units'].lower() == 'mm'):
      target_factor = 0.001
    elif(g.inp['target']['units'].lower() == 'cm'):
      target_factor = 0.01
    else:
      target_factor = 1.0
    
    g.target['height'] = target_factor * g.inp['target']['height']
    g.target['width'] = target_factor * g.inp['target']['width']
    g.target['depth'] = target_factor * g.inp['target']['depth']
    
    
    # LOAD EXPERIMENT
    if('experiment' in g.inp.keys()):
      g.experiment['flux'] = g.inp['experiment']['flux']
      g.experiment['i_time'] = g.inp['experiment']['i_time']
      g.experiment['c_time'] = g.inp['experiment']['c_time']
    
    
    # TOTAL ATOMS
    g.target['atoms_per_m3'] = 0.0
    g.target['atoms'] = 0.0
    for k in g.mat.keys():
      g.target['atoms_per_m3'] = g.target['atoms_per_m3'] + g.mat[k]['atoms_per_m3']
    g.target['atoms'] = g.target['height'] * g.target['width'] * g.target['depth'] * g.target['atoms_per_m3']
      

      
    # REACTION/PRODUCTION RATES
    g.mat_r = {}
    for k in g.mat.keys():
      key = k
      if(key not in g.mat_tally.keys()):
        g.mat_r[key] = {'protons': 0,'nucleons': 0}
      g.mat_r[k]['protons'] = g.mat[k]['protons']
      g.mat_r[k]['nucleons'] = g.mat[k]['nucleons']
      g.mat_r[k]['nd'] = g.mat[k]['atoms_per_m3']


    # MAKE TOTAL TALLY  
    mat_tally = {}
    
    # ALL THE ORIGINAL TARGET ISOTOPES
    for k in g.mat.keys():
      key = k
      protons = g.mat[k]['protons']
      nucleons = g.mat[k]['nucleons']
      percentage = g.mat[k]['mass_percentage']
      
      if(key not in g.mat_tally.keys()):
        g.mat_tally[key] = {
                            'protons': protons,
                            'nucleons': nucleons,
                            'mass_percentage': percentage,
                            'rr': 0.0,
                            'atoms_start': g.mat[k]['atoms_per_m3'] * g.target['height'] * g.target['width'] * g.target['depth'],
                            'atoms_irradiate': 0.0,
                            'atoms_cool': 0.0,
                            'type': 'Target',
                            'nd': g.mat_r[k]['nd'],
                           }    
 
    # RESIDUAL ISOTOPES
    for k in g.mat.keys():
      protons = g.mat[k]['protons']
      nucleons = g.mat[k]['nucleons']
      reaction_list = tendl.read_reactions_list(g.dirs['xs'], projectile_protons, projectile_neutrons, protons, nucleons)
      for r in reaction_list:
        rkey = int(r[3]) * 1000 + int(r[5])
        if(rkey not in g.mat_tally.keys()):
          g.mat_tally[rkey] = {
                              'protons': r[3],
                              'nucleons': r[5],
                              'mass_percentage': 0.0,
                              'rr': 0.0,
                              'atoms_start': 0.0,
                              'atoms_irradiate': 0.0,
                              'atoms_cool': 0.0,
                              'type': 'Residual',
                              'nd': 0.0,
                             }  

    keys = []
    for key in g.mat_tally.keys():
      keys.append(key)

    for k in keys:   
      protons = g.mat_tally[k]['protons']
      nucleons = g.mat_tally[k]['nucleons']                   
      u_iso = isotopes.unique_chain_isotopes(protons, nucleons)
      for u in u_iso:
        ukey = int(u[0]) * 1000 + int(u[1])
        if(ukey not in g.mat_tally.keys()):
          g.mat_tally[ukey] = {
                              'protons': u[0],
                              'nucleons': u[1],
                              'mass_percentage': 0.0,                              
                              'rr': 0.0,
                              'atoms_start': 0.0,
                              'atoms_irradiate': 0.0,
                              'atoms_cool': 0.0,
                              'type': 'Decay',
                              'nd': 0.0,
                              }  
        
    
    # PRE-CACHE REACTIONS
    for k in g.mat.keys():
      target_protons = g.mat[k]['protons']
      target_nucleons = g.mat[k]['nucleons']      
      tendl.cache_reactions(g.dirs['xs'], projectile_protons, projectile_neutrons, target_protons, target_nucleons)
      
      
    # PREPARE STORAGE ARRAYS
    neutrons.prep_arrays()

    g.results_fh.write('\n')
    g.results_fh.write('============================================================================' + '\n')
    g.results_fh.write('TARGET' + '\n')
    g.results_fh.write('============================================================================' + '\n')

    g.results_fh.write(neutrons.pad('Height:', 18) + ' ' + neutrons.pad(g.target['height'], 18) + '\n')
    g.results_fh.write(neutrons.pad('Width:', 18) + ' ' + neutrons.pad(g.target['width'], 18) + '\n')
    g.results_fh.write(neutrons.pad('Depth:', 18) + ' ' + neutrons.pad(g.target['depth'], 18) + '\n')
    g.results_fh.write('\n')
    g.results_fh.write(neutrons.pad('Volume:', 18) + ' ' + neutrons.pad(g.target['height'] * g.target['width'] * g.target['depth'], 18) + '\n')
    g.results_fh.write(neutrons.pad('Density:', 18) + ' ' + neutrons.pad(g.inp['target']['density'], 18) + '\n')
    g.results_fh.write(neutrons.pad('Mass:', 18) + ' ' + neutrons.pad(g.target['height'] * g.target['width'] * g.target['depth'] * g.inp['target']['density'], 18) + '\n')
    g.results_fh.write('\n')
    g.results_fh.write('\n')

    
    # Neutron Spectrum
    flux = g.experiment['flux']
    
    # Default
    x = numpy.linspace(1.0, 10.0, 20)
    y = neutrons.maxwell(x, 2.0)
    
    ne = x[:] * 1000000
    nf = y[:] * (flux / sum(y))
        
    
    if('nspectra' in g.inp.keys() and g.inp['nspectra']['type'][0] == 'point'):
      # nspectra type=point,1000000
      ne = [float(g.inp['nspectra']['type'][1])]
      nf = [flux]
    elif('nspectra' in g.inp.keys() and g.inp['nspectra']['type'][0] == 'mb'):
      # nspectra type=mb,0.001,10.0,2.0,1000000,50
      mb_min = float(g.inp['nspectra']['type'][1])
      mb_max = float(g.inp['nspectra']['type'][2])
      mb_a = float(g.inp['nspectra']['type'][3])
      mb_mult = float(g.inp['nspectra']['type'][4])
      mb_points = int(g.inp['nspectra']['type'][5])
      
      x = numpy.linspace(mb_min, mb_max, mb_points)
      y = neutrons.maxwell(x, mb_a)
      ne = x[:] * mb_mult
      nf = y[:] * (flux / sum(y))
      
      

    g.results_fh.write('\n')
    g.results_fh.write('============================================================================' + '\n')
    g.results_fh.write('Neutron Flux' + '\n')
    g.results_fh.write('============================================================================' + '\n')
    for k in range(len(ne)):
      g.results_fh.write(neutrons.pad(ne[k], 18) + ' ' + neutrons.pad(nf[k], 18))
      g.results_fh.write('\n')
    
    #plt.plot(ne, nf)
    #plt.show()
 
    for i in range(len(ne)):
      neutrons.reaction(ne[i], nf[i])
    
    # SAVE REACTION RATES
    neutrons.save_reaction_rates()
    
    # IRRADIATE
    neutrons.irradiate()
    neutrons.cool()   
    neutrons.save_amounts()
    neutrons.save_activities()
        
    # CALC GAMMAS
    neutrons.calc_gammas()

    g.results_fh.write('\n')
    g.results_fh.write('============================================================================' + '\n')
    g.results_fh.write('GAMMAS' + '\n')
    g.results_fh.write('============================================================================' + '\n')
    for i in range(10):
      tn = int(numpy.floor(i * (len(g.experiment['time_line'])/10)))
      g.results_fh.write(neutrons.pad(g.experiment['time_line'][tn], 18))
      g.results_fh.write(neutrons.pad(g.experiment['gamma_power'][tn], 18))
      g.results_fh.write('\n')
    
    
    
    # CLOSE LOG 
    g.results_fh.write('\n')
    g.results_fh.write('\n')
    g.results_fh.close()
    
    # MAKE PLOTS
    neutrons.make_plots()
    

    
   
  @staticmethod 
  def maxwell(x, a):
    return 0.797884561 * ((x * x * numpy.exp((-x**2)/(2*a**2))) / a**3)
    
    
  @staticmethod 
  def reaction(ne, nf):
  
    projectile_protons = 0
    projectile_neutrons = 1
    
    nt = (g.target['height'] * g.target['width']) / (0.01*0.01)
    
    for iso in g.mat_r.keys():
      target_protons = g.mat_r[iso]['protons']
      target_nucleons = g.mat_r[iso]['nucleons']
      nd = g.mat_r[iso]['nd']
            
      rs = tendl.read_reactions(g.dirs['xs'], projectile_protons, projectile_neutrons, target_protons, target_nucleons)
      for r in rs:
        residual_protons = r['residual_protons']
        residual_nucleons = r['residual_nucleons']
        if(not(residual_protons == target_protons and residual_nucleons == target_nucleons)):
          xs = tendl.get_xs(g.dirs['xs'], projectile_protons, projectile_neutrons, target_protons, target_nucleons, residual_protons, residual_nucleons, ne)
          rr = nt * xs * nd * nf * 1.0e-28 * g.target['depth']
          
          tkey = int(target_protons) * 1000 + int(target_nucleons)
          rkey = int(residual_protons) * 1000 + int(residual_nucleons)
          
          g.mat_tally[rkey]['rr'] = g.mat_tally[rkey]['rr'] + rr
          g.mat_tally[tkey]['rr'] = g.mat_tally[tkey]['rr'] - rr
    

    
  @staticmethod 
  def save_reaction_rates():  
    cols = ['Reaction Rate','Number Density','Stability','Half-Life']
    width = 22
    d = {}
    for key in g.mat_tally.keys():
      protons = g.mat_tally[key]['protons']
      nucleons = g.mat_tally[key]['nucleons']
      if(g.mat_tally[key]['rr'] != 0.0): 
        if(protons not in d.keys()):
          d[protons] = {} 
        stability = 'Unstable'
        if(isotopes.get_stability(protons, nucleons)):
          stability = 'Stable'
          half_life = ''
        else:
          half_life = isotopes.get_half_life(protons, nucleons)          
        d[protons][nucleons] = [g.mat_tally[key]['rr'], g.mat_tally[key]['nd'], stability, half_life]
    neutrons.save_data('REACTION RATES', cols, d)

      
  @staticmethod 
  def save_amounts():      
    cols = ['Atoms Start','Atoms Irradiated', 'Atoms Cooled']
    width = 22
    d = {}
    for key in g.mat_tally.keys():
      protons = g.mat_tally[key]['protons']
      nucleons = g.mat_tally[key]['nucleons']
      if(g.mat_tally[key]['atoms_start'] > 0.0 or g.mat_tally[key]['atoms_irradiate'] > 0.0 or g.mat_tally[key]['atoms_cool'] > 0.0): 
        if(protons not in d.keys()):
          d[protons] = {} 
        stability = 'Unstable'
        if(isotopes.get_stability(protons, nucleons)):
          stability = 'Stable'
          half_life = ''
        else:
          half_life = isotopes.get_half_life(protons, nucleons)          
        d[protons][nucleons] = [g.mat_tally[key]['atoms_start'],g.mat_tally[key]['atoms_irradiate'],g.mat_tally[key]['atoms_cool']]
    neutrons.save_data('ISOTOPE AMOUNTS', cols, d)

      
  @staticmethod 
  def save_activities():      
    i_points = g.experiment['i_points']
    c_points = g.experiment['c_points']
    cols = ['Activity Irradiated', 'Activity Cooled']
    width = 22
    d = {}
    for key in g.mat_tally.keys():
      protons = g.mat_tally[key]['protons']
      nucleons = g.mat_tally[key]['nucleons']
      if(g.mat_tally[key]['activity'][i_points-1, 1] > 0.0 and g.mat_tally[key]['activity'][i_points+c_points-2, 1] > 0.0): 
        if(protons not in d.keys()):
          d[protons] = {} 
        stability = 'Unstable'
        if(isotopes.get_stability(protons, nucleons)):
          stability = 'Stable'
          half_life = ''
        else:
          half_life = isotopes.get_half_life(protons, nucleons)          
        d[protons][nucleons] = [g.mat_tally[key]['activity'][i_points-1, 1],g.mat_tally[key]['activity'][i_points+c_points-2, 1]]
    neutrons.save_data('ISOTOPE ACTIVITIES', cols, d)
      

          
      
  @staticmethod  
  def irradiate():
    i_time = g.experiment['i_time']
    # IRRADIATE
    for k in g.mat_tally.keys():   
      if(g.mat_tally[k]['type'] == 'Target'):
        g.mat_tally[k]['atoms_irradiate'] = g.mat_tally[k]['atoms_irradiate'] + g.mat_tally[k]['atoms_start'] - i_time * g.mat_tally[k]['rr']          
      elif(g.mat_tally[k]['type'] == 'Residual'):
        decay = isotopes.isotope_activities(g.mat_tally[k]['protons'], g.mat_tally[k]['nucleons'], 0.0, g.mat_tally[k]['rr'], i_time)
        for p in decay.keys():
          for n in decay[p].keys():
            key = 1000 * p + n
            g.mat_tally[key]['atoms_irradiate'] = g.mat_tally[key]['atoms_irradiate'] + decay[p][n]
    
    # IRRADIATE OVER TIME
    for k in g.mat_tally.keys():
      time = g.mat_tally[k]['amount'][:, 0]
      for tn in range(len(time)):
        t = time[tn]
        if(t <= i_time):
          if(g.mat_tally[k]['type'] == 'Target'):
            g.mat_tally[k]['amount'][tn, 1] = g.mat_tally[k]['amount'][tn, 1] + g.mat_tally[k]['atoms_start'] - t * g.mat_tally[k]['rr']               
            
          elif(g.mat_tally[k]['type'] == 'Residual'):
            decay = isotopes.isotope_activities(g.mat_tally[k]['protons'], g.mat_tally[k]['nucleons'], 0.0, g.mat_tally[k]['rr'], t)
            for p in decay.keys():
              for n in decay[p].keys():
                key = 1000 * p + n
                g.mat_tally[key]['amount'][tn, 1] = g.mat_tally[key]['amount'][tn, 1] + decay[p][n]
                g.mat_tally[key]['activity'][tn, 1] =  g.mat_tally[key]['activity'][tn, 1] + isotopes.get_activity(p, n, decay[p][n])
                g.experiment['activity'][tn, 1] = g.experiment['activity'][tn, 1] + isotopes.get_activity(p, n, decay[p][n])
    


  @staticmethod 
  def cool():
    i_time = g.experiment['i_time']
    c_time = g.experiment['c_time']
    for k in g.mat_tally.keys():
      if(g.mat_tally[k]['type'] == 'Target'):
        g.mat_tally[k]['atoms_cool'] = g.mat_tally[k]['atoms_irradiate']
      elif(g.mat_tally[k]['type'] == 'Residual'):
        decay = isotopes.isotope_activities(g.mat_tally[k]['protons'], g.mat_tally[k]['nucleons'], g.mat_tally[k]['atoms_irradiate'], 0.0, c_time)
        for p in decay.keys():
          for n in decay[p].keys():
            key = 1000 * p + n
            g.mat_tally[key]['atoms_cool'] = g.mat_tally[key]['atoms_cool'] + decay[p][n]
    
    # COOL OVER TIME
    for k in g.mat_tally.keys():
      time = g.mat_tally[k]['amount'][:, 0]
      for tn in range(len(time)):
        t = time[tn]  - i_time
        if(t > 0):
        
          if(g.mat_tally[k]['type'] == 'Target'):
            g.mat_tally[k]['amount'][tn, 1] = g.mat_tally[k]['atoms_irradiate']
            g.mat_tally[key]['activity'][tn, 1] =  0.0
          elif(g.mat_tally[k]['type'] == 'Residual'):
            decay = isotopes.isotope_activities(g.mat_tally[k]['protons'], g.mat_tally[k]['nucleons'], g.mat_tally[k]['atoms_irradiate'], 0.0, t)
            for p in decay.keys():
              for n in decay[p].keys():
                key = 1000 * p + n
                g.mat_tally[key]['amount'][tn, 1] = g.mat_tally[key]['amount'][tn, 1] + g.mat_tally[key]['amount'][tn, 1] + decay[p][n]
                g.mat_tally[key]['activity'][tn, 1] = g.mat_tally[key]['activity'][tn, 1] + isotopes.get_activity(p, n, decay[p][n])
                g.experiment['activity'][tn, 1] = g.experiment['activity'][tn, 1] + isotopes.get_activity(p, n, decay[p][n])
        
        
        

  @staticmethod
  def prep_arrays():
    # ARRAYS TO STORE TIME/AMOUNT and TIME/ACTIVITY   
    i_time = g.experiment['i_time'] 
    c_time = g.experiment['c_time']
    
    i_points = g.experiment['i_points']
    c_points = g.experiment['c_points']
    m_size = i_points + c_points - 1
    
    a = 0
    b = i_points
    c = i_points
    d = m_size
    
    for key in g.mat_tally.keys():
      g.mat_tally[key]['amount'] = numpy.zeros((m_size, 2,),)
      g.mat_tally[key]['amount'][0:i_points, 0] = numpy.linspace(0.0, i_time, i_points)
      g.mat_tally[key]['amount'][i_points-1:m_size, 0] = numpy.linspace(i_time, i_time + c_time, c_points)
      
      g.mat_tally[key]['activity'] = numpy.zeros((m_size, 2,),)
      g.mat_tally[key]['activity'][0:i_points, 0] = numpy.linspace(0.0, i_time, i_points)
      g.mat_tally[key]['activity'][i_points-1:m_size, 0] = numpy.linspace(i_time, i_time + c_time, c_points)

    g.experiment['activity'] = numpy.zeros((m_size, 2,),)
    g.experiment['activity'][0:i_points, 0] = numpy.linspace(0.0, i_time, i_points)
    g.experiment['activity'][i_points-1:m_size, 0] = numpy.linspace(i_time, i_time + c_time, c_points)
 
    g.experiment['time_line'] = numpy.zeros((m_size,),)   
    g.experiment['time_line'][0:i_points] = numpy.linspace(0.0, i_time, i_points) 
    g.experiment['time_line'][i_points-1:m_size] = numpy.linspace(i_time, i_time + c_time, c_points)
     
    g.experiment['gammas'] = [] 
    g.experiment['gamma_energy'] = [] 
    g.experiment['gamma_power'] = [] 
    g.experiment['gamma_dose'] = [] 
    for i in range(m_size):
      g.experiment['gammas'].append({})
      g.experiment['gamma_energy'].append(0.0)
      g.experiment['gamma_power'].append(0.0)
      g.experiment['gamma_dose'].append(0.0)
    
    
    
  @staticmethod  
  def calc_gammas():    
    for tn in range(len(g.experiment['time_line'])):
      time = g.experiment['time_line'][tn]
      
      for key in g.mat_tally.keys():
        a = g.mat_tally[key]['activity'][tn, 1]
        if(a > 0.0e0):
          protons = g.mat_tally[key]['protons']
          nucleons = g.mat_tally[key]['nucleons']
          gammas = isotopes.get_gammas(protons, nucleons)
          if(gammas is not None and len(gammas) > 0):
            for i in range(len(gammas)):
              # Sum Total
              g.experiment['gamma_energy'][tn] = g.experiment['gamma_energy'][tn] + a * gammas[i,0] * gammas[i,1]
              g.experiment['gamma_power'][tn] = g.experiment['gamma_energy'][tn] * 1.60218e-13
              
              
              # Save individual lines
              if(gammas[i,0] not in g.experiment['gammas'][tn].keys()):
                g.experiment['gammas'][tn][gammas[i,0]] = 0.0
              g.experiment['gammas'][tn][gammas[i,0]] = g.experiment['gammas'][tn][gammas[i,0]] + a * gammas[i,1]
              
    for tn in range(len(g.experiment['time_line'])):
      a = 4.0 * 3.1415927 * g.dose['distance']
      fraction = g.dose['area'] / a
      g.experiment['gamma_dose'][tn] =  g.dose['time'] * (g.experiment['gamma_power'][tn] * fraction) / g.dose['mass']


  @staticmethod  
  def make_plots(): 
  
    # ACTIVITY
    ################################
  
    plt.clf()
    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')

    fig, axs = plt.subplots(1, 1, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    fig.suptitle('Activity vs Time')  
    
    plt.xlabel('Time/s')
    plt.ylabel('Activity/Bq')
    plt.yscale('symlog', linthreshy=1000)
    plt.plot(g.experiment['activity'][:, 0], g.experiment['activity'][:, 1], color='k', ls='solid')
    
    neutrons.plot_output('activity')

  
    # GAMMA ENERGY/J
    ################################
  
    plt.clf()
    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')

    fig, axs = plt.subplots(1, 1, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    fig.suptitle('Gamma Output vs Time')  
    
    plt.xlabel('Time/s')
    plt.ylabel('Power/Micro Watt')
    plt.plot(g.experiment['time_line'][:], g.experiment['gamma_power'][:], color='k', ls='solid')
    
    neutrons.plot_output('gamma_output')
    
    
  
    # GAMMA ENERGY/J
    ################################
  
    plt.clf()
    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')

    fig, axs = plt.subplots(1, 1, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    fig.suptitle('Estimated Dose vs Time (5 minutes exposure)')  
    
    plt.xlabel('Time/s')
    plt.ylabel('Dose/micro Sieverts')
    plt.plot(g.experiment['time_line'][:], g.experiment['gamma_dose'][:], color='k', ls='solid')

    neutrons.plot_output('gamma_dose')
    
    tn = g.experiment['i_points'] - 1
    neutrons.gamma_plot(tn)
    neutrons.plot_output('gamma_lines_' + str(g.experiment['time_line'][tn]))
    
    tn = g.experiment['i_points'] + g.experiment['c_points'] - 2
    neutrons.gamma_plot(tn)
    neutrons.plot_output('gamma_lines_' + str(g.experiment['time_line'][tn]))
    
    
    
  
    # ISOTOPE PLOTS
    ################################
    
    i_points = g.experiment['i_points']
    for key in g.mat_tally.keys():
      if(g.mat_tally[key]['activity'][i_points, 1] > 0.0):
        protons = g.mat_tally[key]['protons']
        nucleons = g.mat_tally[key]['nucleons']
        
        file_name = str(protons) + '_' + isotopes.get_symbol(protons) + '_' + str(nucleons) + '_activity'
    
        plt.clf()
    
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')

        fig, axs = plt.subplots(1, 1, figsize=(12,9))
        fig.tight_layout(pad=5.0)
        fig.suptitle('Activity vs Time ' + isotopes.get_symbol(protons) + ' ' + str(nucleons))  
    
        plt.xlabel('Time/s')
        plt.ylabel('Activity/Bq')
        if(max(g.mat_tally[key]['activity'][:, 1]) > 10000):
          plt.yscale('symlog', linthreshy=1000)
        else:  
          plt.yscale('linear')          
        plt.plot(g.mat_tally[key]['activity'][:, 0], g.mat_tally[key]['activity'][:, 1], color='k', ls='solid')    
        neutrons.plot_output(file_name)

    
    
    
    
    
  @staticmethod
  def gamma_plot(tn):        
    g_max = None
    for gamma in g.experiment['gammas'][tn].keys():
      if(g_max == None or g.experiment['gammas'][tn][gamma] > g_max):
        g_max = g.experiment['gammas'][tn][gamma]
    
    box_x = []
    box_y = []
    for gamma in g.experiment['gammas'][tn].keys():
      if(g.experiment['gammas'][tn][gamma] > 0.001 * g_max):
        box_x.append(gamma)
        box_y.append(g.experiment['gammas'][tn][gamma])
        
    w = (max(box_x) - min(box_x))/200
    
    plt.clf()
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(1, 1, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    fig.suptitle('Gamma Lines at t=' + str(g.experiment['time_line'][tn])) 
    axs.bar(box_x, box_y, color=(0.2, 0.2, 0.2, 0.6), width=w) 
    axs.set_xlabel('Gamma Energy/eV')
    axs.set_ylabel('Activity/Bq')
    
    
    
    
  @staticmethod
  def plot_output(file_name):
    plt.savefig(g.dirs['plots_svg'] + '/' + file_name + '.svg', format='svg')
    plt.savefig(g.dirs['plots_eps'] + '/' + file_name + '.eps', format='eps')
    plt.savefig(g.dirs['plots_png'] + '/' + file_name + '.png', format='png')
    plt.close()

  @staticmethod  
  def pad(inp, length=12):
    if(inp == None):
      inp = ''     
    if(type(inp) == float):
      inp = round(inp, 10)
    inp = str(inp)
    while(len(inp) < length):
      inp = inp + ' '
    return inp
      
      
  @staticmethod  
  def write_line(inp, length=12):
    line = ''
    for i in inp:
      line = line + neutrons.pad(i, length)
    return line + '\n'
    
  @staticmethod
  def save_data(table_name, cols, data, col_width=25):
    d = {}
    for p in data.keys():
      if(int(p) not in d):
        d[int(p)] = {}
      for n in data[p].keys():
        d[int(p)][int(n)] = data[p][n]
        
  
  #  for key in sorted(wordsFreqDict.keys()) :
    g.results_fh.write('\n')
    g.results_fh.write('====================================================================================================================================' + '\n')
    g.results_fh.write(str(table_name) + '\n')
    g.results_fh.write('====================================================================================================================================' + '\n')
    
    pcols = ['Protons','Nucleons'] + cols
    g.results_fh.write(neutrons.write_line(pcols, col_width))

    for p in sorted(d):
      for n in sorted(d[p].keys()):      
        line = [str(p) + ' ' + isotopes.get_symbol(p),str(n)]
        for i in d[p][n]:
          line.append(i)  
        g.results_fh.write(neutrons.write_line(line, col_width))
    g.results_fh.write('\n')
    g.results_fh.write('\n')
        