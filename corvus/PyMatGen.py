from corvus.structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import re
import numpy as np
# Debug: FDV
import pprint

from mp_api.client import MPRester
from pymatgen.io.cif import CifParser,CifWriter
from pymatgen.core import structure
from pymatgen.io.vasp import Vasprun, Outcar, Xdatcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer

pp_debug = pprint.PrettyPrinter(indent=4)


# Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))
subs = lambda L:[{L[j] for j in range(len(L)) if 1<<j&k} for k in range(1,1<<len(L))]
for s in subs(['cell_vectors', 'cell_struct_xyz_red', 'cell_scaling_iso', 'cell_scaling_abc', 'number_density']):
    key = strlistkey(s)
    autodesc = 'Get ' + ', '.join(s) + ' using cif2cell'
    cost = 10
    implemented[key] = {'type':'Exchange','out':list(s),'req':['mp.structure'],
                        'desc':autodesc,'cost':cost}

#implemented['cell_structure'] = {'type':'Exchange','out':['cell_structure'],'cost':0,
#                        'req':['cell_vectors','cell_struct_xyz_red','cell_scaling_iso','cell_scaling_abc'],'desc':'Calculate cell structure from cif file using cif2cell.'}


#implemented['cell_structure'] = {'type':'Exchange','out':['cell_structure','cell_struc_xyz','cell_scaling_abc','cell_scaling_iso'],'cost':0,
#                        'req':['cif_input'],'desc':'Calculate cell structure from cif file using cif2cell.'}

implemented['mp.structure'] = {'type':'Exchange','out':['mp.structure'],'cost':0,
                        'req':['mp_id|cif_input|vasp_xml|vasp_xdatcar'],'desc':'Get pymatgen structure from.'}
implemented['cluster_array'] = {'type':'Exchange','out':['cluster_array'],'cost':0,
                        'req':['mp.structure'],'desc':'Calculate cluster from cif using pymatgen.'}
implemented['supercell'] = {'type':'Exchange','out':['supercell'],'cost':0,
                        'req':['mp.structure'],'desc':'Calculate supercell from cif input using pymatgen.'}



class PyMatGen(Handler):
    def __str__(self):
        return 'PyMatGen Handler'

    @staticmethod
    def canProduce(output):
        if isinstance(output, list) and output and isinstance(output[0], str):
            if len(output) > 1:
                return strlistkey(output) in implemented
            else:
                out = []
                canProduce = False
                for o in output[0].split('|'):
                    if o in implemented:
                       out.append(o)
                       canProduce = True
                       break
                if canProduce: output = out
                #print(output)
                return canProduce 
        elif isinstance(output, str):
            canProduce = False
            for o in output.split('|'):
                if o in implemented:
                   out = o
                   canProduce = True
                   break
            if canProduce: output = out
            #print(output)
            return canProduce 
            # JK - Add support for logical or in requirements, i.e., 'req':['cif_input', 'mp_id|mp_query'
            # Split output by | and check if any are implemented.
        else:
            raise TypeError('Output should be token or list of tokens')

    @staticmethod
    def requiredInputFor(output):
        #print("output:",output)
        if isinstance(output, list) and output and isinstance(output[0], str):
            unresolved = {o for o in output if not PyMatGen.canProduce(o)}
            canProduce = (o for o in output if PyMatGen.canProduce(o))
            additionalInput = (set(implemented[o]['req']) for o in canProduce)
            return list(set.union(unresolved,*additionalInput))
        elif isinstance(output, str):
            if output in implemented:
                return implemented[output]['req']
            else:
                return [output]
        else:
            raise TypeError('Output should be token or list of tokens')

    @staticmethod
    def cost(output):
        if isinstance(output, list) and output and isinstance(output[0], str):
            key = strlistkey(output)
        elif isinstance(output, str):
            key = output
        else:
            raise TypeError('Output should be token or list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using PyMatGen')
        return implemented[key]['cost']

    @staticmethod
    def sequenceFor(output,inp=None):
        if isinstance(output, list) and output and isinstance(output[0], str):
            key = strlistkey(output)
        elif isinstance(output, str):
            key = output
        else:
            raise TypeError('Output should be token of list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using PyMatGen')
        f = lambda subkey : implemented[key][subkey]
        if f('type') == 'Exchange':
            return Exchange(PyMatGen, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_PYMATGEN'
        xcDir = os.path.join(config['cwd'], subdir)
        # Make new output directory if if doesn't exist
        if not os.path.exists(xcDir):
            os.mkdir(xcDir)
        # Store current Exchange directory in configuration
        config['xcDir'] = xcDir

    #@staticmethod
    #def setDefaults(input,target):

    @staticmethod
    def run(config, input, output):
        

        # Loop over targets in output.
        #print("Inside PyMatGen")
        site_tol = 1.0e-6
        if 'cluster_array' in output:

            #print("Number of absorbers:", len(cluster_array))
            #parser = CifParser(input.get("cif_input")[0][0],site_tolerance=site_tol)
            #structure = parser.get_structures()[0]
            structure = input['mp.structure']
            #symprec=input['pymatgen.symprec'][0][0]
            symprec=0.1
            #angle_tolerance=input['pymatgen.angle_tolerance'][0][0]
            #sg_anal = SpacegroupAnalyzer(structure,symprec=symprec, angle_tolerance=angle_tolerance) 
            try:
                is_magnetic = CollinearMagneticStructureAnalyzer(structure).is_magnetic
                #symprec = 1.0e-12
            except:
                is_magnetic = False
            
            #print(input)
            #print('HELP1:', structure.sites[0].label)
            #print('absorbing_atom_by_label' in input)
            if 'absorbing_atom_by_label' in input:
                # Set equivelent indices arrays to single sites.
                structure.equivalent_indices = [[i] for i,x in enumerate(structure.sites)]
            else:
                sg_anal = SpacegroupAnalyzer(structure,symprec=symprec) 
                #conventional_structure = sg_anal.get_conventional_standard_structure()
                #sg_anal2 = SpacegroupAnalyzer(conventional_structure,symprec=symprec) 
                structure = sg_anal.get_symmetrized_structure()

            #print('HELP2:', structure.sites[0].label)
            # Get local denSpacegroupAnalyzersities for possible later use.
            cluster_radius = input['clusterradius'][0][0]
            getLocalDensity(structure,cluster_radius)
            with open('densities.dat', 'w') as fden:
               for site in structure.sites:
                  print(site.label,site.coords.tolist()[0],site.coords.tolist()[1],site.coords.tolist()[2],site.properties["local_density"],file=fden)
            

            #print(input['absorbing_atom_type'])
            if "absorbing_atom_type" in input: # Will set up calculation of all unique absorbers in unit cell.
                absorber_types=[input["absorbing_atom_type"][0][0]]
                absorber_spec = 1
            elif "absorbing_atom_by_label" in input: # Will set up calculation for atoms with label that starts with string.
                absorber_labels = [input["absorbing_atom_by_label"][0][0]]
                absorber_spec = 2
            else:
                # Use all elements in crystal
                absorber_spec = 1
                absorber_types=structure.symbol_set

            #print("Absorber types:", absorber_types)
            disordered_structure = False
            ipot = 1
            n_disord = 1
            cluster_array = []
            abstype = []
            # Get index of all absorbers
            abs_inds = []
            #print('HELP3:', structure.sites[0].label)
            for inds in structure.equivalent_indices:
                xnat = len(inds)
                #print('inds', inds)
                if absorber_spec == 1:
                   for key in structure.sites[inds[0]].species.as_dict().keys():
                        for absorber in absorber_types:
                            if absorber == re.sub('[^a-zA-Z]','',key): abs_inds = abs_inds + [inds[0]]
                            if absorber == re.sub('[^a-zA-Z]','',key): 
                                abstype.append(key)
                elif absorber_spec == 2:
                    for abs_label in absorber_labels:
                        #print(structure.sites[0].label)
                        #print(dir(structure.sites[inds[0]]))
                        #exit()
                        #print(structure.sites[inds[0]].label,abs_label,inds)
                        if structure.sites[inds[0]].label.startswith(abs_label): 
                            abs_inds = abs_inds + inds
                            for k in structure.sites[inds[0]].species.as_dict().keys():
                                abstype.append(k)
                
                    #elif absorber_spec == 2:
                    #    for absorber in absorber_types_by_label:
                    #        absorber_inds = []
                    #        for site in structure.sites
                for ind in inds:
                    mag = 'magmom' in structure.sites[ind].properties
                    #print('checking mag', mag,structure.sites[ind].properties['magmom'])
                    structure.sites[ind].properties['itype'] = []
                    structure.sites[ind].properties['xnat']  = []
                    magmom = []
                    i_spec=0
                    for occ in  structure.sites[ind].species.as_dict().values():
                        #print(structure.sites[ind])
                        #print('occ', occ)
                        structure.sites[ind].properties['itype'] += [ipot+i_spec]
                        structure.sites[ind].properties['xnat'] += [xnat*occ]
                        if not mag: 
                            magmom += [0.0]
                        else:
                            magmom += [structure.sites[ind].properties['magmom']]
                        #print(structure.sites[ind].properties['itype'])
                        if 'magmom_by_label' in input:
                           for elem in input['magmom_by_label']:
                              if structure.sites[ind].label == elem[0]:
                                 magmom = [elem[1]]
                        if occ != 1.0:
                            disordered_structure = True
                            n_disord = input['numberofconfigurations'][0][0]
                        i_spec += 1

                    structure.sites[ind].properties['magmom'] = magmom
                    #print(structure.sites[ind],magmom)
                    #print(structure.sites[ind].properties['magmom'])
                       
                #print(ipot,structure.sites[ind].properties['itype'])
                ipot = ipot + i_spec
                #print('ipot',ipot)
           
            absorber_types = set(abstype)
            #print(abs_inds)
            #exit()
            #print(absorber_types)
            #exit() 
            #print(dir(structure.sites[0].species))
            #print(structure.sites[0].species.as_dict())
            i_disord = 1
            #print('n_disord:', n_disord)
            while i_disord <= n_disord:
                for inds in structure.equivalent_indices:
                    weight = 1 #len(inds)
                    #print(structure.sites[inds[0]].species_string)
                    species = {}
                    #for abs_symbol in absorber_types:
                    iabs=0
                    for abs_ind in abs_inds:
                        abs_symbol = abstype[iabs]
                        #print('abs_symbol',abs_symbol)
                        iabs=iabs+1
                        # remove alphabetical characters from keys in the dictionary.
                        #for key,value in structure.sites[inds[0]].species.as_dict().items():
                        #    species[re.sub('[^a-zA-Z]','',key)] = value
                        #print(structure.sites[inds[0]].species.as_dict())
                        #if any(abs_symbol in spec for spec in structure.sites[inds[0]].species.as_dict().keys()):
                        #print(abs_ind,inds)
                        if abs_ind in inds:
                        #if abs_symbol in structure.sites[inds[0]].species.as_dict():
                            # Set absorber symbol
                            # Make a cluster around this absorber
                            site_cluster = structure.get_neighbors(structure.sites[inds[0]],cluster_radius)            
                            site_cluster = [structure.sites[inds[0]]] + site_cluster
                            cluster = []
                         
                            iclust = 0
                            for site in site_cluster:
                                # Loop over all species at this site
                                nspec = len(site.species.as_dict().keys())
                                # Get total occupancy.
                                tot_occ = 0.0
                                #print(site)
                                #print(site.properties)
                                for occ in site.species.as_dict().values():
                                    # Always put the absorbing atom in the cluster
                                    if iclust == 0:
                                        weight = weight*occ
                                   
                                        tot_occ = tot_occ + occ
                                        if tot_occ > 1.0:
                                            print('Total occupation in cif file for one site is greater than 1')
                              
                                occ_sum = 0.0
                                rnd = np.random.uniform()
                                i_spec = 0
                                for spec_str,occ in site.species.as_dict().items():
                                    #print(iclust,rnd,spec_str,occ_sum)
                                    #print(site)
                                    # Always put the absorbing atom in the cluster
                                    #print(spec_str, absorber_types)
                                    #if spec_str in absorber_types:
                                    #    abs_symbol = spec_str 
                                    
                                    #if absorber_spec == 2: abs_symbol = spec_str 
                                    #print('spec_str,abs_symbol', spec_str,abs_symbol,iclust)
                                    if(iclust == 0):
                                        if spec_str == abs_symbol:
                                            #print('Magnetic moment of absorber:', site.properties.get('magmom'))
                                            #exit()
                                            #print('inside: spec_str,abs_symbol', spec_str,abs_symbol,iclust)
                                            cluster = cluster + [[abs_symbol] + site.coords.tolist()     + 
                                                             [ site.properties['itype'][i_spec] ]    + 
                                                             [site.properties['xnat'][i_spec]]       + 
                                                             [site.properties.get('magmom')[i_spec]] + 
                                                             [site.properties.get('local_density')]  +
                                                             [site.label]]
                                            label = site.label
                                            break
                                        else:
                                            continue
                                    elif occ_sum < rnd <= occ + occ_sum:
                                        #print(site.properties)
                                        cluster = cluster + [[re.sub('[^a-zA-Z]','',spec_str)]       + 
                                                             site.coords.tolist()                    + 
                                                             [ site.properties['itype'][i_spec] ]    + 
                                                             [site.properties['xnat'][i_spec]]       + 
                                                             [site.properties.get('magmom')[i_spec]] + 
                                                             [site.properties.get('local_density')]  +
                                                             [site.label]]
                                        #print(site.properties)
                                        #print(cluster[-1])
                                        #sys.stdin.readline()
                                        
                                    occ_sum = occ_sum + occ
                                    i_spec += 1
                            
                                iclust = iclust + 1
                                
                            # cluster_array is a list of tuples, each with absorbing atom, associated cluster, and 
                            # weighting (stoichiometry for example).
                            cluster_array = cluster_array + [(1,weight,cluster,label)]
                i_disord = i_disord + 1        

            if(len(cluster_array) == 0):
               print("No absorbing atoms found.")
               exit() 
            #print(cluster_array)
            output['cluster_array'] = cluster_array
            #print("CLUSTER ARRAY")
            #print(cluster_array[0][0:1])
            #for line in cluster_array[0][2]:
            #   print(line[4])
        elif 'vasp_md_to_feff' in output:
             print("hello")
        elif set(output.keys()).issubset(set(['supercell', 'cell_vectors', 'cell_struct_xyz_red', 'cell_scaling_iso', 'cell_scaling_abc', 'number_density'])):
        #elif 'supercell' in output:
            structure = input['mp.structure']
            sg_anal = SpacegroupAnalyzer(structure) 
            structure = sg_anal.get_symmetrized_structure()
            if "absorbing_atom_type" in input: # Will set up calculation of all unique absorbers in unit cell.
                absorber_types=[input["absorbing_atom_type"][0][0]]
            else:
                # Use all elements in crystal
                absorber_types=structure.symbol_set

            # Add a type label to the properties that is different for each unique site in the crystal.
            absorbers = []
            itype = 0
            itypes = []
            for inds in structure.equivalent_indices:
                for ind in inds:
                    structure.sites[ind].properties["type"] = itype

                itypes = itypes + [itype]
                itype += 1

            scaling_vector=input['supercell.dimensions'][0]
            structure.make_supercell(scaling_vector)

            # Get equivalent indices of the supercell.
            eq_inds = []
            for itype in itypes:
                inds = []
                for i,s in enumerate(structure.sites):
                    if s.properties['type'] == itype:
                        inds = inds + [i]

                eq_inds = eq_inds + [inds]
                        
            #print(structure.sites)
            #print(eq_inds)
            #print(itypes)
            for abs_symbol in absorber_types:
                for inds in eq_inds:
                    if structure.sites[inds[0]].species.elements[0].value == abs_symbol:
                        absorbers = absorbers + [inds[0]+1]
            #print(absorbers)
            #exit()
            output['absorbing_atom']=[absorbers]
            output['cell_scaling_abc'] = [structure.lattice.lengths]
            output['cell_angles_abg'] = [structure.lattice.angles]
            output['number_of_atoms'] = [[structure.num_sites]]
            output['number_of_species'] = [[structure.ntypesp]]

            spcs=list(set([comp.elements[0] for comp in structure.species_and_occu]))
            species=[]
            for s in spcs:
                species = species + [[s.number,s.value]]
            output['species'] = species

            red_coords = []
            ind = 0
            for s in structure.sites:
                # For now take only the first element if there is disorder. We should augment this later.
                red_coords = red_coords + [[s.species.elements[0].value] + list(s.frac_coords)]
                #red_coords = red_coords + [[s.species] + s.frac_coords]

            #print(red_coords[0])
            output['cell_struc_xyz_red'] = red_coords
            output['supercell']=[scaling_vector]

        elif 'mp.structure' in output: 
            # Note: Right now, we are using mp_id to produce a cif file, while we should 
            # just go directly from the mp structure. However, we need to be able to have
            # or operators available in requiredInput rather than just and operators.
            # Need to add inputs: mp_id, mp_api_key
            if 'mp_id' in input:
                mpr = MPRester(input['mp_apikey'][0][0])
                output['mp.structure'] = mpr.get_structure_by_material_id(input["mp_id"][0][0])
                output['mp.structure'].to(filename=input["mp_id"][0][0] + ".cif")
            elif 'cif_input' in input:
                parser = CifParser(input.get("cif_input"))
                # Only take first structure for now.
                output['mp.structure'] = parser.parse_structures()[0]
            elif 'vasp_xml' in input:
                vr = Vasprun(input['vasp_xml'][0][0])
                struct = vr.structures[input['vasp_snapshot'][0][0]]
                #print(struct)
                if 'vasp_outcar' in input: 
                   oc = Outcar(input['vasp_outcar'][0][0])
                   if len(oc.magnetization) == len(struct.sites): 
                      magmom = [m['tot'] for m in oc.magnetization]
                      struct.add_site_property("magmom",magmom)
                #print(struct)
                output['mp.structure'] = struct
            elif 'vasp_xdatcar' in input:
                xc = Xdatcar(input['vasp_xdatcar'][0][0])
                struct = xc.structures[input['vasp_snapshot'][0][0]]
                if 'vasp_outcar' in input: 
                   oc = Outcar(input['vasp_outcar'][0][0])
                   if len(oc.magnetization) == len(struct.sites): 
                      magmom = [m['tot'] for m in oc.magnetization]
                      struct.add_site_property("magmom",magmom)
               
                output['mp.structure'] = struct
                 
    
            
             
            




    @staticmethod
    def cleanup(config):
        pass


def getLocalDensity(structure,cluster_radius):
   # Define a few simple local structure parameters to define potentials in FEFF. Maybe this should be moved to the feff
   # module along with cluster_array etc.
   sigma = 2.3 # For now use 2.3 angstrom width.
   mindens=1.0e10
   maxdens=0.0
   densities = []
   for site in structure.sites:
       # Create a cluster around the site
       site_cluster = structure.get_neighbors(site,sigma*5)
       coords0 = np.array(site.coords.tolist())
       density = 0.0
       for neighbor in site_cluster:
         coords = np.array(neighbor.coords.tolist()) - coords0
         distance = np.sqrt(np.inner(coords,coords))
          
         #print(neighbor)
         #print(dir(neighbor.species))
         #print(neighbor.species.elements)
         # For now, average over partially occupied sites for the density. May need to change this later.
         i=0
         #print('Elements:')
         #print(neighbor.species.elements)
         for occ in neighbor.species.as_dict().values():   
            #print(occ)
            #print(neighbor.species.elements[i].Z)
            #print(distance)
            density = density + localDensityFunction(distance,sigma)*neighbor.species.elements[i].Z
            #print('Z: ', neighbor.species.elements[i].Z)
            #print(density,distance)
            i=i+1
       
       # Add density to site properties
       densities = densities + [density]
       mindens=min(mindens,density)
       maxdens=max(maxdens,density)
   for i,density in enumerate(densities):
       if maxdens == mindens:
          structure.sites[i].properties["local_density"] = 0.0
       else:
          structure.sites[i].properties["local_density"] = (density - mindens)/(maxdens-mindens)
          #print('Densities: ', density, mindens, maxdens, (density-mindens)/(maxdens-mindens))

   return (mindens,maxdens)
       
def localDensityFunction(r,sigma):
   # Define a local function with a volume integral = 1.
   vol = 4.0*np.pi/3.0*sigma**3
   return r/(6*sigma)*np.exp(-(r/sigma)**2/2.0)/vol

 
       
# def getSymmetryUniqueAtoms(structure):
#     ''' input  - structure of cell 
#         output - modified structure with site properties labeled by number and type of nearest neighbors.   
#     '''
#     # Loop over equivalent sets of atoms
#     for inds in structure.equivalent_indices:
#         # Get number of atoms of given type
#         xnat = len(inds)
                
#         for ind in inds:
#             structure.sites[ind].properties['itype'] = []
#             structure.sites[ind].properties['xnat']  = []

#             # get number and type of nearest neighbors. 
#             #print(structure.sites[ind].species.as_dict().values())
#             i_spec=0
#             for occ in  structure.sites[ind].species.as_dict().values():
#                 structure.sites[ind].properties['itype'] += [ipot+i_spec]
#                 structure.sites[ind].properties['xnat'] += [xnat*occ]
                
#                 if occ != 1.0:
#                     disordered_structure = True
#                     n_disord = input['numberofconfigurations'][0][0]
#                     i_spec += 1
                        
#         ipot += 1
        
#         return 1
