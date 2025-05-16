from corvus.structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import re
import numpy as np
# Debug: FDV
import pprint

from pymatgen.io.cif import CifParser
from pymatgen.core import structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

pp_debug = pprint.PrettyPrinter(indent=4)


# Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))
subs = lambda L:[{L[j] for j in range(len(L)) if 1<<j&k} for k in range(1,1<<len(L))]
for s in subs(['cell_vectors', 'cell_struct_xyz_red', 'cell_scaling_iso', 'cell_scaling_abc', 'number_density']):
    key = strlistkey(s)
    autodesc = 'Get ' + ', '.join(s) + ' using cif2cell'
    cost = 10
    implemented[key] = {'type':'Exchange','out':list(s),'req':['cif_input'],
                        'desc':autodesc,'cost':cost}

#implemented['cell_structure'] = {'type':'Exchange','out':['cell_structure'],'cost':0,
#                        'req':['cell_vectors','cell_struct_xyz_red','cell_scaling_iso','cell_scaling_abc'],'desc':'Calculate cell structure from cif file using cif2cell.'}


#implemented['cell_structure'] = {'type':'Exchange','out':['cell_structure','cell_struc_xyz','cell_scaling_abc','cell_scaling_iso'],'cost':0,
#                        'req':['cif_input'],'desc':'Calculate cell structure from cif file using cif2cell.'}


implemented['cluster_array'] = {'type':'Exchange','out':['cluster_array'],'cost':0,
                        'req':['cif_input'],'desc':'Calculate cluster from cif using pymatgen.'}
implemented['supercell'] = {'type':'Exchange','out':['supercell'],'cost':0,
                        'req':['cif_input'],'desc':'Calculate supercell from cif input using pymatgen.'}



class PyMatGen(Handler):
    def __str__(self):
        return 'PyMatGen Handler'

    @staticmethod
    def canProduce(output):
        if isinstance(output, list) and output and isinstance(output[0], str):
            return strlistkey(output) in implemented
        elif isinstance(output, str):
            return output in implemented
        else:
            raise TypeError('Output should be token or list of tokens')

    @staticmethod
    def requiredInputFor(output):
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
            raise LookupError('Corvus cannot currently produce ' + key + ' using FEFF')
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
            raise LookupError('Corvus cannot currently produce ' + key + ' using FEFF')
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
        if 'cluster_array' in output:
            site_tol = 1.0e-6
            parser = CifParser(input.get("cif_input")[0][0],site_tolerance=site_tol)
            structure = parser.get_structures()[0]
            #symprec=input['pymatgen.symprec'][0][0]
            symprec=1.0e-6
            print('Structure read from cif file.')
            #angle_tolerance=input['pymatgen.angle_tolerance'][0][0]
            #sg_anal = SpacegroupAnalyzer(structure,symprec=symprec, angle_tolerance=angle_tolerance) 
            sg_anal = SpacegroupAnalyzer(structure,symprec=symprec) 
            conventional_structure = sg_anal.get_conventional_standard_structure()
            sg_anal2 = SpacegroupAnalyzer(conventional_structure,symprec=symprec) 
            structure = sg_anal2.get_symmetrized_structure()
            #print(input['absorbing_atom_type'])
            if "absorbing_atom_type" in input: # Will set up calculation of all unique absorbers in unit cell.
                absorber_types=[input["absorbing_atom_type"][0][0]]
            else:
                # Use all elements in crystal
                absorber_types=structure.symbol_set

            #print("Absorber types:", absorber_types)
            disordered_structure = False
            ipot = 1
            n_disord = 1
            cluster_array = []
            abstype = []
            for inds in structure.equivalent_indices:
                xnat = len(inds)
                
                # Here, redefine the absorber_types using regex. Take all species
                # That start with the chemical symbols given in the input.
                for key in structure.sites[inds[0]].species.as_dict().keys():
                    for absorber in absorber_types:
                        if absorber == re.sub('[^a-zA-Z]','',key): abstype.append(key)
                for ind in inds:
                    structure.sites[ind].properties['itype'] = []
                    structure.sites[ind].properties['xnat']  = []
                    
                    i_spec=0
                    for occ in  structure.sites[ind].species.as_dict().values():
                        #print(structure.sites[ind])
                        #print('occ', occ)
                        structure.sites[ind].properties['itype'] += [ipot+i_spec]
                        structure.sites[ind].properties['xnat'] += [xnat*occ]
                        #print(structure.sites[ind].properties['itype'])

                        if occ != 1.0:
                            disordered_structure = True
                            n_disord = input['numberofconfigurations'][0][0]
                        i_spec += 1
                    #print(structure.sites[ind],structure.sites[ind].properties)
                       
                #print(ipot,structure.sites[ind].properties['itype'])
                ipot = ipot + i_spec
                #print('ipot',ipot)

            absorber_types = set(abstype)
            #exit() 
            #print(dir(structure.sites[0].species))
            #print(structure.sites[0].species.as_dict())
            i_disord = 1
            cluster_radius = input['clusterradius'][0][0]
            while i_disord <= n_disord:
                for inds in structure.equivalent_indices:
                    weight = len(inds)
                    #print(structure.sites[inds[0]].species_string)
                    species = {}
                    for abs_symbol in absorber_types:
                        # remove alphabetical characters from keys in the dictionary.
                        #for key,value in structure.sites[inds[0]].species.as_dict().items():
                        #    species[re.sub('[^a-zA-Z]','',key)] = value

                        if any(abs_symbol in spec for spec in structure.sites[inds[0]].species.as_dict().keys()):
                        #if abs_symbol in structure.sites[inds[0]].species.as_dict():
                            # loop over indices in equivalent indices.
                            for ind in inds:
                                # Make a cluster around this absorber
                                site_cluster = structure.get_neighbors(structure.sites[ind],cluster_radius)            
                                site_cluster = [structure.sites[ind]] + site_cluster
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
                                        if spec_str == abs_symbol and iclust == 0:
                                            cluster = cluster + [[spec_str] + site.coords.tolist() + [ site.properties['itype'] ] + [site.properties['xnat']]]
                                            break
                                        elif occ_sum < rnd <= occ + occ_sum:
                                            cluster = cluster + [[spec_str] + site.coords.tolist() + [ site.properties['itype'][i_spec] ] + [site.properties['xnat'][i_spec]]]
                                            #print(site.properties)
                                            #print(cluster[-1])
                                            #sys.stdin.readline()
                                            
                                        occ_sum = occ_sum + occ
                                        i_spec += 1
                                
                                    iclust = iclust + 1
                                    
                                # cluster_array is a list of tuples, each with absorbing atom, associated cluster, and 
                                # weighting (stoichiometry for example).
                                cluster_array = cluster_array + [(1,weight,cluster)]
                i_disord = i_disord + 1        
            #print("Number of absorbers:", len(cluster_array))
            if(len(cluster_array) == 0):
               print("No absorbing atoms of types", absorber_types)
               print("found.")
               exit() 
            output['cluster_array'] = cluster_array

        elif set(output.keys()).issubset(set(['supercell', 'cell_vectors', 'cell_struct_xyz_red', 'cell_scaling_iso', 'cell_scaling_abc', 'number_density'])):
        #elif 'supercell' in output:
            parser = CifParser(input.get("cif_input")[0][0])
            structure = parser.get_structures()[0]
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




    @staticmethod
    def cleanup(config):
        pass



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
