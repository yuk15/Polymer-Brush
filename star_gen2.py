# Star Generator Module

import numpy as np
import pandas as pd
import math
import os
import random
from brush import Brush
import multilayer_base_lower
import multilayer_base_upper

spacing = 0.56123
        
def make_brush(item, mol=1):
    brush = Brush(item['trunk']['lam'], mol=mol,
                  starting_position=item['start'],
                  direction=item['direction'], base_id=item.get('base_id'),
                  graft_type=item.get('graft_type', 1))
    for branch in item['branches']:
        brush.create_branch(branch['site'],branch['lam'])
    return brush

class MaxCalculator():

    def __init__(self, item):
        self.item = item
    
    def atoms(self, system):        
        
        if self.item['molecule'] == 'star':
            kap = self.item['kap']
            lam = self.item['lam']
            if self.item['counterions'] == True:
                
                max_atoms = 2*(kap*lam+1)
            elif self.item['counterions'] == False:
                max_atoms = kap*lam + 1
        elif self.item['molecule'] == 'dummy':
            max_atoms = 0
        elif self.item['molecule'] == 'DNA':
            lam = self.item['lam']
            if self.item['counterions'] == True:
                max_atoms = 2*lam
            elif self.item['counterions'] == False:
                max_atoms = lam
        if self.item['molecule'] == 'salt':
            max_atoms = 2*self.item['concentration']
            if self.item['neutralise'] == True:
                max_atoms += abs(neutraliser(system))
        if self.item['molecule']=='brush':
            brush=make_brush(self.item)
            max_atoms = brush.write_atoms()[1]
        if self.item['molecule']=='multilayer_base_lower':
            max_atoms = multilayer_base_lower.multilayer(self.item['dims'],
                                               spac=self.item['spacing'])[1]
        if self.item['molecule']=='multilayer_base_upper':
            max_atoms = multilayer_base_upper.multilayer(self.item['dims'],
                                               spac=self.item['spacing'])[1]
        
        return max_atoms

    def bonds(self):
        
        if self.item['molecule'] == 'star':
            kap = self.item['kap']
            lam = self.item['lam']
            max_bonds = kap*lam
        elif self.item['molecule'] == 'dummy':
            kap = self.item['kap']
            lam = self.item['lam']
            max_bonds = 0
        elif self.item['molecule'] == 'DNA':
            kap = self.item['kap']
            lam = self.item['lam']
            max_bonds = lam-1
        elif self.item['molecule'] == 'salt':
            max_bonds = 0
        elif self.item['molecule'] == 'brush':
            brush = make_brush(self.item)
            max_bonds = brush.write_bonds()[1]
        elif self.item['molecule'] == 'multilayer_base_lower':
            max_bonds=0
        elif self.item['molecule'] == 'multilayer_base_upper':
            max_bonds=0
        
        return max_bonds

    def angles(self):
        
        if self.item['molecule'] == 'star':
            kap = self.item['kap']
            lam = self.item['lam']
            max_angles = kap*(kap-3+2*lam)/2
        elif self.item['molecule'] == 'dummy':
            max_angles = 0
        elif self.item['molecule'] == 'DNA':
            kap = self.item['kap']
            lam = self.item['lam']
            max_angles = lam-2
        elif self.item['molecule'] == 'salt':
            max_angles = 0
        elif self.item['molecule'] == 'brush':
            max_angles=0
        elif self.item['molecule'] == 'multilayer_base_lower':
            max_angles=0
        elif self.item['molecule'] == 'multilayer_base_upper':
            max_angles=0
        
        return max_angles
    

class FileGenerator():
    """
    Input:

    system: list of dictionaries, with each item generally either a star polymer or DNA molecule:

    """

    def __init__(self, box, fstyle='exp', atom_masses=[1.0]):
        self.box = box
        self.fstyle = fstyle
        self.atom_masses = atom_masses

    def write_comments(self, system):

        """

        Returns string that is formatted as the first few lines of a LAMMPS config data file

        """
        comments = str()
        #kap = system[0]['kap']
        #lam = system[0]['lam']
        #first_line = str('Star Polymer with {} arms which are {} beads in length'.format(kap, lam))
        first_line = 'test'
        second_line = str('secondline')
        comments += str('# {}\n'.format(first_line))
        comments += str('# {}\n'.format(second_line))
        comments += '\n'

        return comments
    
    def write_header(self, system):

        """

        Returns string that is formatted as the header of a LAMMPS config data file

        """
        
        MAX_length = float()
        MAX_atoms = int()
        MAX_bonds = int()
        MAX_angles = int()

        spac = spacing

        atom_type_list = range(1, len(self.atom_masses)+1)
        bond_type_list = [1]
        angle_type_list = [1]
        
            

        for item in system:
            HeadGen = MaxCalculator(item)
            #MAX_length += int(item['lam']) * spac
            MAX_atoms += HeadGen.atoms(system)
            MAX_bonds += HeadGen.bonds()
            MAX_angles += HeadGen.angles()

            try:
                atom_type_list.append(item['atom_type'])
            except:
                None
                
            try:
                bond_type_list.append(item['bond_type'])
            except:
                None
                
            try:
                angle_type_list.append(item['angle_type'])
            except:
                None
        
        n_atoms = MAX_atoms
        m_bonds = MAX_bonds
        l_angles = MAX_angles

        a_atom_types = max(atom_type_list)
        b_bond_types = max(bond_type_list)
        c_angle_types = max(angle_type_list)

        xlo = -self.box
        xhi = self.box
        ylo = -self.box*math.sqrt(3)*0.5
        yhi = self.box*math.sqrt(3)*0.5
        zlo = -self.box - 30
        zhi = self.box + 30
        
        header = str()

        next_line = str()
        next_line += str("{} atoms\n".format(n_atoms))
        next_line += str("{} bonds\n".format(m_bonds))
        next_line += str("{} angles\n".format(l_angles))
        next_line += "\n"
        next_line += str("{} atom types\n".format(a_atom_types))
        next_line += str("{} bond types\n".format(b_bond_types))
        next_line += str("{} angle types\n".format(c_angle_types))
        next_line += "\n"
        next_line += str("{} {} xlo xhi\n".format(xlo, xhi))
        next_line += str("{} {} ylo yhi\n".format(ylo, yhi))
        next_line += str("{} {} zlo zhi\n".format(zlo, zhi))
        header += next_line

        return header

    def write_masses(self):

        
        masses = str()

        for i in range(len(self.atom_masses)):
            next_line = str()
            next_line += str("{} {}\n".format(i+1, self.atom_masses[i]))
            masses += next_line

        return masses
        
    def write_atoms(self, system, system_index):

        """

        Returns string that is formatted as the atoms section of a LAMMPS input data file

        The format is:

        atom-ID mol-ID atom-type q x y z

        """

        atom_list = str()
        item = system[system_index]
        box = self.box
       
        # write function to check that item obeys rules

        CUMU_atoms = int()        
        for i in range(system_index):
            CUMU_atoms += MaxCalculator(system[i]).atoms(system)
        
        atom_ID_shift = CUMU_atoms
        spac = spacing
        
        molecule_id = system_index + 1
                  
        if item['molecule'] == 'brush':
            brush = make_brush(item, mol=molecule_id)
            atom_list = brush.write_atoms(shift=atom_ID_shift)[0]
        
        if item['molecule'] == 'multilayer_base_lower':
            atom_list = multilayer_base_lower.multilayer(item['dims'], z=item['plane'],
                                      spac=item['spacing'], 
                                      mol=molecule_id, 
                                      shift=atom_ID_shift)[0]
            
        if item['molecule'] == 'multilayer_base_upper':
            atom_list = multilayer_base_upper.multilayer(item['dims'], z=item['plane'],
                                      spac=item['spacing'], 
                                      mol=molecule_id, 
                                      shift=atom_ID_shift)[0]
            
                       
        return atom_list
        
    def write_bonds(self, system, system_index):

        """

        Returns string that is formatted as the bonds section of a LAMMPS input data file

        The format is:

        bond-ID bond-type atom-ID_1 atom-ID_2

        """

        bond_list = str()
        item = system[system_index]

        CUMU_atoms = int()
        CUMU_bonds = int()
        for i in range(system_index):
            CUMU_atoms += MaxCalculator(system[i]).atoms(system)
            CUMU_bonds += MaxCalculator(system[i]).bonds()
            
        atom_ID_shift = CUMU_atoms
        bond_ID_shift = CUMU_bonds
        try:
            lam = item['lam'] 
        except:
            None


        if item['molecule'] == 'brush':
            brush = make_brush(item)
            bond_list = brush.write_bonds(bond_shift=bond_ID_shift,
                                          atom_shift=atom_ID_shift)[0]
        
        return bond_list

    def create_filename(self, system):

        """

        Returns string that is the filename for the system

        """
        
        if self.fstyle == 'exp':
            filename = 'exp.dat'
        elif self.fstyle == 'al':
            star = system[0]
            salt = system[1]
            filename = str('al_'+star['kap']+'_'+star['lam']+'_'+salt['conc'])
        elif self.fstyle == 'ssr':
            f_kap = str(system[0]['kap'])
            f_lam = str(system[0]['lam'])
            f_conc = str(system[2]['concentration'])
            filename = 'ssr_'+f_kap+'_'+f_lam+'_'+f_conc+'.dat'
        elif self.fstyle == 'svl':
            f_kap = str(system[0]['kap'])
            f_lam = str(system[0]['lam'])
            filename = 'svl_{}_{}.dat'.format(f_kap, f_lam)
        elif self.fstyle == 'ca':
            f_kap = str(system[0]['kap'])
            f_lam = str(system[0]['lam'])
            f_ang = str(system[0]['central'])
            filename = 'ca_{}_{}_{}.dat'.format(f_kap, f_lam, f_ang)
        elif self.fstyle == 'es':
            f_kap = str(system[0]['kap'])
            f_lam = str(system[0]['lam'])
            filename = 'es_{}_{}.dat'.format(f_kap, f_lam)
        else:
            if len(system) == 1:
                item = system[0]
                filename = item['molecule']+str(item['kap'])+'_'+str(item['lam'])+'.dat'
            elif len(system) == 3:
                f_kap = str(system[0]['kap'])
                f_lam = str(system[0]['lam'])
                f_conc = str(system[2]['concentration'])
                filename = '{}_{}_{}_{}.dat'.format(self.fstyle,
                                                    f_kap, f_lam,
                                                    f_conc)
            elif len(system) == 2 and system[0]['molecule'] == 'star':
                f_kap = str(system[0]['kap'])
                f_lam = str(system[0]['lam'])
                f_conc = str(system[1]['concentration'])
                filename = 'sl_'+f_kap+'_'+f_lam+'_'+f_conc+'.dat'
            else:
                filename = str('exp.dat')

        return filename

    def write_system_to_file(self, system, angles=True):

                   
        with open(self.create_filename(system), 'w') as f:
            f.write(self.write_comments(system))
            f.write(self.write_header(system))
            f.write('\nMasses\n\n')
            f.write(self.write_masses())
            f.write('\nAtoms\n\n')
            
            for i in range(len(system)):
                if i == 0:
                    f.write(self.write_atoms(system, i))
                else:
                    f.write(self.write_atoms(system, i))   
            f.write('\nBonds\n\n')
            
            for i in range(len(system)):
                if i == 0:
                    f.write(self.write_bonds(system, i))
                else:
                    f.write(self.write_bonds(system, i))
            if angles == True:
                f.write('\nAngles\n\n')
                for i in range(len(system)):
                    if i == 0:
                        f.write(self.write_angles(system, i))
                    else:
                        f.write(self.write_angles(system, i))
        print "Writing complete for {}".format(self.create_filename(system))
        return
