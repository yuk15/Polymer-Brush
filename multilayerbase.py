import numpy as np
import math
import pandas as pd

def base_gen(dims, z=-30, spac=0.56123, mol=1, shift=0, out='output'):
    lines = []

    x = dims[0]
    y = dims[1]
    
    xx = np.arange(-x, x, spac)
    yy = np.arange(-y, y, spac*math.sqrt(3))

    mesh = np.meshgrid(xx,yy)

    data=pd.DataFrame()
    
    for i in range(len(mesh[0])):

        # create first mesh in column form
        temp = pd.DataFrame()
        temp['x'] = mesh[0][i]
        temp['y'] = mesh[1][i]
        data = data.append(temp, ignore_index=True)

        # Do the same but accounting for a shift in the {1000} direction
        temp = pd.DataFrame()
        temp['x'] = mesh[0][i] + 0.5*spac
        temp['y'] = mesh[1][i] + 0.5*spac*math.sqrt(3)
        data = data.append(temp, ignore_index=True)
        
    data['z'] = np.zeros(len(data))+z
    data['atom_type']= (np.full((len(data)),3))
    data = data.reset_index(drop=True)
    
    #print data.values

    if out=='df':
        return data
    
    for i in range(len(data)):
        atom = i+1
        charge = 0
        
        x = data['x'][i]
        y = data['y'][i]
        z = data['z'][i]
        atom_type = data['atom_type'][i]
        
        line = '{} {} {} {} {} {} {}\n'.format(atom+shift, mol, 
                                            atom_type, charge, 
                                            x, y, z)
        lines.append(line)

    output = str()
    for line in lines:
        output += line
    return [output, len(lines)]
        

#print base_gen([10, 10], z=10, spac=2**(1/6)*0.5)


def multilayer(dims, thickness=2, spac=0.56123, z=-30, shift=0, mol=1,  out='output'):
    lines = []
    data = pd.DataFrame()
    height = spac*0.5*math.sqrt(3)
    heightz = (spac**2-(spac/2)**2-(height/2)**2)**0.5
           
    for planes in range(thickness):
              
                      
        plane1 = base_gen(dims, z= -30-planes*3*heightz, out='df')    
        plane1['x'] = plane1['x'].values +0.1                    #have to change this for PBC
        plane1['y'] = plane1['y'].values +0.1
        plane1['z'] = plane1['z'].values
        plane1['bulk'] = True
        plane1['boundary']=False
        if planes == 0:
            plane1['surface'] = True
            plane1['atom_type'] = 4
        else:
            plane1['surface'] = False
            plane1['atom_type'] = 3
        data = data.append(plane1, ignore_index=True)
        #data = data.append(plane['y'], ignore_index=True)
        
        plane2 = base_gen(dims, z= -30-heightz-planes*3*heightz, out='df') 
        plane2['x'] = plane2['x'].values +0.1 +spac*0.5 
        plane2['y'] = plane2['y'].values +0.1 +height*0.5
        plane2['z'] = plane2['z'].values
        plane2['atom_type'] = 3
        plane2['bulk'] = True
        plane2['surface'] = False
        plane2['boundary']=False
        data = data.append(plane2, ignore_index=True)
        
        plane3 = base_gen(dims, z= -30-2*heightz-planes*3*heightz, out='df') 
        plane3['x'] = plane3['x'].values +0.1
        plane3['y'] = plane3['y'].values +0.1 +height
        plane3['z'] = plane3['z'].values
        if planes == 1:
            plane3['boundary']=True
        else:
            plane3['boundary']=False
        plane3['bulk'] = True
        plane3['surface'] = False
        plane3['atom_type'] = 3
        data = data.append(plane3, ignore_index=True)
    
    data = data[abs(data['x'])<dims[0]]
    data = data[abs(data['y'])<dims[1]]    
    data = data.reset_index(drop=True)
         
    if out == 'df':
        return data
    
    for i in range(len(data)):
        atom = i+1
        mol = 1
        charge = 0
        
        x=data['x'][i]
        y=data['y'][i]
        z=data['z'][i]
        atom_type = data['atom_type'][i]
        line = '{} {} {} {} {} {} {}\n'.format(atom+shift, mol, 
                                               atom_type, charge,
                                               x, y, z)
        lines.append(line)
    
    output = str()
    for line in lines:
        output+=line
    return [output, len(lines)]

#print multilayer([10,10],thickness=2,spac=0.5)
