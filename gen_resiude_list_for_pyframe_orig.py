#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 10:33:40 2018

@author: cagus
"""


import collections
import numpy
import string

#COLUMNS        DATA  TYPE    FIELD        DEFINITION
#-------------------------------------------------------------------------------------
# 1 -  6        Record name   "ATOM  "
# 7 - 11        Integer       serial       Atom  serial number.
#13 - 16        Atom          name         Atom name.
#17             Character     altLoc       Alternate location indicator.
#18 - 20        Residue name  resClme      Residue name.
#22             Character     chainID      Chain identifier.
#23 - 26        Integer       resSeq       Residue sequence number.
#27             AChar         iCode        Code for insertion of residues.
#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#55 - 60        Real(6.2)     occupancy    Occupancy.
#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#77 - 78        LString(2)    element      Element symbol, right-justified.
#79 - 80        LString(2)    charge       Charge  on the atom.


import os


path_to_script=os.getcwd()
os.chdir(str(path_to_script))

print('Converting files:')
for files in os.listdir(os.getcwd()):
    if not files.startswith('new'):
        if files.startswith('protein_BTV') and files.endswith('.pdb'):
            os.chdir(str(path_to_script))
            print(files)
            #print("--- %s seconds with open file ---" % (time.time() - start_time))
            with open(files) as f: 
                print
                form_inp_data=[]
                PDB_output=[]
                counter=0
                for line in f:
                    counter=counter+1
                    if line[0:4]=='ATOM': #remove all other lines from PDB file
                        str_split_data=line.split()
                        orig_atom_name=str_split_data[2] #split the atom name argument in order to sort on first letter                        
                        if orig_atom_name[0] =='1' or orig_atom_name[0] =='2' or orig_atom_name[0] == '3'or orig_atom_name[0] == '4'or orig_atom_name[0] == '5' or orig_atom_name[0] == '6' or orig_atom_name[0] == '7' or orig_atom_name[0] == '8' or orig_atom_name[0] == '9'or orig_atom_name[0] == '0':#numb_list[numb]:
                            new_atom_name=orig_atom_name[1]
                        elif orig_atom_name[0:2] == 'Cl' or orig_atom_name[0:2] == 'Cl' or orig_atom_name[0:2] == 'CL': # upgrade to include atom names with more than one char with islower() or smt
                            new_atom_name=orig_atom_name[0:2]
                            new_atom_name='Cl'
                        else:
                            new_atom_name=orig_atom_name[0]
                        if len(str_split_data) ==11:
                            residue_type=str_split_data[3]
                            res_nr=str_split_data[4]
                            x_coord= str_split_data[5]
                            y_coord= str_split_data[6]
                            z_coord= str_split_data[7]
                            if residue_type=='SOL':
                                residue_type='WAT'
                        elif len(str_split_data) ==10:
                            if str_split_data[4].startswith('Z'):
                                residue_type=str_split_data[3]
                                res_nr_tmp=str_split_data[4]
                                res_nr=res_nr_tmp[1:]
                                x_coord= str_split_data[5]
                                y_coord= str_split_data[6]
                                z_coord= str_split_data[7]
                                if residue_type=='SOL':
                                    residue_type='WAT'
                        elif len(str_split_data)==12:
                            print('here')
                            print('length str split: ', str_split_data )
                            if str_split_data[4]=='Z':
                                
                                print('length str split: ', str_split_data )
                                residue_type=str_split_data[3]
                                res_nr=str(str_split_data[4])+str(str_split_data[5])
                                x_coord= str_split_data[6]
                                y_coord= str_split_data[7]
                                z_coord= str_split_data[8]
                        else:
                            print('Problem, str_splitdata: ', str_split_data)
                        new_format=[residue_type, str(new_atom_name), float(x_coord), float(y_coord), float(z_coord)]
                        form_inp_data.append(new_format)
                        
                        #create PDB to check output
                        if len(res_nr) == 1:
                            fix=11
                        if len(res_nr) == 2:
                            fix=10
                        if len(res_nr) == 3:
                            fix=9
                        if len(res_nr) == 4:  
                            fix=8
                        #print(new_atom_name)
                        x_mod="%.2f" % float(x_coord)
                        y_mod="%.2f" % float(y_coord)
                        z_mod="%.2f" % float(z_coord)
                        if len(residue_type)==2:
                            #print('na')
                            #PDB_line=['ATOM',str(counter-1).ljust(5), new_atom_name.ljust(2),str(residue_type +'  N'+str(res_nr)).ljust(2), x_mod.ljust(6,'0').rjust(fix), y_mod.ljust(6,'0'),z_mod.ljust(6,'0'),'1.00','0.00','        ',new_atom_name ]
                            continue
                        elif residue_type=='WAT':
                            #print('wat')
                            #PDB_line=['ATOM',str(counter-1).ljust(5), new_atom_name.ljust(2),str(residue_type+ ' W'+str(res_nr)), x_mod.ljust(6,'0').rjust(fix), y_mod.ljust(6,'0'),z_mod.ljust(6,'0'),'1.00','0.00','        ',new_atom_name ]
                            continue
                        elif residue_type=='BTV':
                            continue
                        
                        else:
                            #PDB_line=['ATOM',str(counter-1).ljust(5), new_atom_name.ljust(2),residue_type,'X'+str(res_nr).ljust(1), x_mod.ljust(6,'0').rjust(fix), y_mod.ljust(6,'0'),z_mod.ljust(6,'0'),'1.00','0.00','        ',new_atom_name ]
                            PDB_line=['ATOM',str(counter-1).ljust(5), new_atom_name.ljust(2),residue_type, str(res_nr).ljust(1), x_mod.ljust(6,'0').rjust(fix), y_mod.ljust(6,'0'),z_mod.ljust(6,'0'),'1.00','0.00','        ',new_atom_name ]
                        PDB_output.append(PDB_line)
            f.close()   

            suffix=0
            PDB_output_renum=[]
            i=0
            resnr_list=[]
            for item in PDB_output:
                if i == 0:
                    #item[4]=suffix
                    new_resnr= str('X'+str(suffix))
                    #resnr_list.append(new_resnr+str(item[3].split(' ')))
                    tmp_resnr=item[4]
                    resnr_list.append(str(tmp_resnr[1:])+'_X_'+str(item[3].split(' ')))
                if i > 0:
                    current_atomname=item[2]
                    current_resname_tmp=item[3].split(' ')
                    current_resname=current_resname_tmp[0]
                    current_resnr=item[4]
                    prev_line=PDB_output[i-1]
                    prev_resnr=prev_line[4]
                    prev_resnr_tmp=prev_line[3].split(' ')
                    prev_resname=prev_resnr_tmp[0]
                    if current_resnr == prev_resnr and current_resname == prev_resname :
                        #new_resnr= str('X'+str(suffix))
                        new_resnr= str(suffix)
                    else:
                        suffix=suffix+1
                        new_resnr= str(current_resnr)
                        resnr_list.append(str(current_resnr[1:])+'_X_'+str(current_resname))
                i=i+1

            file_name_pdb=files.split('.')
            os.chdir(str(path_to_script))
            file_name_out='list_unique_resnr_'+file_name_pdb[0]+'.pdb'
            with open(file_name_out,'w')  as f:
                for num in resnr_list:
                    f.write('\''+str(num)+'\'' +',')
                    print(num)
            f.close()
            print('Done')
            os.chdir(path_to_script)
