#!/usr/bin/env python2

from collections import Counter
import os
from operator import itemgetter


#read in all files in the gro_out_files directory with file extension .pdb
path_to_script=os.getcwd()
origin_folder=''
os.chdir(str(path_to_script+origin_folder))

print('Converting files:')
for files in os.listdir(os.getcwd()):
    os.chdir(str(path_to_script+origin_folder))
    if not files.startswith('check'):
        if files.startswith('PFT_ext_frame_') and files.endswith('.pdb'):
            print(files)
            
            #open each file, read every line as a string, split and rename atoms, 
            #extract xyz-coords from pdb file
            with open(files) as f: 
                form_inp_data=[]
                PDB_output=[]
                counter=0
                for line in f:
                    counter=counter+1
                    if line[0:4]=='ATOM': #remove all other lines from PDB file
                        str_split_data=line.split()
                        orig_atom_name=str_split_data[2] #split the atom name argument in order to sort on first letter                        
                        if orig_atom_name[0] =='1' or orig_atom_name[0] =='2' or orig_atom_name[0] == '3'or orig_atom_name[0] == '4'or orig_atom_name[0] == '5' or orig_atom_name[0] == '6' or orig_atom_name[0] == '7' or orig_atom_name[0] == '8' or orig_atom_name[0] == '9'or orig_atom_name[0] == '0':#numb_list[numb]:
                            #print('check')
                            new_atom_name=orig_atom_name[1]
                        elif orig_atom_name[0:2] == 'Na' or orig_atom_name[0:2] == 'Cl' or orig_atom_name[0:2] == 'NA': # upgrade to include atom names with more than one char with islower() or smt
                            #print('check2')
                            new_atom_name=orig_atom_name[0:2]
                            #print(new_atom_name)
                        else:
                            new_atom_name=orig_atom_name[0]
                        if len(str_split_data) ==12:
                            #print(str_split_data)
                            residue_type=str_split_data[3]
                            res_nr=str_split_data[5]
                            x_coord= str_split_data[6]
                            y_coord= str_split_data[7]
                            z_coord= str_split_data[8]
                        elif len(str_split_data) ==10:
                            if str_split_data[4].startswith('X'):
                                res_nr_tmp=str_split_data[3]
                                res_nr=res_nr_tmp[1:]
                                x_coord= str_split_data[4]
                                y_coord= str_split_data[5]
                                z_coord= str_split_data[6]
                            
                        else:
                            print('Probem with PDB file format. Delimiters?')
                        new_format=[residue_type, str(new_atom_name), float(x_coord), float(y_coord), float(z_coord)]
                        form_inp_data.append(new_format)
                        
                        #create PDB to check output
                        if len(str(counter-1)) == 1:
                            space1='  '
                        if len(str(counter-1)) == 2:
                            space1=' '
                        if len(str(counter-1)) == 3:
                            space1=''
                        if len(res_nr) == 1:
                            space2='   '
                        if len(res_nr) == 2:
                            space2='  '
                        if len(res_nr) == 3:
                            space2=' '
                        if len(res_nr) == 4:  
                            space2=''
			if len(res_nr) == 5:  
                            space2=''
			if len(res_nr) == 7:  
                            space2=''

			#print('length of resnr: ',len(res_nr))
                        #print(new_atom_name)
                        x_mod="%.2f" % float(x_coord)
                        y_mod="%.2f" % float(y_coord)
                        z_mod="%.2f" % float(z_coord)
                        PDB_line=['ATOM',space1,str(counter-1), new_atom_name,' '+str(str_split_data[3]),'X',res_nr,space2,x_mod, y_mod,z_mod,'1.00','0.00' ]
                        PDB_output.append(PDB_line)

                        #result= "\n".join("\t".join(map(str,l)) for l in PDB_output)
                        file_name_pdb=files.split('.')
                        file_name_out='check_'+file_name_pdb[0]+'.pdb'
                        output_file_pdb=open(file_name_out,'w')  
                        print >> output_file_pdb,  'CRYST1  148.873  148.873  178.647  90.00  90.00  90.00 P 1           1'
                        for element in PDB_output:
                            str_el=str(element)
                            #str_line=str_el.replace(',','\t').replace('[','').replace(']','').replace('\'','').replace(' ','')
                            str_line=str_el.replace(',',' ').replace('[','').replace(']','').replace('\'','').replace(' X','X')
                            print >> output_file_pdb,  str_line
                            #print >> output_file_pdb,  str_el.replace(',|[|]',"")
                            #print >> output_file_pdb,  element.replace(',|[|]','', PDB_output)
                        print >> output_file_pdb,  'END'
                        
                    #else:
                     #   print('Oxygen in PFT deleted')
                        #print(line)
            f.close()
            #os.chdir('../')
            #os.chdir('../../')
            
    
            #sort atoms after types 
            second_item=itemgetter(1)
            sorted_list=sorted(form_inp_data, key=second_item)
            #sorted_list=sorted(form_inp_data);
            list_atom_names=[item[1] for item in sorted_list]#form_inp_data]
            unique_atom_names=[]
            
           
            #count number of atom types
            for letter in list_atom_names: 
                #print
                if letter not in unique_atom_names:
                    unique_atom_names.append(letter)
            nbr_atom_types=len(unique_atom_names)
            nbr_atoms=Counter(list_atom_names)
            
            # count number of atoms of each type
            type_nbr=[] 
            for atm_tp in unique_atom_names:
                #print(nbr_atoms[atm_tp])
                type_nbr.append(nbr_atoms[atm_tp])
            
            # ------>  currently static number of charges. FIX IF OTHER ATOMS!!    
            atom_charges_alphabetically=[6.0, 1.0, 8.0, 16.0] #['C', 'H', 'O', 'S'] 
            
            # create directory for all Dalton input files
            os.chdir(str(path_to_script)+'/../calculations/')
            
            #write Dalton input formatted file
            str1= 'Atomtypes='+str( nbr_atom_types )+ ' Charge=-4.0 Generators=0 Angstrom'
            name_file=files.split('.')
            new_name=name_file[0]+'.dal'
            output_file=open(new_name,'w')  
            print >> output_file, 'BASIS'
            print >> output_file, 'aug-cc-pVDZ'
            print >> output_file, 'QM region, snapshot from MD, '+ name_file[0]
            print >> output_file, '======================================'
            print >> output_file, str1
            for atom in range(len(unique_atom_names)):
               str2= type_nbr[atom] 
               str3= atom_charges_alphabetically[atom]
               print >> output_file, 'Charge='+str(str3)+' Atoms='+str(str2)
               for line2 in sorted_list:
                   if line2[1] == unique_atom_names[atom]: #  ------> fix formatting with minus sign? test run dalton to see if accepted
                       print >> output_file,' '+ str(line2[1])+'                  '+ str(line2[2]).ljust(7,'0')+ '    '+ str(line2[3]).ljust(7,'0') + '    '+str( line2[4]).ljust(7,'0')
            print >> output_file, '#end of input'
            print >> output_file,  ' '
            print >> output_file, '**DALTON INPUT'
            print >> output_file, '.RUN RESPONSE'
            print >> output_file, '.DIRECT'
            print >> output_file, '.PEQM'
	    print >> output_file, '*PEQM'
            print >> output_file, '.ISOPOL'
            print >> output_file, '**INTEGRALS'
            print >> output_file, '.DIPLEN'
            print >> output_file, '**WAVE FUNCTION'
            print >> output_file, '.DFT'
            print >> output_file, 'CAMB3LYP'
            print >> output_file, '**RESPONSE'
            print >> output_file, '*LINEAR'
            print >> output_file, '.DIPLEN'
            print >> output_file, '.SINGLE RESIDUE'
            print >> output_file, '.ROOTS'
            print >> output_file, ' 10'
            print >> output_file, '*END OF INPUT '
            
            output_file.close()
            os.chdir(path_to_script)
        
    

