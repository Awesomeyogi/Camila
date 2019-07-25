#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 09:50:07 2019

@author: cagus
"""

import os

matching_list=[]
list_of_charges_prot=[]
list_of_charges_water=[]
list_of_pol_prot=[]
list_of_pol_water=[]
list_of_excl_prot=[]
list_of_excl_wat=[]
path_to_script=os.getcwd()
protein_files=os.listdir(os.getcwd())
os.chdir(path_to_script)
list_of_residues=[]
with open('list_unique_resnr_protein22.pdb') as d:
    for line1 in d:
        split_line1=line1.split(',')
        for item in split_line1:
            line2=item.rstrip('\n') 
            name=line2[1:-1]
            list_of_residues.append(name)      
d.close()

with open('charges_prot_small.out') as charges_prot:
	for line1 in charges_prot:
	    split_line1=line1.split(',')
	    for item in split_line1:
		line2=item.rstrip('\n') 
		name=line2[1:-1]
		list_of_charges_prot.append(name)      
charges_prot.close()

with open('polarizabilities_prot_small.out') as pol_prot:
	for line1 in pol_prot:
	    split_line1=line1.split(',')
	    for item in split_line1:
		line2=item.rstrip('\n') 
		name=line2[1:-1]
		list_of_pol_prot.append(name)      
pol_prot.close()

with open('exclist_prot_small.out') as excl_prot:
	for line1 in excl_prot:
	    line2=line1.rstrip('\n')
	    #split_line1=line1.split(',')
	    #for item in split_line1:
	    #    line2=item.rstrip('\n') 
	    #    name=line2[0:-1]
	    list_of_excl_prot.append(line2)      
excl_prot.close()

print('Converting files:')
for files in os.listdir(os.getcwd()):
    if not files.startswith('check'):
        if files.startswith('combined_protein_water_renumbered_frame') and files.endswith('.pdb'):
            print(files)
            file_name=files
            b=file_name.split('_')
            c=file_name.split('.')
            
            #open each file, read every line as a string, split and rename atoms, 
            #extract xyz-coords from pdb file
            with open(files) as f: 
		count_prot_atoms=0
		count_Na_atoms=0
		count_water_atoms=0
                form_inp_data=[]
                PDB_output=[]
                all_atom_lines=[]
		coordinates_section_protein=[]
		coordinates_section_water=[]
                counter=0
                atom_count=0
                for line in f:
                    counter=counter+1
                    if line[0:4]=='ATOM': #remove all other lines from PDB file
                        split_line=line.split()
                        residue_name=split_line[3]
                        residue_number1=split_line[4]
                        residue_count=0
                        for i in list_of_residues: # protein part of PDB file, selected from pyframe list since whole protein is in the PDB file
                            residue_count=residue_count+1
                            items3=i.split('_')
                            
                            residue_number2=residue_number1[1:]
                            #print('residue_number2',residue_number2)
                            if residue_number2=='0000':
                               residue_number=0
                               #print('residue_number',residue_number)
                            else:
                               residue_number=int(residue_number2.lstrip('0'))
                            new_resnumber=residue_number
			    #print('items3(2) resnr ', items3[0],new_resnumber)
			    #print('items3(2) resnr ', items3[2],residue_name)
			    if len(items3)>1:
                            	if items3[2]==residue_name:
                                    if int(items3[0])==new_resnumber:
                                	#print('new number: ',new_resnumber)
                                    	count_prot_atoms=count_prot_atoms+1
                                    	if str(new_resnumber)+'_X_'+str(residue_name) not in matching_list:
                                        	matching_list.append(str(new_resnumber)+'_X_'+str(residue_name))
                                        	#print(matching_list)
                                    	new_name=split_line[2]
                                    	#new_line=[new_name[0],split_line[5],split_line[6],split_line[7], count_prot_atoms ]
                                    	new_lineP=[new_name[0],split_line[5],split_line[6],split_line[7], count_prot_atoms ]
                                    	coordinates_section_protein.append(new_lineP)
                        
                        if residue_name=='WAT' or residue_name=='NA': #the only water in this combined PDB is the 22 Å around PFT
                            count_water_atoms=count_water_atoms+1
                            new_name=split_line[2]
                            if residue_name=='WAT':
                                new_name2=new_name[0]
                            if residue_name=='NA':
                                count_Na_atoms=count_Na_atoms+1
                                new_name2=new_name
                            new_line=[new_name2,split_line[5],split_line[6],split_line[7], count_water_atoms+count_prot_atoms ]
                            coordinates_section_water.append(new_line)    
            print('water atoms',count_water_atoms)

            water_pad_list=[float(0.00000000)]*84
            l = [x * 3 for x in range(1,(count_water_atoms-count_Na_atoms)/3+1)]
            l2=map(lambda x: x + count_prot_atoms, l)
            excl_list=[]
            for i in l2:
                #print i
                line1b=i-2, i-1, i #i, i+1, i+2
                line2b=i-1, i-2, i #i+1, i, i+3
                line3b=i , i-2, i-1#i+2, i, i+1
                #print(line1b)
                #print(line2b)
                #print(line3b)
                water_pad_list=[0]*84
                water_pad_list[:len(line1b)]=line1b
                excl_list.append(water_pad_list)
                water_pad_list=[0]*84
                water_pad_list[:len(line2b)]=line2b
                excl_list.append(water_pad_list)
                water_pad_list=[0]*84
                water_pad_list[:len(line3b)]=line3b
                excl_list.append(water_pad_list)
            Na_pad_list=[0]*83
                

            with open('residue_list.out','w')  as f:
                for num in matching_list:
                    f.write('\''+str(num)+'\'' +',')
                    #print(num)
            f.close()
            print('atom count protein:', count_prot_atoms)  
	    print('length of prot coor list: ',len(coordinates_section_protein))   

            file_name_out='../calculations/'+str(c[0])+'.pot'                    
            output_file=open(file_name_out,'w') 
            print >> output_file, '! potential file for '+str(b[-1])
            print >> output_file, '@COORDINATES'
            print >> output_file, count_water_atoms+count_prot_atoms 
            print >> output_file, 'AA'
            for line3 in coordinates_section_protein:
                if len(str(line3[4]))==1:
                    dist='         '
                elif len(str(line3[4]))==2:
                    dist='        '
                elif len(str(line3[4]))==3:
                    dist='       '
                elif len(str(line3[4]))==4:
                    dist='      '
                else:
                    dist='     '
                print >> output_file, str(line3[0])+'    '+str(line3[1]).ljust(7,'0')+'   '+str(line3[2]).ljust(7,'0')+'   '+str(line3[3]).ljust(7,'0')+dist+str(line3[4])
            for line4 in coordinates_section_water:
                if len(str(line4[4]))==4:
                    dist='      '
                else:
                    dist='     '
                if line4[0]=='O' or line4[0]=='H':
                    dist2='    '
                elif line4[0]=='Na':
                    dist2='   '
                print >> output_file, str(line4[0])+dist2+str(line4[1]).ljust(7,'0')+'   '+str(line4[2]).ljust(7,'0')+'   '+str(line4[3]).ljust(7,'0')+dist+str(line4[4])
            
	    print('length of protein charges: ',len(list_of_charges_prot))
            print >> output_file, '@MULTIPOLES'
            print >> output_file, 'ORDER 0'
            print >> output_file, count_water_atoms+count_prot_atoms 
            count_prot_charges=0
            for line5 in list_of_charges_prot:
                count_prot_charges=count_prot_charges+1
                if len(str(count_prot_charges))==1:
                    dist='         '
                elif len(str(count_prot_charges))==2:
                    dist='        '
                elif len(str(count_prot_charges))==3:
                    dist='       '
                elif len(str(count_prot_charges))==4:
                    dist='      '
                else:
                    dist='     '
                
                charge=line5.split()
                print >> output_file, str(count_prot_charges)+dist+str(charge[-1])
                count_water_charges=0
            for line6 in coordinates_section_water:
                count_water_charges=count_water_charges+1
                if line6[0]=='O':
                    charge_water=str(-0.669000)
                if line6[0]=='H':
                    charge_water=str(0.334500)
                if line6[0]=='Na':
                    charge_water=str(1.000000)
                print >> output_file, str(count_water_charges+count_prot_atoms)+dist+str(charge_water)
            
            print >> output_file, '@POLARIZABILITIES'
            print >> output_file, 'ORDER 1 1'
            print >> output_file, count_water_atoms+count_prot_atoms 
            count_prot_pol=0
            count_water_pol=0
            for line5 in list_of_pol_prot:
                count_prot_pol=count_prot_pol+1
                if len(str(count_prot_pol))==1: 
                    dist='           '
                elif len(str(count_water_pol+count_prot_pol))==2:
                    dist='          '
                elif len(str(count_water_pol+count_prot_pol))==3:
                    dist='         '
                elif len(str(count_water_pol+count_prot_pol))==4:
                    dist='        '
                    fix=7
                else:
                    dist='     '
                space2= '    '
                pol=line5.split()
                print >> output_file, str(count_prot_pol)+dist+str(pol[-6]).rjust(11,' ')+space2+str(pol[-5]).rjust(11,' ')+space2+str(pol[-4]).rjust(11,' ')+space2+str(pol[-3]).rjust(11,' ')+space2+str(pol[-2]).rjust(11,' ')+space2+str(pol[-1])
            for line7 in coordinates_section_water:
                count_water_pol=count_water_pol+1
                if len(str(count_water_pol+count_prot_pol))==4:
                    dist='         '
                    fix=10
                elif len(str(count_water_pol+count_prot_pol))==5:
                    dist='        '
                    fix=10
                else:
                    dist='      '
                    fix=8
                space2= '     '
		space3= '    '
                if line7[0] == 'O':
                    print >> output_file,str(count_prot_pol+count_water_pol)+dist+ str(5.73935000).ljust(fix,'0') +space2 +str(0.00000000).ljust(fix,'0')+ space2 +str(0.00000000).ljust(fix,'0')+space2+str(5.73935000).ljust(fix,'0') +space2 +str(0.00000000).ljust(fix,'0') +space3 +str(5.73935000).ljust(fix,'0') # charge from sep.csv in pyframe
                if line7[0] == 'H':
                    print >> output_file,str(count_prot_pol+count_water_pol)+dist+ str(2.30839000).ljust(fix,'0') +space2 +str(0.00000000).ljust(fix,'0')+ space2 +str(0.00000000).ljust(fix,'0')+space2+str(2.30839000).ljust(fix,'0') +space2 +str(0.00000000).ljust(fix,'0') +space3 +str(2.30839000).ljust(fix,'0')
                if line7[0] == 'Na':
                    print >> output_file,str(count_prot_pol+count_water_pol)+dist+str(0.48203000).ljust(fix,'0') +space2 +str(0.00000000).ljust(fix,'0')+ space2 +str(0.00000000).ljust(fix,'0')+space2+str(0.48203000).ljust(fix,'0') +space2 +str(0.00000000).ljust(fix,'0') +space3 +str(0.48203000).ljust(fix,'0')
            print >> output_file, 'EXCLISTS'
            print >> output_file, str(count_water_atoms+count_prot_atoms)+' 83'
            for line8 in list_of_excl_prot:
                print >> output_file, line8
            count_water_pol=0
            dist2='     '
            dist3='   '
            dist4='      '
            for line9 in excl_list:
                count_water_pol=count_water_pol+1
                a=str(line9).replace(',','     ') 
                b=str(a).strip('[');
                c=str(b).strip(']');
                d=str(c).strip('\'');
                splitD=d.split()
                #if len(splitD[0])==5:
                    
                splitD2=splitD[3:-1]
                splitE=str(splitD2).strip('[')
                splitF=str(splitE).strip(']')
                splitG=str(splitF).replace('\', \'','      ')
                splitH=splitG.strip('\'')
                if len(splitD[0])==4:
                    fix1=9
                elif len(splitD[0])==5:
                    fix1=8
                else:
                    print('problem, length:', len(splitD[0]))
                if len(splitD[1])==4:
                    fix2=7
                elif len(splitD[1])==5:
                    fix2=7
                else:
                    print('problem, length:', len(splitD[1]))
                if len(splitD[2])==4:
                    fix3=10
                elif len(splitD[2])==5:
                    fix3=11
                else:
                    print('problemos')
		    out_name=files.split('.')
		    out_name2=out_name.split('_')
		    with open('error_log_pot_file_'+str(out_name2[-1])) as g:
		   	g.write(files) 			         
		    g.close()
                #print >> output_file,str(splitD[0])+dist2+str(splitD[1])+dist3+str(splitD[2])+dist4+str(splitH)
                print >> output_file,str(splitD[0]).ljust(fix1,' ')+str(splitD[1]).ljust(fix2,' ')+str(splitD[2]).ljust(fix3,' ')+str(splitH)
            number_Na=1

            for line10 in range(count_Na_atoms):
                a=str(Na_pad_list).replace(',','     ') 
                b=str(a).strip('[');
                c=str(b).strip(']');
                e=str(c).strip('\'');
                print >> output_file,str(count_prot_pol+count_water_pol+number_Na) ,str(e)
                number_Na=number_Na+1

            #water_pad_list[:len(x)] = x 
            output_file.close()

                    
            f.close()

    
