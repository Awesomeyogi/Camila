#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#export PATH="/home/x_camgu/anaconda3/bin:$PATH"

import pyframe
import os
system=pyframe.MolecularSystem(input_file='combined_protein_water_renumbered_frame0.pdb')
core=system.get_fragments_by_name(names=['PFT'])
system.set_core_region(fragments=core, program='Dalton', basis ='aug-cc-pVDZ')
protein=system.get_fragments_by_identifier(identifiers=['11_X_VAL','12_X_HIS','13_X_GLN','14_X_LYS','15_X_LEU','16_X_VAL','17_X_PHE','18_X_ALA','19_X_GLU','20_X_ASP','21_X_VAL','24_X_ASN','25_X_LYS','26_X_GLY','27_X_ALA','28_X_ILE','29_X_GLY','30_X_LEU','47_X_GLU','48_X_VAL','49_X_HIS','50_X_GLN','51_X_LYS','52_X_LEU','53_X_VAL','54_X_PHE','55_X_ALA','56_X_GLU','57_X_ASP','58_X_VAL','60_X_SER','61_X_ASN','62_X_LYS','63_X_GLY','64_X_ALA','65_X_ILE','66_X_GLY','67_X_LEU','68_X_MET','71_X_VAL','120_X_TYR','121_X_GLU','122_X_VAL','123_X_HIS','124_X_GLN','125_X_LYS','126_X_LEU','127_X_VAL','128_X_PHE','129_X_ALA','130_X_GLU','131_X_ASP','132_X_VAL','133_X_GLY','134_X_SER','135_X_ASN','136_X_LYS','137_X_GLY','138_X_ALA','139_X_ILE','140_X_GLY','141_X_LEU','142_X_MET','145_X_VAL','190_X_HIS','192_X_SER','194_X_TYR','195_X_GLU','196_X_VAL','197_X_HIS','198_X_GLN','199_X_LYS','200_X_LEU','201_X_VAL','202_X_PHE','203_X_ALA','204_X_GLU','205_X_ASP','206_X_VAL','207_X_GLY','208_X_SER','209_X_ASN','210_X_LYS','211_X_GLY','212_X_ALA','213_X_ILE','214_X_GLY','215_X_LEU','216_X_MET','219_X_VAL','262_X_PHE','264_X_HIS','266_X_SER','268_X_TYR','269_X_GLU','270_X_VAL','271_X_HIS','272_X_GLN','273_X_LYS','274_X_LEU','275_X_VAL','276_X_PHE','277_X_ALA','278_X_GLU','279_X_ASP','280_X_VAL','281_X_GLY','282_X_SER','283_X_ASN','284_X_LYS','285_X_GLY','286_X_ALA','287_X_ILE','288_X_GLY','289_X_LEU','290_X_MET','293_X_VAL','336_X_PHE','338_X_HIS','340_X_SER','342_X_TYR','343_X_GLU','344_X_VAL','345_X_HIS','346_X_GLN','347_X_LYS','348_X_LEU','349_X_VAL','350_X_PHE','351_X_ALA','352_X_GLU','353_X_ASP','354_X_VAL','355_X_GLY','356_X_SER','357_X_ASN','358_X_LYS','359_X_GLY','360_X_ALA','361_X_ILE','362_X_GLY','363_X_LEU','364_X_MET','367_X_VAL','410_X_PHE','412_X_HIS','414_X_SER','416_X_TYR','417_X_GLU','418_X_VAL','419_X_HIS','420_X_GLN','421_X_LYS','422_X_LEU','423_X_VAL','424_X_PHE','425_X_ALA','426_X_GLU','427_X_ASP','428_X_VAL','429_X_GLY','430_X_SER','431_X_ASN','432_X_LYS','433_X_GLY','434_X_ALA','435_X_ILE','436_X_GLY','437_X_LEU','438_X_MET','441_X_VAL','443_X_ALA','484_X_PHE','486_X_HIS','490_X_TYR','491_X_GLU','492_X_VAL','493_X_HIS','494_X_GLN','495_X_LYS','496_X_LEU','497_X_VAL','498_X_PHE','499_X_ALA','500_X_GLU','501_X_ASP','502_X_VAL','503_X_GLY','504_X_SER','505_X_ASN','506_X_LYS','507_X_GLY','508_X_ALA','509_X_ILE','510_X_GLY','511_X_LEU','512_X_MET','515_X_VAL','517_X_ALA','564_X_TYR','565_X_GLU','566_X_VAL','567_X_HIS','568_X_GLN','569_X_LYS','570_X_LEU','571_X_VAL','572_X_PHE','573_X_ALA','574_X_GLU','575_X_ASP','576_X_VAL','577_X_GLY','578_X_SER','579_X_ASN','580_X_LYS','581_X_GLY','582_X_ALA','583_X_ILE','584_X_GLY','585_X_LEU','586_X_MET','589_X_VAL','639_X_GLU','640_X_VAL','641_X_HIS','642_X_GLN','643_X_LYS','644_X_LEU','645_X_VAL','646_X_PHE','647_X_ALA','648_X_GLU','649_X_ASP','650_X_VAL','651_X_GLY','652_X_SER','653_X_ASN','654_X_LYS','655_X_GLY','656_X_ALA','657_X_ILE','658_X_GLY','659_X_LEU','660_X_MET','714_X_VAL','715_X_HIS','716_X_GLN','717_X_LYS','718_X_LEU','719_X_VAL','720_X_PHE','721_X_ALA','722_X_GLU','723_X_ASP','724_X_VAL','730_X_ALA','731_X_ILE','732_X_GLY','733_X_LEU'])
system.add_region(name='protein', fragments=protein, use_mfcc=True, use_multipoles = True, multipole_order=0 , multipole_model= 'LoProp',multipole_method= 'DFT', multipole_xcfun= 'CAMB3LYP',multipole_basis= 'loprop-6-31+G*',use_polarizabilities= True, polarizability_model= 'LoProp',polarizability_method= 'DFT',polarizability_xcfun='CAMB3LYP',polarizability_basis= 'loprop-6-31+G*')

# JMHO: in case you also include water in the PDB you simply add the following line
# JMHO: this will take all water molecules with center-of-mass (COM) that is within 25.0 Å from the
# JMHO: from the COM of the core. If you don't want to use COM, then add "use_center_of_mass=False".
water = system.get_fragments_by_distance_and_name(names=['WAT'], reference=core, distance=22.0, use_center_of_mass=False)
system.add_region(name='solvent', fragments=water, use_standard_potentials=True, standard_potential_model='SEP')
sodions = system.get_fragments_by_distance_and_name(names=['NA'], reference=core, distance=22.0, use_center_of_mass=False)
system.add_region(name='solvent2', fragments=sodions, use_standard_potentials=True, standard_potential_model='SEP')

# JMHO: It is a good idea to check the settings for the calculation because PyFraME will run a calculation
# JMHO: for each protein fragment.
project=pyframe.Project()
print('Imported path:')
print(os.environ['pyframe_tmpdir'])
project.scratch_dir = os.environ['pyframe_tmpdir']
project.print_info()
project.create_embedding_potential(system)
print('Embedding potential created')
project.write_core(system)
project.write_potential(system)
print('Embedding potential written')



# JMHO: What you can do for the following frames is to reuse the output generated by PyFraME, i.e.
# JMHO: the *_dalton_loprop.out files, and then run this same script again but of course with a new
# JMHO: filename for the PDB. By default PyFraME will use a directory named after the PDB file, so
# JMHO: you can create it in advance and put all the *_dalton_loprop.out in there before running the
# JMHO: script again.
#


