#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os


path_to_script= os.getcwd()
list_of_dirs=os.listdir(os.getcwd())
for dirs in list_of_dirs:
	if dirs.endswith('loprop'):
		os.chdir(path_to_script+'/'+dirs)
		list_of_files=os.listdir(os.getcwd())
		for files in list_of_files:
			if files.endswith('.sh'):
				run_file='sbatch '+str(files)
				os.popen(run_file)   
				print('Submitting: '+files )   
