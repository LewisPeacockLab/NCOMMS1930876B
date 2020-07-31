#!/usr/bin/env python
# modified from Mumford's script:
# https://www.dropbox.com/s/xaa974vxxpqvizn/2_prep_bold.py?dl=0
import glob
import os
import sys
import subprocess
import pdfkit

###########################################
# file/directory variables
###########################################
# setup directories
home_dir    = os.path.expanduser('~')
data_home   = os.path.join(home_dir,'lewpealab_dpbx/STUDY/Clearmem/imaging_data')
subject_dir = glob.glob('%s/clearmem_v1_sub*' % (data_home))
                  
for xsubdir in subject_dir:
	print(xsubdir)

	###########################################
	# motion correction
	###########################################   
	bold_files = glob.glob('%s/bold/clearmem*/bold.nii.gz' % (xsubdir))

	if os.path.isdir("%s/motion_assess/" % (data_home))==False:
	  os.system("mkdir %s/motion_assess" % (data_home))

	outhtml = "%s/motion_assess/bold_motion_QA.html" % (data_home)
	out_bad_bold_list = "%s/motion_assess/subs_lose_gt_45_vol_scrub.txt" % (data_home)

	for xbold in list(bold_files):
		print(xbold)
		# Store directory name
		cur_dir = os.path.dirname(xbold)

		# strip off .nii.gz from file name (makes code below easier)
		cur_bold_no_nii = xbold[:-7]

		# Assessing motion.
		if os.path.isdir("%s/motion_assess/" % (cur_dir))==False:
		  os.system("mkdir %s/motion_assess" % (cur_dir))

		# -i: input, -o: output, --fd: use FD (framewise displacement) as metric, -
		# p <filename>: save metric values (e.g. DVARS) as a graphical plot (png format) 
		# -v: verbose mode
		os.system("fsl_motion_outliers -i %s -o %s/motion_assess/confound.txt --fd --thresh=0.9 -p %s/motion_assess/fd_plot -v > %s/motion_assess/outlier_output.txt"
				   % (cur_bold_no_nii, cur_dir, cur_dir, cur_dir))

		# Put confound info into html file for review later on
		os.system("cat %s/motion_assess/outlier_output.txt >> %s"%(cur_dir, outhtml))
		os.system("echo '<p>=============<p>FD plot %s <br><IMG BORDER=0 SRC=%s/motion_assess/fd_plot.png WIDTH=100%s></BODY></HTML>' >> %s"
				   % (cur_dir, cur_dir,'%', outhtml))

		# if we're planning on modeling out scrubbed volumes later
		if os.path.isfile("%s/motion_assess/confound.txt" % (cur_dir))==False:
		  os.system("touch %s/motion_assess/confound.txt" % (cur_dir))

		# Very last, create a list of subjects who exceed a threshold for
		#  number of scrubbed volumes.  This should be taken seriously.  If
		#  most of your scrubbed data are occurring during task, that's
		#  important to consider (e.g. subject with 20 volumes scrubbed
		#  during task is much worse off than subject with 20 volumes
		#  scrubbed during baseline.
		output = subprocess.check_output("grep -o 1 %s/motion_assess/confound.txt | wc -l" 
										  % (cur_dir), shell=True)
		num_scrub = [int(s) for s in output.split() if s.isdigit()]
		if num_scrub[0]>45:
			with open(out_bad_bold_list, "a") as myfile:
			  myfile.write("%s\n"%(xbold))
		
	# Save html to pdf		   	
	pdfkit.from_file(outhtml, "%s/motion_assess/bold_motion_QA_%s.pdf" % (data_home, xsubdir.split('/')[-1]))

	os.system("rm %s" % (outhtml))
	os.system("rm %s" % (out_bad_bold_list))	
	
