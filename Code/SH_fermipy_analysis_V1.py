#!/usr/bin/env python
from time import sleep
from random import randint

#June 16, 25: changed binsperdec to 8 from 10. 


#sleep(randint(1,93600)) # 93600 for full telescope time, for flares, lesser 72000 or something. uncomment for launching at the SLAC farm, 10800: for full times
#sleep(randint(1,10000))
#sleep(randint(1,150800))
#sleep(randint(1,10800))
#sleep(randint(1,7200))


import resource
resource.setrlimit(resource.RLIMIT_CORE, (0,-1))
import random,shutil,yaml
import os,sys
from math import *
import numpy as np
import matplotlib as matplotlib
matplotlib.use('agg')
from fermipy.gtanalysis import GTAnalysis
import matplotlib.pyplot as plt
import pandas as pd
import astropy.io.fits as pyfits
from astropy.table import Table
from pathlib import Path

import xml.etree.ElementTree as ET


#Issue: roiwidth
#consider this as the base code
#**********************************************************************

def fix_parameters(source,xml_file,gta):
    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # Initialize the GTAnalysis object
    #gta = GTAnalysis.create_from_xml(xml_file)

    # Iterate over each source in the XML file
    for source_elem in root.findall('.//source'):
        source_name = source_elem.attrib['name']
        source_type = source_elem.attrib['type']
        
        # Exclude 'galdiff' and 'isodiff' sources
        source = source.replace(" ", "")
        if source_name == 'galdiff' or source_name == 'isodiff' or source_name==source:
            continue

        # Fix the parameters for the source
        for param_elem in source_elem.findall('.//parameter'):
            param_name = param_elem.attrib['name']
            param_free = param_elem.attrib['free']
            
            # If the parameter is free, fix it
            if param_free == '1':
                print(f"Fixing parameter '{param_name}' for source '{source_name}'")
                gta.free_source(source_name,free=False)
                break

    # Perform the fit again with fixed parameters
    gta.fit()




#*********** Setting PFILES Environment *******
def SetPfilesDir(outdir):
	pfile_dir=os.getcwd()+"/PFILES_%s"%(random.random())
	#hea_pdir = os.getenv("HEADAS")+"/syspfiles/"
	gt_pdir  = os.getenv("INST_DIR")+"/syspfiles/"
	if(os.path.isdir(pfile_dir)==False):
		print("\n >> Creating pfile dir "+pfile_dir)
		cmd = "mkdir "+pfile_dir
		os.system(cmd)
	#cmd = "cp %s/*.par %s"%(hea_pdir,pfile_dir)
	#os.system(cmd)
	cmd = "cp %s/*.par %s"%(gt_pdir,pfile_dir)
	os.system(cmd)
	os.putenv("PFILES",pfile_dir)
	return pfile_dir
#*************************************************


def main(cmd_line):

	name1 = cmd_line[1]
	name2 = cmd_line[2]
	#peak=cmd_line[3]
	time1=float(cmd_line[3])
	time2=float(cmd_line[4])
	peak=cmd_line[5]

	# Convert MJD to MET (Fermi Mission Elapsed Time, seconds since 2001-01-01 00:00:00 UTC)
	# MET = (MJD - 51910.0) * 86400.0
	FERMI_MJD_REF = 51910.0
	time1_MET = (time1 - FERMI_MJD_REF) * 86400.0
	time2_MET = (time2 - FERMI_MJD_REF) * 86400.0

	srcname=name1+name2
	print("srcname,time1,time2,window_id:",srcname,time1,time2,peak)
	print("MET range: %.0f to %.0f" % (time1_MET, time2_MET))
	directory='/sdf/data/fermi/u/adinesh'
	print('Data directory: %s' % (directory))
	galdiff = '%s/DATA/Diffuse/gll_iem_v07.fits' % (directory)

	# --- Source-specific data files ---
	# Map each 4FGL source name to its event list and spacecraft file.
	# If you have source-specific pre-filtered data, specify here.
	# Otherwise use the all-sky files: lat_alldata.fits / SC.fits
	DATA_MAP = {
		'4FGLJ0522.9-3628': {
		        'ft1': '%s/DATA/J0522/events_j0522.txt' % directory,
		        'ft2': '%s/DATA/J0522/L260416101627E79FF56C01_SC00.fits' % directory,
		},
		# Add more sources here as needed
		
	}
	if srcname not in DATA_MAP:
	    raise ValueError("No source-specific FT1/FT2 configured for %s" % srcname)
	    
	ft1 = DATA_MAP[srcname]['ft1']
	ft2 = DATA_MAP[srcname]['ft2']
	print('Using source-specific data for %s' % srcname)
	print('ft1: %s' % ft1)
	print('ft2: %s' % ft2)

	# --- Directory structure: source/run_id ---
	file_tag = '%s_%s' % (srcname, peak)
	
	outdir='/sdf/scratch/users/a/adinesh/SH_Alba/%s/%s' % (srcname, peak)
	copy_dir=Path('/sdf/data/fermi/u/adinesh/results/SH_Alba')
	# Get the path of the currently executing script
	current_script_path = Path(__file__).resolve() 
	# Copy the current script to the top-level results directory (once is enough)
	top_copy = copy_dir / srcname
	top_copy.mkdir(parents=True, exist_ok=True)
	if not (top_copy / current_script_path.name).exists():
		shutil.copy2(current_script_path, top_copy / current_script_path.name)

	hdulist = pyfits.open("/sdf/home/a/adinesh/fermi-user/catalogues/table_4LAC_DR3_h.fits")
	tbdata = hdulist[1].data
	tbdata=Table(tbdata)
	df = tbdata.to_pandas()
	glon=df.loc[df['Source_Name'] == name1+" "+name2, 'RAJ2000'].iloc[0]
	glat=df.loc[df['Source_Name'] == name1+" "+name2, 'DEJ2000'].iloc[0]
	print(glon,glat)
	
	#pivot energy
	hdulist1 = pyfits.open("/sdf/home/a/adinesh/fermi-user/catalogues/gll_psc_v31.fit")
	tbdata1 = hdulist1[1].data
	srcname1=name1+" "+name2
	mask = (tbdata1['Source_Name'] == srcname1)
	# Access the Pivot_Energy value for the matching row.
	pivot_energy = tbdata1['Pivot_Energy'][mask][0]
	print(pivot_energy)
	
	#f(os.path.isdir(outdir)==True):
	#       shutil.rmtree(outdir)
	os.makedirs(outdir, exist_ok=True)
	os.chdir(outdir)
	pdir = SetPfilesDir("./")
	
	with open('%s.yaml' % srcname, 'w') as yml:
		yml.write("logging:\n")
		yml.write("  verbosity : 3\n")
		yml.write("#--------#\n")
		yml.write("fileio:\n")
		yml.write("  outdir : output\n")
		yml.write("  logfile : %s\n" % srcname)
		yml.write("  usescratch : False\n")
		yml.write("  scratchdir : scratch\n")
		yml.write("#--------#\n")
		yml.write("data:\n")
		yml.write("  evfile : %s\n" %ft1)
		yml.write("  scfile : %s\n" %ft2)
		yml.write("#--------#\n")
		yml.write("binning:\n")
		yml.write("  roiwidth : 15\n")
		#yml.write("  binsz : 0.05\n")
		yml.write("  binsperdec : 10\n")
		yml.write("#--------#\n")
		yml.write("selection:\n")
		##yml.write("  emin : 80\n")
		##yml.write("  emax : 300000\n")
		##yml.write("  zmax : 90\n")
		yml.write("  target : '%s'\n" % srcname)
		yml.write("  radius : 20\n")
		yml.write("  tmin : %.0f \n" %time1_MET)
		yml.write("  tmax : %.0f \n" %time2_MET) # 574560005
		yml.write("  evclass : 128\n")
		##yml.write("  evtype : 3\n")
		yml.write("  filter : '(DATA_QUAL>0)&&(LAT_CONFIG==1)&&(angsep(%s,%s, RA_SUN, DEC_SUN)>15)'\n"%(glon,glat))
		yml.write("#--------#\n")
		yml.write("gtlike:\n")
		yml.write("  edisp : True\n")
		yml.write("  edisp_bins : -1\n")
		yml.write("  edisp_disable : ['isodiff']\n")
		yml.write("  irfs : 'P8R3_SOURCE_V3'\n")
		yml.write("#--------#\n")
		yml.write("model:\n")
		yml.write("  src_roiwidth : 20\n")
		yml.write("  galdiff : '%s'\n" %galdiff)
		##yml.write("  isodiff : '%s'\n" %isodiff)
		yml.write("  catalogs :\n")
		yml.write("    - '4FGL-DR3'\n")
		yml.write("  extdir : '%s/Extended_archive_v18/'\n" % (directory))
		yml.write("#--------#\n")
		yml.write("components:\n") 
		#yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF3_v1_extrapolated.txt"},\n')
		#yml.write('     selection : { evtype : 32, zmax : 80, emin : 50 , emax : 100 } }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF3_v1.txt"},\n')
		yml.write('     selection : { evtype : 32, zmax : 90, emin : 100, emax : 300 } }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF2_v1.txt"},\n')
		yml.write('     selection : { evtype : 16, zmax : 90, emin : 100, emax : 300} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF3_v1.txt"},\n')
		yml.write('     selection : { evtype : 32, zmax : 100, emin : 300, emax : 1000} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF2_v1.txt"},\n')
		yml.write('     selection : { evtype : 16, zmax : 100, emin : 300, emax : 1000} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF1_v1.txt"},\n')
		yml.write('     selection : { evtype : 8, zmax : 100, emin : 300, emax : 1000} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF3_v1.txt"},\n')
		yml.write('     selection : { evtype : 32, zmax : 105, emin : 1000, emax : 300000} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF2_v1.txt"},\n')
		yml.write('     selection : { evtype : 16, zmax : 105, emin : 1000, emax : 300000} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF1_v1.txt"},\n')
		yml.write('     selection : { evtype : 8, zmax : 105, emin : 1000, emax : 300000} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF0_v1.txt"},\n')
		yml.write('     selection : { evtype : 4, zmax : 105, emin : 1000, emax : 300000} }\n')
		#High energy part
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF3_v1.txt"},\n')
		yml.write('     selection : { evtype : 32, zmax : 105, emin : 300000, emax : 2000000} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF2_v1.txt"},\n')
		yml.write('     selection : { evtype : 16, zmax : 105, emin : 300000, emax : 2000000} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF1_v1.txt"},\n')
		yml.write('     selection : { evtype : 8, zmax : 105, emin : 300000, emax : 2000000} }\n')
		yml.write('  - { model: {isodiff: "/sdf/group/fermi/g/neutrino/shared/software/diffuse/iso_P8R3_SOURCE_V3_PSF0_v1.txt"},\n')
		yml.write('     selection : { evtype : 4, zmax : 105, emin : 300000, emax : 2000000} }\n')
		
		

		yml.write("#--------#\n")
		yml.write("plotting:\n")
		yml.write("  format : png\n")
		yml.write("#--------#\n")
		yml.write("sed:\n")
		yml.write("  use_local_index : True")
	yml.close()
	gta = GTAnalysis('%s.yaml' % srcname,logging={'verbosity' : 3})
	gta.setup()
	#gta.set_energy_range(np.log10(100.),np.log10(300000.)) 
	gta.optimize()
	gta.print_roi()
	gta.free_source('galdiff', pars='norm')
	gta.free_source('isodiff')
	gta.free_sources(minmax_ts=[25,None],pars='norm', distance=9.0)
	gta.free_sources(minmax_ts=[500,None],distance=12.0)
	gta.free_source(srcname)
	gta.fit()
	gta.write_roi('fit_model_1')
	p = np.load('output/fit_model_1.npy',allow_pickle=True).flat[0]
	src = p['sources']['%s %s' %(name1,name2)]
	Spectrum1=src['SpectrumType']
	if src['ts']>=25 and Spectrum1=='PowerLaw':
		Index1=src['param_values'][1]
	elif src['ts']>=25 and Spectrum1=='LogParabola':
		Index1=src['param_values'][1]
	else:
		Index1=2.0
	model1 = {'Index' : np.fabs(Index1), 'SpatialModel' : 'PointSource'}
	finder1=gta.find_sources(prefix='find',model=model1,sqrt_ts_threshold=5.0,
                                 min_separation=0.5,max_iter=10,sources_per_iter=20,tsmap_fitter='tsmap')
	gta.write_roi('fit_model_2')
	gta.optimize()
	gta.print_roi()
	gta.free_source('galdiff', pars='norm')
	gta.free_source('isodiff')
	gta.free_sources(minmax_ts=[25,None],pars='norm',distance=9.0)
	gta.free_sources(minmax_ts=[500,None],distance=12.0)
	gta.free_source(srcname)
	
	
	
	fit3=gta.fit()
	gta.write_roi('fit_model_3')
	#Adding the residualmaps
	gta.residmap('resid_fit_model_3_%s' % file_tag, model=model1, make_plots=True)
	gta.tsmap('tsmap_fit_model_3_%s' % file_tag, model=model1, make_plots=True)
	o = np.load('output/fit_model_3.npy',allow_pickle=True).flat[0]
	src3 = o['sources']['%s %s' %(name1,name2)]
	TS_PL = '%.2f' %src3['ts']
	PL_like=fit3['loglike']
	TS_fit3 = '%.2f' %src3['ts']
	fit3_like=fit3['loglike']
	
	energy_bins=np.array([50, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 2000000])
	log_energy_bins=np.log10(energy_bins)
	# 9 equal log-spaced bins, 100 MeV to 2 TeV (following 2FHL philosophy)
	# Edges at: 100, 300.5, 903, 2714, 8158, 24517, 73681, 221435, 665484, 2000000 MeV
	# Centers (geometric mean): 173, 521, 1566, 4706, 14142, 42502, 127732, 383877, 1153676 MeV
	#energy_bins = np.logspace(np.log10(100), np.log10(2000000), 10)
	#log_energy_bins = np.log10(energy_bins)

	
	sed = gta.sed(srcname,loge_bins=log_energy_bins, free_radius=3, make_plots=True)
	os.rename("%s/output/4fgl_j%s_sed.png"%(outdir,srcname[5:16]), "%s/output/sed_fit_model_3_%s.png" %(outdir,file_tag))
	os.system("rm -rf %s/output/4fgl_j%s_sedlnl.png"%(outdir,srcname[5:16]))
	os.rename("%s/output/4fgl_j%s_sed.fits"%(outdir,srcname[5:16]), "%s/output/sed_fit_model_3_%s.fits" %(outdir,file_tag))
	Spectrum3=src3['SpectrumType']
	
	#fitting powerlaw
	model=gta.set_source_spectrum(srcname, spectrum_type='PowerLaw', spectrum_pars={'Prefactor' : {'value' : 1.0, 'scale' : 1e-12,
                                                         'max' : 1000.0, 'min' : 1e-5, 'free': 1},'Index' : {'value' : -2, 'scale' : 1.0,
                                                         'max' : 5, 'min' : -5, 'free': 1},'Scale' : {'value' : pivot_energy, 'scale' : 1.0,
                                                         'max' : 100000, 'min' : 100, 'free': 0}}, update_source=True)
                     
	fitPL=gta.fit()
	gta.write_roi('fit_model_PL')
	k = np.load('output/fit_model_PL.npy',allow_pickle=True).flat[0]
	src_PL = k['sources']['%s %s' %(name1,name2)]
	TS_PL = '%.2f' %src_PL['ts']
	PL_like=fitPL['loglike']
	PL_fit_quality = fitPL['fit_quality']
	# Extract PL parameters for the summary
	PL_index = src_PL['param_values'][1]        # Photon index
	PL_index_err = src_PL['param_errors'][1]     # Error on index
	PL_flux = src_PL['flux']                      # Photon flux (ph/cm2/s)
	PL_flux_err = src_PL['flux_err']              # Flux error
	PL_eflux = src_PL['eflux']                    # Energy flux (MeV/cm2/s)
	PL_eflux_err = src_PL['eflux_err']            # Energy flux error
	PL_npred = src_PL['npred']                    # Predicted counts
	print("PL fit: Index=%.3f+/-%.3f, Flux=%.2e+/-%.2e, TS=%s, Npred=%.1f, fit_quality=%d" 
	      % (PL_index, PL_index_err, PL_flux, PL_flux_err, TS_PL, PL_npred, PL_fit_quality))
	sed = gta.sed(srcname, loge_bins=log_energy_bins, free_radius=3,make_plots=True)
	os.rename("%s/output/4fgl_j%s_sed.png"%(outdir,srcname[5:16]), "%s/output/sed_PL_%s.png" %(outdir,file_tag))
	os.rename("%s/output/4fgl_j%s_sed.fits"%(outdir,srcname[5:16]), "%s/output/sed_PL_%s.fits" %(outdir,file_tag))
	
	with open('testst_%s.txt' % (file_tag), 'w') as f:
		    f.write(Spectrum3)
		    f.write("\t")
		    f.write(str(fit3_like))
		    f.write("\t")
		    f.write(str(PL_like))
		    f.write("\t")
		    f.write(str(TS_PL))
		    #f.write("\t")
		    #f.write(str(LP_Eb))
	f.close()
	if Spectrum3!='PowerLaw':
	    fit3_TS=2*(fit3_like-PL_like)
	    with open('fit3_testst_%s.txt' % (file_tag), 'w') as f:
		    f.write(Spectrum3)
		    f.write("\t")
		    f.write(str(fit3_like))
		    f.write("\t")
		    f.write(str(PL_like))
		    f.write("\t")
		    f.write(str(fit3_TS))
		    #f.write("\t")
		    #f.write(str(LP_Eb))
	    f.close()
	 
	print("PL_like:"+str(PL_like))
	test_st=[]
	TS2=[]
	BPL_likes=[]
	BPL_Index1=[]
	BPL_Index2=[]
	
	
	#Eb_arraylog=np.logspace(2,4,100)
	Eb_arraylog = np.logspace(np.log10(100000), np.log10(500000), 30) # 100 GeV to 300 GeV, 25 steps
	# --- LEAN BPL LOOP: only fit and record likelihoods ---
	for i in range(len(Eb_arraylog)):
	    model=gta.set_source_spectrum(srcname, spectrum_type='BrokenPowerLaw', spectrum_pars={'Prefactor' : {'value' : 1.0, 'scale' : 1e-12,
                                                         'max' : 1000.0, 'min' : 0.001, 'free': 1},'Index1' : {'value' : -2, 'scale' : 1.0,
                                                         'max' : 5, 'min' : -5, 'free': 1},'BreakValue' : {'value' : Eb_arraylog[i], 'scale' : 1.0,
                                                         'max' : 600000, 'min' : 50000, 'free': 0},'Index2' : {'value' : -1, 'scale' : 1.0,
                                                         'max' : 5, 'min' : -5, 'free': 1}}, update_source=True)
	    BPLfit=gta.fit()
	    gta.write_roi('fit_model_eb_%s'%Eb_arraylog[i])
	    BPL_like=BPLfit['loglike']
	    TS=gta.roi[srcname]['ts']
	    test_st.append(TS)
	    TS2.append(2*(BPL_like-PL_like))
	    BPL_likes.append(BPL_like)
	    # Record the BPL indices at each step
	    bpl_src = gta.roi[srcname]
	    BPL_Index1.append(bpl_src['param_values'][1])
	    BPL_Index2.append(bpl_src['param_values'][2])
	    
	# --- Save TSvsEb profile (with indices for later analysis) ---
	output = '\n'.join('\t'.join(map(str,row)) for row in zip(Eb_arraylog,test_st,TS2,BPL_likes,BPL_Index1,BPL_Index2))
	with open('TSvsEb_%s.txt' % file_tag, 'w') as f:
		f.write('# Eb_MeV\tTS_src\tTS_hardening\tBPL_loglike\tIndex1\tIndex2\n')
		f.write(output)
	f.close()
	
	# --- Find best break energy ---
	ind=np.where(TS2==max(np.array(TS2)))[0]
	best_Eb = Eb_arraylog[ind[0]]
	best_TS_hard = max(TS2)
	print("Maximum value of TS_hardening=%.2f, index=%d, Eb=%.1f MeV"%(best_TS_hard, ind[0], best_Eb))
	
	# --- Estimate E_break confidence interval from TS profile ---
	# 1-sigma: where TS_hardening drops by 1 from maximum (1 d.o.f. on Eb position)
	TS2_arr = np.array(TS2)
	threshold_1sig = best_TS_hard - 1.0
	# Find lower bound
	Eb_lo = Eb_arraylog[0]  # default if no crossing found
	for j in range(ind[0], -1, -1):
	    if TS2_arr[j] < threshold_1sig:
	        Eb_lo = Eb_arraylog[j]
	        break
	# Find upper bound
	Eb_hi = Eb_arraylog[-1]  # default if no crossing found
	for j in range(ind[0], len(TS2_arr)):
	    if TS2_arr[j] < threshold_1sig:
	        Eb_hi = Eb_arraylog[j]
	        break
	Eb_err_lo = best_Eb - Eb_lo   # asymmetric lower error
	Eb_err_hi = Eb_hi - best_Eb   # asymmetric upper error
	print("E_break = %.1f (+%.1f / -%.1f) MeV" % (best_Eb, Eb_err_hi, Eb_err_lo))
	
	# --- Re-fit at best E_break and produce diagnostics ONCE ---
	model=gta.set_source_spectrum(srcname, spectrum_type='BrokenPowerLaw', spectrum_pars={'Prefactor' : {'value' : 1.0, 'scale' : 1e-12,
                                                         'max' : 1000.0, 'min' : 0.001, 'free': 1},'Index1' : {'value' : -2, 'scale' : 1.0,
                                                         'max' : 5, 'min' : -5, 'free': 1},'BreakValue' : {'value' : best_Eb, 'scale' : 1.0,
                                                         'max' : 11000, 'min' : 50, 'free': 0},'Index2' : {'value' : -1, 'scale' : 1.0,
                                                         'max' : 5, 'min' : -5, 'free': 1}}, update_source=True)
	BPLfit_best=gta.fit()
	gta.write_roi('fit_model_BPL_best_%s' % file_tag)
	BPL_best_like = BPLfit_best['loglike']
	BPL_fit_quality = BPLfit_best['fit_quality']
	# Extract BPL best-fit source properties
	src_BPL_best = gta.roi[srcname]
	BPL_npred = src_BPL_best['npred']
	
	# Residual and TS maps only for the best-fit BPL
	gta.residmap('resid_BPL_best_%s' % file_tag, model=model1, make_plots=True)
	gta.tsmap('tsmap_BPL_best_%s' % file_tag, model=model1, make_plots=True)
	
	# SED only for the best-fit BPL
	sed = gta.sed(srcname, loge_bins=log_energy_bins, free_radius=3, make_plots=True)
	os.rename("%s/output/4fgl_j%s_sed.png"%(outdir,srcname[5:16]), "%s/output/sed_BPL_best_%s.png" %(outdir,file_tag))
	os.system("rm -rf %s/output/4fgl_j%s_sedlnl.png"%(outdir,srcname[5:16]))
	os.rename("%s/output/4fgl_j%s_sed.fits"%(outdir,srcname[5:16]), "%s/output/sed_BPL_best_%s.fits" %(outdir,file_tag))
	
	# --- Generate TS_hardening vs E_break profile plot ---
	fig, ax1 = plt.subplots(figsize=(8, 5))
	ax1.plot(Eb_arraylog / 1000.0, TS2_arr, 'b-', linewidth=1.5)
	ax1.axhline(y=12, color='r', linestyle='--', alpha=0.7, label='TS=12 (~3$\\sigma$)')
	ax1.axhline(y=16, color='r', linestyle=':', alpha=0.5, label='TS=16 (~3.5$\\sigma$)')
	ax1.axvline(x=best_Eb / 1000.0, color='k', linestyle='--', alpha=0.5, label='Best E$_{break}$=%.2f GeV' % (best_Eb/1000.0))
	# Shade 1-sigma E_break region
	ax1.axvspan(Eb_lo / 1000.0, Eb_hi / 1000.0, alpha=0.15, color='blue', label='1$\\sigma$ E$_{break}$ range')
	ax1.set_xscale('log')
	ax1.set_xlabel('Break Energy [GeV]', fontsize=13)
	ax1.set_ylabel('TS$_{hardening}$', fontsize=13)
	ax1.set_title('%s   [%s]' % (srcname, peak), fontsize=12)
	ax1.legend(fontsize=9)
	ax1.grid(True, alpha=0.3)
	# Add secondary axis showing BPL indices
	ax2 = ax1.twinx()
	ax2.plot(Eb_arraylog / 1000.0, BPL_Index1, 'g--', alpha=0.5, linewidth=0.8, label='$\\Gamma_1$ (below break)')
	ax2.plot(Eb_arraylog / 1000.0, BPL_Index2, 'm--', alpha=0.5, linewidth=0.8, label='$\\Gamma_2$ (above break)')
	ax2.set_ylabel('Photon Index', fontsize=12, color='gray')
	ax2.legend(fontsize=8, loc='lower right')
	plt.tight_layout()
	plt.savefig('%s/output/TSprofile_%s.png' % (outdir, file_tag), dpi=150)
	plt.savefig('%s/output/TSprofile_%s.pdf' % (outdir, file_tag))
	plt.close()
	
	# --- Determine if this is hardening or softening ---
	best_idx1 = BPL_Index1[ind[0]]
	best_idx2 = BPL_Index2[ind[0]]
	delta_gamma = best_idx1 - best_idx2
	# Hardening: spectrum gets harder (flatter) above the break -> Gamma2 < Gamma1
	# In fermipy convention with negative indices: Index2 > Index1 (less negative = harder)
	is_hardening = (best_idx2 > best_idx1)
	
	# --- Write comprehensive summary file ---
	with open('summary_%s.txt' % file_tag, 'w') as f:
		f.write('# ============================================\n')
		f.write('# Summary for %s\n' % file_tag)
		f.write('# ============================================\n')
		f.write('#\n')
		f.write('# --- Run definition ---\n')
		f.write('# run_id: %s\n' % peak)
		f.write('# MJD_start: %.2f\n' % time1)
		f.write('# MJD_stop: %.2f\n' % time2)
		f.write('# duration_days: %.2f\n' % (time2 - time1))
		f.write('#\n')
		f.write('# --- Power Law fit ---\n')
		f.write('# PL_TS: %s\n' % TS_PL)
		f.write('# PL_Index: %.4f\n' % PL_index)
		f.write('# PL_Index_err: %.4f\n' % PL_index_err)
		f.write('# PL_flux_ph_cm2_s: %.4e\n' % PL_flux)
		f.write('# PL_flux_err: %.4e\n' % PL_flux_err)
		f.write('# PL_eflux_MeV_cm2_s: %.4e\n' % PL_eflux)
		f.write('# PL_eflux_err: %.4e\n' % PL_eflux_err)
		f.write('# PL_Npred: %.1f\n' % PL_npred)
		f.write('# PL_loglike: %.4f\n' % PL_like)
		f.write('# PL_fit_quality: %d\n' % PL_fit_quality)
		f.write('#\n')
		f.write('# --- Broken Power Law (best E_break) ---\n')
		f.write('# Best_Eb_MeV: %.1f\n' % best_Eb)
		f.write('# Best_Eb_GeV: %.3f\n' % (best_Eb / 1000.0))
		f.write('# Eb_err_lo_MeV: %.1f\n' % Eb_err_lo)
		f.write('# Eb_err_hi_MeV: %.1f\n' % Eb_err_hi)
		f.write('# BPL_Index1: %.4f\n' % best_idx1)
		f.write('# BPL_Index2: %.4f\n' % best_idx2)
		f.write('# Delta_Gamma: %.4f\n' % delta_gamma)
		f.write('# BPL_loglike: %.4f\n' % BPL_best_like)
		f.write('# BPL_Npred: %.1f\n' % BPL_npred)
		f.write('# BPL_fit_quality: %d\n' % BPL_fit_quality)
		f.write('#\n')
		f.write('# --- Spectral hardening test ---\n')
		f.write('# TS_hardening: %.4f\n' % best_TS_hard)
		f.write('# is_hardening: %s\n' % str(is_hardening))
		f.write('# TS_hardening_gt12: %s\n' % str(best_TS_hard >= 12.0))
		f.write('# TS_hardening_gt16: %s\n' % str(best_TS_hard >= 16.0))
		f.write('#\n')
		f.write('# --- Quality flags ---\n')
		f.write('# PL_TS_gt25: %s\n' % str(float(TS_PL) >= 25.0))
		f.write('# PL_converged: %s\n' % str(PL_fit_quality >= 3))
		f.write('# BPL_converged: %s\n' % str(BPL_fit_quality >= 3))
	f.close()
	print("Summary written to summary_%s.txt" % file_tag)
	
	# --- Cleanup: remove per-step BPL roi files (keep only best) ---
	for k, val in enumerate(Eb_arraylog):
	    if val!=best_Eb:
	        os.system('rm -rf %s/output/fit_model_eb_%s*'%(outdir,val))
	  
	
	
	#os.system('rm -rf lt.fits output/*.par PFILES* galdiff.fits isodiff.txt output/srcmap_00.fits output/srcmdl_00.xml output/ltcube_00.fits output/ft1_00.fits output/bexpmap_00.fits output/bexpmap_roi_00.fits output/ccube.fits output/ccube_00.fits output/bexpmap_roi_00.fits output/find_sourcefind_00_pointsource* *.log output/fit_model_*.fits output/srcmap*.fits output/find_sourcefind_00_pointsource* *.log output/srcmap*.fits  output/srcmap_*.fits output/fit_model_1.npy output/fit_model_2.npy output/fit_model_3.npy')
	os.system('rm -rf lt.fits output/*.par PFILES* galdiff.fits isodiff.txt output/srcmap_*.fits output/srcmdl_*.xml output/ltcube_*.fits output/ft1_*.fits output/ccube.fits output/ccube_*.fits output/bexpmap_*.fits')
	
	os.chdir('../')
	
	# Copy to permanent storage: results/SH_Alba/source/run_id/
	final_dest = '%s/%s/%s' % (copy_dir, srcname, peak)
	os.makedirs(os.path.dirname(final_dest), exist_ok=True)
	if os.path.exists(final_dest):
	    shutil.rmtree(final_dest)
	shutil.copytree(outdir, final_dest)
	#shutil.rmtree(outdir)
	return
################################################
if __name__=="__main__":
	main(sys.argv)
