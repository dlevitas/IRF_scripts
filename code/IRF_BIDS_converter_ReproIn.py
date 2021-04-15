#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:52:42 2021

IRF BIDS Conversion tool, with the assumption being
that all protocol acquisitions follow the ReproIn naming convention.

Beginning 2/2021, all new IRF protocols are being renamed to ReproIn.

Following conversion, QA is performed to ensure things are good. 

@author: dlevitas
"""

import os, re, glob, json, shutil, smtplib
import pandas as pd
from datetime import datetime
from os.path import basename
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

prisma_dir = '/N/dcwan/projects/irf/Prisma'
base_dir = '/N/dcwan/projects/irf/BIDS'

# Begin:

# Import/download additional python packages
try:
    import quickshear
except:
    os.system('pip install quickshear --user')
    import quickshear

try:
    import pydicom
except ImportError:
    os.system('pip install -U pydicom --user')
    import pydicom
    
try:
    import heudiconv
except:
    print('Installing heudiconv')
    print('')
    os.system('pip install heudiconv --user')
    import heudiconv
    
try:
    import nipype
except:
    print('Installing nipype')
    print('')
    os.system('pip install nipype')
    
if nipype.__version__ != '1.6.0':
    os.system('pip install --upgrade nipype --user')
    

if shutil.which('dcm2niix') is None:
    os.environ["PATH"] += os.pathsep + os.pathsep.join(['/N/dcwan/projects/irf/software_packages/dcm2niix/build/bin'])
if shutil.which('pigz') is None:
    os.environ["PATH"] += os.pathsep + os.pathsep.join(['/N/dcwan/projects/irf/software_packages/pigz-2.4'])

#Get date/time info
date = datetime.now()
year = str(date.year)

if date.month < 10:
    month = '0' + str(date.month)
else:
    month = str(date.month)
if date.day < 10:
    day = '0' + str(date.day)
else:
    day = str(date.day)

date = year + month + day

# Find scans by date (year-month-day)
today_scans = [x for x in os.listdir(prisma_dir) if date in x and 'ReproIn' in x]

# Complete BIDS conversion for each scan today
for scan in today_scans:
    # Determine study name, subjectID, and sessionID (when applicable)
    study_name = re.sub('[^A-Za-z0-9]+', '', scan.split(date)[-1].split('_sub')[0])
    subID = scan.split('sub-')[-1].split('_ses')[0]
    if 'ses' in scan:
        sesID = scan.split('ses-')[-1].split('.')[0]
    else:
        sesID = 'N/A'
    
    bids_root_dir = '{}/{}'.format(base_dir, study_name)

    if not os.path.isdir(bids_root_dir):
        os.mkdir(bids_root_dir)
    
    # Create error log file
    if not os.path.isfile('{}/error_log.txt'.format(bids_root_dir)):
        os.system('touch {}/error_log.txt'.format(bids_root_dir))
        log = open('{}/error_log.txt'.format(bids_root_dir), "r+")
        log.write('###################################################\n')
        if sesID == 'N/A':
            scan_info = 'study {}, subject {}'.format(study_name, subID)
            log.write('Error log for {}\n'.format(scan_info))
        else:
            scan_info = 'study {}, subject {}, session {}'.format(study_name, subID, sesID)
            log.write('Error log for {}\n'.format(scan_info))
        log.write('###################################################\n')
        log_num = 1

    
    # Perform BIDS conversion via Heudiconv/ReproIn combo
    if sesID == 'N/A':
        if not os.path.isdir('{}/sub-{}'.format(bids_root_dir, subID)):
            convert = 'yes'
            print('')
            print('Performing BIDS conversion using Heudiconv/ReproIn pairing for {}'.format(scan_info))
            print('-------------------------------------------------------------------------------------')
            os.system('heudiconv -f reproin --bids --subjects {} --locator {} -o {} --files {}/{}'.format(subID, bids_root_dir, bids_root_dir, prisma_dir, scan))
        else:
            convert = 'no'
    else:
        if not os.path.isdir('{}/sub-{}/ses-{}'.format(bids_root_dir, subID, sesID)):
            convert = 'yes'
            print('')
            print('Performing BIDS conversion using Heudiconv/ReproIn pairing for'.format(scan_info))
            print('-------------------------------------------------------------------------------------')
            os.system('heudiconv -f reproin --bids --subjects {} --ses {} --locator {} -o {} --files {}/{}'.format(subID, sesID, bids_root_dir, bids_root_dir, prisma_dir, scan))
        else:
            convert = 'no'
    
    if convert == 'yes':
    
        if sesID == 'N/A':
            path = '{}/sub-{}'.format(bids_root_dir, subID)
            converted_files = pd.read_csv('{}/sub-{}_scans.tsv'.format(path, subID), sep='\t')['filename'].tolist()
            converted_files = ['{}/'.format(path) + x for x in converted_files]
    
        else:
            path = '{}/sub-{}/ses-{}'.format(bids_root_dir, subID, sesID)
            converted_files = pd.read_csv('{}/sub-{}_ses-{}_scans.tsv'.format(path, subID, sesID), sep='\t')['filename'].tolist()
            converted_files = ['{}/'.format(path) + x for x in converted_files]
            
        os.system('chmod -R 777 {}/fmap/*.json'.format(path))
        os.system('chmod -R 777 {}/anat/*.nii.gz'.format(path))
            
        # Get ordered protocol list
        unique_dicom_series = glob.glob('{}/{}/*_*_000001.dcm'.format(prisma_dir, scan))
        unique_dicom_series = sorted(unique_dicom_series)
        protocol_list = []
        for acquisition in unique_dicom_series:
            protocol_list.append([int(pydicom.dcmread(acquisition).SeriesNumber), pydicom.dcmread(acquisition).SeriesDescription])
    
        json_files = [x.split('.nii.gz')[0] + '.json' for x in converted_files]
        data_check = [[json.load(open(x))['SeriesNumber'],json.load(open(x))['SeriesDescription']] for x in json_files]
        [value.append(converted_files[index].split('.nii.gz')[0] + '.json') for index,value in enumerate(data_check)]
        
        for x in range(len(protocol_list)):
            if protocol_list[x][0] not in [y[0] for y in data_check]:
                data_check.insert(protocol_list[x][0]-1, protocol_list[x])
        data_check.sort()
        
        section_indices = [index for index, value in enumerate(data_check) if 'scout' in value[1] or 'localizer' in value[1]]
        
        if len(section_indices) == 0:
            section_indices = [0]
        
        if 0 not in section_indices:
            section_indices.insert(0,0)
            
        for j in range(len(section_indices)):
            section_start = section_indices[j]
            
            try:
                section_end = section_indices[j+1]
            except:
                section_end = len(protocol_list)
        
            section = data_check[section_start:section_end]
            
            fmap_se_dwi = [value[-1] for index,value in enumerate(section) if 'fmap' in value[1] and 'dwi' in value[1] and '_dup' not in value[1]]
            fmap_se_func = [value[-1] for index,value in enumerate(section) if 'fmap_dir' in value[1] and 'dwi' not in value[1] and '_dup' not in value[1]]
            fmap_magphase = [value[-1] for index,value in enumerate(section) if 'fmap' in value[1] and 'dwi' not in value[1] and 'dir' not in value[1] and '_dup' not in value[1]]
    
            if fmap_se_dwi:
                fmap_se_dwi_ped = [json.load(open(x))['PhaseEncodingDirection'] for x in fmap_se_dwi]
                dwi = [value[-1] for index,value in enumerate(section) if 'fmap' not in value[1] and 'dwi' in value[1] and '_dup' not in value[1]]
                dwi_ped = [json.load(open(x))['PhaseEncodingDirection'] for x in dwi]
                if not len(dwi):
                    log.write('\n')
                    log.write('{}. There is no DWI acquisition associated with the DWI fmap. Will remove DWI fmap acquisition from BIDS output\n'.format(log_num))
                    log_num += 1
                    os.system('rm -rf {}*'.format(fmap_se_dwi[0].split('.json')[0]))
                elif fmap_se_dwi_ped[0][0] != dwi_ped[0][0] or fmap_se_dwi_ped[0] == dwi_ped[0]:
                    log.write('\n')
                    log.write('{}. Improper DWI fmap PED for associated DWI. Will remove DWI fmap acquisition from BIDS output\n'.format(log_num))
                    log_num += 1
                    os.system('rm -rf {}*'.format(fmap_se_dwi[0].split('.json')[0]))
                else:                
                    intended_for = [x.split('/sub-{}/'.format(subID))[-1].split('.json')[0] + '.nii.gz' for x in dwi]
                    fmap_json = json.load(open(fmap_se_dwi[0]))
                    fmap_json.update({"IntendedFor": intended_for})
                    with open(fmap_se_dwi[0], 'w') as fp:
                        json.dump(fmap_json, fp, indent=3)
    
            if fmap_se_func:
                fmap_se_func_ped = [json.load(open(x))['PhaseEncodingDirection'] for x in fmap_se_func]
                func = [value[-1] for index,value in enumerate(section) if 'func-bold' in value[1] and '_dup' not in value[1]]
                if not len(func):
                    log.write('\n')
                    log.write('{}. There are no func-bold acquisition associated with the spin-echo fmaps. Will remove spin-echo fmap acquisitions from BIDS output\n'.format(log_num))
                    log_num += 1
                    for bad in fmap_se_func:
                        os.system('rm -rf {}*'.format(bad.split('.json')[0]))
                        
                elif len(fmap_se_func) == 1:
                    log.write('\n')
                    log.write('{}. Only one spin-echo fmap acquisition found, must have two. Will remove spin-echo fmap acquisition from BIDS output\n'.format(log_num))
                    log_num += 1
                    os.system('rm -rf {}*'.format(fmap_se_func[0].split('.json')[0]))
                elif fmap_se_func_ped[0][0] != fmap_se_func_ped[1][0] or fmap_se_func_ped[0] == fmap_se_func_ped[1]:
                    log.write('\n')
                    log.write('{}. PEDs of spin-echo fmaps are not flipped. Will remove spin-echo fmap acquisitions from BIDS output\n'.format(log_num))
                    log_num += 1
                    for bad in fmap_se_func:
                        os.system('rm -rf {}*'.format(bad.split('.json')[0]))
                else:
                    intended_for = [x.split('/sub-{}/'.format(subID))[-1].split('.json')[0] + '.nii.gz' for x in func]
                    for i in fmap_se_func:
                        fmap_json = json.load(open(i))
                        fmap_json.update({"IntendedFor": intended_for})
                        with open(i, 'w') as fp:
                            json.dump(fmap_json, fp, indent=3)
                        
            if fmap_magphase:
                if len(fmap_magphase) == 1:
                    log.write('\n')
                    log.write('{}. Only one magnitude or phasediff image found, must have 2 or 3. Will remove magnitude/phasediff acquisition from BIDS output\n'.format(log_num))
                    log_num += 1
                    os.system('rm -rf {}*'.format(fmap_magphase[0].split('.json')[0]))
                elif 'magnitude' not in fmap_magphase[0] and 'phasediff' not in fmap_magphase[-1]:
                    log.write('\n')
                    log.write('{}. Magntude & Phasediff pair not found. Will remove other acquisition(s) from BIDS output\n'.format(log_num))
                    log_num += 1
                    for bad in fmap_magphase:
                        os.system('rm -rf {}*'.format(bad.split('.json')[0]))                
                else:
                    intended_for = [x.split('/sub-{}/'.format(subID))[-1].split('.json')[0] + '.nii.gz' for x in func]
                    for i in fmap_se_func:
                        fmap_json = json.load(open(i))
                        fmap_json.update({"IntendedFor": intended_for})
                        with open(i, 'w') as fp:
                            json.dump(fmap_json, fp, indent=3)    
                                   
        log.close()
        # Send email to Hu if fmap issues have occurred. These errors suggest a hardware issue
        if any('1.' in x for x in open('{}/error_log.txt'.format(bids_root_dir), "r+").readlines()):
            def send_mail(send_from, send_to, subject, text, file=None, server="127.0.0.1"):
                msg = MIMEMultipart()
                msg['From'] = send_from
                msg['To'] = send_to
                msg['Subject'] = subject
            
                msg.attach(MIMEText(text))
            
                with open(file, "r+") as fil:
                    part = MIMEApplication(fil.read(), Name=basename(file))
                # After the file is closed
                part['Content-Disposition'] = 'attachment; filename="%s"' % basename(file)
                msg.attach(part)
            
                smtp = smtplib.SMTP(server)
                smtp.sendmail(send_from, send_to, msg.as_string())
                smtp.close()
                
            send_mail('dlevitas@iu.edu', 'dlevitas@iu.edu', 'IRF BIDS conversion issue for {}'.format(scan_info), 'Hi Hu,\nThere appears to be a field map issue, see attachment', file='{}/error_log.txt'.format(bids_root_dir))
        os.remove('{}/error_log.txt'.format(bids_root_dir))
        
        # Remove duplicate acquisition files
        os.system('rm -rf {}/*/*_dup*'.format(path))
    
        # Remove files with _part- label, which is likely an indication that ReproIn thinks there's something weird
        os.system('rm -rf {}/*/*_part-*'.format(path))
        
        # Remove miscellaneous files generated by Heudiconv/ReproIn that aren't BIDS-compliant
        remove_files = glob.glob('{}/task-*.json'.format(bids_root_dir))
        remove_files = remove_files + glob.glob('{}/*scans*'.format(path))
        
        for remove in remove_files:
            os.remove(remove)  
        
        # Remove stuff in sourcedata
        os.system('rm -rf {}/sourcedata/*'.format(bids_root_dir))
        
        # Deface anatomical images
        anats = glob.glob('{}/anat/*.nii.gz'.format(path))
        for anat in anats:
            print('Defacing {}'.format(anat))
            print('')
            os.system('/N/dcwan/projects/irf/software_packages/ROBEX/runROBEX.sh {} {}'.format(anat, anat.split('.nii.gz')[0] + '_mask.nii.gz'))
            os.system('quickshear {} {} {}'.format(anat, anat.split('.nii.gz')[0] + '_mask.nii.gz', anat))
            os.system('rm -rf {}'.format(anat.split('.nii.gz')[0] + '_mask.nii.gz'))
                 
        # Adjust participants.tsv file
        os.remove('{}/participants.tsv'.format(bids_root_dir))
        subjects = [x for x in os.listdir(bids_root_dir) if 'sub-' in x]
        df = pd.DataFrame(subjects ,columns=['participant_id'])
        df.to_csv('{}/participants.tsv'.format(bids_root_dir), sep='\t', index=False, header=True)
        

print("Finished converting today's {} scans: {}".format(date, today_scans))
    
        





