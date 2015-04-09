#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os
import subprocess
import numpy   as np
import cPickle as pickle
import scipy   as sp
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import stats 
import module_running_jobs as m_jobs

class Map_From_sra(dict):
   
   def __init__(self,samp_info, genome_file,anno_file, dir_name,sftw_name ):
      self['sam_info'] = {}
      self['sample'] = []
      self['infile'] = {'info_file':samp_info, 'genome_file':genome_file,'anno_file':anno_file }
      self['stage'] = { 'name':[] }
      self['dir_name']  = dir_name
      self['sftw_name'] = sftw_name
      
   def load_samp(self):
      
      self['sam_info']['samp_brief'] = {}
      self['sam_info']['type']     = {}
      self['sam_info']['stage']    = {}
      self['sam_info']['dilute']   = {}
      self['sam_info']['stage_sam']= {}
      self['sam_info']['RFP_mols'] = {}
      self['sam_info']['GFP_mols'] = {}
      self['sam_info']['CRE_mols'] = {}
      self['sam_info']['stage_sam']= {}
      self['sam_info']['data_type']= {}
      
      info_file = self['infile']['info_file']
      file = open(info_file,"r")
      in_h = file.readline()
      for line in file:
         line = line.strip('\n')
         f = line.split()
         samp       = f[0]
         brief_name = f[1]
         stage      = f[2]
         ltype      = f[3]
         ERCC_dilute= float( f[4] )
         RFP_mols   = float( f[5] )
         GFP_mols   = float( f[6] )
         CRE_mols   = float( f[7] )
         data_type  = f[8]  # PE or SE

         self['sample'].append( samp )
         self['sam_info']['samp_brief'][ samp ] = brief_name
         self['sam_info']['type'][  samp ]      = ltype
         self['sam_info']['stage'][ samp ]      = stage
         self['sam_info']['dilute'][samp ]      = ERCC_dilute
         self['sam_info']['RFP_mols'][ samp ]   = RFP_mols
         self['sam_info']['GFP_mols'][ samp ]   = GFP_mols
         self['sam_info']['CRE_mols'][ samp ]   = CRE_mols
         self['sam_info']['data_type'][samp ]   = data_type
         
         if stage not in self['stage']['name']:
            self['stage']['name'].append( stage )
         if stage not in self['sam_info']['stage_sam']:
            self['sam_info']['stage_sam'][stage] = []
            
         self['sam_info']['stage_sam'][stage].append( samp )
      file.close()
 
   def SRA2fastq(self):
      home_dir    = os.path.abspath('./')
      raw_dir     = self['dir_name']['raw_data']
      fq_dir      = self['dir_name']['fastq_data']
      
      if not os.path.isdir( fq_dir ):
         os.mkdir( fq_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      
      fqDump       = self['sftw_name'].fastqDump
      
      sh_file       = "%s/scripts/s01.SRA2Fastq.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/s01.SRA2Fastq_work.sh" % (home_dir)
      
      sh_info = """
samp_name=$1
fqDump=$2
raw_dir=$3
fq_dir=$4

$fqDump --split-files --gzip --outdir $fq_dir/${samp_name} $raw_dir/${samp_name}.sra 

mv $fq_dir/${samp_name}/${samp_name}_1.fastq.gz $fq_dir/${samp_name}/${samp_name}.1.fq.gz && \\
mv $fq_dir/${samp_name}/${samp_name}_2.fastq.gz $fq_dir/${samp_name}/${samp_name}.2.fq.gz
      """
      
      sh_work = ""
      for samp_name in self['sample']:
         sh_work  += " sh %s  %s %s %s %s\n" % ( sh_file,  samp_name,fqDump,  raw_dir,fq_dir  )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )
      
      
   
   def run_tophat_tran(self):
      home_dir     = os.path.abspath('./')
      
      fq_dir       = self['dir_name']['fastq_data']
      tophat_dir   = self['dir_name']['tophat_dir']
      
      if not os.path.isdir(  tophat_dir):
         os.mkdir( tophat_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      
      tophat_py    = self['sftw_name'].tophat
      samtools_exe = self['sftw_name'].samtools
      
      sh_file      = "%s/s02.1.tophatTrans.sh"      % (script_dir)
      sh_work_file = "%s/s02.1.tophatTrans_work.sh" % (script_dir)
      
      sh_info = """
tophat_py=$1
fq_dir=$2
samp_name=$3
brief_name=$4
tophat_dir=$5
genome=$6
gtf_file=$7
samtools_exe=$8

$tophat_py  \\
   -p 8 -G $gtf_file                                                 \\
   --library-type fr-unstranded    --phred64-quals                   \\
   -o $tophat_dir/$brief_name                                        \\
   $genome                                                           \\
   $fq_dir/$samp_name/$samp_name.1.fq.gz $fq_dir/$samp_name/$samp_name.2.fq.gz
      """ 
      sh_work = ""
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         sh_work += "sh %s  %s %s %s %s  %s %s %s %s\n" % ( sh_file, tophat_py, fq_dir, samp, brief_name, tophat_dir, self['infile']['genome_file'],self['infile']['anno_file'],samtools_exe  )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=6 )
#      my_job.running_SGE( vf="7g",maxjob=100 )

      
   def run_tophat_mannual(self):
      home_dir     = os.path.abspath('./')
      
      fq_dir       = self['dir_name']['fastq_data']
      tophat_dir   = self['dir_name']['tophat_mannual_dir']
      
      if not os.path.isdir(  tophat_dir):
         os.mkdir( tophat_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      
      tophat_py    = self['sftw_name'].tophat
      samtools_exe = self['sftw_name'].samtools
      
      sh_file      = "%s/s02.2.tophatMannual.sh"      % (script_dir)
      sh_work_file = "%s/s02.2.tophatMannual_work.sh" % (script_dir)
      
      sh_info = """
tophat_py=$1
fq_dir=$2
samp_name=$3
brief_name=$4
tophat_dir=$5
genome=$6
gtf_file=$7
samtools_exe=$8

$tophat_py  \\
   -p 8                                                              \\
   --read-edit-dist 3                                                \\
   --read-realign-edit-dist 3                                        \\
   --phred64-quals                                                   \\
   -o $tophat_dir/$brief_name                                        \\
   $genome                                                           \\
   $fq_dir/$samp_name/$samp_name.1.fq.gz $fq_dir/$samp_name/$samp_name.2.fq.gz
      """ 
      sh_work = ""
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         sh_work += "sh %s  %s %s %s %s  %s %s %s %s\n" % ( sh_file, tophat_py, fq_dir, samp, brief_name, tophat_dir, self['infile']['genome_file'],self['infile']['anno_file'],samtools_exe  )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=6 )
#      my_job.running_SGE( vf="7g",maxjob=100 )

      


   def run_hisat_tran(self):
      home_dir     = os.path.abspath('./')
      
      fq_dir       = self['dir_name']['fastq_data']
      hisat_dir    = self['dir_name']['hisat_dir']
      
      if not os.path.isdir(  hisat_dir):
         os.mkdir( hisat_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      
      hisat        = self['sftw_name'].hisat
      samtools_exe = self['sftw_name'].samtools
      
      sh_file      = "%s/s02.3.hisat.sh"      % (script_dir)
      sh_work_file = "%s/s02.3.hisat_work.sh" % (script_dir)
      
      sh_info = """
hisat=$1
fq_dir=$2
samp_name=$3
brief_name=$4
hisat_dir=$5
genome=$6
splice_file=$7
samtools_exe=$8

$hisat         -p 8 -x $genome --phred64        \\
   -1 $fq_dir/$samp_name/$samp_name.1.fq.gz     \\
   -2 $fq_dir/$samp_name/$samp_name.2.fq.gz     \\
   -S /dev/stdout                               \\
   --known-splicesite-infile $splice_file       \\
   2>$hisat_dir/$brief_name/log                |\\
awk '{if($1 ~ /^@/) print $0; else{ for(i=1;i<=NF;i++) if($i!~/^XS/) printf("%s\\t",$i);else XS0=$i;  XS1=((and($2, 0x10) && and($2, 0x40)) || (and($2,0x80) && !and($2,0x10)))?"XS:A:+":"XS:A:-"; print XS1 } }' | $samtools_exe view -Sb -q 1 - >$hisat_dir/$brief_name/accepted_hits.raw.bam &&\\ 
$samtools_exe sort -m 2000000000 $hisat_dir/$brief_name/accepted_hits.raw.bam $hisat_dir/$brief_name/accepted_hits
      """ 
      sh_work = ""
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         if not os.path.isdir( "%s/%s" % (hisat_dir,brief_name) ):
            os.mkdir( "%s/%s" % (hisat_dir,brief_name) )
         genome     = "%s.hisat"       % (self['infile']['genome_file'])
         splice_file= "%s.spliceSites" % ( ".".join( self['infile']['anno_file'].split(".")[:-1] )  )
         sh_work += "sh %s  %s %s %s %s  %s %s %s %s\n" % ( sh_file, hisat, fq_dir, samp, brief_name, hisat_dir, genome, splice_file, samtools_exe  )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=6 )
#      my_job.running_SGE( vf="7g",maxjob=100 )


   def run_hisat_mannual(self):
      home_dir     = os.path.abspath('./')
      
      fq_dir       = self['dir_name']['fastq_data']
      hisat_dir    = self['dir_name']['hisat_mannual_dir']
      
      if not os.path.isdir(  hisat_dir):
         os.mkdir( hisat_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      
      hisat        = self['sftw_name'].hisat
      samtools_exe = self['sftw_name'].samtools
      
      sh_file      = "%s/s02.4.hisatMannual.sh"      % (script_dir)
      sh_work_file = "%s/s02.4.hisatMannual_work.sh" % (script_dir)
      
      sh_info = """
hisat=$1
fq_dir=$2
samp_name=$3
brief_name=$4
hisat_dir=$5
genome=$6
splice_file=$7
samtools_exe=$8

$hisat         -p 8 -x $genome --phred64        \\
   -1 $fq_dir/$samp_name/$samp_name.1.fq.gz     \\
   -2 $fq_dir/$samp_name/$samp_name.2.fq.gz     \\
   -S /dev/null                                 \\
   --novel-splicesite-outfile $splice_file      \\
   2>$hisat_dir/$brief_name/log              && \\
$hisat         -p 8 -x $genome --phred64        \\
   -1 $fq_dir/$samp_name/$samp_name.1.fq.gz     \\
   -2 $fq_dir/$samp_name/$samp_name.2.fq.gz     \\
   -S /dev/stdout                               \\
   --novel-splicesite-infile  $splice_file      \\
   2>$hisat_dir/$brief_name/log.2              |\\
awk '{if($1 ~ /^@/) print $0; else{ for(i=1;i<=NF;i++) if($i!~/^XS/) printf("%s\\t",$i);else XS0=$i;  XS1=((and($2, 0x10) && and($2, 0x40)) || (and($2,0x80) && !and($2,0x10)))?"XS:A:+":"XS:A:-"; print XS1 } }' | $samtools_exe view -Sb -q 1 - >$hisat_dir/$brief_name/accepted_hits.raw.bam &&\\ 
$samtools_exe sort -m 2000000000 $hisat_dir/$brief_name/accepted_hits.raw.bam $hisat_dir/$brief_name/accepted_hits
      """ 
      sh_work = ""
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         if not os.path.isdir( "%s/%s" % (hisat_dir,brief_name) ):
            os.mkdir( "%s/%s" % (hisat_dir,brief_name) )
         genome     = "%s.hisat"          % ( self['infile']['genome_file'] )
         splice_file= "%s/%s.spliceSites" % ( hisat_dir,  brief_name        )
         sh_work += "sh %s  %s %s %s %s  %s %s %s %s\n" % ( sh_file, hisat, fq_dir, samp, brief_name, hisat_dir, genome, splice_file, samtools_exe  )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=6 )
#      my_job.running_SGE( vf="7g",maxjob=100 )
