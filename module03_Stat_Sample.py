from __future__ import division
import re,sys,os
import module_StatInfo as Stat
import scipy
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import module_CountInfo  as m_cnt
import module_GTFFeature   as m_gtf
import module_Matrix       as m_mat
import module_running_jobs as m_jobs
import module_RepeatCount   as m_repcnt

class SampStat(dict):
   def __init__(self,samp_info,ercc_info,genome_gtf,dir_name):
      self['sam_info'] = {}
      self['sample'] = []
      self['infile'] = {'info_file':samp_info,'ercc_info':ercc_info,'anno_file':genome_gtf}
      self['stage'] = { 'name':[] }
      self['dir_name'] = dir_name
      self.ERCC_info   = { 'len':{}, 'mol':{} }
      self.__load_name()
      self.__load_ERCC_info()
      
   def __load_ERCC_info(self):
      f_file = open( self['infile']['ercc_info'],"r")
      f_file.readline()
      for line in f_file:
         line = line.strip('\n')
         f    = line.split()
         ERCC_id = f[0]
         self.ERCC_info['len'][ERCC_id] = int(f[1])
         self.ERCC_info['mol'][ERCC_id] = float(f[2])
      f_file.close()
      
   def __load_name(self):
      self.tophat       = self['dir_name']['bam_dir']
      self.HTS          = self['dir_name']['HTSeq_result_dir']
      self.HTS_k        = self['dir_name']['HTSeq_known_dir']
      self.HTS_u        = self['dir_name']['HTSeq_unknown_dir']
      self.cufflink_u   = self['dir_name']['cufflinks_unknown_dir']
      self.cuffquant    = self['dir_name']['cuffquant_dir']
      self.cuffnorm     = self['dir_name']['cuffnorm_dir']
      self.cuffquant_ercc=self['dir_name']['cuffquant_ERCC_dir']
      self.cuffnorm_ercc =self['dir_name']['cuffnorm_ERCC_dir']
      self.repeatCount  = self['dir_name']['repeat_counts_dir']
      self.repeat_mrg   = self['dir_name']['repeat_mrg_dir']
      
      home_dir          = os.path.abspath('./')
      self.script_dir   = "%s/scripts" % (home_dir)
      self.bin_dir      = "%s/bin"     % (home_dir)
      self.data_dir     = "%s/Database"% (home_dir)
      self.statInfo     = "%s/StatInfo"% (home_dir)
      if not os.path.isdir( self.statInfo ):
         os.mkdir( self.statInfo )

   def load_samp(self):      
      self['sam_info']['samp_brief'] = {}
      self['sam_info']['brief_samp'] = {}
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
         self['sam_info']['brief_samp'][ brief_name ] = samp
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
   
   def __get_HTS_reads(self, HTS_info, samp):
      idx = self['sample'].index( samp )
      return HTS_info.sam_tot_reads[idx]
      
   def __get_HTS_clean_split(self):
      sh_file        =  "%s/p.HTSeq_split.sh"       %  (self.script_dir)
      sh_work_file   =  "%s/p.HTSeq_split_work.sh"  %  (self.script_dir)
      sh_info = """
infile=$1
out_Refseq=$2
out_NONCODE=$3
out_NSMB=$4
inNeo=$5
outNeo=$6

grep -v -P '^NONHSAG|XLOC_' $infile >$out_Refseq
head -n 1 $infile >$out_NONCODE && grep -P '^NONHSAG' $infile >>$out_NONCODE
head -n 1 $infile >$out_NSMB    && grep -P '^XLOC'    $infile >>$out_NSMB

head -n 1 $inNeo >$outNeo
for i in `cut -f 1 %s/novo_lnc_raw.combined.FPKM0.5_rep0.25.multiExon.genlen | uniq`;do grep -w $i $inNeo ;done >>$outNeo
      """ % ( self.data_dir )
      
      infile      = "%s/merge.dexseq_clean.gene.xls"          % ( self.HTS )
      out_Refseq  = "%s/merge.dexseq_clean_refseq.gene.xls"   % ( self.HTS )
      out_NONCODE = "%s/merge.dexseq_clean_NONCODE.gene.xls"  % ( self.HTS )
      out_NSMB    = "%s/merge.dexseq_clean_NSMB.gene.xls"     % ( self.HTS )
      inNeo       = "%s/merge.dexseq_NeoRaw.gene.xls"         % ( self.HTS )
      outNeo      = "%s/merge.dexseq_NeoPass.gene.xls"        % ( self.HTS )
      sh_work = "sh %s  %s %s %s %s %s %s" % ( sh_file,infile,out_Refseq,out_NONCODE,out_NSMB,inNeo,outNeo )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=1 )
      
   def Basic_Stat(self):
      """
      Stat for QC, Tophat mapping, ERCC RGC count
      """
      out_file   = "%s/01.BasicInfo_QC_map_SpikeIn.xls" % (self.statInfo)
      f_out_file = open( out_file,"w" )
      out_info   = "Sample\tBrief_samp\t"    +\
                   "Pre_Map_Reads\tAligned_Reads\tHTSseq_Known_Reads\t"    +\
                   "HTSeq_Refseq_Reads\tHTSeq_NONCODE_V4_Reads\tHTSeq_Nsmb_Reads\tHTSeq_Neo_Reads\t" +\
                   "RFP_Reads\tGFP_Reads\tCRE_Reads\tERCC_Reads\t"   +\
                   "RFP_Mols\tGFP_Mols\tCRE_Mols\tERCC_Mols"
                   
      print >>f_out_file, out_info

      l_brief_samp = [ self['sam_info']['samp_brief'][samp] for samp in self['sample'] ]

      '''
         Load refseq reads
      '''
      HTS_info = m_cnt.CountInfo( self.HTS_k,l_brief_samp,"dexseq_clean",self.HTS )
      HTS_info.load_mat()
      HTS_info.sam_tot_reads()
      
      self.__get_HTS_clean_split()
      
      Refseq_info = m_cnt.CountInfo( self.HTS_k,l_brief_samp,"dexseq_clean_refseq",self.HTS )
      Refseq_info.load_mat()
      Refseq_info.sam_tot_reads()

      NONCODE_info = m_cnt.CountInfo( self.HTS_k,l_brief_samp,"dexseq_clean_NONCODE",self.HTS )
      NONCODE_info.load_mat()
      NONCODE_info.sam_tot_reads()
      
      NSMB_info = m_cnt.CountInfo( self.HTS_k,l_brief_samp,"dexseq_clean_NSMB",self.HTS )
      NSMB_info.load_mat()
      NSMB_info.sam_tot_reads()

      NeoPass_info = m_cnt.CountInfo( self.HTS_u,l_brief_samp,"dexseq_NeoPass",self.HTS )
      NeoPass_info.load_mat()
      NeoPass_info.sam_tot_reads()
      
      Cnt_Info = m_cnt.CountInfo(  self.HTS_k, l_brief_samp, "dexseq_ERCC_RGCPloyA", self.HTS  )
      Cnt_Info.generate_mat()
      
      '''
         Load other information
      '''
      for idx,samp in enumerate(self['sample']):
         brief_name  = self['sam_info']['samp_brief'][samp]
         Tophat_log  =  "%s/%s/align_summary.txt"           % ( self.tophat,brief_name )
         HTSeq_SpikeIn= "%s/%s/%s.dexseq_ERCC_RGCPloyA.txt" % ( self.HTS_k ,brief_name, brief_name )
         
         print brief_name
         
         SpikeIn_info= Stat.SpikeIn( HTSeq_SpikeIn, self['infile']['ercc_info'] )
         SpikeIn_info.load_HTS_file()
                  
         if os.path.isfile( Tophat_log ):
            MapStat_info= Stat.TophatStat( Tophat_log )
            MapStat_info.read_infile()
         else:
            Hisat_log   = "%s/%s/log" % ( self.tophat,brief_name )
            MapStat_info= Stat.HisatStat( Hisat_log )
            MapStat_info.read_infile()

         
         pre_map_read = MapStat_info['statInfo']['totalRead']
         aligned_read = MapStat_info['statInfo']['mappedRead']
         
         if self['sam_info']['data_type'][samp ] == "PE":         
            HTSseq_read  = self.__get_HTS_reads( HTS_info    ,samp ) * 2
            Refseq_read  = self.__get_HTS_reads( Refseq_info ,samp ) * 2
            NONCODE_read = self.__get_HTS_reads( NONCODE_info,samp ) * 2
            NSMB_read    = self.__get_HTS_reads( NSMB_info   ,samp ) * 2
            NeoPass_read = self.__get_HTS_reads( NeoPass_info,samp ) * 2
            read_RFP   =  SpikeIn_info.RGC_count['RGC-mRFP'] * 2
            read_GFP   =  SpikeIn_info.RGC_count['RGC-GFP' ] * 2
            read_CRE   =  SpikeIn_info.RGC_count['RGC-CRE' ] * 2
            read_ERCC  =  SpikeIn_info.ERCC_total            * 2
         else:
            HTSseq_read  = self.__get_HTS_reads( HTS_info    ,samp )
            Refseq_read  = self.__get_HTS_reads( Refseq_info ,samp )
            NONCODE_read = self.__get_HTS_reads( NONCODE_info,samp )
            NSMB_read    = self.__get_HTS_reads( NSMB_info   ,samp )
            NeoPass_read = self.__get_HTS_reads( NeoPass_info,samp )             
            read_RFP   =  SpikeIn_info.RGC_count['RGC-mRFP']
            read_GFP   =  SpikeIn_info.RGC_count['RGC-GFP' ]
            read_CRE   =  SpikeIn_info.RGC_count['RGC-CRE' ]
            read_ERCC  =  SpikeIn_info.ERCC_total
            
         
         mol_RFP  = self['sam_info']['RFP_mols'][ samp ]
         mol_GFP  = self['sam_info']['GFP_mols'][ samp ]
         mol_CRE  = self['sam_info']['CRE_mols'][ samp ]
         mol_ERCC = self['sam_info']['dilute'][   samp ] * 6.023*10**10
         
         out_info =  "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%1.2e\t%1.2e\t%1.2e\t%1.2e"   \
            % ( samp, brief_name,            \
                pre_map_read, aligned_read, HTSseq_read,                         \
                Refseq_read,NONCODE_read,NSMB_read,NeoPass_read,                 \
                read_RFP, read_GFP,read_CRE, read_ERCC,                          \
                mol_RFP , mol_GFP ,mol_CRE , mol_ERCC )
         print  >>f_out_file, out_info
      f_out_file.close()
      
      
   def Repeat_Stat(self):

      l_brief_samp = [ self['sam_info']['samp_brief'][samp] for samp in self['sample'] ]
     
      l_align_reads = []
      for brief_samp in l_brief_samp:
         Tophat_log  =  "%s/%s/align_summary.txt"           % ( self.tophat,brief_samp )

         if os.path.isfile( Tophat_log ):
            MapStat_info= Stat.TophatStat( Tophat_log )
            MapStat_info.read_infile()
         else:
            Hisat_log   = "%s/%s/log" % ( self.tophat,brief_samp )
            MapStat_info= Stat.HisatStat( Hisat_log )
            MapStat_info.read_infile()

         l_align_reads.append( MapStat_info['statInfo']['mappedRead'] )
         
      np_align_reads = np.array( l_align_reads )
     
      RepCnt = m_repcnt.RepeatCount( self.repeatCount,l_brief_samp, self.repeat_mrg )
      RepCnt.generate_mat()
      RepCnt.element_group_subgroup_FPKM_sum(np_align_reads)