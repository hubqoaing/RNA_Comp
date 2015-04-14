#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os
import RNA_Comp.module_Matrix   as m_mat
import cPickle         as pickle
import numpy           as np
from optparse import OptionParser
import pysam
import logging
logging.basicConfig(level=logging.INFO,format='%(levelname)-5s @ %(asctime)s: %(message)s ',stream=sys.stderr)

def cigar_len(cigar):
   l0 = 0
   for pairs in cigar:
      ltype, leng = pairs
      if ltype == 0 or ltype == 3:
         l0 += leng
   return l0

class ReadFile( object ):
   def __init__( self, infile ):
      self.infile = infile
      logging.info( 'Reading file ==> %s' % ( infile ) )
      self.f_bam = pysam.Samfile( infile,"rb" )
      self.id    = ""
      self.cache = {}
      
   def loadID( self, id_ref ):
      comp_idx   = cmp( self.id,id_ref )
      if comp_idx == -1:
         try:
            '''  self.id == ""  '''
            read = self.f_bam.next()
            self.__load_cache( read )
            self.id = read.qname
         except:
            '''  The last self.id, rest refID were all above that.   '''
            a = 1
            
      else:
         while comp_idx == 0:
            try:
               read = self.f_bam.next()
               self.__load_cache( read )
               self.id  = read.qname
            except:  
               '''all reads were loaded.''' 
               comp_idx = 1
               break
            
            comp_idx = cmp( self.id,id_ref )
   
   def __load_cache( self, read ):
      chrom  = self.f_bam.header['SQ'][read.rname]['SN']
      chrpos = "*"
      qname  = read.qname
      if chrom != "*":
         beg = read.pos + 1
         cigar = read.cigar
         end   = read.pos + 1 + cigar_len(cigar)
         chrpos = "%s:%d-%d" % (chrom,beg,end)
      
      gene = "-"
      tags = dict(read.tags)
      if 'XF' in tags:  # not fusion junctions
         gene = tags['XF']
      
      if qname not in self.cache:
         self.cache[qname] = { 'seq':"", 'chrpos':[], 'gene':[] }
      
      self.cache[qname]['seq']    = read.seq
      self.cache[qname]['chrpos'].append( chrpos )
      self.cache[qname]['gene'  ].append(  gene  )

class ReadsMapping(ReadFile):
   def __init__( self, l_files, id_list ):
      self.l_files  = l_files
      self.id_list  = id_list
      self.bam_files= []
      self.all_gene = []
      self.all_samp = []
      
   def load_idList(self):
      f_id_list = open( self.id_list,"r" )
      self.id_list = f_id_list.readlines()
      self.id_list = [ ID.strip('\n') for ID in self.id_list ]
      f_id_list.close()
   
   def load_bamInfo(self):
      for file in self.l_files:
         m_file = ReadFile( file )
         self.bam_files.append( m_file )
         logging.info( 'File %s were initiated.' % ( file ) )
      
      for id_ref in self.id_list:
         l_NH     = []
         l_chrpos = []
         l_gene   = []
         seqs     = "-"
         for i,file in enumerate(self.l_files):
            self.bam_files[i].loadID( id_ref )
            if id_ref in self.bam_files[i].cache:
               NH  = len( self.bam_files[i].cache[ id_ref ]['chrpos'] )
               l_NH.append( str(NH) )
               if NH > 0:
                  l_chrpos.append( ",".join( [ self.bam_files[i].cache[ id_ref ]['chrpos'][j] for j in range(NH)  ]) )
                  l_gene.append(   ",".join( [ self.bam_files[i].cache[ id_ref ]['gene'][j]   for j in range(NH)  ]) )
                  seqs = self.bam_files[i].cache[ id_ref ]['seq']
            else:
               l_NH.append( "0" )
               l_chrpos.append( "-" )
               l_gene.append(   "-" )

         dif = len( set(l_chrpos) )
         if dif > 1:
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ( id_ref.split('.')[0],id_ref.split('.')[1],id_ref.split('.')[1], seqs, ",".join(l_NH),  "\t".join(l_chrpos), "\t".join(l_gene) )
         
         for i,file in enumerate(self.l_files):
            if id_ref in self.bam_files[i].cache:
               del self.bam_files[i].cache[ id_ref ]
   

def prepare_optparser():
   usage ="""usage: %s [options] 

   Using -h or --help for more information

Example:
   python %s /datd/huboqiang/test_hisat_tophat/test/list_1000 /datd/huboqiang/test_hisat_tophat/test/test.bam
   
   """ % (sys.argv[0],sys.argv[0])

   description = " The methylation analysis standard pipeline. "

   optparser = OptionParser(version="%s v0.1 20150317" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser
   
def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   
   if len(args)<2:
      prepare_optparser().print_help()
      sys.exit(1)
      
   try:
      listFile= args[0]
      l_files = args[1:]

   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)

   m_comp = ReadsMapping( l_files,listFile )
   m_comp.load_idList()
   m_comp.load_bamInfo()

if __name__ == '__main__':
   main()
