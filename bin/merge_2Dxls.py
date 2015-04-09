#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os
import RNA_Comp.module_Matrix   as m_mat
import cPickle         as pickle
import numpy           as np
from optparse import OptionParser
import logging
logging.basicConfig(level=logging.INFO,format='%(levelname)-5s @ %(asctime)s: %(message)s ',stream=sys.stderr)

class ReadFile( m_mat.Matrix_info ):
   def __init__( self, infile,inf_column=1, in_dtype="int",header=1 ):
      super( ReadFile,self ).__init__( infile, inf_column, in_dtype,header )
      logging.info( 'Reading file ==> %s' % ( infile ) )

class CompFPKM(ReadFile):
   def __init__( self, l_files ):
      self.l_files  = l_files
      self.m_files = []
      self.all_gene = []
      self.all_samp = []
      for file in self.l_files:
         m_file = ReadFile( file ,inf_column=1, in_dtype="float",header=1  )
         m_file.load_mat()
         self.m_files.append( m_file )
         self.all_gene +=  m_file.rowname
         self.all_samp +=  m_file.colname
         logging.info( 'Reading %s Done. %d genes and %d samples.' % ( file,  len( m_file.rowname ) ,len( m_file.colname )  ) )


   def common_gene( self ):
      self.s_gene_common = sorted(list( set( self.all_gene ) ))
      self.s_gene_common = np.array( list(self.s_gene_common),dtype="string" )
   
      logging.info( '%d shared genes detected.' % ( len(self.s_gene_common) ) )
   
      print  "Gene\t%s" % "\t".join( self.all_samp )
   
      for gene in sorted(self.s_gene_common):
         l_val = []
         for i,file in enumerate(self.l_files):
            if gene in self.m_files[i].rowname: 
               j      = list( self.m_files[i].rowname ).index( gene )
               l_val += list( self.m_files[i].matrix[j,:] )
            else:
               l_val += list( np.zeros( len(self.all_samp) ) )
         
         print "%s\t%s" % ( gene, "\t".join( np.array(l_val,dtype="string") ) )
   

def prepare_optparser():
   usage ="""usage: %s [options] 

   Using -h or --help for more information

Example:
   python %s a/genes.fpkm_table b/genes.fpkm_table >out.xls
   
   """ % (sys.argv[0],sys.argv[0])

   description = " The methylation analysis standard pipeline. "

   optparser = OptionParser(version="%s v0.1 20150317" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser
   
def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   
   if len(args)<1:
      prepare_optparser().print_help()
      sys.exit(1)
      
   try:
      l_files = args

   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)

   m_comp = CompFPKM( l_files )
   m_comp.common_gene()

if __name__ == '__main__':
   main()
