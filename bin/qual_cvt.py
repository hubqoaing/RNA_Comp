#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os,gzip
import RNA_Comp.module_Matrix   as m_mat
import cPickle         as pickle
import numpy           as np
from optparse import OptionParser
import logging
logging.basicConfig(level=logging.INFO,format='%(levelname)-5s @ %(asctime)s: %(message)s ',stream=sys.stderr)


class ReadFastq( object ):
   def __init__( self, infile,outfile, in_min, out_min ):
      self.infile = infile
      self.outfile= outfile
      logging.info( 'Reading file ==> %s' % ( self.infile ) )
      
      self.f_in   = gzip.open( infile ,"rb" )
      logging.info( 'Writing file ==> %s' % ( self.outfile ) )
      self.f_out  = gzip.open( outfile,"wb" )
      self.in_min = in_min
      self.out_min= out_min
      
   def loadID( self ):
      while 1:
         self.__get_one_record()
         if not self.line1:
            break
      self.f_in.close()
      logging.info( 'Reading file ==> %s done.' % ( self.infile ) )
      self.f_out.close()
      logging.info( 'Writing file ==> %s done.' % ( self.outfile ) )
   
   def __get_one_record(self):
      line1 = self.f_in.readline()
      line2 = self.f_in.readline()
      line3 = self.f_in.readline()
      line4 = self.f_in.readline()
      self.line1 = line1
      line1 = line1.strip()
      line2 = line2.strip()
      line3 = line3.strip()
      line4 = line4.strip()
      
      if len(line1)>0:
         line4 = self.__change_phred( line4 )
         print >>self.f_out, line1
         print >>self.f_out, line2
         print >>self.f_out, line3
         print >>self.f_out, line4
      
   def __change_phred(self,line4):
      l_line4 = np.array( list(line4), dtype="string" )
      idx     = (l_line4=="!") + (l_line4=='"')
      l_line4[ idx ] = ";"
      
      l_char = []
      for i in l_line4:
         val = ord(i)-self.in_min+self.out_min
         if val <64:
            val=64
         l_char.append( chr(val) )

#      print line4,idx,char
      return "".join( l_char )


def prepare_optparser():
   usage ="""usage: %s [options] 

   Using -h or --help for more information

Example:
   python %s /datd/huboqiang/test_hisat_tophat.nature07509/SRR015286_1.fastq.gz /datd/huboqiang/test_hisat_tophat.nature07509/test.fq.gz 59 64
   
   """ % (sys.argv[0],sys.argv[0])

   description = " The methylation analysis standard pipeline. "

   optparser = OptionParser(version="%s v0.1 20150317" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser
   
def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      in_gz  = args[0]
      out_gz = args[1]
      in_min  = int( args[2] )
      out_min = int( args[3] )

   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)

   m_fq_cvt = ReadFastq( in_gz,out_gz, in_min, out_min )
   m_fq_cvt.loadID()
   

if __name__ == '__main__':
   main()
