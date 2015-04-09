from __future__ import division
import re,sys,os
import numpy as np
from optparse   import OptionParser

def prepare_optparser():
   usage ="""usage: %s [options] 

   Using -h or --help for more information

Example:
   python %s sample_v2_20150228.xls merge.FPKM.xls reorder_sample.xls
   
   """ % (sys.argv[0],sys.argv[0])
   
   description = " Task 1.3"

   optparser = OptionParser(version="%s v0.2 20150305" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser

def get_convert_idx( seq1,seq2 ):
   seq1 = np.array( seq1,dtype="string" )
   seq2 = np.array( seq2,dtype="string" )
   cons = []
   for word in seq1:
      if word in seq2:
         cons.append(  np.where( seq2 == word )[0][0] )
   return np.array( cons,dtype="int" )

def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      file_name_in = args[0]
      file_FPKM_in = args[1]
      file_order   = args[2]
      
   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)
      
   f_file_name_in  = open( file_name_in ,"r" )
   f_file_FPKM_in  = open( file_FPKM_in ,"r" )
   f_file_order    = open( file_order   ,"r" )
   
   file_name_out   =  "%s.%s.xls" % ( '.'.join( file_name_in.split('.')[:-1] ),  file_order.split('.')[0] )
   file_FPKM_out   =  "%s.%s.xls" % ( '.'.join( file_FPKM_in.split('.')[:-1] ),  file_order.split('.')[0] )
   
   f_file_name_out = open( file_name_out,"w" )
   f_file_FPKM_out = open( file_FPKM_out,"w" )
   
   l_order = []
   l_samp  = []
   l_brief = []
   
   '''
      Reading the order
   '''
   for line in f_file_order:
      line = line.strip()
      l_order.append( line )
   f_file_order.close()

   '''
      Reading names.
   '''
   head = f_file_name_in.readline()
   head_n = head.strip()
   l_lines = []
   for line in f_file_name_in:
      line = line.strip()
      f    = line.split()
      l_lines.append( line )
      l_samp.append( f[1] )
      l_brief.append(f[9] )
   f_file_name_in.close()
      
   '''
      Reading RPKM matrix
   '''
   head   = f_file_FPKM_in.readline()
   l_head = head.split()

   l_samout= []
   l_order = np.array( l_order,dtype="string" )
   l_brief = np.array( l_brief,dtype="string" )
   l_samp  = np.array( l_samp ,dtype="string" )
   l_cons  = get_convert_idx( l_order,l_brief )
   
#   print >>f_file_FPKM_out, "%s\t%s" %  ( l_head[0],"\t".join( l_samp[l_cons] ) )
   print >>f_file_FPKM_out, "%s\t%s" %  ( l_head[0],"\t".join( l_brief[l_cons] ) )

   for line in f_file_FPKM_in:
      line    = line.strip()
      f       = line.split()
      l_FPKM  = f[1:]
      np_FPKM = np.array( l_FPKM,dtype="string" )
      print >>f_file_FPKM_out, "%s\t%s" %  ( f[0],"\t".join( np_FPKM[l_cons] ) )
   
   print >>f_file_name_out, head_n
   np_lines= np.array( l_lines,dtype="string" )

   for line in np_lines[l_cons]:
      line = line.strip()
      print >>f_file_name_out, line
   
   f_file_name_in.close()
   f_file_FPKM_in.close()

   f_file_name_out.close()
   f_file_FPKM_out.close()
   
      
if __name__ == '__main__':
   main()