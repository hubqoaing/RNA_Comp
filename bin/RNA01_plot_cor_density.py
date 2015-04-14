from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import scipy.stats
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch

class matrix_cor(dict):
   def __init__(self,infile):
      self['infile'] = infile
   
   def load_matrix(self):
      f_infile = open( self['infile'],"r" )
      in_h = f_infile.readline()
      in_h = in_h.strip('\n')
      f_h  = in_h.split()
      self['sample'] = f_h[1:]
      
      l_mat = []
      for line in f_infile:
         line = line.strip('\n')
         f    = line.split()
         l_mat.append(  np.array(f[1:],dtype=float)  )
      f_infile.close()
         
      self['raw_matrix'] = np.array( l_mat,dtype=float )
   
   def get_samp_tag(self):
      self['l_stage'] = []
      self['l_tissue']= []
      self['l_ltype'] = []

      dict_stage = {}
      dict_tissue = {}
      dict_ltype = {}
      
      id_stage = 0
      id_tissue= 0
      id_ltype = 0
      
      for samp in self['sample']:
         fs = samp.split('_')
         print fs
         stage  = fs[0]
         tissue = fs[0]
         ltype  = fs[1]
         if stage  not in self['l_stage']:
            dict_stage[  stage ]  = id_stage
            id_stage             += 1
         if tissue not in self['l_tissue']:
            dict_tissue[ tissue]  = id_tissue
            id_tissue            += 1
         if ltype  not in dict_ltype:
            dict_ltype[ ltype  ]  = id_ltype
            id_ltype             += 1
      
         self['l_stage'].append(  stage  )
         self['l_tissue'].append( tissue )
         self['l_ltype'].append(  ltype  )
         
      ## Initiate index for each tissue 
   
      l_stage_idx  = []
      l_tissue_idx = []
      l_ltype_idx  = []
   
      y_pos = []
      lab_p = []
      y_lab = []
   
      for stage  in self['l_stage']:
         l_stage_idx.append(  dict_stage[  stage  ] )
      for tissue in self['l_tissue']:
         l_tissue_idx.append( dict_tissue[ tissue ] )   
      for i,ltype  in enumerate( self['l_ltype'] ):
         if dict_ltype[  ltype  ] not in l_ltype_idx:
            y_pos.append( i-0.5 )
            y_lab.append( ltype )
         l_ltype_idx.append(  dict_ltype[  ltype  ] )
   
      np_ltype = np.array( self['l_ltype'],dtype='string' )
   
      self['l_ltype_less'] = []
      for ltype in self['l_ltype']:
         if ltype not in self['l_ltype_less']:
            self['l_ltype_less'].append( ltype )
   
      for i,ltype  in enumerate(self['l_ltype_less']):
         lab_p.append( y_pos[i]-0.5+len( np_ltype[ np_ltype == ltype ] )/2 )
      print lab_p 
      print self['l_ltype_less']
      
      self['stage_idx']  = np.array(l_stage_idx ,dtype=float)
      self['tissue_idx'] = np.array(l_tissue_idx,dtype=float)
      self['ltype_idx']  = np.array(l_ltype_idx ,dtype=float)

      self['lab_p'] = lab_p
      self['y_pos'] = y_pos
      self['y_lab'] = y_lab

   def get_cor_matrix(self):
      np_corr = np.ones(  [len(self['sample']),len(self['sample'])]   )
      for i in range( 0,len(self['sample']) ):
         for j in range( i+1,len(self['sample']) ):
            idx1 = (self['raw_matrix'][:,i] > 0.1)
            idx2 = (self['raw_matrix'][:,j] > 0.1)
            idx  = idx1*idx2
#            print >>sys.stderr, "Calculating %d and %d" % ( i,j )
            val_cor,p    = scipy.stats.pearsonr( np.log10(self['raw_matrix'][idx,i]),np.log10(self['raw_matrix'][idx,j]) )
            print val_cor
            np_corr[i,j] = val_cor
            np_corr[j,i] = val_cor
      self['cor_mat'] = np_corr

   def output_cor_matrix(self,out_cor_file):
      f_out_cor_file = open(out_cor_file,"w")
      out = "Sample\t%s" % ( "\t".join(self['sample']) )
      print    >>f_out_cor_file, out
      for i,samp in enumerate(self['sample']):
         out = "%s\t%s" % ( samp, "\t".join( np.array( self['cor_mat'][i,:] ,dtype="string") ) )
         print >>f_out_cor_file, out
      f_out_cor_file.close()
      
   def plot_cor_density(self,out_cor_file):
      """docstring for plot_cor_density"""
      l_fpkm_sample = np.array( self['sample'],dtype="string" )
      n = len( l_fpkm_sample )
      
      dpi = 150
      
      fig = plt.figure(figsize=(0.95*5*n,0.95*5*n))
      
      for i,sam1    in enumerate( l_fpkm_sample ):
         for j,sam2 in enumerate( l_fpkm_sample ):
            
            print >>sys.stderr, "Calculating %s and %s" % ( sam1,sam2)
            
            ax = fig.add_subplot( len(l_fpkm_sample),len(l_fpkm_sample),i*len(l_fpkm_sample)+j+1 )
            if i==j:
               ax.plot( [0,10],[0,10],visible=False )
#               ax.text( 0.5,5,"%s" % (sam1) )
               for sp in ax.spines.values():
                  sp.set_visible(False)
               ax.get_xaxis().set_ticks( [] )
               ax.get_yaxis().set_ticks( [] )
               ax.get_yaxis().set_ticklabels( [] )
               ax.get_xaxis().set_ticklabels( [] )
            elif i<j:
#               dat1 = self['raw_matrix'][:,i]
#               dat2 = self['raw_matrix'][:,j]
#               data1 = dat1[ (dat1>0.0001)+(dat2>0.0001) ]
#               data2 = dat2[ (dat1>0.0001)+(dat2>0.0001) ]
#               ax.hist2d( np.log10( data1+0.01), np.log10( data2+0.01), bins=50,cmap="Purples")
               ax.plot( np.log10(self['raw_matrix'][:,i]+1), np.log10(self['raw_matrix'][:,j]+1),'.',color="#3857A2")
               ax.set_xlim(-1,7)
               ax.set_ylim(-1,7)
            elif i>j:
               ax.plot( [0,10],[0,10],visible=False )
               ax.text( 2,5,"r = %1.3f" % ( self['cor_mat'][j,i] ),fontsize=30 )
               ax.get_xaxis().set_ticks( [] )
               ax.get_yaxis().set_ticks( [] )
               ax.get_yaxis().set_ticklabels( [] )
               ax.get_xaxis().set_ticklabels( [] )
               

      plt.savefig('%s.correlation.density.png' % (out_cor_file),dpi = dpi)
         

   def plot_cor_matrix(self,out_cor_file):
      def __load_cor_matrix(self,out_cor_file):
         f_infile = open( out_cor_file,"r" )
         in_h = f_infile.readline()
         in_h = in_h.strip('\n')
         f_h  = in_h.split()
         self['sample'] = f_h[1:]
      
         l_mat = []
         for line in f_infile:
            line = line.strip('\n')
            f    = line.split()
            l_mat.append(  np.array(f[1:],dtype=float)  )
         f_infile.close()
         self['cor_mat'] = np.array( l_mat,dtype=float )

      def __plot_ltype(self, figure):
         """docstring for plot_ltype"""
         ax01 = figure.add_axes([ 0.25, 0.91, 0.5, 0.02])   
         ax01.get_xaxis().set_ticks( [] )
         ax01.get_yaxis().set_ticks( [] )
         ax01.get_yaxis().set_ticklabels( [] )
         ax01.get_xaxis().set_ticklabels( [] )   
         for i,pos in enumerate( self['lab_p'] ):
            ax01.text( pos,-2,"%s" % ( self['l_ltype_less'][i] ),rotation=90,ha = "right",size=15 )
         ax01.tick_params(colors="white")
         ax01.grid(True,color='white',linestyle='-',linewidth=3)
         for sp in ax01.spines.values():
            sp.set_visible(False)
         cax1 = ax01.imshow( np.array( [ self['ltype_idx'] ] ), aspect='auto', cmap='rainbow',interpolation='nearest')
   
      def __plot_stage(self, figure):
         """docstring for plot_stage"""
         ax1 = figure.add_axes([ 0.22, 0.3  ,0.02, 0.6 ]  )
         ax1.get_xaxis().set_ticks( [] )
         ax1.get_yaxis().set_ticks( [] )
         ax1.get_yaxis().set_ticklabels( [] )
         ax1.get_xaxis().set_ticklabels( [] )
         ax1.tick_params(colors="white")
         ax1.grid(True,color='white',linestyle='-',linewidth=3)
         for sp in ax1.spines.values():
            sp.set_visible(False)
         cax1 = ax1.imshow( np.array( [ self['stage_idx'] ] ).T, aspect='auto', cmap='rainbow',interpolation='nearest')
   
      def __plot_tissue(self,figure):
         """docstring for plot_tissue"""
         ax2 = figure.add_axes([ 0.25, 0.275, 0.5, 0.02])
         matrix_lab2 = []
         ax2.get_xaxis().set_ticks( [] )
         ax2.get_yaxis().set_ticks( [] )
         ax2.get_yaxis().set_ticklabels( [] )
         ax2.get_xaxis().set_ticklabels( [] )
         ax2.tick_params(colors="white")
         ax2.grid(True,color='white',linestyle='-',linewidth=3)
         l_fpkm_sample = np.array( self['sample'],dtype='string')
         for i in  range( 0,len(l_fpkm_sample) ) :
            ax2.text( i+0.5,2,'%s' % ( l_fpkm_sample[i]),rotation=270,ha = "right",size=15)
         for sp in ax2.spines.values():
            sp.set_visible(False)
         cax2 = ax2.imshow( np.array( [  self['tissue_idx']  ] ), aspect='auto', cmap='rainbow',interpolation='nearest')
   
      def __plot_cor_matrix(self, figure):
         """docstring for plot_cor_matrix"""
         axmatrix = figure.add_axes([ 0.25, 0.3, 0.5, 0.6])
         axmatrix.grid(True,color='black',linestyle='-',linewidth=3)
#         cax1 = axmatrix.imshow(dataMatrixOrdered[ den['leaves'] ][ den['leaves'] ], aspect='auto', cmap='bwr',interpolation='nearest')
         cax1 = axmatrix.imshow(self['cor_mat'], aspect='auto', cmap='jet',interpolation='nearest')
         axmatrix.get_xaxis().set_ticks( self['y_pos'] )
         axmatrix.get_yaxis().set_ticks( self['y_pos'] )
         axmatrix.get_yaxis().set_ticklabels([])
         axmatrix.get_xaxis().set_ticklabels([])
         l_fpkm_sample = np.array( self['sample'],dtype='string')
         for i in range( 0,len(l_fpkm_sample) ):
            axmatrix.text( len(l_fpkm_sample)+0.5,i,'%s' % (l_fpkm_sample[i]),ha = "left",size=15)
         cbaxes = figure.add_axes([0.05,0.87,0.08,0.08])
         cbar = plt.colorbar(cax1,cax=cbaxes)
         cbar.ax.yaxis.set_ticks_position('left')
         cbar.ax.yaxis.set_label_position('left')
         cbar.set_label('Correlation', size=17,rotation=90)
         cbar.ax.tick_params(labelsize=10)

      """  sub main function """
      if 'cor_mat' not in self:
         __load_cor_matrix(self,out_cor_file)
         self.get_samp_tag()
         
      fig = plt.figure(figsize=(24,24))
      __plot_ltype(      self,fig )
      __plot_stage(      self,fig )
      __plot_tissue(     self,fig )
      __plot_cor_matrix( self,fig ) 
      plt.savefig('%s.pdf' % (out_cor_file),format='pdf')
         
 

      
   



def show_help():
   print >>sys.stderr,"\n\tpython",sys.argv[0],"merge.chip-level.1kb.xls "
   
def main():
   try:
      infile = sys.argv[1]
   except IndexError:
      show_help()
      sys.exit(1)

   out_file = "%s.correlation.pearson.py.xls" % (infile)

   samp_info = matrix_cor( infile )

   samp_info.load_matrix()
   samp_info.get_samp_tag()
   samp_info.get_cor_matrix()
   samp_info.output_cor_matrix(out_file)
   samp_info.plot_cor_density(out_file)
   
   samp_info.plot_cor_matrix(out_file)






if __name__ == '__main__':
   main()