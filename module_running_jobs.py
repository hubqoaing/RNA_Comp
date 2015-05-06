from __future__ import division
import re,sys,os
import subprocess,time

class running_jobs(object):
   def __init__( self, sh_file,sh_work_file ):
      self.sh_file      = sh_file
      self.sh_work_file = sh_work_file
      path              = os.path.realpath(__file__)
      self.sge          = "%s/bin/qsub-sge.pl"      % ( "/".join( path.split('/')[:-1] )  )
      self.multi        = "%s/bin/multi-process.pl" % ( "/".join( path.split('/')[:-1] )  )
   
   def load_sh_file(self,sh_scripts):
      f_out = open( self.sh_file,"w" )
      print >>f_out, sh_scripts
      f_out.close()
   
   def load_sh_work_file(self,sh_scripts):
      f_out = open( self.sh_work_file,"w" )
      print >>f_out, sh_scripts
      f_out.close()
   
   def running_SGE(self,vf,maxjob=100):
      shell_work = 'perl %s --resource vf=%s  --maxjob %d %s ' % ( self.sge, vf, maxjob,  self.sh_work_file)
      p = subprocess.Popen(shell_work,shell='True')
      while 1:
         run_cnt = 0
         if p.poll() is None:
            run_cnt += 1
            time.sleep(3)
         if run_cnt == 0:
            break
      
   def running_multi(self,cpu=8):
      log_file1 = "%s.%s.out" % (  self.sh_work_file, time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime())  )
      log_file2 = "%s.%s.err" % (  self.sh_work_file, time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime())  )
      shell_work = 'perl %s -cpu %d  %s >%s 2>%s' % ( self.multi, cpu,  self.sh_work_file, log_file1, log_file2 )
      p = subprocess.Popen(shell_work,shell='True')
      while 1:
         run_cnt = 0
         if p.poll() is None:
            run_cnt += 1
            time.sleep(3)
         if run_cnt == 0:
            break
      
