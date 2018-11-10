from AUTOclui import *
import subprocess
def calculate_initialsolution(system,nconf,nbody,datafile):
  "calculate a initial condition from the initial solution given on file" 
  file=open("nbp_constants.h","w")
  file.write("  INTEGER,PARAMETER::nconf="+`nconf`+"\n")
  file.write("  INTEGER,PARAMETER::nbody="+`nbody`+"\n")
  file.close()
  subprocess.check_call(["cp",str(datafile),"initial_condition.h"]) 
  subprocess.check_call(["gfortran","-c","nbp.f90"])
  subprocess.check_call(["gfortran","-o","orbcalc.exe","orbcalc.f90","nbp.f90"])
  subprocess.check_call(["./orbcalc.exe"])
  subprocess.check_call(["cp","nbp.dat",system+".dat"])
  userdata(system)
  dl(system)
  subprocess.check_call(["mv","s.dat","s."+system])
  return
