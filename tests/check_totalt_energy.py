#!/usr/bin/env python

if __name__ == '__main__':
  fname = 'input.00.log'
  fref  = 'input.ref.log'
  
  for line in open(fname):
    if "final tot" in line:
      print line.split("=")[1]
      Etot = float(line.split("=")[1].split()[0])
  for line in open(fref):
    if "final tot" in line:
      Eref = float(line.split("=")[1].split()[0])

  delta = Etot - Eref
  print "Reference TOTAL ENERGY: %15.8f"%Eref
  print "Current   TOTAL ENERGY: %15.8f"%Etot
  print "deviation             : %15.8f"%(Etot-Eref)

  if(abs(Etot-Eref) < 1.0e-6):
    print "Test status: pass"
  else:
    print "Test status: fail"

