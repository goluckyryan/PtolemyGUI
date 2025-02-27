# To run Fresco, 

./fresco < infile > outfile

I have fresco in other folder, I set the PATH to include the execute of fresco

# haha.in

## Elastic scattering

60Ni(p,p)60Ni@30.00 MeV/u Elastic
NAMELIST
&FRESCO  hcm=0.1  rmatch=60 
         jtmin=0.0  jtmax=50 absend=0.001 
         thmin=0.0  thmax=180.0 thinc=1.0 
         chans=1  smats=2 xstabl=1 
         elab=30.0  /
&PARTITION  namep='p' massp=1.00 zp=1 namet='60Ni' masst=60 zt=28 qval=0.000  nex=1 /
&STATES   jp=0.5  ptyp=1  ep=0.000 jt=0.0  ptyt=1  et=0.000 cpot=1/
&partition /
&POT kp=1        ap=0  at=60  rc=1.258 /
&POT kp=1 type=1 p1=47.937 p2=1.20 p3=0.669 p4=2.853  p5=1.20 p6=0.669 /
&POT kp=1 type=2 p1=0.0    p2=1.20 p3=0.669 p4=6.878  p5=1.28 p6=0.550 /
&POT kp=1 type=3 p1=5.250  p2=1.02 p3=0.590 p4=-0.162 p5=1.02 p6=0.590 /
&pot /
&overlap /
&coupling /

In the first line of &POT, ap and at are used to calculate the radius. use ap = 0 to match Ptolemy r0target for zero-range approximation, all potential inputs are same as Ptolemy

in &POT, type=1 WS, type=2 surface, type=3 SO.