#!/usr/bin/python3

import pandas as pd
import urllib.request
import re
import numpy

me = 0.51099895000 # +- 15
mp = 938.27208816 # +- 29
mn = 939.56542052 # +- 54
ma = 3727.37915 

class Isotope:
  def __int__(self):

    self.livechart = "https://nds.iaea.org/relnsd/v0/data?"
    self.data = None

    # Read the saved CSV file back into a DataFrame
    try :
      self.data = pd.read_csv('IAEA_NuclearData.csv')
    except FileNotFoundError:
      # the service URL
      url = self.livechart + "fields=ground_states&nuclides=all"
      self.data = self.lc_read_csv(url)
      self.data.insert(0, 'A', self.data['z'] + self.data['n'])
      self.data.to_csv('IAEA_NuclearData.csv', index=False)
      self.data = pd.read_csv('IAEA_NuclearData.csv')

  def lc_read_csv(self, url):    
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
    return pd.read_csv(urllib.request.urlopen(req))


  # print(haha.columns)
  ## ['A', 'z', 'n', 'symbol', 'radius', 'unc_r', 'abundance', 'unc_a',
  ##        'energy_shift', 'energy', 'unc_e', 'ripl_shift', 'jp', 'half_life',
  ##        'operator_hl', 'unc_hl', 'unit_hl', 'half_life_sec', 'unc_hls',
  ##        'decay_1', 'decay_1_%', 'unc_1', 'decay_2', 'decay_2_%', 'unc_2',
  ##        'decay_3', 'decay_3_%', 'unc_3', 'isospin', 'magnetic_dipole', 'unc_md',
  ##        'electric_quadrupole', 'unc_eq', 'qbm', 'unc_qb', 'qbm_n', 'unc_qbmn',
  ##        'qa', 'unc_qa', 'qec', 'unc_qec', 'sn', 'unc_sn', 'sp', 'unc_sp',
  ##        'binding', 'unc_ba', 'atomic_mass', 'unc_am', 'massexcess', 'unc_me',
  ##        'me_systematics', 'discovery', 'ENSDFpublicationcut-off',
  ##        'ENSDFauthors', 'Extraction_date']


  def GetExList(self, ASym : str, maxEx : float):
    try:
      exList = self.lc_read_csv(self.livechart + "fields=levels&nuclides=" + ASym)
      exJpi = exList[['energy', 'jp']]
      exJpi = exJpi[exJpi['energy'] < (maxEx * 1000)]
      return exJpi
    except:
      return pd.DataFrame()

  def BreakDownName(self, ASym : str):
    match =  re.match(r'(\d+)(\D+)', ASym)
    return [int(match.group(1)), match.group(2) ]

  def GetAZ(self, ASym : str):
    [A, sym] = self.BreakDownName(ASym)
    try:
      dudu = self.data[(self.data['symbol']==sym) & (self.data['A']==A)]
      Z  = int(dudu['z'].iloc[0])
      return [A, Z]
    except:
      return [A, numpy.nan]

  def GetBindingPerA(self, ASym : str) -> float: 
    [A, sym] = self.BreakDownName(ASym)
    try:
      dudu = self.data[(self.data['symbol']==sym) & (self.data['A']==A)]
      Z  = int(dudu['z'].iloc[0])
      N  = A - Z
      return dudu['binding'].iloc[0]/1000
    except:
      return numpy.nan

  def GetMassFromSym(self, ASym : str) -> float:
    [A, sym] = self.BreakDownName(ASym)
    try:
      dudu = self.data[(self.data['symbol']==sym) & (self.data['A']==A)]
      Z  = int(dudu['z'].iloc[0])
      N  = A - Z
      binding = dudu['binding'].iloc[0]/1000
      return Z*mp + N*mn - binding*A
    except:
      return numpy.nan

  def GetMassFromAZ(self, A : int, Z : int) -> float:
    try:
      dudu = self.data[(self.data['z']==Z) & (self.data['A']==A)]
      Z  = int(dudu['z'].iloc[0])
      N  = A - Z
      binding = dudu['binding'].iloc[0]/1000
      return Z*mp + N*mn - binding*A
    except:
      return numpy.nan

  def GetSymbol(self, A : int, Z : int) -> str:
    try:
      dudu = self.data[(self.data['z']==Z) & (self.data['A']==A)]
      return "%d%s" % (A , dudu['symbol'].iloc[0])
    except:
      return "0x"

  def GetJpi(self, ASym : str):
    [A, sym] = self.BreakDownName(ASym)
    try:
      dudu = self.data[(self.data['symbol']==sym) & (self.data['A']==A)]
      return dudu['jp'].iloc[0]
    except:
      return "unknown"

  def GetHalfLife(self, ASym : str) -> float:
    [A, sym] = self.BreakDownName(ASym)
    try:
      dudu = self.data[(self.data['symbol']==sym) & (self.data['A']==A)]
      return dudu['half_life_sec'].iloc[0]
    except:
      return "unknown"

  def GetSn(self, ASym : str) -> float:
    [A, Z] = self.GetAZ(ASym)
    if numpy.isnan(Z) :
      return numpy.nan
    else:
      mass0 = self.GetMassFromAZ(A, Z)
      mass1 = self.GetMassFromAZ(A-1, Z)
      return mass1 + mn - mass0  

  def GetSp(self, ASym : str):
    [A, Z] = self.GetAZ(ASym)
    if numpy.isnan(Z) :
      return numpy.nan
    else:
      mass0 = self.GetMassFromAZ(A, Z)
      mass1 = self.GetMassFromAZ(A-1, Z-1)
      return mass1 + mp - mass0

  def GetSa(self, ASym : str):
    [A, Z] = self.GetAZ(ASym)
    if numpy.isnan(Z) :
      return numpy.nan
    else:
      mass0 = self.GetMassFromAZ(A, Z)
      mass1 = self.GetMassFromAZ(A-4, Z-2)
      return mass1 + ma - mass0
    
  def PrintIso(self, ASym: str):
    [A, Z] = self.GetAZ(ASym)
    print("========================= ", ASym)
    print("A : %d, Z : %d, N : %d" % (A, Z, A-Z))
    print("      Jpi : ", self.GetJpi(ASym))
    print("half-live : %.2f sec" % (self.GetHalfLife(ASym)))
    print("     Mass : %9.2f MeV" % (self.GetMassFromSym(ASym) ))
    print("  Binding : %9.2f MeV/A" % (self.GetBindingPerA(ASym)))
    print("  Binding : %9.2f MeV" % (self.GetBindingPerA(ASym) * A))
    print("        Sn: %9.2f MeV" % self.GetSn(ASym))
    print("        Sp: %9.2f MeV" % self.GetSp(ASym))
    print("        Sa: %9.2f MeV" % self.GetSa(ASym))
    print("=============================")

  def PrintIsoWeb(self, ASym : str):
    [A, Z] = self.GetAZ(ASym)
    print("<br>========================= ", ASym)
    print("<br>A : %d, Z : %d, N : %d" % (A, Z, A-Z))
    print("<br>      Jpi : ", self.GetJpi(ASym))
    print("<br>half-live : %.2f sec" % (self.GetHalfLife(ASym)))
    print("<br>     Mass : %9.2f MeV" % (self.GetMassFromSym(ASym) ))
    print("<br>  Binding : %9.2f MeV/A" % (self.GetBindingPerA(ASym)))
    print("<br>  Binding : %9.2f MeV" % (self.GetBindingPerA(ASym) * A))
    print("<br>        Sn: %9.2f MeV" % self.GetSn(ASym))
    print("<br>        Sp: %9.2f MeV" % self.GetSp(ASym))
    print("<br>        Sa: %9.2f MeV" % self.GetSa(ASym))
    print("<br>=============================")


  def PrintIsoExWeb(self, ASym : str, maxEx : float):
    exList = self.GetExList(ASym, maxEx)
    if exList.empty:
      print("<br> cannot find Ex data")
    else:
      print("<table>")
      for i in range(0,len(exList)):
        print("<tr><td style=\"text-align:right\" width=80>%9.3f</td><td style=\"text-align:right\" width=100>%s</td></tr>" % (exList['energy'].iloc[i]/1000, exList['jp'].iloc[i].replace(' ', ',')))
      print("</table>")

