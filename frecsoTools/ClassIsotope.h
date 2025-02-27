/***********************************************************************
 * 
 *  This is ClassIsotope.h, To extract the isotope mass from nubase_4.mas20.txt
 * 
 *-------------------------------------------------------
 *  created by Ryan (Tsz Leung) Tang, Sept-6, 2024
 *  email: goluckyryan@gmail.com
 * ********************************************************************/


#ifndef ISOTOPE_H
#define ISOTOPE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "constant.h" // amu
#include <stdlib.h>  //atoi
#include <algorithm>

std::string massData="nubase_4.mas20.txt";

class Isotope {
public:
  int A, Z;
  double Mass, MassError, BEA; // BEA in keV
  short twoSpin; // 2*J
  short parity;
  std::string Name, Symbol; // Name = 6Li, Symbol = Li
  std::string dataSource;
  
  Isotope(){dataSource = massData;};
  Isotope(int a, int z){ dataSource = massData; SetIso(a,z);  };
  Isotope(std::string name){ dataSource = massData; SetIsoByName(name); };

  void SetMassTablePath(std::string path){ dataSource = path;}

  void SetIso(int a, int z);
  void SetIsoByName(std::string name);

  double CalSp(int Np, int Nn); // this for the Np-proton, Nn-neutron removal 
  double CalSp2(int a, int z); // this is for (a,z) nucleus removal

  double CalBeta(double T){
    double Etot = Mass + T;
    double gamma = 1 + T/Mass;
    double beta = sqrt(1 - 1 / gamma / gamma ) ;
    return beta;
  }

  void Print();
  void ListShell();
   
private:

  void InvalidIso();
  void BreakName(std::string nameStr);
  void BreakJpi(std::string jPiStr);
  void DigestLine(std::string line, int A, int Z);

  void FindMassByAZ(int a, int z); // give mass, massError, BEA, Name, Symbol;
  void FindMassByName(std::string name); // give Z, mass, massError, BEA;

  int TwoJ(int nShell);
  std::string Orbital(int nShell);
  int magic(int i){
    switch (i){
      case 0: return   2; break;
      case 1: return   8; break;
      case 2: return  20; break;
      case 3: return  28; break;
      case 4: return  40; break;
      case 5: return  50; break;
      case 6: return  82; break;
      case 7: return 128; break;
    }
    return 0;
  }

  int magicShellID(int i){
    switch (i){
      case 0: return   0; break;
      case 1: return   2; break;
      case 2: return   5; break;
      case 3: return   6; break;
      case 4: return   9; break;
      case 5: return  10; break;
      case 6: return  15; break;
      case 7: return  21; break;
    }
    return 0;
  }
  
  int fileStartLine;
  int fileEndLine;
  int lineMass050_099;
  int lineMass100_149;
  int lineMass150_199;
  int lineMass200_249;
  int lineMass250;

  int colA[2];
  int colZ[2];
  int colName[2];
  int colMassExcess[2];
  int colMassError[2];
  int colSpin[2];
  
  void setFileLinesAndCol(){
    fileStartLine = 26;
    fileEndLine = 5869;
    
    lineMass050_099 = 638;
    lineMass100_149 = 1718;
    lineMass150_199 = 3123;
    lineMass200_249 = 4591;
    lineMass250     = 5495;

    colA[0] = 0;
    colA[1] = 3;

    colZ[0] = 4;
    colZ[1] = 8;

    colName[0] = 11;
    colName[1] = 16;

    colMassExcess[0] = 18;
    colMassExcess[1] = 31;

    colMassError[0] = 31;
    colMassError[1] = 42;

    colSpin[0] = 88;
    colSpin[1] = 97;
  }
  
  bool isFindOnce;
  
};

void Isotope::SetIso(int a, int z){
  FindMassByAZ(a,z); 
}

void Isotope::SetIsoByName(std::string name){
  FindMassByName(name); 
}

void Isotope::InvalidIso(){
  this->BEA  = std::nan("");
  this->Mass = std::nan("");
  this->MassError = std::nan("");
  this->Symbol = "non";
  this->Name   = "non";
  this->twoSpin = 0;
  this->parity = 0;
}

void Isotope::BreakName(std::string nameStr){
  std::string::size_type lastDigit = nameStr.find_last_of("0123456789");

  this->Symbol = nameStr.substr(lastDigit + 1);
  if (this->Symbol.length() == 1) {
    this->Symbol += " ";
  }

  this->A = std::stoi(nameStr.substr(0, lastDigit + 1));

}

void Isotope::BreakJpi(std::string jPiStr){

  jPiStr.erase(std::remove(jPiStr.begin(), jPiStr.end(), '('), jPiStr.end());
  jPiStr.erase(std::remove(jPiStr.begin(), jPiStr.end(), ')'), jPiStr.end());

  std::string::size_type lastDigit = jPiStr.find_last_of("0123456789");

  std::string spinStr = jPiStr.substr(0, lastDigit + 1);

  std::string::size_type lastSlash = spinStr.find_last_of('/');

  if (lastSlash != std::string::npos) {
    twoSpin = std::stoi(spinStr.substr(0, lastSlash + 1));
  } else {
    twoSpin = 2 * std::stoi(spinStr);
  }

  this->parity = (jPiStr.substr(lastDigit + 1) == "+") ? 1 : -1;

}

void Isotope::DigestLine(std::string line, int A, int Z){

  double massExcess = atof((line.substr(colMassExcess[0],colMassExcess[1]-colMassExcess[0])).c_str());

  Mass = massExcess/1000. + A * amu - Z * me;
  MassError =  atof((line.substr(colMassError[0],colMassError[1]-colMassError[0])).c_str()) / 1000.;

  BEA = (Z * mp + (A-Z) * mn - Mass)/A *1000;

  std::string str = line.substr(colName[0],colName[1]-colName[0]);
  str.erase(remove(str.begin(), str.end(), ' '), str.end());
  Name = str;

  BreakName(Name);

  if( this->Name == "1H" ) this->Name = "p";
  if( this->Name == "2H" ) this->Name = "d";
  if( this->Name == "3H" ) this->Name = "t";
  if( this->Name == "4He" ) this->Name = "a";

  str = line.substr(colSpin[0], colSpin[1]-colSpin[0]);
  str.erase(remove(str.begin(), str.end(), ' '), str.end());
  str.erase(std::remove(str.begin(), str.end(), '*'), str.end());
  str.erase(std::remove(str.begin(), str.end(), '#'), str.end());
  std::string Jpi = str;

  BreakJpi(Jpi);

  // printf("A %d  Z: %d | %s | %f | %f |%s|\n", A, Z, Name.c_str(), massExcess, Mass, Jpi.c_str());
}


void Isotope::FindMassByAZ(int A, int Z){
  std::string line;
  int    lineNum=0;

  std::ifstream myfile;
  int    flag=0;

  setFileLinesAndCol();

  int numLineStart = fileStartLine;
  int numLineEnd  = fileEndLine;

  this->A = A;
  this->Z = Z;

  if ( A >= 50 && A < 100) numLineStart = lineMass050_099;
  if ( A >=100 && A < 150) numLineStart = lineMass100_149;
  if ( A >=150 && A < 200) numLineStart = lineMass150_199;
  if ( A >=200 && A < 250) numLineStart = lineMass200_249;
  if ( A >=250           ) numLineStart = lineMass250;

  myfile.open(dataSource.c_str());

  if (myfile.is_open()) {
    while (/*! myfile.eof() &&*/ flag == 0 && lineNum <numLineEnd){
      lineNum ++ ;
      // printf("%3d  ",lineNum);
      getline (myfile,line);

      if (lineNum >= numLineStart ){
        int list_Z     = std::stoi((line.substr(colZ[0],colZ[1]-colZ[0])));
      	int list_A     = std::stoi((line.substr(colA[0],colA[1]-colA[0])));

      	if ( A == list_A && Z * 10 == list_Z) {
    
          DigestLine(line, A, Z);          

     		  flag = 1;
      	}else if ( list_A > A) {
          InvalidIso();
          break;
        }

      }
    }
    
    myfile.close();
  }else {
    printf("Unable to open %s\n", dataSource.c_str());
  }
}

void Isotope::FindMassByName(std::string name){

  // done seperate the Mass number and the name 
  if( name == "n" ) {
    this->Name = "1n";
    this->BEA       = 0;
    this->Mass      = mn;
    this->MassError = 0;
    this->Name      = "n";
    this->A         = 1;
    this->Z         = 0;
    this->twoSpin   = 1;
    this->parity    = 1;
    return;
  }

  if( name == "p" ) name = "1H";
  if( name == "d" ) name = "2H";
  if( name == "t" ) name = "3H";
  if( name == "a" ) name = "4He";

  this->Name = name;
  BreakName(name);

  // find the nucleus in the data
  std::string line;
  int    lineNum=0;

  std::ifstream myfile;
  int    flag=0;

  setFileLinesAndCol();

  int numLineStart = fileStartLine;
  int numLineEnd  = fileEndLine;

  if ( A >= 50 && A < 100) numLineStart = lineMass050_099;
  if ( A >=100 && A < 150) numLineStart = lineMass100_149;
  if ( A >=150 && A < 200) numLineStart = lineMass150_199;
  if ( A >=200 && A < 250) numLineStart = lineMass200_249;
  if ( A >=250           ) numLineStart = lineMass250;

  myfile.open(dataSource.c_str());

  if (myfile.is_open()) {
    while (/*! myfile.eof() &&*/ flag == 0 && lineNum <numLineEnd){
      lineNum ++ ;
      //printf("%3d  ",lineNum);
      getline (myfile,line);

      if (lineNum >= numLineStart ){

        std::string nameStr = line.substr(colName[0],colName[1]-colName[0]);
        nameStr.erase(remove(nameStr.begin(), nameStr.end(), ' '), nameStr.end());

        int list_A     = std::stoi((line.substr(colA[0],colA[1]-colA[0])));
        int list_Z     = std::stoi((line.substr(colZ[0],colZ[1]-colZ[0])));

        // printf("A %d Z %d | %s \n", list_A, list_Z, nameStr.c_str());

        if ( list_A == A && list_Z % 10 == 0  &&  Name == nameStr) {

          Z = list_Z/10;

          DigestLine(line, list_A, Z);

          flag = 1;
        }else if ( list_A > this->A) {
          InvalidIso();
          break;
        }
      }
    }
    myfile.close();
  }else {
    printf("Unable to open %s\n", dataSource.c_str());
  }
}

double Isotope::CalSp(int Np, int Nn){
  Isotope nucleusD(A - Np - Nn, Z - Np);  

  if( !std::isnan(nucleusD.Mass) ){
    return nucleusD.Mass + Nn*mn + Np*mp - this->Mass;
  }else{
    return std::nan("");
  }
}

double Isotope::CalSp2(int a, int z){
  Isotope nucleusD(A - a , Z - z);
  Isotope nucleusS(a,z);  

  if( !std::isnan(nucleusD.Mass) && !std::isnan(nucleusS.Mass) ){
    return nucleusD.Mass + nucleusS.Mass - this->Mass;
  }else{
    return std::nan("");
  }
}

int Isotope::TwoJ(int nShell){

  switch(nShell){
    case  0: return  1; break; // 0s1/2
    case  1: return  3; break; // 0p3/2
    case  2: return  1; break; // 0p1/2 -- 8
    case  3: return  5; break; // 0d5/2
    case  4: return  1; break; // 1s1/2
    case  5: return  3; break; // 0d3/2 -- 20
    case  6: return  7; break; // 0f7/2 -- 28
    case  7: return  3; break; // 1p3/2
    case  8: return  1; break; // 1p1/2
    case  9: return  5; break; // 0f5/2 -- 40
    case 10: return  9; break; // 0g9/2 -- 50
    case 11: return  7; break; // 0g7/2
    case 12: return  5; break; // 1d5/2
    case 13: return 11; break; // 0h11/2
    case 14: return  3; break; // 1d3/2
    case 15: return  1; break; // 2s1/2 -- 82
    case 16: return  9; break; // 0h9/2
    case 17: return  7; break; // 1f7/2
    case 18: return 13; break; // 0i13/2
    case 19: return  3; break; // 2p3/2
    case 20: return  5; break; // 1f5/2
    case 21: return  1; break; // 1p1/2 -- 126
    case 22: return  9; break; // 1g9/2
    case 23: return 11; break; // 0i11/2
    case 24: return 15; break; // 0j15/2
    case 25: return  5; break; // 2d5/2
    case 26: return  1; break; // 3s1/2
    case 27: return  3; break; // 2d3/2
    case 28: return  7; break; // 1g7/2
  }

  return 0;
}

std::string Isotope::Orbital(int nShell){

  switch(nShell){
    case  0: return  "0s1 "; break;  //
    case  1: return  "0p3 "; break;  //
    case  2: return  "0p1 "; break;  //-- 8
    case  3: return  "0d5 "; break;  //
    case  4: return  "1s1 "; break;  //
    case  5: return  "0d3 "; break;  //-- 20
    case  6: return  "0f7 "; break;  //-- 28
    case  7: return  "1p3 "; break;  //
    case  8: return  "1p1 "; break;  //
    case  9: return  "0f5 "; break;  //-- 40
    case 10: return  "0g9 "; break;  //-- 50
    case 11: return  "0g7 "; break;  //
    case 12: return  "1d5 "; break;  //
    case 13: return  "0h11"; break;  //
    case 14: return  "1d3 "; break;  //
    case 15: return  "2s1 "; break;  //-- 82
    case 16: return  "0h9 "; break;  //
    case 17: return  "1f7 "; break;  //
    case 18: return  "0i13"; break;  //
    case 19: return  "2p3 "; break;  //
    case 20: return  "1f5 "; break;  //
    case 21: return  "1p1 "; break;  //-- 126
    case 22: return  "1g9 "; break;  //
    case 23: return  "0i11"; break;  //
    case 24: return  "0j15"; break;  //
    case 25: return  "2d5 "; break;  //
    case 26: return  "3s1 "; break;  //
    case 27: return  "2d3 "; break;  //
    case 28: return  "1g7 "; break;  //
  }

  return "nan";
}

void Isotope::ListShell(){

  if( Mass < 0 ) return;

  int n = A-Z;
  int p = Z;

  int k = std::min(n,p);
  int nMagic = 0;
  for( int i = 0; i < 7; i++){
    if( magic(i) < k && k <= magic(i+1) ){
      nMagic = i;
      break;
    }
  }

  int coreShell = magicShellID(nMagic-1);
  int coreZ1 = magic(nMagic-1);
  int coreZ2 = magic(nMagic);

  Isotope core1( 2*coreZ1, coreZ1);
  Isotope core2( 2*coreZ2, coreZ2);

  printf("------------------ Core:%3s, inner Core:%3s \n", (core2.Name).c_str(), (core1.Name).c_str());
  printf("         || ");
  int t = std::max(n,p);
  int nShell = 0;
  do{
    int occ = TwoJ(nShell)+1;
    if( nShell > coreShell) {
      printf("%4s", Orbital(nShell).c_str());
         if( nShell == 0 || nShell == 2 || nShell == 5 || nShell ==6 || nShell == 9  || nShell == 10 || nShell == 15 || nShell == 21){
        printf("|");
      }else{
        printf(",");
      } 
    }
    t = t - occ;
    nShell++;
  }while( t > 0  && nShell < 29);
  for( int i = 1; i <= 6; i++){
    if (nShell < 28) {
      printf("%4s,", Orbital(nShell).c_str());
    }else if( nShell == 28 ) {
      printf("%4s", Orbital(nShell).c_str());
    }
    nShell ++;
  }
  if (nShell < 29) printf("%4s", Orbital(nShell).c_str());
  printf("\n");


  printf(" Z = %3d || ", p);
  nShell = 0;
  do{
    int occ = TwoJ(nShell)+1;
    if( nShell > coreShell ){
      if( p > occ ) {
         printf("%-4d", occ);
         if( nShell == 0 || nShell == 2 || nShell == 5 || nShell ==6 || nShell == 9  || nShell == 10 || nShell == 15 || nShell == 21){
          printf("|");
        }else{
          printf(",");
        } 
      }else{
        printf("%-4d", p);
      }
    }
    p = p - occ;
    nShell++;
  }while( p > 0  && nShell < 29);
  printf("\n");

  printf(" N = %3d || ", n);
  nShell = 0;  
  do{
    int occ = TwoJ(nShell)+1;
    if ( nShell > coreShell ){
      if( n > occ ) {
         printf("%-4d", occ);
         if( nShell == 0 || nShell == 2 || nShell == 5 || nShell ==6 || nShell == 9  || nShell == 10 || nShell == 15 || nShell == 21){
          printf("|");
        }else{
          printf(",");
        } 
      }else{
        printf("%-4d", n);
      }
    }
    n = n - occ;
    nShell++;
  }while( n > 0  && nShell < 29);
  printf("\n");

  printf("------------------ \n");
}



void Isotope::Print(){

  if (Mass > 0){  
    
    printf(" using mass data : %s \n", dataSource.c_str());

    std::string jPiStr = "";
    if( twoSpin % 2  == 0) {
      jPiStr = std::to_string(twoSpin/2).c_str();
    }else{
      jPiStr = std::to_string(twoSpin) + "/2";
    }

    printf(" mass of \e[47m\e[31m%s\e[m nucleus (Z,A)=(%3d,%3d) is \e[47m\e[31m%12.5f\e[m MeV, J-pi \e[31m%s%s\e[m\n", 
            Name.c_str(), 
            Z, 
            A, 
            Mass,
            jPiStr.c_str(),
            parity == 1 ? "+" : "-");     

    printf(" total BE    : %12.5f MeV\n",BEA*A/1000.);
    printf(" total BE/A  : %12.5f MeV\n",BEA/1000.);
    printf(" mass in amu : %12.5f u\n",Mass/amu);
    printf(" mass excess : %12.5f MeV +- %12.5f\n", Mass + Z*0.510998950 - A*amu, MassError);
    printf("-------------- Seperation energy \n");
    printf(" S1p: %8.4f| S1n: %8.4f| S(2H ): %8.4f| S1p1n : %8.4f\n", CalSp(1, 0), CalSp(0, 1), CalSp2(2, 1), CalSp(1, 1));
    printf(" S2p: %8.4f| S2n: %8.4f| S(3He): %8.4f| S(3H) : %8.4f\n", CalSp(2, 0), CalSp(0, 2), CalSp2(3, 2), CalSp2(3, 1));
    printf(" S3p: %8.4f| S3n: %8.4f| S(4He): %8.4f|\n",               CalSp(3, 0), CalSp(0, 3), CalSp2(4, 2));
    printf(" S4p: %8.4f| S4n: %8.4f| \n",                             CalSp(4, 0), CalSp(0, 4));

  }else{
    printf("Error %6.0f, no nucleus with (Z,A) = (%3d,%3d). \n", Mass, Z, A);
  }


}

#endif
