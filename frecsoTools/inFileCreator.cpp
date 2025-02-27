#include <fstream>
#include <stdlib.h>     /* atof */
#include <cmath>
#include <vector>
#include <string>
#include <stdlib.h>     /* atof */
#include "ClassIsotope.h" // for geting Z
#include "potentials.h"

#define RED   "\033[31m"
#define GREEN "\033[32m"
#define BLUE  "\033[34m"
#define RESET "\033[0m"

int GetLValue(std::string spdf){
  
  if( spdf == "s" ) return 0;
  if( spdf == "p" ) return 1;
  if( spdf == "d" ) return 2;
  if( spdf == "f" ) return 3;
  if( spdf == "g" ) return 4;
  if( spdf == "h" ) return 5;
  if( spdf == "i" ) return 6;
  if( spdf == "j" ) return 7;  
  return -1;
}

FILE * file_out = nullptr;

std::vector<std::string> SplitStr(std::string tempLine, std::string splitter, int shift = 0){

  std::vector<std::string> output;

  size_t pos;
  do{
    pos = tempLine.find(splitter); /// fine splitter
    if( pos == 0 ){ ///check if it is splitter again
      tempLine = tempLine.substr(pos+1);
      continue;
    }

    std::string secStr;
    if( pos == std::string::npos ){
      secStr = tempLine;
    }else{
      secStr = tempLine.substr(0, pos+shift);
      tempLine = tempLine.substr(pos+shift);
    }

    ///check if secStr is begin with space
    while( secStr.substr(0, 1) == " ") secStr = secStr.substr(1);
    
    ///check if secStr is end with space
    while( secStr.back() == ' ') secStr = secStr.substr(0, secStr.size()-1);

    output.push_back(secStr);
    ///printf(" |%s---\n", secStr.c_str());
    
  }while(pos != std::string::npos );

  return output;
};

std::vector<std::string> Breakdown_AabB( std::string reactStr ){

  std::vector<std::string> str1 = SplitStr(reactStr, "(", 0);
  std::vector<std::string> str2 = SplitStr(str1[1], ")", 1);
  std::vector<std::string> str3 = SplitStr(str2[0], ",", 0);

  std::vector<std::string> output;
  output.push_back ( str1[0] );
  output.push_back ( str3[0] );
  str3[1].pop_back();
  output.push_back ( str3[1] );
  output.push_back ( str2[1] );

  return output;

}

float GetTotalEnergy( std::string ELabStr, float beamA ){
  //get Beam energy, distingusih MeV or MeV/u
  int pos = ELabStr.length() - 1;
  for( int i = pos; i >= 0 ; i--){
    if( isdigit(ELabStr[i]) ) {
      pos = i; 
      break;
    }
  }
  std::string unit = ELabStr.substr(pos+1);
  int factor = 1;
  if( unit == "MeV/u") factor = beamA;
  return atof(ELabStr.substr(0, pos+1).c_str()) * factor;
}

int InFileCreator(std::string readFile, std::string infile, double angMin, double angMax, double angStep) {
   
  //================= read infile. extract the reactions, write pptolemy infile for each reaction
  std::ifstream file_in;
  file_in.open(readFile.c_str(), std::ios::in);

  if( !file_in ){
    printf(" cannot read file. \n");
    return 0 ; 
  }
   
  printf("Save to infile : %s \n", infile.c_str()); 
  file_out = fopen (infile.c_str(), "w+");

  printf("Angle setting (%5.2f, %5.2f) deg | Step : %5.2f deg\n", angMin, angMax, angStep);
  printf("---------------------------\n");

  //^============================== extract reaction line by line
  int numOfReaction = 0;
  while( file_in.good() ) {
    std::string tempLine;
    std::getline(file_in, tempLine );

    if( tempLine.substr(0, 1) == "#" ) continue;
    if( tempLine.size() < 5 ) continue;

    //split line using space
    std::vector<std::string> strList = SplitStr(tempLine, " ");
    if ( strList.size() == 0 ) continue;

    printf("  %s\n", tempLine.c_str());

    std::vector<std::string> AabB = Breakdown_AabB( strList[0] );

    if( AabB.size() != 4 ) {
      printf("     %sUnable to do 3-body reactions.%s\n", RED, RESET);
      continue;
    }
    
    Isotope iso_A(AabB[0]);
    Isotope iso_a(AabB[1]);
    Isotope iso_b(AabB[2]);
    Isotope iso_B(AabB[3]);

    float spin_A = atof(strList[1].c_str());

    short twoSpin = 0, parity = 0;
    if( strList[3] != "none"){ /// extrac spin and parity
      std::string Jpi = strList[3];
      std::string::size_type lastDigit = Jpi.find_last_of("0123456789");
      std::string spinStr = Jpi.substr(0, lastDigit + 1);
      std::string::size_type lastSlash = spinStr.find_last_of('/');
      if (lastSlash != std::string::npos) {
        twoSpin = std::stoi(spinStr.substr(0, lastSlash + 1));
      } else {
        twoSpin = 2 * std::stoi(spinStr);
      }
      parity = (Jpi.substr(lastDigit + 1) == "+") ? 1 : -1;
    }    


    float Ex = atof(strList[4].c_str());
    float totalEnergy = GetTotalEnergy(strList[5], iso_a.A );

    float extra = 0;
    if( strList.size() == 8 ){ 
      extra = atof( strList[7].c_str() );
    }

    // printf(" Ex : %.2f, totalEnergy : %.2f \n", Ex, totalEnergy);

    fprintf(file_out, "##############################################\n");

    //^===================================== elastics + inelastics
    if( iso_a.Mass == iso_b.Mass ) {

      if( Ex == 0.0 ){
 
        printf("=========== elatsic \n");

        fprintf(file_out, "%s@%.2f MeV/u Elastic\n" , strList[0].c_str(), totalEnergy/iso_a.A );
        fprintf(file_out, "NAMELIST\n");
        
        fprintf(file_out, "&FRESCO  hcm=0.1  rmatch=60 \n");
        fprintf(file_out, "         jtmin=0.0  jtmax=50 absend=0.001 \n");
        fprintf(file_out, "         thmin=%.1f  thmax=%.1f thinc=%.1f \n", angMin, angMax, angStep);
        fprintf(file_out, "         chans=1  smats=2 xstabl=1 \n");
        fprintf(file_out, "         elab=%.1f  /\n", totalEnergy);

        fprintf(file_out, "&PARTITION  namep='%s' massp=%.2f zp=%d namet='%s' masst=%.2f zt=%d qval=0.000  nex=1 /\n", iso_a.Name.c_str(), iso_a.Mass/amu, iso_a.Z, iso_A.Name.c_str(), iso_A.Mass/amu, iso_A.Z);

        fprintf(file_out, "&STATES   jp=%.1f  ptyp=%d  ep=0.000 jt=%.1f  ptyt=%d  et=0.000 cpot=1 /\n", iso_a.twoSpin/2., iso_a.parity, spin_A, iso_A.parity);

        fprintf(file_out, "&partition /\n");
        
        fprintf(file_out, "&POT kp=1        ap=%d  at=%d  rc=1.2 /\n", iso_a.A, iso_A.A); // Coulomb
        fprintf(file_out, "&POT kp=1 type=1 p1=40.00 p2=1.2 p3=0.65 p4=10.0 p5=1.2 p6=0.500 /\n"); // Woods-Saxon (V0, r0, a0, Wi, ri, ai)

        fprintf(file_out, "&pot /\n");
        
        fprintf(file_out, "&overlap /\n");
        fprintf(file_out, "&coupling /\n");
        
      }else{

        printf("=========== inelatsic \n");
        
        fprintf(file_out, "%s@%.2f MeV/u Inelastic\n" , strList[0].c_str(), totalEnergy/iso_a.A );
        fprintf(file_out, "NAMELIST\n");
        
        fprintf(file_out, "&FRESCO  hcm=0.1  rmatch=60 \n");
        fprintf(file_out, "         jtmin=0.0  jtmax=50 absend=0.001 \n");
        fprintf(file_out, "         thmin=%.1f  thmax=%.1f thinc=%.1f \n", angMin, angMax, angStep);
        fprintf(file_out, "         chans=1  smats=2 xstabl=1 \n");
        fprintf(file_out, "         elab=%.1f  /\n", totalEnergy);

        fprintf(file_out, "&PARTITION  namep='%s' massp=%.2f zp=%d namet='%s' masst=%.2f zt=%d qval=0.000  nex=2 /\n", iso_a.Name.c_str(), iso_a.Mass/amu, iso_a.Z, iso_A.Name.c_str(), iso_A.Mass/amu, iso_A.Z);

        fprintf(file_out, "&STATES   jp=%.1f  ptyp=%d  ep=0.000 jt=%.1f  ptyt=%d  et=0.000 cpot=1 /\n", iso_a.twoSpin/2., iso_a.parity, spin_A, iso_A.parity);
        fprintf(file_out, "&STATES   jp=%.1f  ptyp=%d  ep=0.000 jt=%.1f  ptyt=%d  et=%.3f cpot=1 /\n", iso_a.twoSpin/2., iso_a.parity, twoSpin/2., parity, Ex);

        fprintf(file_out, "&partition /\n");
        
        fprintf(file_out, "&POT kp=1         ap=%d  at=%d  rc=1.2 /\n", iso_a.A, iso_A.A); // Coulomb
        fprintf(file_out, "&POT kp=1 type=1  p1=40.00 p2=1.2 p3=0.65 p4=10.0 p5=1.2 p6=0.500 /\n"); // volume, Woods-Saxon (V0, r0, a0, Wi, ri, ai)
        fprintf(file_out, "&POT kp=1 type=11 p1= 0.00 p2=%.2f/\n", extra); // target deformation add to the volume, deformation length delta_2 = 1.3 fm
        fprintf(file_out, "&POT kp=1 type=2  p1= 0.00 p2=1.2 p3=0.65 p4=6.0 p5=1.2 p6=0.500 /\n"); // surface, d(Woods-Saxon)/dr (V0, r0, a0, Wi, ri, ai)
        fprintf(file_out, "&POT kp=1 type=11 p1= 0.00 p2=%.2f/\n", extra); // target deformation add to the suface

        fprintf(file_out, "&pot /\n");
        
        fprintf(file_out, "&overlap /\n");
        fprintf(file_out, "&coupling /\n");
      }

    }

    //^====================================== single neutron transfer
    if( abs( iso_a.A - iso_b.A ) == 1  ){

      printf("=========== single nucleon transfer \n");

      fprintf(file_out, "%s@%.2f MeV/u Inelastic\n" , strList[0].c_str(), totalEnergy/iso_a.A );
      fprintf(file_out, "NAMELIST\n");
      
      fprintf(file_out, "&FRESCO  hcm=0.1  rmatch=60 \n");
      fprintf(file_out, "         rintp=0.2  hnl=0.1   rnl=5.00 centre=0.0 \n");
      fprintf(file_out, "         iter=1 nnu=36 \n");
      fprintf(file_out, "         jtmin=0.0  jtmax=50 absend=0.001 \n");
      fprintf(file_out, "         thmin=%.1f  thmax=%.1f thinc=%.1f \n", angMin, angMax, angStep);
      fprintf(file_out, "         chans=1  smats=2 xstabl=1 \n");
      fprintf(file_out, "         elab=%.1f  /\n", totalEnergy);

      fprintf(file_out, "&PARTITION  namep='%s' massp=%.2f zp=%d namet='%s' masst=%.2f zt=%d qval=0.000  nex=1 /\n", iso_a.Name.c_str(), iso_a.Mass/amu, iso_a.Z, iso_A.Name.c_str(), iso_A.Mass/amu, iso_A.Z);
      fprintf(file_out, "&STATES   jp=%.1f  ptyp=%d  ep=0.000 jt=%.1f  ptyt=%d  et=0.000 cpot=1 /\n", iso_a.twoSpin/2., iso_a.parity, spin_A, iso_A.parity);

      float QValue = iso_a.Mass + iso_A.Mass  - iso_b.Mass - iso_B.Mass - Ex;
      fprintf(file_out, "&PARTITION  namep='%s' massp=%.2f zp=%d namet='%s' masst=%.2f zt=%d qval=%.3f  nex=1 /\n", iso_b.Name.c_str(), iso_b.Mass/amu, iso_b.Z, iso_B.Name.c_str(), iso_B.Mass/amu, iso_B.Z, QValue);
      fprintf(file_out, "&STATES   jp=%.1f  ptyp=%d  ep=0.000 jt=%.1f  ptyt=%d  et=%.2f cpot=1 /\n", iso_b.twoSpin/2., iso_b.parity, twoSpin/2., parity, Ex);

      fprintf(file_out, "&partition /\n");
      
      //*================ a + A
      fprintf(file_out, "&POT kp=1         ap=%d  at=%d  rc=1.2 /\n", iso_a.A, iso_A.A); // Coulomb
      fprintf(file_out, "&POT kp=1 type=1  p1=40.00 p2=1.2 p3=0.65 p4=10.0 p5=1.2 p6=0.500 /\n"); // volume, Woods-Saxon (V0, r0, a0, Wi, ri, ai)

      //*================ b + B      
      fprintf(file_out, "&POT kp=2         ap=%d  at=%d  rc=1.2 /\n", iso_b.A, iso_B.A); // Coulomb
      fprintf(file_out, "&POT kp=2 type=1  p1=40.00 p2=1.2 p3=0.65 p4=10.0 p5=1.2 p6=0.500 /\n"); // volume, Woods-Saxon (V0, r0, a0, Wi, ri, ai)

      //*================ a
      fprintf(file_out, "&POT kp=3         at=%d  rc=1.2 /\n", iso_a.A); // Coulomb
      fprintf(file_out, "&POT kp=3 type=1  p1=50.00 p2=1.2 p3=0.65 /\n"); // volume, Woods-Saxon (V0, r0, a0)
      fprintf(file_out, "&POT kp=3 type=3  p1= 6.00 p2=1.2 p3=0.65 /\n"); // S-O, projectile (V0, r0, a0)

      //*================ a - N

      //*================ A + (a-N)
           

      fprintf(file_out, "&pot /\n");
      
      fprintf(file_out, "&overlap /\n");
      fprintf(file_out, "&coupling /\n");

    }

    //^====================================== two neutron transfer
    if( abs( iso_a.A - iso_b.A ) == 2  ){


    }

  }

  fprintf(file_out, "#>>>>>############################################# End of File\n");

  file_in.close();
  fclose(file_out);

  return 0;

}