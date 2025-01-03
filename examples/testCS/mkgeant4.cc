#include <math.h>
#include "globals.hh"
#include "DarkMatter.hh"

#define nPointsE0 8
#define nPointsMass 62
#define RELATIVE_DELTA 0.005 // 0.5% delta to check if calculated tot cs is 'equal' to reference cs

#include "DarkPhotons.hh"
#include "DarkPhotonsSigmaTotETL.inc" // include file with ETL reference values

using std::cout;
using std::endl;
using std::cerr;

int main() {

/*GeV*/double EThresh0 = 0.01; // for total CS testing applicable to mass range mA<1MeV
/*GeV*/double EThresh = 0.01;
/*GeV*/double E0, massTested;
/*GeV*/double valuesE0[nPointsE0]={20., 40., 60., 80., 100., 120., 140., 160.};
/*MeV*/double testMassValues[nPointsMass]={0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 , 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000};

       double csCalcResult, csRefResult;
       int countNotEqual = 0;
       int arrayNotEqual[nPointsMass][nPointsE0] = {};

       for(int ii=0; ii<nPointsE0; ii++){
         E0 = valuesE0[ii];
         double* pETLii = pdataETL[ii];
         for(int jj=0; jj<nPointsMass; jj++) {
           massTested = testMassValues[jj]/1000.;                             //convert mass to GeV    
           EThresh = EThresh0;
           if(EThresh < 2.*massTested) EThresh = 2.*massTested;
           DarkMatter* myDarkMatter = new DarkPhotons(massTested, EThresh);  //Initialize DM by default for Pb with eps=0.0001
           csCalcResult = myDarkMatter->TotalCrossSectionCalc(E0);            //get calculated ETL cs from DMG4       
           csRefResult = (*pETLii); pETLii++;                                 //get reference ETL cs from include file

           if (fabs(csCalcResult-csRefResult)/csRefResult>RELATIVE_DELTA) {countNotEqual++; arrayNotEqual[jj][ii]++;}  //check if calculated tot cs is 'equal' to reference cs
         }
       }
       if(countNotEqual) {
         for(int jj=0; jj<nPointsMass; jj++) {cerr<<testMassValues[jj]<<"\t"; 
           for(int ii=0; ii<nPointsE0; ii++){cerr<<arrayNotEqual[jj][ii]<<" ";}
           cerr<<endl;
         }
         cerr<<" Error: Cross section test for DarkPhotons failed!"<<endl; 
         cerr<<"        Found "<<countNotEqual<<" deviation(s) from "<<nPointsMass*nPointsE0<<" reference points, exiting..."<<endl;
         exit(1);
       }
       else cout<<" All OK!:"<<endl;
       return 0;
}
