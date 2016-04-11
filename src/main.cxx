/**
 * @file      main.cxx
 * @brief     Main body of the system.
 * @author    Wooyoung Jang (wyjang)
 * @date      2015. 04. 01
 */
#include <iostream>
#include <cstring>

#ifndef __AMSINC__
#define __AMSINC__
#include "amschain.h"
#include "TrdKCluster.h"
#endif

#ifndef __ACSOFTINC__
#define __ACSOFTINC__
#include "FileManager.hh"
#include "AMSRootSupport.hh"
#include "AMSRootParticleFactory.hh"
#include "Utilities.hh"
#include "SlowControlLookup.hh"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStopwatch.h>
#include <TTree.h>
#endif

#include "selector.h"

// Some global variables
char releaseName[16];
char softwareName[10] = "ACCTOFInt";
char versionNumber[5] = "0.01";
bool debugMode = false;

bool RegisterFileToChain(int argc, char* argv[], AMSChain* amsChain);

/**
 * @brief This is main() function.
 * @return 0 : The program terminated properly. Otherwise : The program unintentionally terminated by error.
 */

using namespace std;

int main(int argc, char* argv[])
{
  sprintf(releaseName, "%s%s", softwareName, versionNumber);

  AMSChain *amsChain = new AMSChain();
  unsigned int  nEntries = 0;            // Number of entries to be analyzed.
  char outputFileName[256];     // The name of the output file to be stored in the disk.

  bool isDataLoadedOk = RegisterFileToChain(argc, argv, amsChain);

  if( isDataLoadedOk != true )
  {
    std::cout << "[" << releaseName <<"] FATAL ERROR: File Open Failed! " << std::endl;
    return -1;
  }

  nEntries = amsChain->GetEntries();
  if( argc == 1 ) strcpy(outputFileName, "testrun.root");
  else if( argc == 2 ) strcpy(outputFileName, "testrun.root");
  else if( argc == 3 ) strcpy(outputFileName, argv[2]);
  else if( argc == 4 ) strcpy(outputFileName, argv[2]);

  std::cout << "[" << releaseName << "] TOTAL EVENTS   : " << nEntries << std::endl;
  std::cout << "[" << releaseName << "] OUTPUT FILE NAME : " << outputFileName << std::endl;

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //
   //  ANALYSIS PHASE
   //
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  AMSSetupR::RTI::UseLatest();
  TkDBc::UseFinal();

  const bool setAMSRootDefaults = true;
  AMSRootSupport amsRootSupport(AC::ISSRun, setAMSRootDefaults);
  Analysis::EventFactory& eventFactory = amsRootSupport.EventFactory();
  Analysis::AMSRootParticleFactory& particleFactory = amsRootSupport.ParticleFactory();

  // Variable declaration which will be used for tree construction
  unsigned int  nCuts = 7;/*{{{*/           // Number of cuts
  unsigned int  nProcessed = 0;             // Number of processed events
  unsigned int  nRun;                       // Run number
  unsigned int  nEvent;                     // Event number
  unsigned int  nProcessedNumber;           // Number of processed events
  unsigned int  nLevel1;                    // Number of Level1 triggers
  unsigned int  nParticle;                  // Number of particles
  unsigned int  nCharge;                    // Number of charges
  unsigned int  nTrTrack;                   // Number of the Tracker tracks which are successfully reconstructed
  unsigned int  nTrdTrack;                  // Number of the TRD tracks which are successfully reconstructed
  unsigned int  nAntiCluster;               // Number of clusters on the ACC
  unsigned int  nTofClustersInTime;         // Number of in-time clusters on the TOF
  unsigned int  nRichRing;                  // Number of successfully reconstructed the RICH rings
  unsigned int  nRichRingB;                 // Number of successfully reconstructed the RICH rings with algorithm B
  unsigned int  nBeta;                      // Number of successfully estimated beta(v/c) values
  unsigned int  nBetaB;                     // Number of successfully estimated beta(v/c) values with algorithm B
  unsigned int  nBetaH;                     // Number of successfully estimated beta(v/c) values with algorithm H
  unsigned int  nShower;                    // Number of the ECAL shower objects
  unsigned int  nVertex;                    // Number of vertices in this event
  unsigned int  particleType;               // Type of ParticleR
  float         livetime;                   // Livetime fraction
  float         utcTime;                    // UTC time
  float         orbitAltitude;              // (cm) in GTOD coordinates system.
  float         orbitLatitude;              // (rad) in GTOD coordinates system.
  float         orbitLongitude;             // (rad) in GTOD coordinates system.
  float         orbitLatitudeM;             // (rad) in eccentric dipole coordinate system.
  float         orbitLongitudeM;            // (rad) in eccentric dipole coordinate system.
  float         velR;                       // Speed of the ISS in radial direction
  float         velTheta;                   // Angular speed of the ISS in polar angle direction
  float         velPhi;                     // Angular speed of the ISS in azimuthal angle direction
  float         yaw;                        // A parameter describing ISS attitude (tilted angle with respect to the x-axis)
  float         pitch;                      // A parameter describing ISS attitude (tilted angle with respect to the y-axis)
  float         roll;                       // A parameter describing ISS attitude (tilted angle with respect to the z-axis)
  float         gLongitude;                 // Galactic longitude of the incoming particle.
  float         gLatitude;                  // Galactic latitude of the incoming particle.
  int           gCoordCalcResult;           // Return value for galactic coordinate calculation.
  int           gCoordCalcResultBit;        // Result of GetGalCoo written in bit form
  float         sunPosAzimuth;              // Azimuthal angle of the position of the Sun.
  float         sunPosElevation;            // Elevation angle of the position of the Sun.
  int           sunPosCalcResult;           // Return value for the Sun's position calculation.
  unsigned int  unixTime;                   // UNIX time
  float         solarArrayCoord[3];
  int           isInShadow;                 // Value for check whether the AMS is in ISS solar panel shadow or not.
  float         zenithAngle;
  int           isInSAA;
  unsigned int  ptlCharge;                  // ParticleR::Charge value
  float         ptlMomentum;                // ParticleR::Momentum value
  float         ptlTheta;                   // Direction of the incoming particle (polar angle)
  float         ptlPhi;                     // Direction of the incoming particle (azimuthal angle)
  float         ptlCoo[3];                  // Coo (1st[last tracker plane) cm
  float         ptlTofCoo[4][3];            // Track extrapolated point [TOF Layer][x/y/z]
  float         ptlTofTrLength[4];          // Track length till TOF planes crossing
  float         ptlAntiCoo[2][5];           // Track extrapolated in ACC [Direction][x/y/z/theta/phi]
  float         ptlEcalCoo[3][3];           // Track extrapolated in ECAL [Entering position/ Center of gravity/ Exiting position][x/y/z]
  float         ptlTrCoo[9][3];             // Track extrapolated in Tracker [Layer][x/y/z]
  float         ptlTrdCoo[2][3];            // Track extrapolated in TRD [at center / at top][x/y/z]
  float         ptlRichCoo[2][3];           // Track extrapolated in RICH [at radiator / at PMT][x/y/z]
  float         ptlRichPath[2];
  float         ptlRichPathBeta[2];
  float         ptlRichLength;              // Estimated pathlength of particle within RICH radiator (cm).
  int           ptlRichParticles;
  float         ptlCutOffStoermer;
  float         ptlCutOffDipole;
  float         ptlCutOffMax[2];
  float         showerEnergyD;
  float         showerEnergyDL[18];
  float         showerEnergyE;
  float         showerEnergyCorrected;
  float         showerBDT;
  float         showerCofG[3];
  float         showerCofGDist;
  float         showerCofGdX;
  float         showerCofGdY;
  int           tofNCluster;
  int           tofNClusterH;
  int           tofNUsedHits;
  int           tofNUnusedHits;
  int           tofNUsedLayersForQ;
  float         tofBeta;
  float         tofInvBetaErr;
  float         tofMass;
  float         tofMassError;
  int           isGoodBeta;
  int           isTkTofMatch;
  float         tofReducedChisqT;
  float         tofReducedChisqC;
  float         tofDepositedEnergyOnLayer[4];
  float         tofEstimatedChargeOnLayer[4];
  float         tofCharge;
  float         tofUpperCharge;
  float         tofLowerCharge;
  float         tofChargeOnLayer[4];
  int           trkFitCodeFS;
  float         trkRigidityFS;
  float         trkRigidityInverseErrorFS;
  float         trkReducedChisquareFSX;
  float         trkReducedChisquareFSY;
  int           trkFitCodeMS;
  float         trkRigidityMS;
  float         trkRigidityInverseErrorMS;
  float         trkReducedChisquareMSX;
  float         trkReducedChisquareMSY;
  int           trkFitCodeInner;
  float         trkRigidityInner;
  float         trkRigidityInverseErrorInner;
  float         trkReducedChisquareInnerX;
  float         trkReducedChisquareInnerY;
  int           trkEdepLayerJXSideOK[9];
  int           trkEdepLayerJYSideOK[9];
  float         trkEdepLayerJ[9];
  float         trkCharge;
  float         trkInnerCharge;
  int           trkHasExtLayers;
  int           richRebuild;
  int           richIsGood;
  int           richIsClean;
  int           richIsNaF;
  int           richUsedHits;
  float         richRingWidth;
  int           richNHits;
  int           richNPMTsOnRing;
  float         richBeta;
  float         richBetaError;
  float         richChargeSquared;
  float         richKolmogorovProbability;
  float         richPhotoelectrons;                         // Return value from RichRingR::getPhotoelectrons(bool corr=true)
  float         richPhotoElectrons;                         // Return value from RichRingR::getPhotoElectrons(int pmt, bool corr)
  float         richExpectedPhotoelectrons;                 // Return value from RichRingR::getExpectedPhotoElectrons(bool corr=true)
  float         richExpectedPhotoElectrons;                 // Return value from RichRingR::getExpectedPhotoElectrons(int pmt, bool corr)
  int           richNUsedHits;
  float         richTheta;
  float         richPhi;
  int           trdNClusters;
  int           trdNUnusedHits;
  int           trdNTracks;
  float         trdTrackPhi;
  float         trdTrackTheta;
  float         trdTrackChi2;
  int           trdTrackPattern;
  float         trdTrackCharge;
  float         trdTrackEdepL[20];
  float         trdDepositedEnergyOnLayer[20];
  float         trdTrackMeanDepositedEnergy;
  float         trdTrackTotalDepositedEnergy;
  int           trdQtNActiveLayer;
  int           trdQtIsCalibrationGood;
  int           trdQtIsSlowControlDataGood;
  int           trdQtIsInsideTrdGeometricalAcceptance;
  int           trdQtIsValid;
  int           trdQtNActiveStraws;
  int           trdQtNActiveLayers;
  int           trdQtNTRDVertex;
  float         trdQtElectronToProtonLogLikelihoodRatio;
  float         trdQtHeliumToElectronLogLikelihoodRatio;
  float         trdQtHeliumToProtonLogLikelihoodRatio;
  int           trdKNRawHits;
  int           trdKIsReadAlignmentOK;
  int           trdKIsReadCalibOK;
  int           trdKNHits;
  int           trdKIsValid;
  float         trdKElectronToProtonLogLikelihoodRatio;
  float         trdKHeliumToElectronLogLikelihoodRatio;
  float         trdKHeliumToProtonLogLikelihoodRatio;
  float         trdKCharge;
  float         trdKChargeError;
  int           trdKNUsedHitsForCharge;
  float         trdKAmpLayer[20];
  float         trdKTotalPathLength;
  float         trdKTotalAmp;
  float         trdKElectronLikelihood;
  float         trdKProtonLikelihood;
  float         trdKHeliumLikelihood;
  float         trdPElectronToProtonLogLikelihoodRatio;
  float         trdPHeliumToElectronLogLikelihoodRatio;
  float         trdPHeliumToProtonLogLikelihoodRatio;
  int           accNHits;
  int           accNRecoClPG;
  vector<int>  accSector;
  vector<float> accTime;
  vector<float> accHitPosZ;
  vector<float> accChi2;
  vector<int>   accNPairs;
  vector<float> accUnfoldedHitPosZ;
  vector<float> accRawCharge;
  vector<float> accHitPosZfromADC;
  vector<float> accUnfoldedHitPosZfromADC;
  float         accCrossingPhiAngle;
  vector<float>  accEdep;
  vector<float>  accTimePG;
  /*}}}*/

  // Declare event counter
  TH1D* hEvtCounter = new TH1D("hEvtCounter", "Event Collecting Status", nCuts, 0., (float)nCuts);
  hEvtCounter->GetXaxis()->SetBinLabel(1, "No AMSEventR");
  hEvtCounter->GetXaxis()->SetBinLabel(1, "Science-run check");
  hEvtCounter->GetXaxis()->SetBinLabel(2, "DAQ H/W status check");
  hEvtCounter->GetXaxis()->SetBinLabel(3, "Unbiased physics trigger check");
  hEvtCounter->GetXaxis()->SetBinLabel(4, "Multiple particles");
  hEvtCounter->GetXaxis()->SetBinLabel(5, "Tracker alignment test");
  hEvtCounter->GetXaxis()->SetBinLabel(6, "Good track test");
  hEvtCounter->GetXaxis()->SetBinLabel(7, "SAA rejection");

  // Make a TFile and TTree
  TFile* resultFile = new TFile(outputFileName, "RECREATE");
  TTree* tree  = new TTree(softwareName, releaseName);

  // Assign branches
  tree->Branch("nRun",                                    &nRun,                                    "nRun/i");/*{{{*/
  tree->Branch("nEvent",                                  &nEvent,                                  "nEvent/i");
  tree->Branch("nProcessedNumber",                        &nProcessedNumber,                        "nProcessedNumber/i");
  tree->Branch("nLevel1",                                 &nLevel1,                                 "nLevel1/i");
  tree->Branch("nParticle",                               &nParticle,                               "nParticle/i");
  tree->Branch("nCharge",                                 &nCharge,                                 "nCharge/i");
  tree->Branch("nTrTrack",                                &nTrTrack,                                "nTrTrack/i");
  tree->Branch("nTrdTrack",                               &nTrdTrack,                               "nTrdTrack/i");
  tree->Branch("nAntiCluster",                            &nAntiCluster,                            "nAntiCluster/i");
  tree->Branch("nTofClustersInTime",                      &nTofClustersInTime,                      "nTofClustersInTime/I");
  tree->Branch("nRichRing",                               &nRichRing,                               "nRichRing/i");
  tree->Branch("nRichRingB",                              &nRichRingB,                              "nRichRingB/i");
  tree->Branch("nBeta",                                   &nBeta,                                   "nBeta/i");
  tree->Branch("nBetaB",                                  &nBetaB,                                  "nBetaB/i");
  tree->Branch("nBetaH",                                  &nBetaH,                                  "nBetaH/i");
  tree->Branch("nShower",                                 &nShower,                                 "nShower/i");
  tree->Branch("nVertex",                                 &nVertex,                                 "nVertex/i");
  tree->Branch("particleType",                            &particleType,                            "particleType/i");
  tree->Branch("livetime",                                &livetime,                                "livetime/F");
  tree->Branch("utcTime",                                 &utcTime,                                 "utcTime/F");
  tree->Branch("orbitAltitude",                           &orbitAltitude,                           "orbitAltitude/F");
  tree->Branch("orbitLatitude",                           &orbitLatitude,                           "orbitLatitude/F");
  tree->Branch("orbitLongitude",                          &orbitLongitude,                          "orbitLongitude/F");
  tree->Branch("orbitLatitudeM",                          &orbitLatitudeM,                          "orbitLatitudeM/F");
  tree->Branch("orbitLongitudeM",                         &orbitLongitudeM,                         "orbitLongitudeM/F");
  tree->Branch("velR",                                    &velR,                                    "velR/F");
  tree->Branch("velTheta",                                &velTheta,                                "velTheta/F");
  tree->Branch("velPhi",                                  &velPhi,                                  "velPhi/F");
  tree->Branch("yaw",                                     &yaw,                                     "yaw/F");
  tree->Branch("pitch",                                   &pitch,                                   "pitch/F");
  tree->Branch("roll",                                    &roll,                                    "roll/F");
  tree->Branch("gLongitude",                              &gLongitude,                              "gLongitude/F");
  tree->Branch("gLatitude",                               &gLatitude,                               "gLatitude/F");
  tree->Branch("gCoordCalcResult",                        &gCoordCalcResult,                        "gCoordCalcResult/I");
  tree->Branch("gCoordCalcResultBit",                     &gCoordCalcResultBit,                     "gCoordCalcResultBit/I");
  tree->Branch("sunPosAzimuth",                           &sunPosAzimuth,                           "sunPosAzimuth/F");
  tree->Branch("sunPosElevation",                         &sunPosElevation,                         "sunPosElevation/F");
  tree->Branch("sunPosCalcResult",                        &sunPosCalcResult,                        "sunPosCalcResult/I");
  tree->Branch("unixTime",                                &unixTime,                                "unixTime/i");
  tree->Branch("solarArrayCoord",                         &solarArrayCoord,                         "solarArrayCoord[3]/F");
  tree->Branch("isInShadow",                              &isInShadow,                              "isInShadow/I");
  tree->Branch("zenithAngle",                             &zenithAngle,                             "zenithAngle/F");
  tree->Branch("isInSAA",                                 &isInSAA,                                 "isInSAA/I");
  tree->Branch("ptlCharge",                               &ptlCharge,                               "ptlCharge/i");
  tree->Branch("ptlMomentum",                             &ptlMomentum,                             "ptlMomentum/F");
  tree->Branch("ptlTheta",                                &ptlTheta,                                "ptlTheta/F");
  tree->Branch("ptlPhi",                                  &ptlPhi,                                  "ptlPhi/F");
  tree->Branch("ptlCoo",                                  &ptlCoo,                                  "ptlCoo[3]/F");
  tree->Branch("ptlTofCoo",                               &ptlTofCoo,                               "ptlTofCoo[4][3]/F");
  tree->Branch("ptlTofTrLength",                          &ptlTofTrLength,                          "ptlTofTrLength[4]/F");
  tree->Branch("ptlAntiCoo",                              &ptlAntiCoo,                              "ptlAntiCoo[2][5]/F");
  tree->Branch("ptlEcalCoo",                              &ptlEcalCoo,                              "ptlEcalCoo[3][3]/F");
  tree->Branch("ptlTrCoo",                                &ptlTrCoo,                                "ptlTrCoo[9][3]/F");
  tree->Branch("ptlTrdCoo",                               &ptlTrdCoo,                               "ptlTrdCoo[2][3]/F");
  tree->Branch("ptlRichCoo",                              &ptlRichCoo,                              "ptlRichCoo[2][3]/F");
  tree->Branch("ptlRichPath",                             &ptlRichPath,                             "ptlRichPath[2]/F");
  tree->Branch("ptlRichPathBeta",                         &ptlRichPathBeta,                         "ptlRichPathBeta[2]/F");
  tree->Branch("ptlRichLength",                           &ptlRichLength,                           "ptlRichLength/F");
  tree->Branch("ptlRichParticles",                        &ptlRichParticles,                        "ptlRichParticles/I");
  tree->Branch("ptlCutOffStoermer",                       &ptlCutOffStoermer,                       "ptlCutOffStoermer/F");
  tree->Branch("ptlCutOffDipole",                         &ptlCutOffDipole,                         "ptlCutOffDipole/F");
  tree->Branch("ptlCutOffMax",                            &ptlCutOffMax,                            "ptlCutOffMax[2]/F");
  tree->Branch("showerEnergyD",                           &showerEnergyD,                           "showerEnergyD/F");
  tree->Branch("showerEnergyDL",                          &showerEnergyDL,                          "showerEnergyDL[18]/F");
  tree->Branch("showerEnergyE",                           &showerEnergyE,                           "showerEnergyE/F");
  tree->Branch("showerEnergyCorrected",                   &showerEnergyCorrected,                   "showerEnergyCorrected/F");
  tree->Branch("showerBDT",                               &showerBDT,                               "showerBDT/F");
  tree->Branch("showerCofG",                              &showerCofG,                              "showerCofG[3]/F");
  tree->Branch("showerCofGDist",                          &showerCofGDist,                          "showerCofGDist/F");
  tree->Branch("showerCofGdX",                            &showerCofGdX,                            "showerCofGdX/F");
  tree->Branch("showerCofGdY",                            &showerCofGdY,                            "showerCofGdY/F");
  tree->Branch("tofNCluster",                             &tofNCluster,                             "tofNCluster/I");
  tree->Branch("tofNClusterH",                            &tofNClusterH,                            "tofNClusterH/I");
  tree->Branch("tofNUsedHits",                            &tofNUsedHits,                            "tofNUsedHits/I");
  tree->Branch("tofNUnusedHits",                          &tofNUnusedHits,                          "tofNUnusedHits/I");
  tree->Branch("tofNUsedLayersForQ",                      &tofNUsedLayersForQ,                      "tofNUsedLayersForQ/I");
  tree->Branch("tofBeta",                                 &tofBeta,                                 "tofBeta/F");
  tree->Branch("tofInvBetaErr",                           &tofInvBetaErr,                           "tofInvBetaErr/F");
  tree->Branch("tofMass",                                 &tofMass,                                 "tofMass/F");
  tree->Branch("tofMassError",                            &tofMassError,                            "tofMassError/F");
  tree->Branch("isGoodBeta",                              &isGoodBeta,                              "isGoodBeta/I");
  tree->Branch("isTkTofMatch",                            &isTkTofMatch,                            "isTkTofMatch/I");
  tree->Branch("tofReducedChisqT",                        &tofReducedChisqT,                        "tofReducedChisqT/F");
  tree->Branch("tofReducedChisqC",                        &tofReducedChisqC,                        "tofReducedChisqC/F");
  tree->Branch("tofDepositedEnergyOnLayer",               &tofDepositedEnergyOnLayer,               "tofDepositedEnergyOnLayer[4]/F");
  tree->Branch("tofEstimatedChargeOnLayer",               &tofEstimatedChargeOnLayer,               "tofEstimatedChargeOnLayer[4]/F");
  tree->Branch("tofCharge",                               &tofCharge,                               "tofCharge/F");
  tree->Branch("tofUpperCharge",                          &tofUpperCharge,                          "tofUpperCharge/F");
  tree->Branch("tofLowerCharge",                          &tofLowerCharge,                          "tofLowerCharge/F");
  tree->Branch("tofChargeOnLayer",                        &tofChargeOnLayer,                        "tofChargeOnLayer[4]/F");
  tree->Branch("trkFitCodeFS",                            &trkFitCodeFS,                            "trkFitCodeFS/I");
  tree->Branch("trkRigidityFS",                           &trkRigidityFS,                           "trkRigidityFS/F");
  tree->Branch("trkRigidityInverseErrorFS",               &trkRigidityInverseErrorFS,               "trkRigidityInverseErrorFS/F");
  tree->Branch("trkReducedChisquareFSX",                  &trkReducedChisquareFSX,                  "trkReducedChisquareFSX/F");
  tree->Branch("trkReducedChisquareFSY",                  &trkReducedChisquareFSY,                  "trkReducedChisquareFSY/F");
  tree->Branch("trkFitCodeMS",                            &trkFitCodeMS,                            "trkFitCodeMS/I");
  tree->Branch("trkRigidityMS",                           &trkRigidityMS,                           "trkRigidityMS/F");
  tree->Branch("trkRigidityInverseErrorMS",               &trkRigidityInverseErrorMS,               "trkRigidityInverseErrorMS/F");
  tree->Branch("trkReducedChisquareMSX",                  &trkReducedChisquareMSX,                  "trkReducedChisquareMSX/F");
  tree->Branch("trkReducedChisquareMSY",                  &trkReducedChisquareMSY,                  "trkReducedChisquareMSY/F");
  tree->Branch("trkFitCodeInner",                         &trkFitCodeInner,                         "trkFitCodeInner/I");
  tree->Branch("trkRigidityInner",                        &trkRigidityInner,                        "trkRigidityInner/F");
  tree->Branch("trkRigidityInverseErrorInner",            &trkRigidityInverseErrorInner,            "trkRigidityInverseErrorInner/F");
  tree->Branch("trkReducedChisquareInnerX",               &trkReducedChisquareInnerX,               "trkReducedChisquareInnerX/F");
  tree->Branch("trkReducedChisquareInnerY",               &trkReducedChisquareInnerY,               "trkReducedChisquareInnerY/F");
  tree->Branch("trkEdepLayerJXSideOK",                    &trkEdepLayerJXSideOK,                    "trkEdepLayerJXSideOK[9]/I");
  tree->Branch("trkEdepLayerJYSideOK",                    &trkEdepLayerJYSideOK,                    "trkEdepLayerJYSideOK[9]/I");
  tree->Branch("trkEdepLayerJ",                           &trkEdepLayerJ,                           "trkEdepLayerJ[9]/F");
  tree->Branch("trkCharge",                               &trkCharge,                               "trkCharge/F");
  tree->Branch("trkInnerCharge",                          &trkInnerCharge,                          "trkInnerCharge/F");
  tree->Branch("trkHasExtLayers",                         &trkHasExtLayers,                         "trkHasExtLayers/F");
  tree->Branch("richRebuild",                             &richRebuild,                             "richRebuild/I");
  tree->Branch("richIsGood",                              &richIsGood,                              "richIsGood/I");
  tree->Branch("richIsClean",                             &richIsClean,                             "richIsClean/I");
  tree->Branch("richIsNaF",                               &richIsNaF,                               "richIsNaF/I");
  tree->Branch("richUsedHits",                            &richUsedHits,                            "richUsedHits/I");
  tree->Branch("richRingWidth",                           &richRingWidth,                           "richRingWidth/F");
  tree->Branch("richNHits",                               &richNHits,                               "richNHits/I");
  tree->Branch("richNPMTsOnRing",                         &richNPMTsOnRing,                         "richNPMTsOnRing/I");
  tree->Branch("richBeta",                                &richBeta,                                "richBeta/F");
  tree->Branch("richBetaError",                           &richBetaError,                           "richBetaError/F");
  tree->Branch("richChargeSquared",                       &richChargeSquared,                       "richChargeSquared/F");
  tree->Branch("richKolmogorovProbability",               &richKolmogorovProbability,               "richKolmogorovProbability/F");
  tree->Branch("richPhotoelectrons",                      &richPhotoelectrons,                      "richPhotoelectrons/F");
  tree->Branch("richExpectedPhotoelectrons",              &richExpectedPhotoelectrons,              "richExpectedPhotoelectrons/F");
  tree->Branch("richTheta",                               &richTheta,                               "richTheta/F");
  tree->Branch("richPhi",                                 &richPhi,                                 "richPhi/F");
  tree->Branch("trdNClusters",                            &trdNClusters,                            "trdNClusters/I");
  tree->Branch("trdNUnusedHits",                          &trdNUnusedHits,                          "trdNUnusedHits/I");
  tree->Branch("trdNTracks",                              &trdNTracks,                              "trdNTracks/I");
  tree->Branch("trdTrackTheta",                           &trdTrackTheta,                           "trdTrackTheta/F");
  tree->Branch("trdTrackPhi",                             &trdTrackPhi,                             "trdTrackPhi/F");
  tree->Branch("trdTrackChi2",                            &trdTrackChi2,                            "trdTrackChi2/F");
  tree->Branch("trdTrackPattern",                         &trdTrackPattern,                         "trdTrackPattern/I");
  tree->Branch("trdTrackCharge",                          &trdTrackCharge,                          "trdTrackCharge/F");
  tree->Branch("trdTrackEdepL",                           &trdTrackEdepL,                           "trdTrackEdepL[20]/F");
  tree->Branch("trdDepositedEnergyOnLayer",               &trdDepositedEnergyOnLayer,               "trdDepositedEnergyOnLayer[2]/F");
  tree->Branch("trdTrackMeanDepositedEnergy",             &trdTrackMeanDepositedEnergy,             "trdTrackMeanDepositedEnergy/F");
  tree->Branch("trdTrackTotalDepositedEnergy",            &trdTrackTotalDepositedEnergy,            "trdTrackTotalDepositedEnergy/F");
  /*
  tree->Branch("trdQtNActiveLayer",                       &trdQtNActiveLayer,                       "trdQtNActiveLayer/I");
  tree->Branch("trdQtIsCalibrationGood",                  &trdQtIsCalibrationGood,                  "trdQtIsCalibrationGood/I");
  tree->Branch("trdQtIsSlowControlDataGood",              &trdQtIsSlowControlDataGood,              "trdQtIsSlowControlDataGood/I");
  tree->Branch("trdQtIsInsideTrdGeometricalAcceptance",   &trdQtIsInsideTrdGeometricalAcceptance,   "trdQtIsInsideTrdGeometricalAcceptance/I");
  tree->Branch("trdQtIsValid",                            &trdQtIsValid,                            "trdQtIsValid/I");
  tree->Branch("trdQtNActiveStraws",                      &trdQtNActiveStraws,                      "trdQtNActiveStraws/I");
  tree->Branch("trdQtNActiveLayers",                      &trdQtNActiveLayers,                      "trdQtNActiveLayers/I");
  tree->Branch("trdQtNTRDVertex",                         &trdQtNTRDVertex,                         "trdQtNTRDVertex/I");
  tree->Branch("trdQtElectronToProtonLogLikelihoodRatio", &trdQtElectronToProtonLogLikelihoodRatio, "trdQtElectronToProtonLogLikelihoodRatio/F");
  tree->Branch("trdQtHeliumToElectronLogLikelihoodRatio", &trdQtHeliumToElectronLogLikelihoodRatio, "trdQtHeliumToElectronLogLikelihoodRatio/F");
  tree->Branch("trdQtHeliumToProtonLogLikelihoodRatio",   &trdQtHeliumToProtonLogLikelihoodRatio,   "trdQtHeliumToProtonLogLikelihoodRatio/F");
  */
  tree->Branch("trdKNRawHits",                            &trdKNRawHits,                            "trdKNRawHits/I");
  tree->Branch("trdKIsReadAlignmentOK",                   &trdKIsReadAlignmentOK,                   "trdKIsReadAlignmentOK/I");
  tree->Branch("trdKIsReadCalibOK",                       &trdKIsReadCalibOK,                       "trdKIsReadCalibOK/I");
  tree->Branch("trdKNHits",                               &trdKNHits,                               "trdKNHits/I");
  tree->Branch("trdKIsValid",                             &trdKIsValid,                             "trdKIsValid/I");
  tree->Branch("trdKElectronToProtonLogLikelihoodRatio",  &trdKElectronToProtonLogLikelihoodRatio,  "trdKElectronToProtonLogLikelihoodRatio/F");
  tree->Branch("trdKHeliumToElectronLogLikelihoodRatio",  &trdKHeliumToElectronLogLikelihoodRatio,  "trdKHeliumToElectronLogLikelihoodRatio/F");
  tree->Branch("trdKHeliumToProtonLogLikelihoodRatio",    &trdKHeliumToProtonLogLikelihoodRatio,    "trdKHeliumToProtonLogLikelihoodRatio/F");
  tree->Branch("trdKCharge",                              &trdKCharge,                              "trdCharge/F");
  tree->Branch("trdKChargeError",                         &trdKChargeError,                         "trdKChargeError/F");
  tree->Branch("trdKNUsedHitsForCharge",                  &trdKNUsedHitsForCharge,                  "trdKNUsedHitsForCharge/I");
  tree->Branch("trdKAmpLayer",                            &trdKAmpLayer,                            "trdKAmpLayer[20]/F");
  tree->Branch("trdKTotalPathLength",                     &trdKTotalPathLength,                     "trdKTotalPathLength/F");
  tree->Branch("trdKTotalAmp",                            &trdKTotalAmp,                            "trdKTotalAmp/F");
  tree->Branch("trdKElectronLikelihood",                  &trdKElectronLikelihood,                  "trdKElectronLikelihood/F");
  tree->Branch("trdKProtonLikelihood",                    &trdKProtonLikelihood,                    "trdKProtonLikelihood/F");
  tree->Branch("trdKHeliumLikelihood",                    &trdKHeliumLikelihood,                    "trdKHeliumLikelihood/F");
  tree->Branch("trdPElectronToProtonLogLikelihoodRatio",  &trdPElectronToProtonLogLikelihoodRatio,  "trdPElectronToProtonLogLikelihoodRatio/F");
  tree->Branch("trdPHeliumToElectronLogLikelihoodRatio",  &trdPHeliumToElectronLogLikelihoodRatio,  "trdPHeliumToElectronLogLikelihoodRatio/F");
  tree->Branch("trdPHeliumToProtonLogLikelihoodRatio",    &trdPHeliumToProtonLogLikelihoodRatio,    "trdPHeliumToProtonLogLikelihoodRatio/F");
  tree->Branch("accNHits",                                &accNHits,                                "accNHits/I");
//  tree->Branch("accNRecoClPG", &accNRecoClPG, "accNRecoClPG/I");
//  tree->Branch("accSector", &accSector);
  tree->Branch("accTime", &accTime);
//  tree->Branch("accHitPosZ", &accHitPosZ);
//  tree->Branch("accChi2", &accChi2);
//  tree->Branch("accNPairs", &accNPairs);
//  tree->Branch("accUnfoldedHitPosZ", &accUnfoldedHitPosZ);
//  tree->Branch("accRawCharge", &accRawCharge);
//  tree->Branch("accHitPosZfromADC", &accHitPosZfromADC);
//  tree->Branch("accUnfoldedHitPosZfromADC", &accUnfoldedHitPosZfromADC);
//  tree->Branch("accCrossingPhiAngle", &accCrossingPhiAngle, "accCrossingPhiAngle/F");
  tree->Branch("accEdep", &accEdep);
//  tree->Branch("accTimePG", &accTimePG);
  /*}}}*/

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //
   // Begin of the event loop !!
   //
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  unsigned int nProcessCheck = 10000;

  for(unsigned int e = 0; e < nEntries; e++)
  {
    AMSEventR* pev   = NULL;
    HeaderR* header = NULL;
    pev = amsChain->GetEvent(e);
    header = &(pev->fHeader);

    // Basic cut processes
    if( IsBadRun(pev) ) { printf("Run number [%d] is bad run. Skip this run.", pev->Run()); break; }/*{{{*/
    if( !pev ) continue;
    hEvtCounter->Fill(0);
    if( !IsScienceRun(pev) ) continue;
    hEvtCounter->Fill(0);
    if( !IsHardwareStatusGood(pev) ) continue;
    hEvtCounter->Fill(1);
    if( IsUnbiasedPhysicsTriggerEvent(pev) ) continue;
    hEvtCounter->Fill(2);
    // Here we need to think about how to select a good particle among several particles.
    // In the study of deuteron flux, a particle need to be defined by the following several variables.
    //
    // 1. Momentum
    // 2. Charge
    //
    // To measure particle momentum correctly, track should be measured correctly.
    //
    // 04-02-15 : Currently, just use single particle case.
    //
    // If ParticleR->pTrTrack is empty, discard the event.
    if( pev->nParticle() != 1 ) continue;
    hEvtCounter->Fill(3);
    if( !IsTrkAlignmentGood(pev) ) continue;
    hEvtCounter->Fill(4);
    int iParticle = GetGoodParticleIndex(pev);
    if( iParticle == -1 ) continue;
    hEvtCounter->Fill(6);/*}}}*/

    // Save data
    nRun            = pev->Run();
    nEvent          = pev->Event();
    nLevel1         = pev->nLevel1();
    nParticle       = pev->nParticle();
    nCharge         = pev->nCharge();
    nTrTrack        = pev->nTrTrack();
    nTrdTrack       = pev->nTrdTrack();
    nAntiCluster    = pev->nAntiCluster();
    nRichRing       = pev->nRichRing();
    nRichRingB      = pev->nRichRingB();
    nBeta           = pev->nBeta();
    nBetaB          = pev->nBetaB();
    nBetaH          = pev->nBetaH();
    nShower         = pev->nEcalShower();
    nVertex         = pev->nVertex();
    livetime        = pev->LiveTime();
    utcTime         = (float)header->UTCTime(0);
    orbitAltitude   = header->RadS;
    orbitLatitude   = header->ThetaS;
    orbitLongitude  = header->PhiS;
    orbitLatitudeM  = header->ThetaM;
    orbitLongitudeM = header->PhiM;
    velR            = header->VelocityS;
    velTheta        = header->VelTheta;
    velPhi          = header->VelPhi;
    yaw             = header->Yaw;
    pitch           = header->Pitch;
    roll            = header->Roll;
    zenithAngle     = header->Zenith();
    isInSAA         = pev->IsInSAA();

    double tmpglong;
    double tmpglat;
    gCoordCalcResult  = pev->GetGalCoo(
        gCoordCalcResultBit,
        tmpglong,
        tmpglat);
    gLongitude  = tmpglong;
    gLatitude   = tmpglat;

    double tmpSunPosAzimuth;
    double tmpSunPosElevation;
    sunPosCalcResult   = header->getSunAMS(tmpSunPosAzimuth, tmpSunPosElevation);
    sunPosAzimuth      = tmpSunPosAzimuth;
    sunPosElevation    = tmpSunPosElevation;

    AMSPoint solarArray;
    isInShadow         = pev->isInShadow(solarArray, 0);
    solarArrayCoord[0] = solarArray.x();
    solarArrayCoord[1] = solarArray.y();
    solarArrayCoord[2] = solarArray.z();

    ParticleR* pParticle = pev->pParticle(iParticle);
    if(!pParticle) continue;
    ptlCharge           = (unsigned int)pParticle->Charge;
    // Unfold below to see the particle coordinates
    ptlMomentum         = pParticle->Momentum;/*{{{*/
    ptlTheta            = pParticle->Theta;
    ptlPhi              = pParticle->Phi;
    ptlCoo[0]           = pParticle->Coo[0];
    ptlCoo[1]           = pParticle->Coo[1];
    ptlCoo[2]           = pParticle->Coo[2];
    ptlTofCoo[0][0]     = pParticle->TOFCoo[0][0];
    ptlTofCoo[0][1]     = pParticle->TOFCoo[0][1];
    ptlTofCoo[0][2]     = pParticle->TOFCoo[0][2];
    ptlTofCoo[1][0]     = pParticle->TOFCoo[1][0];
    ptlTofCoo[1][1]     = pParticle->TOFCoo[1][1];
    ptlTofCoo[1][2]     = pParticle->TOFCoo[1][2];
    ptlTofCoo[2][0]     = pParticle->TOFCoo[2][0];
    ptlTofCoo[2][1]     = pParticle->TOFCoo[2][1];
    ptlTofCoo[2][2]     = pParticle->TOFCoo[2][2];
    ptlTofCoo[3][0]     = pParticle->TOFCoo[3][0];
    ptlTofCoo[3][1]     = pParticle->TOFCoo[3][1];
    ptlTofCoo[3][2]     = pParticle->TOFCoo[3][2];
    ptlTofTrLength[0]   = pParticle->TOFTLength[0];
    ptlTofTrLength[1]   = pParticle->TOFTLength[1];
    ptlTofTrLength[2]   = pParticle->TOFTLength[2];
    ptlTofTrLength[3]   = pParticle->TOFTLength[3];
    ptlAntiCoo[0][0]    = pParticle->AntiCoo[0][0];
    ptlAntiCoo[0][1]    = pParticle->AntiCoo[0][1];
    ptlAntiCoo[0][2]    = pParticle->AntiCoo[0][2];
    ptlAntiCoo[0][3]    = pParticle->AntiCoo[0][3];
    ptlAntiCoo[0][4]    = pParticle->AntiCoo[0][4];
    ptlAntiCoo[1][0]    = pParticle->AntiCoo[1][0];
    ptlAntiCoo[1][1]    = pParticle->AntiCoo[1][1];
    ptlAntiCoo[1][2]    = pParticle->AntiCoo[1][2];
    ptlAntiCoo[1][3]    = pParticle->AntiCoo[1][3];
    ptlAntiCoo[1][4]    = pParticle->AntiCoo[1][4];
    ptlEcalCoo[0][0]    = pParticle->EcalCoo[0][0];
    ptlEcalCoo[0][1]    = pParticle->EcalCoo[0][1];
    ptlEcalCoo[0][2]    = pParticle->EcalCoo[0][2];
    ptlEcalCoo[1][0]    = pParticle->EcalCoo[1][0];
    ptlEcalCoo[1][1]    = pParticle->EcalCoo[1][1];
    ptlEcalCoo[1][2]    = pParticle->EcalCoo[1][2];
    ptlEcalCoo[2][0]    = pParticle->EcalCoo[2][0];
    ptlEcalCoo[2][1]    = pParticle->EcalCoo[2][1];
    ptlEcalCoo[2][2]    = pParticle->EcalCoo[2][2];
    ptlTrCoo[0][0]      = pParticle->TrCoo[0][0];
    ptlTrCoo[0][1]      = pParticle->TrCoo[0][1];
    ptlTrCoo[0][2]      = pParticle->TrCoo[0][2];
    ptlTrCoo[1][0]      = pParticle->TrCoo[1][0];
    ptlTrCoo[1][1]      = pParticle->TrCoo[1][1];
    ptlTrCoo[1][2]      = pParticle->TrCoo[1][2];
    ptlTrCoo[2][0]      = pParticle->TrCoo[2][0];
    ptlTrCoo[2][1]      = pParticle->TrCoo[2][1];
    ptlTrCoo[2][2]      = pParticle->TrCoo[2][2];
    ptlTrCoo[3][0]      = pParticle->TrCoo[3][0];
    ptlTrCoo[3][1]      = pParticle->TrCoo[3][1];
    ptlTrCoo[3][2]      = pParticle->TrCoo[3][2];
    ptlTrCoo[4][0]      = pParticle->TrCoo[4][0];
    ptlTrCoo[4][1]      = pParticle->TrCoo[4][1];
    ptlTrCoo[4][2]      = pParticle->TrCoo[4][2];
    ptlTrCoo[5][0]      = pParticle->TrCoo[5][0];
    ptlTrCoo[5][1]      = pParticle->TrCoo[5][1];
    ptlTrCoo[5][2]      = pParticle->TrCoo[5][2];
    ptlTrCoo[6][0]      = pParticle->TrCoo[6][0];
    ptlTrCoo[6][1]      = pParticle->TrCoo[6][1];
    ptlTrCoo[6][2]      = pParticle->TrCoo[6][2];
    ptlTrCoo[7][0]      = pParticle->TrCoo[7][0];
    ptlTrCoo[7][1]      = pParticle->TrCoo[7][1];
    ptlTrCoo[7][2]      = pParticle->TrCoo[7][2];
    ptlTrCoo[8][0]      = pParticle->TrCoo[8][0];
    ptlTrCoo[8][1]      = pParticle->TrCoo[8][1];
    ptlTrCoo[8][2]      = pParticle->TrCoo[8][2];
    ptlTrdCoo[0][0]     = pParticle->TRDCoo[0][0];
    ptlTrdCoo[0][1]     = pParticle->TRDCoo[0][1];
    ptlTrdCoo[0][2]     = pParticle->TRDCoo[0][2];
    ptlTrdCoo[1][0]     = pParticle->TRDCoo[1][0];
    ptlTrdCoo[1][1]     = pParticle->TRDCoo[1][1];
    ptlTrdCoo[1][2]     = pParticle->TRDCoo[1][2];
    ptlRichCoo[0][0]    = pParticle->RichCoo[0][0];
    ptlRichCoo[0][1]    = pParticle->RichCoo[0][1];
    ptlRichCoo[0][2]    = pParticle->RichCoo[0][2];
    ptlRichCoo[1][0]    = pParticle->RichCoo[1][0];
    ptlRichCoo[1][1]    = pParticle->RichCoo[1][1];
    ptlRichCoo[1][2]    = pParticle->RichCoo[1][2];
    ptlRichPath[0]      = pParticle->RichPath[0];
    ptlRichPath[1]      = pParticle->RichPath[1];
    ptlRichPathBeta[0]  = pParticle->RichPathBeta[0];
    ptlRichPathBeta[1]  = pParticle->RichPathBeta[1];
    ptlRichLength       = pParticle->RichLength;
    ptlRichParticles    = pParticle->RichParticles;
    ptlCutOffStoermer   = pParticle->CutoffS;
    ptlCutOffDipole     = pParticle->Cutoff;/*}}}*/

    tofNCluster         = pev->nTofCluster();
    tofNClusterH        = pev->nTofClusterH();
    BetaHR* pBeta       = pParticle->pBetaH();
    if(!pBeta) continue;
    tofBeta             = pBeta->GetBeta();
    tofInvBetaErr       = pBeta->GetEBetaV();
    tofMass             = pBeta->GetMass();
    tofMassError        = pBeta->GetEMass();
    tofNUsedHits        = pBeta->NTofClusterH();
    tofNUnusedHits      = tofNCluster - tofNUsedHits;
    if( pBeta->IsGoodBeta()   == true ) isGoodBeta = 1;
    else isGoodBeta     = 0;
    if( pBeta->IsTkTofMatch() == true ) isTkTofMatch = 1;
    else isTkTofMatch   = 0;
    tofReducedChisqT    = pBeta->GetNormChi2T();
    tofReducedChisqC    = pBeta->GetNormChi2C();
    pBeta->GetQ(tofNUsedLayersForQ, tofCharge);
    tofChargeOnLayer[0] = pBeta->GetQL(0);
    tofChargeOnLayer[1] = pBeta->GetQL(1);
    tofChargeOnLayer[2] = pBeta->GetQL(2);
    tofChargeOnLayer[3] = pBeta->GetQL(3);

    for(int i = 0; i < pBeta->NTofClusterH(); ++i)
    {
      if( pBeta->TestExistHL(i) == true )
        tofDepositedEnergyOnLayer[i] = pBeta->GetClusterHL(i)->GetEdep();
      else
        tofDepositedEnergyOnLayer[i] = -1;
    }

    int ncls[4] = {0, 0, 0, 0};
    nTofClustersInTime  = pev->GetNTofClustersInTime(pBeta, ncls);

    TrTrackR* pTrTrack  = pParticle->pTrTrack();
    if(!pTrTrack) continue;
    int id_fullspan     = pTrTrack->iTrTrackPar(1, 0, 0);
    int id_inner        = pTrTrack->iTrTrackPar(1, 3, 20);
    int id_maxspan      = pTrTrack->iTrTrackPar(1, 0, 20);

    if( id_fullspan < 0 ) { cerr << "Warning! - Event : " << e << " / id_fullspan = " << id_fullspan << endl; }

    // Tracker variables from full span setting
    if( id_fullspan >= 0)
    {
      trkFitCodeFS              = id_fullspan;
      trkRigidityFS             = pTrTrack->GetRigidity(id_fullspan);
      trkRigidityInverseErrorFS = pTrTrack->GetErrRinv(id_fullspan);
      trkReducedChisquareFSX    = pTrTrack->GetNormChisqX(id_fullspan);
      trkReducedChisquareFSY    = pTrTrack->GetNormChisqY(id_fullspan);
    }
    else
    {
      trkFitCodeFS              = id_fullspan;
      trkRigidityFS             = -99999.;
      trkRigidityInverseErrorFS = -99999.;
      trkReducedChisquareFSX    = -99999.;
      trkReducedChisquareFSY    = -99999.;
    }

    // Tracker variables from maximum span setting
    if( id_maxspan >= 0 )
    {
      trkFitCodeMS              = id_maxspan;
      trkRigidityMS             = pTrTrack->GetRigidity(id_maxspan);
      trkRigidityInverseErrorMS = pTrTrack->GetErrRinv(id_maxspan);
      trkReducedChisquareMSX    = pTrTrack->GetNormChisqX(id_maxspan);
      trkReducedChisquareMSY    = pTrTrack->GetNormChisqY(id_maxspan);
    }
    else
    {
      trkFitCodeMS              = id_maxspan;
      trkRigidityMS             = -99999.;
      trkRigidityInverseErrorMS = -99999.;
      trkReducedChisquareMSX    = -99999.;
      trkReducedChisquareMSY    = -99999.;
    }

    // Tracker variables from inner tracker only setting
    if( id_inner >= 0 )
    {
      trkFitCodeInner              = id_inner;
      trkRigidityInner             = pTrTrack->GetRigidity(id_inner);
      trkRigidityInverseErrorInner = pTrTrack->GetErrRinv(id_inner);
      trkReducedChisquareInnerX    = pTrTrack->GetNormChisqX(id_inner);
      trkReducedChisquareInnerY    = pTrTrack->GetNormChisqY(id_inner);
    }
    else
    {
      trkFitCodeInner              = id_inner;
      trkRigidityInner             = -99999.;
      trkRigidityInverseErrorInner = -99999.;
      trkReducedChisquareInnerX    = -99999.;
      trkReducedChisquareInnerY    = -99999.;
    }

    TrRecHitR* pTrRecHit = NULL;          // This should be ParticleR associated hit.
    for(int ii = 0; ii < 9; ii++) trkEdepLayerJ[ii] = 0.;

    for(int ilayer = 0; ilayer < 9; ilayer++)
    {
      pTrRecHit = pTrTrack->GetHitLJ(ilayer);
      if(!pTrRecHit) continue;
      if(pTrRecHit->GetEdep(0) != 0) trkEdepLayerJXSideOK[ilayer] = 1;
      else trkEdepLayerJXSideOK[ilayer] = 0;
      if(pTrRecHit->GetEdep(1) != 0) trkEdepLayerJYSideOK[ilayer] = 1;
      else trkEdepLayerJYSideOK[ilayer] = 0;
      trkEdepLayerJ[ilayer] += pTrRecHit->GetEdep(0) + pTrRecHit->GetEdep(1);
    }
    trkCharge       = pTrTrack->GetQ();
    trkInnerCharge  = pTrTrack->GetInnerQ();
    trkHasExtLayers = pTrTrack->HasExtLayers();

    TrdTrackR* pTrdTrack = pParticle->pTrdTrack();
    if( !pTrdTrack ) continue;
    /*
    int nTrdSegment;
    nTrdSegment = pTrdTrack->NTrdSegment();
    */

    //for(int i =  0; i < nTrdSegment; i++)
    //  TrdSegmentR* pTrdSegment = pTrdTrack->pTrdSegment(i);

    // Initialize ECAL shower related variables
    showerEnergyD         = -1;
    for(int i = 0; i < 18; ++i) showerEnergyDL[i] = 0;
    showerEnergyE         = -1;
    showerEnergyCorrected = -1;
    showerBDT             = -1;
    showerCofG[0]         = -99;
    showerCofG[1]         = -99;
    showerCofG[2]         = -99;
    showerCofGDist        = -99;
    showerCofGdX          = -99;
    showerCofGdY          = -99;

    // Save ECAL shower related variables
    EcalShowerR* amsEcalShower = pParticle->pEcalShower();
    if( amsEcalShower )
    {
      showerEnergyCorrected = amsEcalShower->GetCorrectedEnergy(2, 2);
      if( showerEnergyCorrected < 0.5 ) continue;
      showerEnergyD = (amsEcalShower->EnergyD)/1000.;
      for(int i = 0; i < amsEcalShower->NEcal2DCluster(); ++i)
      {
        for(int j = 0; j < amsEcalShower->pEcal2DCluster(i)->NEcalCluster(); ++j)
        {
          if (amsEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Edep > 0)
            showerEnergyDL[amsEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Plane] += amsEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Edep;
          else
            showerEnergyDL[amsEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Plane] = -1;
        }
      }
      showerEnergyE = amsEcalShower->EnergyE;
      showerBDT     = amsEcalShower->GetEcalBDT();
      showerCofG[0] = amsEcalShower->CofG[0];
      showerCofG[1] = amsEcalShower->CofG[1];
      showerCofG[2] = amsEcalShower->CofG[2];
      AMSPoint trackPoint;
      AMSDir   trackDirection;
      pTrTrack->Interpolate(amsEcalShower->CofG[2], trackPoint, trackDirection, id_fullspan);
      amsEcalShower->NormaliseVariableLAPP();
      showerCofGDist = std::sqrt(std::pow(trackPoint[0] - showerCofG[0], 2) + std::pow( trackPoint[1] - showerCofG[1], 2));
      showerCofGdX   = std::fabs(trackPoint[0] - showerCofG[0]);
      showerCofGdY   = std::fabs(trackPoint[1] - showerCofG[1]);
    }

    // ACSoft related lines
    particleFactory.SetAMSTrTrackR(pTrTrack);
    particleFactory.SetAMSEcalShowerR(amsEcalShower);
    particleFactory.SetAMSBetaHR(pBeta);

    if( !amsRootSupport.SwitchToSpecificTrackFitById(id_maxspan) ) continue;
    Analysis::Event& event = amsRootSupport.BuildEvent(amsChain, pev);

    // Only do this if you need access to TRD segments/tracks and vertices
    eventFactory.PerformTrdTracking(event);
    eventFactory.PerformTrdVertexFinding(event);

    // If you want to use TrdQt, you should always add Analysis::CreateSplineTrack, so it can use the
    // Tracker track extrapolation to pick up the TRD hits in the tubes and calculate the path length
    // and also Analysis::FillTrdQt so that the likelihoods are calculated.
    // Additionally you shold pass Analysis::CreateTrdTrack if you want to examinge the TRD tracks
    // found by the previous PerformTrdTracking() call, via the Analysis::Particle interface.

    int productionSteps = 0;
    productionSteps |= Analysis::CreateSplineTrack;
    productionSteps |= Analysis::CreateTrdTrack;
//    productionSteps |= Analysis::FillTrdQt;
    eventFactory.FillParticles(event, productionSteps);

    const Analysis::Particle* particle = event.PrimaryParticle();
    assert(particle);

    tofUpperCharge = particle->UpperTofCharge();
    tofLowerCharge = particle->LowerTofCharge();

    trdPElectronToProtonLogLikelihoodRatio = particle->CalculateElectronProtonLikelihood();
    trdPHeliumToProtonLogLikelihoodRatio   = particle->CalculateHeliumProtonLikelihood();
    trdPHeliumToElectronLogLikelihoodRatio = particle->CalculateHeliumElectronLikelihood();

    RichRingR* richRing = pParticle->pRichRing();
    if(!richRing)
    {
      richRebuild               = -1;
      richIsGood                = -1;
      richIsClean               = -1;
      richIsNaF                 = -1;
      richRingWidth             = -1;
      richUsedHits              = -1;
      richNHits                 = -1;
      richNPMTsOnRing           = -1;
      richBeta                  = -1;
      richBetaError             = -1;
      richChargeSquared         = -1;
      richKolmogorovProbability = -1;
      richPhotoelectrons        = -9999;
      richExpectedPhotoelectrons = -9999;
      richTheta                 = -9;
      richPhi                   = -9;
    }
    else
    {
      richRebuild                = (int)richRing->Rebuild();
      richIsGood                 = (int)richRing->IsGood();
      richIsClean                = (int)richRing->IsClean();
      richIsNaF                  = (int)richRing->IsNaF();
      richUsedHits               = (int)richRing->getUsedHits();
      richRingWidth              = (float)richRing->RingWidth();
      richNHits                  = richRing->getHits();
      richNPMTsOnRing            = richRing->getPMTs();
      richBeta                   = richRing->getBeta();
      richBetaError              = richRing->getBetaError();
      richChargeSquared          = richRing->getCharge2Estimate();
      richKolmogorovProbability  = richRing->getProb();
      richPhotoelectrons         = richRing->getPhotoElectrons();
      richExpectedPhotoelectrons = richRing->getExpectedPhotoelectrons();
      richTheta                  = richRing->getTrackTheta();
      richPhi                    = richRing->getTrackPhi();
    }

    TrdTrackR* trdTrack = pParticle->pTrdTrack();
    if(!trdTrack) continue;
    trdNClusters = pev->nTrdCluster();
    int trdNUsedHits = 0;
    int trdNUsedSegment = trdTrack->NTrdSegment();
    for(int i = 0; i < trdNUsedSegment; i++)
    {
      TrdSegmentR* pTrdSegment = trdTrack->pTrdSegment(i);
      trdNUsedHits += pTrdSegment->NTrdCluster();
    }
    trdNUnusedHits = trdNClusters - trdNUsedHits;

    trdNTracks  = pev->nTrdTrack();
    if(!trdTrack)
    {
      trdTrackTheta     = -9.;
      trdTrackPhi       = -9.;
      trdTrackPattern   = -9;
      trdTrackCharge    = -9;
      for(int i = 0; i < 20; ++i) trdTrackEdepL[i] = -9.;
    }
    else
    {
      trdTrackTheta   = trdTrack->Theta;
      trdTrackPhi     = trdTrack->Phi;
      trdTrackPattern = trdTrack->Pattern;
      trdTrackCharge  = trdTrack->Q;
      trdTrackChi2    = trdTrack->Chi2;
      trdTrackMeanDepositedEnergy = 0.;
      int ntrdlayers = 0;
      for(int i = 0; i < trdTrack->NTrdSegment(); i++)
      {
        for(int j = 0; j < trdTrack->pTrdSegment(i)->NTrdCluster(); j++)
        {
          int trdLayer = trdTrack->pTrdSegment(i)->pTrdCluster(j)->Layer;
          float trdEdep = trdTrack->pTrdSegment(i)->pTrdCluster(j)->EDep;
          trdTrackEdepL[trdLayer] = trdEdep;
          trdTrackTotalDepositedEnergy += trdEdep;
          ntrdlayers++;
        }
      }
      trdTrackMeanDepositedEnergy = trdTrackTotalDepositedEnergy/(float)ntrdlayers;
    }

    /*
    const Analysis::TrdQt* trdQtFromTrackerTrack        = particle->GetTrdQtInfo();
    if(!trdQtFromTrackerTrack) continue;
    trdQtIsCalibrationGood                              = trdQtFromTrackerTrack->IsCalibrationGood();
    trdQtIsSlowControlDataGood                          = trdQtFromTrackerTrack->IsSlowControlDataGood();
    trdQtIsInsideTrdGeometricalAcceptance               = kTRUE;
    trdQtIsValid                                        = 1;
    trdQtNActiveStraws                                  = 1;
    trdQtNActiveLayers                                  = 1;
    trdQtNTRDVertex                                     = 0;
    trdQtElectronToProtonLogLikelihoodRatio             = -1.;
    trdQtHeliumToElectronLogLikelihoodRatio             = -1.;
    trdQtHeliumToProtonLogLikelihoodRatio               = -1.;

    const std::vector<Analysis::TrdHit>& TrdHit         = trdQtFromTrackerTrack->GetAssignedHits();
    const std::vector<Analysis::TrdVertex>& verticesXZ  = event.TrdVerticesXZ();
    const std::vector<Analysis::TrdVertex>& verticesYZ  = event.TrdVerticesYZ();

    for( std::vector<Analysis::TrdVertex>::const_iterator xzIter = verticesXZ.begin(); xzIter != verticesXZ.end(); ++xzIter)
    {
      const Analysis::TrdVertex& xzVertex = *xzIter;
      for( std::vector<Analysis::TrdVertex>::const_iterator yzIter = verticesYZ.begin(); yzIter != verticesYZ.end(); ++yzIter)
      {
        const Analysis::TrdVertex& yzVertex = *yzIter;
        if( std::max(xzVertex.NumberOfSegments(), yzVertex.NumberOfSegments() ) < 3)
          continue;
        if( std::fabs( xzVertex.Z() - yzVertex.Z() ) < std::fabs( xzVertex.ErrorZ() + yzVertex.ErrorZ() ) )
        {
          trdQtNTRDVertex++;
        }
      }
    }

    // Below of this for TrdQt
    if ( const Analysis::TrdQt* trdQtFromTrackerTrack = particle->GetTrdQtInfo() )
    {
      // TrdQt likelihood variables based on ECAL energy scale.
      const_cast<Analysis::TrdQt*>(trdQtFromTrackerTrack)->SetRigidity(showerEnergyCorrected);
      trdQtElectronToProtonLogLikelihoodRatio = trdQtFromTrackerTrack->LogLikelihoodRatioElectronProton();
      trdQtHeliumToElectronLogLikelihoodRatio = trdQtFromTrackerTrack->LogLikelihoodRatioHeliumElectron();
      trdQtHeliumToProtonLogLikelihoodRatio   = trdQtFromTrackerTrack->LogLikelihoodRatioHeliumProton();
    }
    */

    // Below of this for TrdK
    trdKNRawHits                                = -1;
    trdKIsReadAlignmentOK                       = -1;
    trdKIsReadCalibOK                           = -1;
    trdKNHits                                   = -1;
    trdKIsValid                                 = -1;
    trdKElectronToProtonLogLikelihoodRatio      = -1;
    trdKHeliumToElectronLogLikelihoodRatio      = -1;
    trdKHeliumToProtonLogLikelihoodRatio        = -1;
    trdKCharge                                  = -1;
    trdKChargeError                             = -1;
    trdKNUsedHitsForCharge                      = -1;
    for(int k = 0; k < 20; k++) trdKAmpLayer[k] = -1;
    trdKTotalPathLength                         = -1;
    trdKElectronLikelihood                      = -1;
    trdKProtonLikelihood                        = -1;
    trdKHeliumLikelihood                        = -1;

    trdKNRawHits = pev->NTrdRawHit();
    if(trdKNRawHits > 0)
    {
      double trdKLikelihoodRatio[3] = {-1., -1., -1.};
      double trdKLikelihood[3]      = {-1., -1., -1.};

      TrdKCluster* trdK = new TrdKCluster(pev, pTrTrack, id_fullspan);
      if(!trdK) continue;

      trdKIsReadAlignmentOK = trdK->IsReadAlignmentOK;
      trdKIsReadCalibOK     = trdK->IsReadCalibOK;
      trdKIsValid           = trdK->GetLikelihoodRatio_TrTrack(15, trdKLikelihoodRatio, trdKLikelihood, trdKNHits, trdKTotalPathLength, trdKTotalAmp, -1, 0);
      if(trdKIsValid != 0 && trdKNHits != 0)
      {
        trdK->CalculateTRDCharge();
        trdKCharge                              = trdK->GetTRDCharge();
        trdKChargeError                         = trdK->GetTRDChargeError();
        trdKNUsedHitsForCharge                  = trdK->GetQTRDHitCollectionNuclei().size();
        trdKElectronToProtonLogLikelihoodRatio  = trdKLikelihoodRatio[0];
        trdKHeliumToElectronLogLikelihoodRatio  = trdKLikelihoodRatio[1];
        trdKHeliumToProtonLogLikelihoodRatio    = trdKLikelihoodRatio[2];
        trdKElectronLikelihood                  = trdKLikelihood[0];
        trdKProtonLikelihood                    = trdKLikelihood[1];
        trdKHeliumLikelihood                    = trdKLikelihood[2];

        AMSPoint trExtraP0;
        AMSDir   trExtraDir;
        trdK->GetTrTrackExtrapolation(trExtraP0, trExtraDir);

        for(int l = 0; l < trdK->NHits(); l++)
        {
          TrdKHit* trdKHit = trdK->GetHit(l);
          int tmpL = 0;
          tmpL = trdKHit->TRDHit_Layer;
          float tmpAmp = 0;
          tmpAmp = trdKHit->TRDHit_Amp;
          float tmpPathLength = 0;
          tmpPathLength = trdKHit->Tube_Track_3DLength(&trExtraP0, &trExtraDir);

          if( tmpAmp < 15 || tmpPathLength <= 0 ) continue;
          if( trdKHit->IsAligned == 0 || trdKHit->IsCalibrated == 0 ) continue;
          trdKAmpLayer[tmpL] += tmpAmp;
        }
      }
    }

    // The following code lines are about ACC

    pev->RebuildAntiClusters();
    accNHits = pev->nAntiCluster();

    if( accNHits != 0 )
    {
      for( int i = 0; i < accNHits; ++i)
      {
        AntiClusterR* anti = pev->pAntiCluster(i);
        if( !anti ) continue;
        if( anti->Edep != 0 )
          accEdep.push_back(anti->Edep);
        else
          accEdep.push_back(-1);
        accTime.push_back(anti->time);
      }
    }


    if( e % nProcessCheck == 0 || e == nEntries - 1 )
      cout << "[" << releaseName << "] Processed " << e << " out of " << nEntries << " (" << (float)e/nEntries*100. << "%)" << endl;

    nProcessedNumber = nProcessed;
    tree->Fill();

    accSector.clear();
    accTime.clear();
    accHitPosZ.clear();
    accChi2.clear();
    accNPairs.clear();
    accUnfoldedHitPosZ.clear();
    accRawCharge.clear();
    accHitPosZfromADC.clear();
    accUnfoldedHitPosZfromADC.clear();
    accEdep.clear();
    accTimePG.clear();

    nProcessed++;
  }

  if( resultFile->Write() ) cout << "[" << releaseName << "] The result file [" << resultFile->GetName() << "] is successfully written." << endl;
  resultFile->cd("/");
  if( hEvtCounter->Write() ) cout << "[" << releaseName << "] The counter histogram is successfully written." << endl;
  resultFile->Close();

  cout << "[" << releaseName << "] The program is terminated successfully. " << nProcessed << " events are stored." << endl;
  return 0;
}

bool RegisterFileToChain(int argc, char* argv[], AMSChain* amsChain)
{
  // Unfold below to see the details of file opening
  /*{{{*/
  /*************************************************************************************************************
   *
   * FILE OPENING PHASE
   *
   *************************************************************************************************************/

  /*
   * 1. Test mode
   *   1-1. Skirmish test
   *     Use the designated run to test the code. No arguments are required. ( argc == 1)
   *   1-2. Skirmish test within given number of events
   *     User should pass the number of events to be analyzed as the program argument.
   *     ./a.out <number of events to be analyzed>  ( argc == 2 )
   *   1-3. Test using user-given run file within given number of events
   *     ./a.out <run file to be analyzed> <output file> <number of events to be analyzed> ( argc == 4 )
   *
   * 2. Job submission mode
   *   2-1. Standard job submission mode
   *     ./a.out <list file which contains list of files to be analyzed> <output file> ( argc == 3 )
   */

  //char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/Data/AMS02/2011B/ISS.B620/pass4/1323051106.00000001.root";   // Path of test run
  char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1323051106.00000001.root";   // Path of test run
  //char skirmishRunPath[] = "1323051106.00000001.root";
  //char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/MC/AMS02/2014/d.B1030/d.pl1.0_520_GG_Blic/872728226.00000001.root";
  char inputFileName[256];      // File path for single run. (Cat 3.)

  if( argc == 1 )
  {
    std::cout << "[" << releaseName << "] RUN MODE : Single Test Run (Cat. 1)" << endl;

    if(amsChain->Add(skirmishRunPath) != 1)
    {
      std::cerr << "[" << releaseName << "] ERROR    : File open error, [" << skirmishRunPath << "] can not be found!" << endl;
      return false;
    }
    else return true;
  }
  else if( argc == 2 ) // Process as many events as user requested. (second argument is number of events to be processed)
  {
    std::cout << "[" << releaseName << "] RUN MODE : Single Test Run (Cat. 2)" << endl;

    if(amsChain->Add(skirmishRunPath) != 1)
    {
      std::cerr << "[" << releaseName << "] ERROR    : File open error, [" << skirmishRunPath << "] can not be found!" << endl;
      return false;
    }
    else return true;
  }
  else if( argc == 3 )
  {
    std::cout << "[" << releaseName << "] RUN MODE : Batch-job Mode" << endl;

    char listFileName[256];       // The name of the list file.
    char inputFileName[256];

    FILE* fp;                     // File pointer to read list file.

    strcpy(listFileName, argv[1]);     // First argument is list file

    if( ( fp = fopen( listFileName, "r") ) == NULL )
    {
      std::cerr << "[" << releaseName << "] ERROR     : Failed to open file [" << listFileName << "]!" << endl;
      return false;
    }

    char* line_p;                 // Character pointer to filter out CRLF at the end of each line.
    while( fgets( inputFileName, 256, fp ) != NULL )
    {
      if( ( line_p = strchr(inputFileName, '\n') ) != NULL) *line_p = 0;  // Filter out \n at the end of lines.

      if(amsChain->Add( inputFileName ) != 1)
      {
        std::cerr << "[" << releaseName << "] ERROR     : Failed to open file [" << inputFileName << "]!" << endl;
        return false;
      }
      else
      {
        std::cout << "[" << releaseName << "] The file [" << inputFileName << "] is added to the chain." << std::endl;
        std::cout << "[" << releaseName << "] Currently loaded events : " << amsChain->GetEntries() << std::endl;
      }
    }

    fclose(fp);
    return true;

    //nEntries = amsChain.GetEntries();
  }
  else if( argc >= 4 )
  {
    std::cout << "[" << releaseName << "] RUN MODE : Single Test Run (Cat. 3)" << endl;

    strcpy(inputFileName, argv[1]);

    if(amsChain->Add(inputFileName) != 1)
    {
      std::cerr << "[" << releaseName << "] ERROR   : File open error, [" << inputFileName << "] can not found!" << endl;
      return false;
    }
    else
    {
      std::cout << "[" << releaseName << "] The file [" << inputFileName << "] is added to the chain." << std::endl;
      return true;
    }
  }/*}}}*/

  return false;
}
