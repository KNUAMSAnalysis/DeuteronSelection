/**
 * @file    main.h
 * @brief   This is a header file for main.cxx
 * @author  Wooyoung Jang (wyjang)
 * @date    23. Oct. 2015
 */

char releaseName[16];
char softwareName[10] = "ACCTOFInt";
char versionNumber[5] = "0.01";
bool debugMode = false;

struct Event
{
  unsigned int  nRun;
  unsigned int  nEvent;
  unsigned int  nLevel1;
  unsigned int  nParticle;
  unsigned int  nCharge;
  unsigned int  nTrTrack;
  unsigned int  nTrdTrack;
  unsigned int  nAntiCluster;
  unsigned int  nTofClustersInTime;
  unsigned int  nRichRing;
  unsigned int  nRichRingB;
  unsigned int  nBeta;
  unsigned int  nBetaB;
  unsigned int  nBbetaH;
  unsigned int  nVertex;
  unsigned int  particleType;
  float         livetime;
  float         utcTime;
  float         orbitAltitude;
  float         orbitLatitude;
  float         orbitLongitude;
  float         orbitLatitudeM;
  float         orbitLongitudeM;
  float         velR;
  float         velTheta;
  float         velPhi;
  float         yaw;
  float         pitch;
  float         roll;
  float         gLongitude;
  float         gLatitude;
  int           gCoordCalcResult;
  float         sunPosAzimuth;
  float         sunPosElevation;
  int           sunPosCalcResult;
  unsigned int  unixTime;
  int           isInShadow;
  unsigned int  ptlCharge;
  float         ptlMomentum;
  float         ptlTheta;
  float         ptlPhi;
  float         ptlCoo[3];
  float         ptlCutOffStoermer;
  float         ptlCutOffDipole;
  float         ptlCutOffMax[2];
};
struct EcalShower
{
  float         showerEnergyD;
  float         showerEnergyE;
  float         showerEnergyCorrected;
  float         showerBDT;
  float         showerCofG[3];
  float         showerCofGDist;
  float         showerCofGdX;
  float         showerCofGdY;
};
struct TofHit
{
  int           layer;
  float         edep;
  float         charge;
};
struct Tof
{
  int           tofNClusters;
  int           tofNUsedHits;
  int           isGoodBeta;
  int           isTkTofMatch;
  float         tofReducedChisqT;
  float         tofReducedChisqC;
  std::vector<TofHit> tofHitInfo;
  float         tofCharge;
};
struct Tracker
{
  std::vector<Track>       trackFitInfo;
  std::vector<trackerHits> trackerHits;
};
