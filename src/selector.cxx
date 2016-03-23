#ifndef __AMSINC__
#define __AMSINC__
#include "amschain.h"
#include "selector.h"
#endif

extern bool debugMode;

/**
 * @brief   This function checks whether current event contained in bad-run.
 * @return  true : The event is contained in bad-run. (Therefore, it should be rejected.) / false : The event is good to analyze further.
 */
bool IsBadRun(AMSEventR* thisEvent)
{/*{{{*/
  unsigned int runNumber = (unsigned int)thisEvent->Run();
  if( thisEvent->fHeader.Run == 1306219312 || thisEvent->fHeader.Run == 1306219522 || thisEvent->fHeader.Run == 1306233745 )
  {
    fprintf(stderr, "The program is aborted due to input run: %d is a bad run(1).\n", runNumber);
    exit(0);
  }
  else if( thisEvent->fHeader.Run >= 1307125541 && thisEvent->fHeader.Run <= 1307218054 )
  {
    fprintf(stderr, "The program is aborted due to input run: %d is a bad run(2).\n", runNumber);
    exit(0);
  }
  else if( thisEvent->fHeader.Run == 132119816 )
  {
    fprintf(stderr, "The program is aborted due to input run: %d is a bad run(3).\n", runNumber);
    exit(0);
  }
  else if( thisEvent->isBadRun(thisEvent->fHeader.Run) )
  {
    fprintf(stderr, "The program is aborted due to input run: %d is a bad run(4).\n", runNumber);
    exit(0);
  }
  else
    return false;
}/*}}}*/



/**
 * @brief This function checks whether the run containing input event is science-run or not.
 * @return true : The event is a science-run / false : The run containing input event is not a science-run
 */
bool IsScienceRun(AMSEventR* thisEvent)
{/*{{{*/
  HeaderR* header = &(thisEvent->fHeader);
  if( !header ) return false;
  if( ( header->RunType >> 12 ) != 0xf ) return false;
  return true;
}/*}}}*/



/**
 * @brief This function returns hardware status.
 * @return true : Hardware status is good / false : Hardware status is not good
 */
bool IsHardwareStatusGood(AMSEventR* thisEvent)
{/*{{{*/
  int nDaqEvent = thisEvent->nDaqEvent();
  bool goodHW = true;
  bool error = false;

  for(int i = 0; i < nDaqEvent; i++)
  {
    DaqEventR* daqEvent = thisEvent->pDaqEvent(i);
    for( int iJINJ = 0; iJINJ < 4; iJINJ++ ) error |= (bool)( ( daqEvent->JINJStatus[iJINJ]>>8 ) & 0x7F );
    for( int iJErr = 0; iJErr < 24;iJErr++ ) error |= (bool)(   daqEvent->JError[iJErr] & 0x7F );
    if(error) goodHW &= false;
  }

  return goodHW;
}/*}}}*/



/**
 * @brief This function checks whether the event triggered by unbiased physics trigger or not.
 * @return true : The trigger is unbiased physics trigger / false : The trigger is physics trigger
 */
bool IsUnbiasedPhysicsTriggerEvent(AMSEventR* thisEvent)
{/*{{{*/
  int nLevel1 = thisEvent->nLevel1();
  bool unbiased = false;

  for(int i = 0; i < nLevel1; i++)
  {
    Level1R* pLevel1 = thisEvent->pLevel1(i);
    bitset<8> physicsBitPattern(pLevel1->PhysBPatt);

    if(physicsBitPattern.test(1) || physicsBitPattern.test(2) || physicsBitPattern.test(3) || physicsBitPattern.test(4) || physicsBitPattern.test(5))
      unbiased |= false;
    else
      unbiased |= true;
  }

  return unbiased;
}/*}}}*/


/**
 * @brief This function checks whether number of ACC hit exceed 5 or not.
 * @return true : Number of ACC hits smaller than or equal to 5. false : Number of ACC hits greater than 5.
 */
bool IsACCPatternGood(AMSEventR* thisEvent)
{/*{{{*/
  int nACCHit = 0;
  int nLevel1 = thisEvent->nLevel1();

  for(int i = 0; i < nLevel1; i++)
  {
    Level1R* pLevel1 = thisEvent->pLevel1(i);
    for( int j = 0; j < 8; j++)
    {
      if(((pLevel1->AntiPatt>>i)&1)==1) nACCHit++;
    }
  }
  if( nACCHit > 5 ) return false;
  return true;
}/*}}}*/



/**
 * @brief  This function checks two aspects from TOF measured information.
 *         First is hit configuration. BetaR::Pattern variable contains information about the configuration of TOF hit pattern.
 *         If BetaR::Pattern is greater than 5, it means there are at least three hits on TOF layers out of 4 layers.
 *         Second is it's velocity. We restricted our interest of particle velocity to the range of 0.3 < beta <~ 1.
 *         So we discarded slow particles to ensure the accuracy of measurement.
 * @return true : The BetaR object has good pattern and the speed of the particle is enough to be analyzed.
 *         false : The BetaR object doesn't have good hit pattern or the speed of the particle is not enough to be analyzed.
 */
bool IsGoodBeta(BetaR* thisBeta)
{/*{{{*/
  if( !thisBeta ) return false;

  if( thisBeta->Pattern > 5 )
    return false;
  if( thisBeta->Beta < 0.3 )
    return false;

  return true;
}/*}}}*/



/**
 * @brief
 * @return
 */
bool IsGoodLiveTime(AMSEventR* thisEvent)
{/*{{{*/
  return true;
}/*}}}*/



/**
 * @brief
 * @return
 */
bool IsInSouthAtlanticAnomaly(AMSEventR* thisEvent)
{/*{{{*/
  return false;
}/*}}}*/



/**
 * @brief
 * @return
 */
bool IsInSolarArrays(AMSEventR* thisEvent)
{/*{{{*/
  return false;
}/*}}}*/



/**
 * @brief
 * @return
 */
bool IsGoodTrTrack(TrTrackR* thisTrack)
{/*{{{*/
  bool debugMode = false;

  if( !thisTrack )
  {
    if( debugMode ) cerr << "No track pointer. Reject it." << endl;
    return false;
  }
  if( thisTrack->IsFake() )
  {
    if( debugMode ) cerr << "Found a fake track. Reject it." << endl;
    return false;
  }

  int id_fullspan;
  id_fullspan = thisTrack->iTrTrackPar(1, 7, 0);
  int id_maxspan;
  id_maxspan = thisTrack->iTrTrackPar(1, 0, 20);
  //int id_inner;

  /*
  if( id_fullspan < 0 )
  {
    if( debugMode == true )
    {
      switch(id_fullspan)
      {
        case -1: cerr << "[FS]The requested fit cannot be performed on this track." << endl; break;
        case -2: cerr << "[FS]The requested fit it is not available without refitting." << endl; break;
        case -3: cerr << "[FS]The refit failed." << endl; break;
        case -4: cerr << "[FS]Should not happen!! Contact the developers!." << endl; break;
        case -5: cerr << "[FS]The refit failed because there was a problem retrieving dynamic alignment for Ext Layers" << endl; break;
        default: break;
      }
    }

    return false;
  }
  */

  if( id_maxspan < 0 )
  {
    if( debugMode == true )
    {
      switch( id_maxspan )
      {
        case -1: cerr << "[MS]The requested fit cannot be performed on this track." << endl; break;
        case -2: cerr << "[MS]The requested fit it is not available without refitting." << endl; break;
        case -3: cerr << "[MS]The refit failed." << endl; break;
        case -4: cerr << "[MS]Should not happen!! Contact the developers!." << endl; break;
        case -5: cerr << "[MS]The refit failed because there was a problem retrieving dynamic alignment for Ext Layers" << endl; break;
        default: break;
      }
    }

    return false;
  }

  /*
  if( id_inner < 0 )
  {
    if( debugMode )
    {
      switch( id_inner )
      {
        case -1: cerr << "[IN]The requested fit cannot be performed on this track." << endl; break;
        case -2: cerr << "[IN]The requested fit it is not available without refitting." << endl; break;
        case -3: cerr << "[IN]The refit failed." << endl; break;
        case -4: cerr << "[IN]Should not happen!! Contact the developers!." << endl; break;
        case -5: cerr << "[IN]The refit failed because there was a problem retrieving dynamic alignment for Ext Layers" << endl; break;
        default: break;
      }
    }

    return false;
  }
  */

  bool hitOnLayerJ[9];
  for(int i = 0; i < 9; i++) hitOnLayerJ[i] = thisTrack->TestHitBitsJ(i, id_fullspan);

  if( !(hitOnLayerJ[2] || hitOnLayerJ[3]) ) { if( debugMode ){cerr << "No hit on layer 2/3 " << endl;} return false; }
  if( !(hitOnLayerJ[4] || hitOnLayerJ[5]) ) { if( debugMode ){cerr << "No hit on layer 4/5 " << endl;} return false; }
  if( !(hitOnLayerJ[6] || hitOnLayerJ[7]) ) { if( debugMode ){cerr << "No hit on layer 6/7 " << endl;} return false; }

  // Recject events with rigidity less than 0.5 GV
  if( thisTrack->GetRigidity(id_maxspan) < 0.5 ) return false;

  return true;
}/*}}}*/


/**
 * @brief
 * @return
 */
bool IsShowerTrackMatched(EcalShowerR* thisShower, TrTrackR* thisTrack)
{/*{{{*/
  if( !thisShower || !thisTrack ) return false;

  int      id_inner = thisTrack->iTrTrackPar(1, 3, 0);
  AMSPoint tkPoint;
  AMSDir   tkDir;

  thisTrack->Interpolate(thisShower->CofG[2], tkPoint, tkDir, id_inner);
  double dX = TMath::Abs(tkPoint[0] - thisShower->CofG[0] );
  double dY = TMath::Abs(tkPoint[1] - thisShower->CofG[1] );

  if( dX < 3.6 && dY < 7.2 ) return true;
  else return false;
}/*}}}*/


/**
 * @brief
 * @return
 */
bool IsTrackInsideEcalFiducialVolume(TrTrackR* thisTrack)
{/*{{{*/
  if( !thisTrack) return false;

  int id_fullspan = thisTrack->iTrTrackPar(1, 7, 0);
  AMSPoint entryPoint, exitPoint;
  AMSDir   entryDir,   exitDir;

  float zEntry = -142.792;
  float zExit  = -158.457;

  thisTrack->Interpolate(zEntry, entryPoint, entryDir, id_fullspan);
  thisTrack->Interpolate(zExit,  exitPoint,  exitDir,  id_fullspan);

  bool Entry_in_32_4 = (TMath::Abs( entryPoint.x() )<32.4)  && (TMath::Abs( entryPoint.y() )<32.4);
  bool Exit_in_32_4  = (TMath::Abs( exitPoint.x() ) <32.4)  && (TMath::Abs( exitPoint.y() ) <32.4);
  bool Entry_in_31_5 = (TMath::Abs( entryPoint.x() )<31.5)  && (TMath::Abs( entryPoint.y() )<31.5);
  bool Exit_in_31_5  = (TMath::Abs( exitPoint.x() ) <31.5)  && (TMath::Abs( exitPoint.y() ) <31.5);

  // Request: Shower axis in ECAL volume (Entry&Exit<32.4), at least Entry||Exit within 1 cell (0.5 Moliere radius) from the border
  bool inacc = (Exit_in_32_4 && Entry_in_32_4) && ( Exit_in_31_5 || Entry_in_31_5 );
  if( inacc ) return true;
  return false;
}/*}}}*/



/**
 * @brief
 * @return
 */
bool IsTrkAlignmentGood(AMSEventR* thisEvent)
{/*{{{*/
  AMSPoint pn1, pn9, pd1, pd9;
  thisEvent->GetRTIdL1L9(0, pn1, pd1, thisEvent->UTime(), 60);
  thisEvent->GetRTIdL1L9(1, pn9, pd9, thisEvent->UTime(), 60);
  if(pd1.y() > 35 || pd9.y() > 45)
    return false;
  else
    return true;
}/*}}}*/

/**
 * @brief
 * @return
 */
int GetGoodParticleIndex(AMSEventR* thisEvent)
{/*{{{*/
  bool debugMode = false;
  int nParticle = 0;
  int nGoodParticle = 0;
  int iGoodParticle = -1;

  ParticleR*   pParticle       = NULL;
  BetaHR*      pBetaH          = NULL;
  TrTrackR*    pTrTrack        = NULL;
  EcalShowerR* pEcalShower     = NULL;
  RichRingR*   pRichRing       = NULL;
  bool         GoodBeta        = false;
  bool         TkTofMatch      = false;
  bool         GoodTrTrack     = false;
  bool         ShowerTkMatch   = false;
  bool         FidVolumeTest   = false;
  bool         GoodRing        = false;
  bool         result          = false;

  nParticle = thisEvent->nParticle();

  if( nParticle < 1 )
    return -1;

  if( debugMode ) cerr << "Event : " << thisEvent->Event() << " / ";
  for( int i = 0; i < nParticle; ++i )
  {
    pParticle   = thisEvent->pParticle(i);
    if( debugMode ) cerr << "ParticleR : ";
    if( !pParticle   )
    {
      if( debugMode ) cerr << "0";
      return -2;
    }
    else
      if( debugMode ) cerr << "1";
    pBetaH      = pParticle->pBetaH();
    if( debugMode ) cerr << " / BetaHR : ";
    if( !pBetaH      )
    {
      if( debugMode ) cerr << "0";
      return -3;
    }
    else
      if( debugMode ) cerr << "1";
    pTrTrack    = pParticle->pTrTrack();
    if( debugMode ) cerr << " / TrTrackR : ";
    if( !pTrTrack    )
    {
      if( debugMode ) cerr << "0";
      return -4;
    }
    else
      if( debugMode ) cerr << "1";
    pEcalShower = pParticle->pEcalShower();
    if( debugMode ) cerr << " / EcalShowerR : ";
    if( !pEcalShower )
    {
      if( debugMode ) cerr << "0";
      return -5;
    }
    else
      if( debugMode ) cerr << "1";
//    pRichRing   = pParticle->pRichRing();

    // BetaHR::IsGoodBeta() will return true when track passes 4 layers of TOF measured hits.
    // BetaHR::IsTkTofMatch() will return true when Tracker track and TOF geometry are match.
    GoodBeta      = pBetaH->IsGoodBeta();
    TkTofMatch    = pBetaH->IsTkTofMatch();
    GoodTrTrack   = IsGoodTrTrack(pTrTrack);
    ShowerTkMatch = IsShowerTrackMatched(pEcalShower, pTrTrack);
    FidVolumeTest = IsTrackInsideEcalFiducialVolume(pTrTrack);

    if( debugMode ) cerr << " / Result: " << GoodBeta << TkTofMatch << GoodTrTrack << ShowerTkMatch << FidVolumeTest;
    result = GoodBeta && TkTofMatch && GoodTrTrack && ShowerTkMatch && FidVolumeTest;
    if( debugMode ) cerr << " / Final : " << result;
    nGoodParticle++;
    if( nGoodParticle > 1 )
    {
      cout << "More than two good particles" << endl;
      return -1; // return -1 if there are more than one good particle
    }
    else if( result == false ) return -1;
    else
      iGoodParticle = i;
  }
  if( debugMode ) cerr << " / iGoodParticle : " << iGoodParticle << endl;
  return iGoodParticle;
}/*}}}*/
