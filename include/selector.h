bool IsBadRun(AMSEventR*);
bool IsScienceRun(AMSEventR*);
bool IsHardwareStatusGood(AMSEventR*);
bool IsUnbiasedPhysicsTriggerEvent(AMSEventR*);
bool IsACCPatternGood(AMSEventR*);
bool IsGoodBeta(AMSEventR*);
bool IsGoodLiveTime(AMSEventR*);
bool IsInSouthAtlanticAnomaly(AMSEventR*);
bool IsInSolarArrays(AMSEventR*);
bool IsGoodTrTrack(TrTrackR*);
bool IsShowerTrackMatched(EcalShowerR*, TrTrackR*);
bool IsTrackInsideEcalFiducialVolume(TrTrackR*);
bool IsTrkAlignmentGood(AMSEventR*);
int  GetGoodParticleIndex(AMSEventR*);

class Event : public TObject
{
  private:
    char fName[25];

};
