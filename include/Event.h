class Event : public TObject
{
  private:
    char fName[20];
    Int_t fNtrack;
    Float_t fRigidity;

    ClassDef(Event,1)
};
