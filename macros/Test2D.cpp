TMultiDim* Build_dsdO()
{
  Char_t  TextLine[4096];
  Char_t* Token;
  FILE* File_dsdO;
  Double_t Value;

  Int_t EnergyBins = 157;
  Int_t ThetaBins = 61;

  TMultiDim* dsdO = new TMultiDim(2);
  dsdO->SetNBins(EnergyBins, ThetaBins);

  //Read 2-D distribution from matrix text file
  File_dsdO = fopen("/home/schumans/Projects/Phobos/data/MAID/PPi0_dsdO.txt", "r");
  for(Int_t th=0; th<ThetaBins; th++)
  {
    fgets(TextLine, 4096, File_dsdO);
    Token = strtok(TextLine, " ");
    sscanf(Token, "%lf", &Value);
    Value*=TMath::Sin(3.0*th*TMath::DegToRad());
    dsdO->SetValue(0, th, Value);
    for(Int_t e=1; e<EnergyBins; e++)
    {
      Token = strtok(0, " ");
      sscanf(Token, "%lf", &Value);
      Value*=TMath::Sin(3.0*th*TMath::DegToRad());
      dsdO->SetValue(e, th, Value);
    }
  }
  fclose(File_dsdO);

  for(Int_t e=0; e<EnergyBins; e++)
   dsdO->SetCoord(0, e, 10.0*e);
  for(Int_t th=0; th<ThetaBins; th++)
    dsdO->SetCoord(1, th, 3.0*th);

  TH2F* TwoD = new TH2F("TwoD", "TwoD", 1570, 0.0, 1570.0, 180, 0.0, 180.0);
  for(Int_t e=0; e<1570; e++)
    for(Int_t th=0; th<180; th++)
      TwoD->Fill(e, th, dsdO->Interpolate(1.0*e, 1.0*th));
  TwoD->Draw("col");

  return dsdO;
}

//------------------------------------------------------------------------------

TMultiDim* BuildFiveFold(Int_t Obs)
{
  Int_t EtagBins = 7;
  Int_t ThetagBins = 7;
  Int_t Thetapi0Bins = 7;
  Int_t Phipi0Bins = 7;
  Int_t EgBins = 40;

  Int_t EtagBinWidth = 25;
  Int_t EtagBinOffset = 325;
  Int_t ThetagBinWidth = 30;
  Int_t Thetapi0BinWidth = 30;
  Int_t Phipi0BinWidth = 60;

  Char_t egaXXX[16];
  Char_t thgaXXX[16];
  Char_t thpiXXX[16];
  Char_t phiXXX[16];

  Char_t DotParName[4096];
  FILE* DotParFile;
  Int_t ScanRes;

  Double_t Eg;
  Double_t DontKnow;
  Double_t CrossSect;

  Char_t PPi0G_FiveFoldPath[1024];
  sprintf(PPi0G_FiveFoldPath, "/data/disk0/a2daphne/home/schumans/3v17/ppi0g.3p0");
  Char_t PPi0G_ParamSet[16];
  sprintf(PPi0G_ParamSet, "3p0");

  TMultiDim* dsdO = new TMultiDim(5);
  dsdO->SetNBins(EtagBins, ThetagBins, Thetapi0Bins, Phipi0Bins, EgBins);

  //Read .par files
  printf("Starting to read .par files from\n %s\n", PPi0G_FiveFoldPath);
  for(Int_t a=0; a<EtagBins; a++)
    for(Int_t b=0; b<ThetagBins; b++)
      for(Int_t c=0; c<Thetapi0Bins; c++)
        for(Int_t d=0; d<Phipi0Bins; d++)
        {
          sprintf(egaXXX, "ega%d", a*EtagBinWidth + EtagBinOffset);
          sprintf(thgaXXX, "thga%d", b*ThetagBinWidth);
          sprintf(thpiXXX, "thpi%d", c*Thetapi0BinWidth);
          sprintf(phiXXX, "phi%d", d*Phipi0BinWidth);
          sprintf(DotParName, "%s/%d/gap_gadelp_5f_%s_%s_%s_%s_delelkappa%s_pole_omega_born.par",
                  PPi0G_FiveFoldPath, a*EtagBinWidth + EtagBinOffset, egaXXX, thgaXXX, thpiXXX, phiXXX, PPi0G_ParamSet);
          DotParFile = fopen(DotParName, "r");
          if(!DotParFile) printf("Could not open \n%s\n", DotParName);
          for(Int_t e=0; e<EgBins; e++)
          {
            ScanRes = fscanf(DotParFile, "%lf %lf %lf", &Eg, &DontKnow, &CrossSect);
            if(ScanRes!=3) CrossSect = 0.0; //'nan' in .par file read
            dsdO->SetValue(a, b, c, d, e, CrossSect);
          }
          fclose(DotParFile);
        }
  printf("Finished to read .par files from\n %s\n", PPi0G_FiveFoldPath);

  //Values for all bins in all dimensions
  for(Int_t t=0; t<7; t++) //E_tag
    dsdO->SetCoord(0, t, t*25.0 + 325.0);
  for(Int_t t=0; t<7; t++) //Theta_g'
    dsdO->SetCoord(1, t, t*30.0);
  for(Int_t t=0; t<7; t++) //Theta_pi0
    dsdO->SetCoord(2, t, t*30.0);
  for(Int_t t=0; t<7; t++) //Phi_pi0
    dsdO->SetCoord(3, t, t*60.0);
  for(Int_t t=0; t<40; t++) //E_gamma'
    dsdO->SetCoord(4, t, t*5.0 + 20.0);

  if(Obs==0)
  {
    TH1F* OneD1 = new TH1F("OneD1", "OneD1", 200, 0.0, 200.0);
    for(Int_t n=0; n<10000; n++)
    {
      Double_t th_g   = gRandom->Rndm()*180.0;
      Double_t th_pi  = gRandom->Rndm()*180.0;
      Double_t phi_pi = gRandom->Rndm()*360.0;
      for(Int_t e=0; e<200; e++)
        OneD1->Fill(1.0*e, dsdO->Interpolate(400.0, th_g, th_pi, phi_pi, 1.0*e));
    }
    OneD1->Draw();
  }

  if(Obs==1)
  {
    TH1F* OneD2 = new TH1F("OneD2", "OneD2", 360, 0.0, 360.0);
    for(Int_t n=0; n<10000; n++)
    {
      Double_t th_g   = gRandom->Rndm()*180.0;
      Double_t th_pi  = gRandom->Rndm()*180.0;
      Double_t e_g    = gRandom->Rndm()*200.0;
      for(Int_t ph=0; ph<360; ph++)
        OneD2->Fill(1.0*ph, dsdO->Interpolate(400.0, th_g, th_pi, 1.0*ph, e_g));
    OneD2->Draw();
    }
  }

  if(Obs==2)
  {
    TH1F* OneD2 = new TH1F("OneD3", "OneD3", 180, 0.0, 180.0);
    for(Int_t n=0; n<10000; n++)
    {
      Double_t th_pi  = gRandom->Rndm()*180.0;
      Double_t ph_pi  = gRandom->Rndm()*180.0;
      Double_t e_g    = gRandom->Rndm()*200.0;
      for(Int_t th=0; th<360; th++)
        OneD3->Fill(1.0*th, dsdO->Interpolate(400.0, 1.0*th, th_pi, ph_pi, e_g));
    OneD3->Draw();
    }
  }

  if(Obs==3)
  {
    TH1F* OneD2 = new TH1F("OneD4", "OneD4", 200, 300.0, 500.0);
    for(Int_t n=0; n<10000; n++)
    {
      Double_t th_g   = gRandom->Rndm()*180.0;
      Double_t th_pi  = gRandom->Rndm()*180.0;
      Double_t ph_pi  = gRandom->Rndm()*180.0;
      Double_t e_g    = gRandom->Rndm()*200.0;
      for(Int_t b=325; b<475; b++)
        OneD4->Fill(1.0*b, dsdO->Interpolate(1.0*b, th_g, th_pi, ph_pi, e_g));
    OneD4->Draw();
    }
  }

  return dsdO;
}

