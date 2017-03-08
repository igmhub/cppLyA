bool append2d(TString filename, TH2D*& h_sum_wx, TH2D*& h_sum_w, TH2D*& h_sum_w2var) {
  
  TFile *file = new TFile(filename);
  if(!file || !file->IsOpen()) {
    cout << "error, file " << filename << " does not exist" << endl;
    return false;
  }
  
  TH2D* h_sum_wdd_2d = (TH2D*) file->Get("h_sum_wdd_2d");
  TH2D* h_sum_w_2d = (TH2D*) file->Get("h_sum_w_2d");
  TH2D* h_sum_w2var_2d = (TH2D*) file->Get("h_sum_w2var_2d");
  
  if(h_sum_w_2d==0) {
    cout << "error, file " << filename << " has NO w histogram" << endl;
    return false;
  }
  
  if(h_sum_w_2d->GetEntries()==0) {
    cout << "error, file " << filename << " has empty w histogram" << endl;
    return false;
  }
  if(h_sum_w2var_2d->GetEntries()==0) {
    cout << "warning, file " << filename << " has empty w2var histogram" << endl;
  }
  


  if(h_sum_wx==0) {
    h_sum_wx=h_sum_wdd_2d ;
    h_sum_w=h_sum_w_2d;
    h_sum_w2var=h_sum_w2var_2d;
    return true;
  }
  h_sum_wx->Add(h_sum_wdd_2d);
  h_sum_w->Add(h_sum_w_2d);
  h_sum_w2var->Add(h_sum_w2var_2d);
  file->Close();
  return true;
}
bool append2d_plates(TString filename, TH2D*& h_sum_wx, TH2D*& h_sum_w) {
  
  TFile *file = new TFile(filename);
  if(!file || !file->IsOpen()) {
    cout << "error, file " << filename << " does not exist" << endl;
    return false;
  }
  TH2D* h_sum_wdd_2d = (TH2D*) file->Get("h_sum_wdd_2d_plates");
  if(h_sum_wdd_2d==0) return false;
  TH2D* h_sum_w_2d = (TH2D*) file->Get("h_sum_w_2d_plates");
  
  if(h_sum_wx==0) {
    h_sum_wx=h_sum_wdd_2d ;
    h_sum_w=h_sum_w_2d;
    return true;
  }
  if(h_sum_wx->GetXaxis()->GetNbins() != h_sum_wdd_2d->GetXaxis()->GetNbins()) {
    cout << "problem here incompatible X axis with file " << filename << endl;
    return false;
  }
  if(h_sum_wx->GetYaxis()->GetNbins() != h_sum_wdd_2d->GetYaxis()->GetNbins()) {
    cout << "problem here incompatible Y axis with file " << filename << endl;
    return false;
  }
  h_sum_wx->Add(h_sum_wdd_2d);
  h_sum_w->Add(h_sum_w_2d);
  file->Close();
  return true;
}
bool append1d(TString filename, TH1D*& h_sum_wx, TH1D*& h_sum_w, TH1D*& h_sum_w2var, TString what) {
  
  TFile *file = new TFile(filename);
  if(!file || !file->IsOpen()) {
    cout << "error, file " << filename << " does not exist" << endl;
    return false;
  }
  TString toto;
  toto="h_sum_wdd_"; toto+=what;
  TH1D* h_sum_wdd_tmp = (TH1D*) file->Get(toto);
  toto="h_sum_w_"; toto+=what;
  TH1D* h_sum_w_tmp = (TH1D*) file->Get(toto);
  toto="h_sum_w2var_"; toto+=what;
  TH1D* h_sum_w2var_tmp = (TH1D*) file->Get(toto);
  
  if(h_sum_w_tmp==0) {
    cout << "error, file " << filename << " has NO w histogram" << endl;
    return false;
  }
  
  if(h_sum_w_tmp->GetEntries()==0) {
    cout << "error, file " << filename << " has empty w histogram" << endl;
    return false;
  }
  if(h_sum_w2var_tmp->GetEntries()==0) {
    cout << "warning, file " << filename << " has empty w2var histogram" << endl;
  }
  


  if(h_sum_wx==0) {
    h_sum_wx=h_sum_wdd_tmp ;
    h_sum_w=h_sum_w_tmp;
    h_sum_w2var=h_sum_w2var_tmp;
    return true;
  }
  h_sum_wx->Add(h_sum_wdd_tmp);
  h_sum_w->Add(h_sum_w_tmp);
  h_sum_w2var->Add(h_sum_w2var_tmp);
  file->Close();
  return true;
}

bool append_nspec_per_plates(TString filename,TH1D*& hsum1,TH1D*& hsum2) {
  TFile *file = new TFile(filename);
  if(!file || !file->IsOpen()) {
    cout << "error, file " << filename << " does not exist" << endl;
    return false;
  }
  TH1D* hsum1_tmp = (TH1D*) file->Get("h_nspec1_per_plate");
  TH1D* hsum2_tmp = (TH1D*) file->Get("h_nspec2_per_plate");
  
  if(hsum1_tmp==0) {
    cout << "error, file " << filename << " has NO h_nspec1_per_plate" << endl;
    return false;
  }
  
  
  
  


  if(hsum1==0) {
    hsum1=hsum1_tmp ;
    hsum2=hsum2_tmp ;
    //file->Close();
    return true;
  }
  hsum1->Add(hsum1_tmp);
  hsum2->Add(hsum2_tmp);
  file->Close();
  return true;
}

void merge_res(TString ofilename="res-merged.root",TString dirname=".") {

  gSystem->Unlink(ofilename.Data());

  /* find a way to get list of files */
  vector<TString> filenames;
  {
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (! files) {
      cout << "empty directory" << endl;
      return;
    }
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
	filenames.push_back(fname);
      }
    }
    if(1)
      for(int i=0;i<filenames.size();i++) {
	cout << i << " " << filenames[i] << endl;
      }
  }

  
  
  TFile *nfile=new TFile(ofilename.Data(),"recreate");
  
  for(int wedge=0; wedge<3; wedge++) {
    TH1D* h_sum_wx=0;
    TH1D* h_sum_w=0;
    TH1D* h_sum_w2var=0;
    
    TString what;
    if(wedge==0) what="para";
    if(wedge==1) what="mid";
    if(wedge==2) what="perp";

    for(int i=0;i<filenames.size();i++) {
      append1d(filenames[i],h_sum_wx,h_sum_w,h_sum_w2var,what);
    }
    nfile->cd();
    h_sum_wx->Write();
    h_sum_w->Write();
    h_sum_w2var->Write();
  }

  {
    TH2D* h2_sum_wx=0;
    TH2D* h2_sum_w=0;
    TH2D* h2_sum_w2var=0;
    for(int i=0;i<filenames.size();i++) {
      append2d(filenames[i],h2_sum_wx,h2_sum_w,h2_sum_w2var);
      if(h2_sum_wx==0) {cout << "error for append2d" << endl; break;}
    }
    nfile->cd();
    if(h2_sum_wx) {
      h2_sum_wx->Write();
      h2_sum_w->Write();
      h2_sum_w2var->Write();
    }
  }
  
  
   {
    TH2D* h2_sum_wx=0;
    TH2D* h2_sum_w=0;
    for(int i=0;i<filenames.size();i++) {
      append2d_plates(filenames[i],h2_sum_wx,h2_sum_w);
      if(h2_sum_wx==0) {cout << "error for append2d_plates" << endl; break;}
    }
    nfile->cd();
    if(h2_sum_wx) {
      h2_sum_wx->Write();
      h2_sum_w->Write();
    }
   }
   
   {
    TH1D* hsum1=0;
    TH1D* hsum2=0;
    // only the first one
    append_nspec_per_plates(filenames[0],hsum1,hsum2);
    nfile->cd();
    if(hsum1) {
      hsum1->Write();
      hsum2->Write();
    }
   }
   
  nfile->Close();
  cout << "wrote " << ofilename << endl;
}
