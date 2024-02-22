#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <errno.h>
#include <assert.h>

//#ifdef __CLING__
/////R__ADD_INCLUDE_PATH(/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/Delphes);
/////R__ADD_INCLUDE_PATH(/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/Delphes/external/);
/////R__ADD_LIBRARY_PATH(/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/Delphes/);
/////#include "classes/DelphesClasses.h"
/////#include "ExRootAnalysis/ExRootTreeReader.h"
/////R__LOAD_LIBRARY(Delphes);

#include "/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/ExRootAnalysis/ExRootAnalysis/ExRootClasses.h"
#include "/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/ExRootAnalysis/ExRootAnalysis/ExRootLHEFReader.h"
#include "/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/ExRootAnalysis/ExRootAnalysis/ExRootTreeReader.h"
R__ADD_INCLUDE_PATH(/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/ExRootAnalysis/ExRootAnalysis/)
R__ADD_LIBRARY_PATH(/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/ExRootAnalysis/)
R__LOAD_LIBRARY(/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/ExRootAnalysis/libExRootAnalysis.so)
//#endif



bool AnalyseEvents(ExRootTreeReader *treeReader, TString filename, int RWnstep=20)
{
	double a_matrix[8][8] = {0};
	for(int aindex=1; aindex<=8; aindex++){
		for(int aindex2=1; aindex2<=8; aindex2++){
			a_matrix[aindex-1][aindex2-1] = 0;
		}
	}
	double gL_SQ = (-0.2766)*(-0.2766);
	double gR_SQ = (0.2234)*(0.2234);
	a_matrix[1-1][1-1] = gR_SQ;
	a_matrix[2-1][2-1] = gR_SQ;
	a_matrix[3-1][3-1] = gR_SQ - 0.5*gL_SQ;
	a_matrix[4-1][4-1] = gR_SQ - gL_SQ;
	a_matrix[5-1][5-1] = gR_SQ - gL_SQ;
	a_matrix[6-1][6-1] = gR_SQ;
	a_matrix[7-1][7-1] = gR_SQ;
	a_matrix[8-1][8-1] = 0.5*gL_SQ - gR_SQ;
	a_matrix[1-1][6-1] = gL_SQ;
	a_matrix[2-1][7-1] = gL_SQ;
	a_matrix[3-1][8-1] = (sqrt(3)/2.0)*gL_SQ;
	a_matrix[6-1][1-1] = gL_SQ;
	a_matrix[7-1][2-1] = gL_SQ;
	a_matrix[8-1][3-1] = (sqrt(3)/2.0)*gL_SQ;
	for(int aindex=1; aindex<=8; aindex++){
		for(int aindex2=1; aindex2<=8; aindex2++){
			a_matrix[aindex-1][aindex2-1] /= (gL_SQ-gR_SQ);
		}
	}

	if( !(access(filename,0)) ) { 
		TFile file(filename, "read");
		TTree *tree = (TTree*)(file.Get("tree"));
		bool event_exist = ( (tree->GetEntries() > 0) ? true : false );
		if(event_exist) cout << filename << " contains " << tree->GetEntries() << " entries " << endl;
		delete tree;
		file.Close();
		return event_exist;
	}


	Long64_t numberOfEntries = treeReader->GetEntries();

	// create a root file to store the variables
	TFile file(filename, "recreate");
	TTree *tree = new TTree("tree", Form("tree_from_%s", filename.Data()));

	double sum_weights = 0;
	double CX = 0;
	double SF = 1.0;

	int inputNum=0;
	int event;
	string lep_fst="";
	double M4l; 
	double Cos_Scatter_Angle;
	double deltaM; 
	double deltaM1234; 
	double deltaM1432; 
	double Ptll1; 
	double Ptll2; 
	double Etall1; 
	double Etall2; 
	double Mll1; 
	double Mll2; 
	double Pt4l;
	double P4l;
	double E4l;
	double met; 
	double met_byhand; 
	double meta; 
	double Mjj; 
	double Ptjj; 
	double deltaRjj; 
	double Mrecoil; 
	double ptMu_ptlead; 
	double ptMu_pttail; 
	double etaMu_ptlead; 
	double etaMu_pttail; 
	double phiMu_ptlead; 
	double phiMu_pttail;
	double ptEl_ptlead; 
	double ptEl_pttail; 
	double etaEl_ptlead; 
	double etaEl_pttail; 
	double phiEl_ptlead; 
	double phiEl_pttail;
	double ptLep_ptlead; 
	double ptLep_pttail; 
	double etaLep_ptlead; 
	double etaLep_pttail; 
	double phiLep_ptlead; 
	double phiLep_pttail;
	double ptJ_ptlead; 
	double ptJ_pttail; 
	double etaJ_ptlead; 
	double etaJ_pttail; 
	double phiJ_ptlead; 
	double phiJ_pttail;
	double MassJ_ptlead;
	double MassJ_pttail;
	double deltaR12; 
	double deltaR34;
	double deltaR13;
	double deltaR24; 
	double deltaR_ll1;
	double deltaR_ll2; 
	int numMu; 
	int numLep; 
	int numGen; 
	int numEl;
	int numJ;
	std::vector<double>* weights = 0;
	double event_weight = 0;
	double  ptLepp_Boson1; 
	double  ptLepm_Boson1; 
	double etaLepp_Boson1; 
	double etaLepm_Boson1; 
	double phiLepp_Boson1; 
	double phiLepm_Boson1;
	double massLepp_Boson1; 
	double massLepm_Boson1;
	double  ptLepp_Boson2; 
	double  ptLepm_Boson2; 
	double etaLepp_Boson2; 
	double etaLepm_Boson2; 
	double phiLepp_Boson2; 
	double phiLepm_Boson2;
	double massLepp_Boson2; 
	double massLepm_Boson2;
	double THETA_p;
	double THETA_m;
	double PHI_p;
	double PHI_m;
	double P1_plus_Boson1;
	double P2_plus_Boson1;
	double P3_plus_Boson1;
	double P4_plus_Boson1;
	double P5_plus_Boson1;
	double P6_plus_Boson1;
	double P7_plus_Boson1;
	double P8_plus_Boson1;
	double P1_plus_Boson2;
	double P2_plus_Boson2;
	double P3_plus_Boson2;
	double P4_plus_Boson2;
	double P5_plus_Boson2;
	double P6_plus_Boson2;
	double P7_plus_Boson2;
	double P8_plus_Boson2;
	double P1_tilde_Boson1;
	double P2_tilde_Boson1;
	double P3_tilde_Boson1;
	double P4_tilde_Boson1;
	double P5_tilde_Boson1;
	double P6_tilde_Boson1;
	double P7_tilde_Boson1;
	double P8_tilde_Boson1;
	double P1_tilde_Boson2;
	double P2_tilde_Boson2;
	double P3_tilde_Boson2;
	double P4_tilde_Boson2;
	double P5_tilde_Boson2;
	double P6_tilde_Boson2;
	double P7_tilde_Boson2;
	double P8_tilde_Boson2;


	//std::vector<double> ptMu;
	//std::vector<double> etaMu;
	//std::vector<double> phiMu;
	//std::vector<double> chargeMu;
	//std::vector<double> px;
	//std::vector<double> py;
	//std::vector<double> pz;
	//std::vector<double> E;
	//std::vector<double> ptEl;
	//std::vector<double> etaEl;
	//std::vector<double> phiEl;
	//std::vector<double> chargeEl;

	tree->Branch("event", &event);
	tree->Branch("event_weight", &event_weight);
	tree->Branch("weights", &weights);
	tree->Branch("lep_fst", &lep_fst);
	tree->Branch("SF", &SF);
	tree->Branch("numMu", &numMu);
	tree->Branch("numEl", &numEl);
	tree->Branch("numJ", &numJ);
	tree->Branch("numLep", &numLep);
	tree->Branch("met",&met);
	tree->Branch("met_byhand",&met_byhand);
	tree->Branch("meta",&meta);
	tree->Branch("Mjj", &Mjj);
	tree->Branch("Ptjj",&Ptjj);
	tree->Branch("deltaRjj",&deltaRjj);
	tree->Branch("M4l", &M4l);
	tree->Branch("Cos_Scatter_Angle", &Cos_Scatter_Angle);
	tree->Branch("Pt4l",&Pt4l);
	tree->Branch("P4l",&P4l);
	tree->Branch("deltaM", &deltaM);
	tree->Branch("deltaM1234", &deltaM1234);
	tree->Branch("deltaM1432", &deltaM1432);
	tree->Branch("Mll1",&Mll1);
	tree->Branch("Ptll1",&Ptll1);
	tree->Branch("Etall1",&Etall1);
	tree->Branch("Mll2",&Mll2);
	tree->Branch("Ptll2",&Ptll2);
	tree->Branch("Etall2",&Etall2);
	tree->Branch("Mrecoil",&Mrecoil);
	tree->Branch("ptMu_ptlead",&ptMu_ptlead); 
	tree->Branch("ptMu_pttail",&ptMu_pttail); 
	tree->Branch("etaMu_ptlead",&etaMu_ptlead); 
	tree->Branch("etaMu_pttail",&etaMu_pttail); 
	tree->Branch("phiMu_ptlead",&phiMu_ptlead); 
	tree->Branch("phiMu_pttail",&phiMu_pttail);
	tree->Branch("ptEl_ptlead",&ptEl_ptlead); 
	tree->Branch("ptEl_pttail",&ptEl_pttail); 
	tree->Branch("etaEl_ptlead",&etaEl_ptlead); 
	tree->Branch("etaEl_pttail",&etaEl_pttail); 
	tree->Branch("phiEl_ptlead",&phiEl_ptlead); 
	tree->Branch("phiEl_pttail",&phiEl_pttail);
	tree->Branch("ptLep_ptlead",&ptLep_ptlead); 
	tree->Branch("ptLep_pttail",&ptLep_pttail); 
	tree->Branch("etaLep_ptlead",&etaLep_ptlead); 
	tree->Branch("etaLep_pttail",&etaLep_pttail); 
	tree->Branch("phiLep_ptlead",&phiLep_ptlead); 
	tree->Branch("phiLep_pttail",&phiLep_pttail);
	tree->Branch("ptJ_ptlead",&ptJ_ptlead); 
	tree->Branch("ptJ_pttail",&ptJ_pttail); 
	tree->Branch("etaJ_ptlead",&etaJ_ptlead); 
	tree->Branch("etaJ_pttail",&etaJ_pttail); 
	tree->Branch("phiJ_ptlead",&phiJ_ptlead); 
	tree->Branch("phiJ_pttail",&phiJ_pttail);
	tree->Branch("MassJ_ptlead",&MassJ_ptlead);
	tree->Branch("MassJ_pttail",&MassJ_pttail);
	tree->Branch("deltaR12",&deltaR12); 
	tree->Branch("deltaR34",&deltaR34);
	tree->Branch("deltaR13",&deltaR13);
	tree->Branch("deltaR24",&deltaR24); 
	tree->Branch("deltaR_ll1",&deltaR_ll1);
	tree->Branch("deltaR_ll2",&deltaR_ll2); 
	tree->Branch("ptLepp_Boson1",&ptLepp_Boson1); 
	tree->Branch("ptLepm_Boson1",&ptLepm_Boson1); 
	tree->Branch("etaLepp_Boson1",&etaLepp_Boson1); 
	tree->Branch("etaLepm_Boson1",&etaLepm_Boson1); 
	tree->Branch("phiLepp_Boson1",&phiLepp_Boson1); 
	tree->Branch("phiLepm_Boson1",&phiLepm_Boson1);
	tree->Branch("massLepp_Boson1",&massLepp_Boson1); 
	tree->Branch("massLepm_Boson1",&massLepm_Boson1);
	tree->Branch("ptLepp_Boson2",&ptLepp_Boson2); 
	tree->Branch("ptLepm_Boson2",&ptLepm_Boson2); 
	tree->Branch("etaLepp_Boson2",&etaLepp_Boson2); 
	tree->Branch("etaLepm_Boson2",&etaLepm_Boson2); 
	tree->Branch("phiLepp_Boson2",&phiLepp_Boson2); 
	tree->Branch("phiLepm_Boson2",&phiLepm_Boson2);
	tree->Branch("massLepp_Boson2",&massLepp_Boson2); 
	tree->Branch("massLepm_Boson2",&massLepm_Boson2);
	tree->Branch("THETA_p",&THETA_p);
	tree->Branch("THETA_m",&THETA_m);
	tree->Branch("PHI_p",&PHI_p);
	tree->Branch("PHI_m",&PHI_m);
	tree->Branch("P1_plus_Boson1",&P1_plus_Boson1);
	tree->Branch("P2_plus_Boson1",&P2_plus_Boson1);
	tree->Branch("P3_plus_Boson1",&P3_plus_Boson1);
	tree->Branch("P4_plus_Boson1",&P4_plus_Boson1);
	tree->Branch("P5_plus_Boson1",&P5_plus_Boson1);
	tree->Branch("P6_plus_Boson1",&P6_plus_Boson1);
	tree->Branch("P7_plus_Boson1",&P7_plus_Boson1);
	tree->Branch("P8_plus_Boson1",&P8_plus_Boson1);
	tree->Branch("P1_plus_Boson2",&P1_plus_Boson2);
	tree->Branch("P2_plus_Boson2",&P2_plus_Boson2);
	tree->Branch("P3_plus_Boson2",&P3_plus_Boson2);
	tree->Branch("P4_plus_Boson2",&P4_plus_Boson2);
	tree->Branch("P5_plus_Boson2",&P5_plus_Boson2);
	tree->Branch("P6_plus_Boson2",&P6_plus_Boson2);
	tree->Branch("P7_plus_Boson2",&P7_plus_Boson2);
	tree->Branch("P8_plus_Boson2",&P8_plus_Boson2);

	tree->Branch("P1_tilde_Boson1",&P1_tilde_Boson1);
	tree->Branch("P2_tilde_Boson1",&P2_tilde_Boson1);
	tree->Branch("P3_tilde_Boson1",&P3_tilde_Boson1);
	tree->Branch("P4_tilde_Boson1",&P4_tilde_Boson1);
	tree->Branch("P5_tilde_Boson1",&P5_tilde_Boson1);
	tree->Branch("P6_tilde_Boson1",&P6_tilde_Boson1);
	tree->Branch("P7_tilde_Boson1",&P7_tilde_Boson1);
	tree->Branch("P8_tilde_Boson1",&P8_tilde_Boson1);

	tree->Branch("P1_tilde_Boson2",&P1_tilde_Boson2);
	tree->Branch("P2_tilde_Boson2",&P2_tilde_Boson2);
	tree->Branch("P3_tilde_Boson2",&P3_tilde_Boson2);
	tree->Branch("P4_tilde_Boson2",&P4_tilde_Boson2);
	tree->Branch("P5_tilde_Boson2",&P5_tilde_Boson2);
	tree->Branch("P6_tilde_Boson2",&P6_tilde_Boson2);
	tree->Branch("P7_tilde_Boson2",&P7_tilde_Boson2);
	tree->Branch("P8_tilde_Boson2",&P8_tilde_Boson2);

	SF = 1.0/((double)numberOfEntries);
	event_weight = 0;

	TClonesArray *branchEvent = treeReader->UseBranch("Event");
	//TClonesArray *branchWeight = treeReader->UseBranch("Weight");
	TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");

	inputNum = numberOfEntries;

	cout << filename << "'s input file contains " << numberOfEntries << " events" << endl;
	TLorentzVector l1;
	TLorentzVector l2;
	TLorentzVector l3;
	TLorentzVector l4;
	int mother1_l1;
	int mother1_l2;
	int mother1_l3;
	int mother1_l4;
	TLorentzVector lep_lv;
	TLorentzVector jet_lv;
	TLorentzVector Z_p;
	TLorentzVector Z_m;
	TLorentzVector Lepp_in_Z_p;
	TVector3 vBoson_Z_p;
	TVector3 boost_Z_p;
	TVector3 vLepp_Z_p;
	TLorentzVector Lepm_in_Z_m;
	TVector3 vBoson_Z_m;
	TVector3 boost_Z_m;
	TVector3 vLepm_Z_m;
	TVector3 v3_zaxis_lab(0,0,1);

	double Mz = 91.0 ;
	double MuMass=1.056600e-01; 
	double ElMass=5.110000e-04;
	int mup_ind=0;
	int mum_ind=0;
	int elp_ind=0;
	int elm_ind=0;

	for(int count=0; count < numberOfEntries; count++) {
		cout<<"Processing event #"<<count<<endl;
		event=-1;
		event_weight=1;
		if(weights && !(weights->empty())) weights->clear();
		lep_fst="";
		met=-1000.0;
		meta=-1000.0;
		Mjj=-1000.0;
		Ptjj=-1000.0;
		deltaRjj=-1000.0;
		numMu=0;
		numEl=0;
		numJ=0;

		ptMu_ptlead=-1000.0; 
		ptMu_pttail=-1000.0; 
		etaMu_ptlead=-1000.0; 
		etaMu_pttail=-1000.0; 
		phiMu_ptlead=-1000.0; 
		phiMu_pttail=-1000.0;
		ptEl_ptlead=-1000.0; 
		ptEl_pttail=-1000.0; 
		etaEl_ptlead=-1000.0; 
		etaEl_pttail=-1000.0; 
		phiEl_ptlead=-1000.0; 
		phiEl_pttail=-1000.0;
		ptJ_ptlead=-1000.0; 
		ptJ_pttail=-1000.0; 
		etaJ_ptlead=-1000.0; 
		etaJ_pttail=-1000.0; 
		phiJ_ptlead=-1000.0; 
		phiJ_pttail=-1000.0;
		MassJ_ptlead=-1000.0; 
		MassJ_pttail=-1000.0; 
		deltaR12=-1000.0; 
		deltaR34=-1000.0;
		deltaR13=-1000.0;
		deltaR24=-1000.0; 
		deltaR_ll1=-1000.0;
		deltaR_ll2=-1000.0; 

		numGen=0;
		numLep=0;
		M4l=-1000.0;
		Cos_Scatter_Angle=0;
		Pt4l=-1000.0;
		P4l=-1000.0;
		deltaM=-1000.0;
		deltaM1234=-1000.0;
		deltaM1432=-1000.0;
		Ptll1=-1000.0;
		Ptll2=-1000.0;
		Etall1=-1000.0;
		Etall2=-1000.0;
		Mll1=-1000.0;
		Mll2=-1000.0;
		Mrecoil=-1000.0;
		E4l=-1000.0; 

		ptLepp_Boson1=0; 
		ptLepm_Boson1=0; 
		etaLepp_Boson1=0; 
		etaLepm_Boson1=0; 
		phiLepp_Boson1=0; 
		phiLepm_Boson1=0;
		massLepp_Boson1=0; 
		massLepm_Boson1=0;
		ptLepp_Boson2=0; 
		ptLepm_Boson2=0; 
		etaLepp_Boson2=0; 
		etaLepm_Boson2=0; 
		phiLepp_Boson2=0; 
		phiLepm_Boson2=0;
		massLepp_Boson2=0; 
		massLepm_Boson2=0;
		THETA_p=0;
		THETA_m=0;
		PHI_p=0;
		PHI_m=0;
		P1_plus_Boson1=0;
		P2_plus_Boson1=0;
		P3_plus_Boson1=0;
		P4_plus_Boson1=0;
		P5_plus_Boson1=0;
		P6_plus_Boson1=0;
		P7_plus_Boson1=0;
		P8_plus_Boson1=0;

		P1_plus_Boson2=0;
		P2_plus_Boson2=0;
		P3_plus_Boson2=0;
		P4_plus_Boson2=0;
		P5_plus_Boson2=0;
		P6_plus_Boson2=0;
		P7_plus_Boson2=0;
		P8_plus_Boson2=0;

		P1_tilde_Boson1=0;
		P2_tilde_Boson1=0;
		P3_tilde_Boson1=0;
		P4_tilde_Boson1=0;
		P5_tilde_Boson1=0;
		P6_tilde_Boson1=0;
		P7_tilde_Boson1=0;
		P8_tilde_Boson1=0;

		P1_tilde_Boson2=0;
		P2_tilde_Boson2=0;
		P3_tilde_Boson2=0;
		P4_tilde_Boson2=0;
		P5_tilde_Boson2=0;
		P6_tilde_Boson2=0;
		P7_tilde_Boson2=0;
		P8_tilde_Boson2=0;

		//if(!ptMu.empty())ptMu.clear(); 
		//if(!etaMu.empty())etaMu.clear(); 
		//if(!phiMu.empty())phiMu.clear(); 
		//if(!chargeMu.empty())chargeMu.clear(); 
		//if(!ptEl.empty())ptEl.clear(); 
		//if(!etaEl.empty())etaEl.clear(); 
		//if(!phiEl.empty())phiEl.clear(); 
		//if(!chargeEl.empty())chargeEl.clear(); 
		//if(!px.empty())px.clear(); 
		//if(!py.empty())py.clear(); 
		//if(!pz.empty())pz.clear(); 
		//if(!E.empty())E.clear();           

		mother1_l1 = -1000;
		mother1_l2 = -1000;
		mother1_l3 = -1000;
		mother1_l4 = -1000;

		//********************************************************************
		treeReader->ReadEntry(count) ;

		TRootLHEFEvent* ev = (TRootLHEFEvent*)branchEvent->At(0);
		event=ev->Number;
		TRootWeight* ew = 0;

		// !!!!!!!!!!!! CAUTION !!!!!!!!!!!!
		// The last 3 weights in the weights vector in branchWeight of a reweighted sample are merged weights?
		// The last third one is the merged weight of weights series named by paramerters setting
		// The last second is the merged weight of weights series named by rwgt name
		// The last first is the merged weight of original weight from default initial paramerters setting
		// !!!!!!!!!!!! CAUTION !!!!!!!!!!!!
		//if( ((filename).Contains("reweightscan")) ) {
		//	for(int k=(RWnstep); k<=((2*RWnstep)); k++) 
		//	{
		//		ew = (TRootWeight*)branchWeight->At(k);
		//		weights->push_back(ew->Weight);
		//	}
		//	ew = (TRootWeight*)branchWeight->At(42);
		//	event_weight=ew->Weight;
		//}

		if( !((filename).Contains("reweightscan")) ) {
			weights->push_back(ev->Weight);
			event_weight=ev->Weight;
		}


		TRootLHEFParticle* particle = 0;

		numGen = branchGenParticle->GetEntries();
//---cout<<"numGen = "<<numGen<<endl;

		numMu=0; numEl=0;
		for(int pn=0; pn<numGen; pn++){
			particle = (TRootLHEFParticle*) branchGenParticle->At(pn);
			if(abs(particle->PID)==11 && particle->Status==1) numEl+=1;
			if(abs(particle->PID)==13 && particle->Status==1) numMu+=1;
		}
		numLep = numMu + numEl;

		if( !( ((numMu==2)&&(numEl==2)) || ((numMu==4)&&(numEl==0)) || ((numMu==0)&&(numEl==4)) ) )
		//if( !( ((numMu==2)&&(numEl==2)) ) ) 
		{ 
			continue; 
		}


		lep_fst = string( Form("%iMu_%iEl", numMu, numEl) );
		mup_ind=0;
		mum_ind=0;
		elp_ind=0;
		elm_ind=0;
		if((numMu==2)&&(numEl==2)){
			for(int j=0; j<numGen; j++){
				particle = (TRootLHEFParticle*) branchGenParticle->At(j);
				if(mup_ind==0 && particle->Status==1 && particle->PID == -13) { l1.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l1=particle->Mother1; mup_ind+=1; continue; }
				if(mum_ind==0 && particle->Status==1 && particle->PID == 13)  { l2.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l2=particle->Mother1; mum_ind+=1; continue; }
				if(elp_ind==0 && particle->Status==1 && particle->PID == -11) { l3.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l3=particle->Mother1; elp_ind+=1; continue; }
				if(elm_ind==0 && particle->Status==1 && particle->PID == 11)  { l4.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l4=particle->Mother1; elm_ind+=1; continue; }

			}
			if( (elp_ind-elm_ind)!=0 || (mup_ind-mum_ind)!=0 ) 
			{ 
				//cout<<"Charges do not conserve!"<<endl;
				//cout<<"count = "<<count<<endl;
				continue;
			}
//---cout<<"mother1_l1 ="<<mother1_l1<<endl;
//---cout<<"mother1_l2 ="<<mother1_l2<<endl;
//---cout<<"mother1_l3 ="<<mother1_l3<<endl;
//---cout<<"mother1_l4 ="<<mother1_l4<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l1-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l1-1))->PID<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l2-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l2-1))->PID<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l3-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l3-1))->PID<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l4-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l4-1))->PID<<endl;
			if(mother1_l1 == -1000 || mother1_l2 == -1000 || mother1_l3 == -1000 || mother1_l4 == -1000) continue;
			if(mother1_l1 != mother1_l2 || mother1_l3 != mother1_l4) continue;
			///---if(
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l1-1))->PID != 23 ) ||
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l2-1))->PID != 23 ) ||
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l3-1))->PID != 23 ) ||
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l4-1))->PID != 23 )
			///---) continue;
			deltaM = fabs((l1+l2).M()-Mz) + fabs((l3+l4).M()-Mz);
			deltaM1234=deltaM;
			deltaM1432=deltaM;
			Mll1=(l1+l2).M();
			Ptll1=(l1+l2).Pt();
			Etall1=(l1+l2).Eta();
			deltaR_ll1=l1.DeltaR(l2);
			Mll2=(l3+l4).M();
			Ptll2=(l3+l4).Pt();
			Etall2=(l3+l4).Eta();
			deltaR_ll2=l3.DeltaR(l4);
			ptMu_ptlead = l1.Pt();
			etaMu_ptlead = (l1.Eta());
			phiMu_ptlead = l1.Phi();
			ptMu_pttail = l3.Pt();
			etaMu_pttail = (l3.Eta());
			phiMu_pttail = l3.Phi();
			ptEl_ptlead = l2.Pt();
			etaEl_ptlead = (l2.Eta());
			phiEl_ptlead = l2.Phi();
			ptEl_pttail = l4.Pt();
			etaEl_pttail = (l4.Eta());
			phiEl_pttail = l4.Phi();
			//--------------
			ptLepp_Boson1 = l1.Pt();
			etaLepp_Boson1 = l1.Eta();
			phiLepp_Boson1 = l1.Phi();
			massLepp_Boson1 = MuMass;
			ptLepm_Boson1 = l2.Pt();
			etaLepm_Boson1 = l2.Eta();
			phiLepm_Boson1 = l2.Phi();
			massLepm_Boson1 = MuMass;
			ptLepp_Boson2 = l3.Pt();
			etaLepp_Boson2 = l3.Eta();
			phiLepp_Boson2 = l3.Phi();
			massLepp_Boson2 = ElMass;
			ptLepm_Boson2 = l4.Pt();
			etaLepm_Boson2 = l4.Eta();
			phiLepm_Boson2 = l4.Phi();
			massLepm_Boson2 = ElMass;
			// THETA_p
			Z_p = l1 + l2;
			TVector3 vBoson_Z_p = Z_p.Vect();
			Cos_Scatter_Angle = v3_zaxis_lab.Dot(vBoson_Z_p.Unit());
			TLorentzVector Lepp_in_Z_p = l1;
			TVector3 boost_Z_p = -(Z_p.BoostVector());
			Lepp_in_Z_p.Boost(boost_Z_p);
			TVector3 vLepp_Z_p = Lepp_in_Z_p.Vect();
			THETA_p = vBoson_Z_p.Angle(vLepp_Z_p);
			//THETA_p = Lepp_in_Z_p.Theta();
			//PHI_p = Lepp_in_Z_p.Phi();
			TVector3 v3_n_p = v3_zaxis_lab.Cross(vBoson_Z_p);
			double angle_zaxis_to_Z_p = v3_zaxis_lab.Angle(vBoson_Z_p);
			v3_n_p *= (1.0/TMath::Sin(angle_zaxis_to_Z_p));
			TVector3 v3_r_p = (1.0/TMath::Sin(angle_zaxis_to_Z_p))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_p))*vBoson_Z_p);
			TVector3 unit_vBoson_Z_p = vBoson_Z_p.Unit();
			TVector3 ProjectTovBoson_vLepp_Z_p = unit_vBoson_Z_p;
			ProjectTovBoson_vLepp_Z_p *= (unit_vBoson_Z_p.Dot(vLepp_Z_p));
			TVector3 PerpendicularComponent_vLepp_Z_p = ( vLepp_Z_p - ProjectTovBoson_vLepp_Z_p );
			PHI_p = v3_n_p.Angle(PerpendicularComponent_vLepp_Z_p);
			if( (v3_r_p.Dot(PerpendicularComponent_vLepp_Z_p) < 0) ) PHI_p *= -1.0;
			//cout<<"angle_zaxis_to_Z_p = "<<angle_zaxis_to_Z_p<<endl;
			//cout<<"v3_n_p.Mag() = "<<v3_n_p.Mag()<<endl;
			//cout<<"(v3_n_p * v3_zaxis_lab) = "<<(v3_n_p * v3_zaxis_lab)<<endl;
			//cout<<"(v3_n_p * vBoson_Z_p) = "<<(v3_n_p * vBoson_Z_p)<<endl;
			//cout<<"unit_vBoson_Z_p.Mag() = "<<unit_vBoson_Z_p.Mag()<<endl;
			//cout<<"(PerpendicularComponent_vLepp_Z_p * vBoson_Z_p) = "<<(PerpendicularComponent_vLepp_Z_p * vBoson_Z_p)<<endl;
			// THETA_m
			Z_m = l3 + l4;
			TVector3 vBoson_Z_m = Z_m.Vect();
			TLorentzVector Lepm_in_Z_m = l3;
			TVector3 boost_Z_m = -(Z_m.BoostVector());
			Lepm_in_Z_m.Boost(boost_Z_m);
			TVector3 vLepm_Z_m = Lepm_in_Z_m.Vect();
			THETA_m = vBoson_Z_m.Angle(vLepm_Z_m);
			//THETA_m = Lepm_in_Z_m.Theta();
			//PHI_m = Lepm_in_Z_m.Phi();
			TVector3 v3_n_m = v3_zaxis_lab.Cross(vBoson_Z_m);
			double angle_zaxis_to_Z_m = v3_zaxis_lab.Angle(vBoson_Z_m);
			v3_n_m *= (1.0/TMath::Sin(angle_zaxis_to_Z_m));
			TVector3 v3_r_m = (1.0/TMath::Sin(angle_zaxis_to_Z_m))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_m))*vBoson_Z_m);
			TVector3 unit_vBoson_Z_m = vBoson_Z_m.Unit();
			TVector3 ProjectTovBoson_vLepm_Z_m = unit_vBoson_Z_m;
			ProjectTovBoson_vLepm_Z_m *= (unit_vBoson_Z_m.Dot(vLepm_Z_m));
			TVector3 PerpendicularComponent_vLepm_Z_m = ( vLepm_Z_m - ProjectTovBoson_vLepm_Z_m );
			PHI_m = v3_n_m.Angle(PerpendicularComponent_vLepm_Z_m);
			if( (v3_r_m.Dot(PerpendicularComponent_vLepm_Z_m) < 0) ) PHI_m *= -1.0;
			if(l1.Pt() < l3.Pt()) {
				ptMu_ptlead = l3.Pt();
				etaMu_ptlead = (l3.Eta());
				phiMu_ptlead = l3.Phi();
				ptMu_pttail = l1.Pt();
				etaMu_pttail = (l1.Eta());
				phiMu_pttail = l1.Phi();
			}
			if(l2.Pt() < l4.Pt()) {
				ptEl_ptlead = l4.Pt();
				etaEl_ptlead = (l4.Eta());
				phiEl_ptlead = l4.Phi();
				ptEl_pttail = l2.Pt();
				etaEl_pttail = (l2.Eta());
				phiEl_pttail = l2.Phi();
			}

			if(ptEl_ptlead > ptMu_ptlead) {
				ptLep_ptlead  = ptEl_ptlead ; 
				etaLep_ptlead = etaEl_ptlead;
				phiLep_ptlead = phiEl_ptlead;
			}
			else {
				ptLep_ptlead  = ptMu_ptlead ; 
				etaLep_ptlead = etaMu_ptlead;
				phiLep_ptlead = phiMu_ptlead;
			}

			if(ptEl_pttail < ptMu_pttail) {
				ptLep_pttail  = ptEl_pttail ;
				etaLep_pttail = etaEl_pttail;
				phiLep_pttail = phiEl_pttail;
			}
			else {
				ptLep_ptlead  = ptMu_pttail ; 
				etaLep_ptlead = etaMu_pttail;
				phiLep_ptlead = phiMu_pttail;
			}

		}

		if((numMu==4)&&(numEl==0)){
			for(int j=0; j<numGen; j++){
				particle = (TRootLHEFParticle*) branchGenParticle->At(j);
				if(mup_ind==0 && particle->Status==1 && particle->PID == -13) { l1.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l1=particle->Mother1; mup_ind+=1; continue; }
				if(mum_ind==0 && particle->Status==1 && particle->PID == 13)  { l2.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l2=particle->Mother1; mum_ind+=1; continue; }
				if(mup_ind==1 && particle->Status==1 && particle->PID == -13) { l3.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l3=particle->Mother1; elp_ind+=1; continue; }
				if(mum_ind==1 && particle->Status==1 && particle->PID == 13)  { l4.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l4=particle->Mother1; elm_ind+=1; continue; }

			}
			if( (mup_ind-mum_ind)!=0 ) 
			{ 
				continue;
			}
//---cout<<"mother1_l1 ="<<mother1_l1<<endl;
//---cout<<"mother1_l2 ="<<mother1_l2<<endl;
//---cout<<"mother1_l3 ="<<mother1_l3<<endl;
//---cout<<"mother1_l4 ="<<mother1_l4<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l1-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l1-1))->PID<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l2-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l2-1))->PID<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l3-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l3-1))->PID<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l4-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l4-1))->PID<<endl;
			if(mother1_l1 == -1000 || mother1_l2 == -1000 || mother1_l3 == -1000 || mother1_l4 == -1000) continue;
			if( (mother1_l1 != mother1_l2 && mother1_l1 != mother1_l4) || (mother1_l3 != mother1_l4 && mother1_l3 != mother1_l2) ) continue;
			///---if(
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l1-1))->PID != 23 ) ||
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l2-1))->PID != 23 ) ||
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l3-1))->PID != 23 ) ||
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l4-1))->PID != 23 )
			///---) continue;
			deltaM1234 = fabs((l1+l2).M()-Mz) + fabs((l3+l4).M()-Mz);
			deltaM1432 = fabs((l1+l4).M()-Mz) + fabs((l2+l3).M()-Mz);
			//if(deltaM1234>=deltaM1432){ 
			if(mother1_l1 == mother1_l4 && mother1_l3 == mother1_l2){
				deltaM=deltaM1432;
				Mll1=(l1+l4).M();
				Ptll1=(l1+l4).Pt();
				Etall1=(l1+l4).Eta();
				deltaR_ll1=l1.DeltaR(l4);
				Mll2=(l2+l3).M();
				Ptll2=(l2+l3).Pt();
				Etall2=(l2+l3).Eta();
				deltaR_ll2=l2.DeltaR(l3);
				//--------------
				ptLepp_Boson1 = l1.Pt();
				etaLepp_Boson1 = l1.Eta();
				phiLepp_Boson1 = l1.Phi();
				massLepp_Boson1 = MuMass;
				ptLepm_Boson1 = l4.Pt();
				etaLepm_Boson1 = l4.Eta();
				phiLepm_Boson1 = l4.Phi();
				massLepm_Boson1 = MuMass;
				ptLepp_Boson2 = l3.Pt();
				etaLepp_Boson2 = l3.Eta();
				phiLepp_Boson2 = l3.Phi();
				massLepp_Boson2 = MuMass;
				ptLepm_Boson2 = l2.Pt();
				etaLepm_Boson2 = l2.Eta();
				phiLepm_Boson2 = l2.Phi();
				massLepm_Boson2 = MuMass;
				// THETA_p
				Z_p = l1 + l4;
				TVector3 vBoson_Z_p = Z_p.Vect();
				Cos_Scatter_Angle = v3_zaxis_lab.Dot(vBoson_Z_p.Unit());
				TLorentzVector Lepp_in_Z_p = l1;
				TVector3 boost_Z_p = -(Z_p.BoostVector());
				Lepp_in_Z_p.Boost(boost_Z_p);
				TVector3 vLepp_Z_p = Lepp_in_Z_p.Vect();
				THETA_p = vBoson_Z_p.Angle(vLepp_Z_p);
				//THETA_p = Lepp_in_Z_p.Theta();
				//PHI_p = Lepp_in_Z_p.Phi();
				TVector3 v3_n_p = v3_zaxis_lab.Cross(vBoson_Z_p);
				double angle_zaxis_to_Z_p = v3_zaxis_lab.Angle(vBoson_Z_p);
				v3_n_p *= (1.0/TMath::Sin(angle_zaxis_to_Z_p));
				TVector3 v3_r_p = (1.0/TMath::Sin(angle_zaxis_to_Z_p))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_p))*vBoson_Z_p);
				TVector3 unit_vBoson_Z_p = vBoson_Z_p.Unit();
				TVector3 ProjectTovBoson_vLepp_Z_p = unit_vBoson_Z_p;
				ProjectTovBoson_vLepp_Z_p *= (unit_vBoson_Z_p.Dot(vLepp_Z_p));
				TVector3 PerpendicularComponent_vLepp_Z_p = ( vLepp_Z_p - ProjectTovBoson_vLepp_Z_p );
				PHI_p = v3_n_p.Angle(PerpendicularComponent_vLepp_Z_p);
				if( (v3_r_p.Dot(PerpendicularComponent_vLepp_Z_p) < 0) ) PHI_p *= -1.0;
				// THETA_m
				Z_m = l3 + l2;
				TVector3 vBoson_Z_m = Z_m.Vect();
				TLorentzVector Lepm_in_Z_m = l3;
				TVector3 boost_Z_m = -(Z_m.BoostVector());
				Lepm_in_Z_m.Boost(boost_Z_m);
				TVector3 vLepm_Z_m = Lepm_in_Z_m.Vect();
				THETA_m = vBoson_Z_m.Angle(vLepm_Z_m);
				//THETA_m = Lepm_in_Z_m.Theta();
				//PHI_m = Lepm_in_Z_m.Phi();
				TVector3 v3_n_m = v3_zaxis_lab.Cross(vBoson_Z_m);
				double angle_zaxis_to_Z_m = v3_zaxis_lab.Angle(vBoson_Z_m);
				v3_n_m *= (1.0/TMath::Sin(angle_zaxis_to_Z_m));
				TVector3 v3_r_m = (1.0/TMath::Sin(angle_zaxis_to_Z_m))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_m))*vBoson_Z_m);
				TVector3 unit_vBoson_Z_m = vBoson_Z_m.Unit();
				TVector3 ProjectTovBoson_vLepm_Z_m = unit_vBoson_Z_m;
				ProjectTovBoson_vLepm_Z_m *= (unit_vBoson_Z_m.Dot(vLepm_Z_m));
				TVector3 PerpendicularComponent_vLepm_Z_m = ( vLepm_Z_m - ProjectTovBoson_vLepm_Z_m );
				PHI_m = v3_n_m.Angle(PerpendicularComponent_vLepm_Z_m);
				if( (v3_r_m.Dot(PerpendicularComponent_vLepm_Z_m) < 0) ) PHI_m *= -1.0;
			}
			//if(deltaM1234<deltaM1432){ 
			if(mother1_l1 == mother1_l2 && mother1_l3 == mother1_l4){
				deltaM=deltaM1234;
				Mll1=(l1+l2).M();
				Ptll1=(l1+l2).Pt();
				Etall1=(l1+l2).Eta();
				deltaR_ll1=l1.DeltaR(l2);
				Mll2=(l3+l4).M();
				Ptll2=(l3+l4).Pt();
				Etall2=(l3+l4).Eta();
				deltaR_ll2=l3.DeltaR(l4);
				//--------------
				ptLepp_Boson1 = l1.Pt();
				etaLepp_Boson1 = l1.Eta();
				phiLepp_Boson1 = l1.Phi();
				massLepp_Boson1 = MuMass;
				ptLepm_Boson1 = l2.Pt();
				etaLepm_Boson1 = l2.Eta();
				phiLepm_Boson1 = l2.Phi();
				massLepm_Boson1 = MuMass;
				ptLepp_Boson2 = l3.Pt();
				etaLepp_Boson2 = l3.Eta();
				phiLepp_Boson2 = l3.Phi();
				massLepp_Boson2 = MuMass;
				ptLepm_Boson2 = l4.Pt();
				etaLepm_Boson2 = l4.Eta();
				phiLepm_Boson2 = l4.Phi();
				massLepm_Boson2 = MuMass;
				// THETA_p
				Z_p = l1 + l2;
				TVector3 vBoson_Z_p = Z_p.Vect();
				Cos_Scatter_Angle = v3_zaxis_lab.Dot(vBoson_Z_p.Unit());
				TLorentzVector Lepp_in_Z_p = l1;
				TVector3 boost_Z_p = -(Z_p.BoostVector());
				Lepp_in_Z_p.Boost(boost_Z_p);
				TVector3 vLepp_Z_p = Lepp_in_Z_p.Vect();
				THETA_p = vBoson_Z_p.Angle(vLepp_Z_p);
				//THETA_p = Lepp_in_Z_p.Theta();
				//PHI_p = Lepp_in_Z_p.Phi();
				TVector3 v3_n_p = v3_zaxis_lab.Cross(vBoson_Z_p);
				double angle_zaxis_to_Z_p = v3_zaxis_lab.Angle(vBoson_Z_p);
				v3_n_p *= (1.0/TMath::Sin(angle_zaxis_to_Z_p));
				TVector3 v3_r_p = (1.0/TMath::Sin(angle_zaxis_to_Z_p))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_p))*vBoson_Z_p);
				TVector3 unit_vBoson_Z_p = vBoson_Z_p.Unit();
				TVector3 ProjectTovBoson_vLepp_Z_p = unit_vBoson_Z_p;
				ProjectTovBoson_vLepp_Z_p *= (unit_vBoson_Z_p.Dot(vLepp_Z_p));
				TVector3 PerpendicularComponent_vLepp_Z_p = ( vLepp_Z_p - ProjectTovBoson_vLepp_Z_p );
				PHI_p = v3_n_p.Angle(PerpendicularComponent_vLepp_Z_p);
				if( (v3_r_p.Dot(PerpendicularComponent_vLepp_Z_p) < 0) ) PHI_p *= -1.0;
				// THETA_m
				Z_m = l3 + l4;
				TVector3 vBoson_Z_m = Z_m.Vect();
				TLorentzVector Lepm_in_Z_m = l3;
				TVector3 boost_Z_m = -(Z_m.BoostVector());
				Lepm_in_Z_m.Boost(boost_Z_m);
				TVector3 vLepm_Z_m = Lepm_in_Z_m.Vect();
				THETA_m = vBoson_Z_m.Angle(vLepm_Z_m);
				//THETA_m = Lepm_in_Z_m.Theta();
				//PHI_m = Lepm_in_Z_m.Phi();
				TVector3 v3_n_m = v3_zaxis_lab.Cross(vBoson_Z_m);
				double angle_zaxis_to_Z_m = v3_zaxis_lab.Angle(vBoson_Z_m);
				v3_n_m *= (1.0/TMath::Sin(angle_zaxis_to_Z_m));
				TVector3 v3_r_m = (1.0/TMath::Sin(angle_zaxis_to_Z_m))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_m))*vBoson_Z_m);
				TVector3 unit_vBoson_Z_m = vBoson_Z_m.Unit();
				TVector3 ProjectTovBoson_vLepm_Z_m = unit_vBoson_Z_m;
				ProjectTovBoson_vLepm_Z_m *= (unit_vBoson_Z_m.Dot(vLepm_Z_m));
				TVector3 PerpendicularComponent_vLepm_Z_m = ( vLepm_Z_m - ProjectTovBoson_vLepm_Z_m );
				PHI_m = v3_n_m.Angle(PerpendicularComponent_vLepm_Z_m);
				if( (v3_r_m.Dot(PerpendicularComponent_vLepm_Z_m) < 0) ) PHI_m *= -1.0;
			}
			ptMu_ptlead = l1.Pt();
			etaMu_ptlead = (l1.Eta());
			phiMu_ptlead = l1.Phi();
			ptMu_pttail = l2.Pt();
			etaMu_pttail = (l2.Eta());
			phiMu_pttail = l2.Phi();
			if(l2.Pt() > l1.Pt()) {
				ptMu_ptlead = l2.Pt();
				etaMu_ptlead = (l2.Eta());
				phiMu_ptlead = l2.Phi();
				ptMu_pttail = l1.Pt();
				etaMu_pttail = (l1.Eta());
				phiMu_pttail = l1.Phi();
			}
			if(l3.Pt() > ptMu_ptlead) {
				ptMu_ptlead = l3.Pt();
				etaMu_ptlead = (l3.Eta());
				phiMu_ptlead = l3.Phi();
			}
			if(l3.Pt() < ptMu_pttail) {
				ptMu_pttail = l3.Pt();
				etaMu_pttail = (l3.Eta());
				phiMu_pttail = l3.Phi();
			}
			if(l4.Pt() > ptMu_ptlead) {
				ptMu_ptlead = l4.Pt();
				etaMu_ptlead = (l4.Eta());
				phiMu_ptlead = l4.Phi();
			}
			if(l4.Pt() < ptMu_pttail) {
				ptMu_pttail = l4.Pt();
				etaMu_pttail = (l4.Eta());
				phiMu_pttail = l4.Phi();
			}

			{
				ptLep_ptlead  = ptMu_ptlead ; 
				etaLep_ptlead = etaMu_ptlead;
				phiLep_ptlead = phiMu_ptlead;
				ptLep_ptlead  = ptMu_pttail ; 
				etaLep_ptlead = etaMu_pttail;
				phiLep_ptlead = phiMu_pttail;
			}
		}

		if((numMu==0)&&(numEl==4)){
			for(int j=0; j<numGen; j++){
				particle = (TRootLHEFParticle*) branchGenParticle->At(j);
				if(elp_ind==0 && particle->Status==1 && particle->PID == -11) { l1.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l1=particle->Mother1;; elp_ind+=1; continue; }
				if(elm_ind==0 && particle->Status==1 && particle->PID == 11)  { l2.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l2=particle->Mother1;; elm_ind+=1; continue; }
				if(elp_ind==1 && particle->Status==1 && particle->PID == -11) { l3.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l3=particle->Mother1;; elp_ind+=1; continue; }
				if(elm_ind==1 && particle->Status==1 && particle->PID == 11)  { l4.SetPtEtaPhiM(particle->PT,particle->Eta,particle->Phi,particle->M); mother1_l4=particle->Mother1;; elm_ind+=1; continue; }

			}
			if( (elp_ind-elm_ind)!=0 ) 
			{ 
				continue;
			}
//---cout<<"mother1_l1 ="<<mother1_l1<<endl;
//---cout<<"mother1_l2 ="<<mother1_l2<<endl;
//---cout<<"mother1_l3 ="<<mother1_l3<<endl;
//---cout<<"mother1_l4 ="<<mother1_l4<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l1-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l1-1))->PID<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l2-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l2-1))->PID<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l3-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l3-1))->PID<<endl;
//---cout<<"((TRootLHEFParticle*) branchGenParticle->At(mother1_l4-1))->PID ="<<((TRootLHEFParticle*) branchGenParticle->At(mother1_l4-1))->PID<<endl;
			if(mother1_l1 == -1000 || mother1_l2 == -1000 || mother1_l3 == -1000 || mother1_l4 == -1000) continue;
			if( (mother1_l1 != mother1_l2 && mother1_l1 != mother1_l4) || (mother1_l3 != mother1_l4 && mother1_l3 != mother1_l2) ) continue;
			///---if(
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l1-1))->PID != 23 ) ||
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l2-1))->PID != 23 ) ||
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l3-1))->PID != 23 ) ||
			///---	( ((TRootLHEFParticle*) branchGenParticle->At(mother1_l4-1))->PID != 23 )
			///---) continue;
			deltaM1234 = fabs((l1+l2).M()-Mz) + fabs((l3+l4).M()-Mz);
			deltaM1432 = fabs((l1+l4).M()-Mz) + fabs((l2+l3).M()-Mz);
			//if(deltaM1234>=deltaM1432){ 
			if(mother1_l1 == mother1_l4 && mother1_l3 == mother1_l2){
				deltaM=deltaM1432;
				Mll1=(l1+l4).M();
				Ptll1=(l1+l4).Pt();
				Etall1=(l1+l4).Eta();
				deltaR_ll1=l1.DeltaR(l4);
				Mll2=(l2+l3).M();
				Ptll2=(l2+l3).Pt();
				Etall2=(l2+l3).Eta();
				deltaR_ll2=l2.DeltaR(l3);
				//--------------
				ptLepp_Boson1 = l1.Pt();
				etaLepp_Boson1 = l1.Eta();
				phiLepp_Boson1 = l1.Phi();
				massLepp_Boson1 = ElMass;
				ptLepm_Boson1 = l4.Pt();
				etaLepm_Boson1 = l4.Eta();
				phiLepm_Boson1 = l4.Phi();
				massLepm_Boson1 = ElMass;
				ptLepp_Boson2 = l3.Pt();
				etaLepp_Boson2 = l3.Eta();
				phiLepp_Boson2 = l3.Phi();
				massLepp_Boson2 = ElMass;
				ptLepm_Boson2 = l2.Pt();
				etaLepm_Boson2 = l2.Eta();
				phiLepm_Boson2 = l2.Phi();
				massLepm_Boson2 = ElMass;
				// THETA_p
				Z_p = l1 + l4;
				TVector3 vBoson_Z_p = Z_p.Vect();
				Cos_Scatter_Angle = v3_zaxis_lab.Dot(vBoson_Z_p.Unit());
				TLorentzVector Lepp_in_Z_p = l1;
				TVector3 boost_Z_p = -(Z_p.BoostVector());
				Lepp_in_Z_p.Boost(boost_Z_p);
				TVector3 vLepp_Z_p = Lepp_in_Z_p.Vect();
				THETA_p = vBoson_Z_p.Angle(vLepp_Z_p);
				//THETA_p = Lepp_in_Z_p.Theta();
				//PHI_p = Lepp_in_Z_p.Phi();
				TVector3 v3_n_p = v3_zaxis_lab.Cross(vBoson_Z_p);
				double angle_zaxis_to_Z_p = v3_zaxis_lab.Angle(vBoson_Z_p);
				v3_n_p *= (1.0/TMath::Sin(angle_zaxis_to_Z_p));
				TVector3 v3_r_p = (1.0/TMath::Sin(angle_zaxis_to_Z_p))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_p))*vBoson_Z_p);
				TVector3 unit_vBoson_Z_p = vBoson_Z_p.Unit();
				TVector3 ProjectTovBoson_vLepp_Z_p = unit_vBoson_Z_p;
				ProjectTovBoson_vLepp_Z_p *= (unit_vBoson_Z_p.Dot(vLepp_Z_p));
				TVector3 PerpendicularComponent_vLepp_Z_p = ( vLepp_Z_p - ProjectTovBoson_vLepp_Z_p );
				PHI_p = v3_n_p.Angle(PerpendicularComponent_vLepp_Z_p);
				if( (v3_r_p.Dot(PerpendicularComponent_vLepp_Z_p) < 0) ) PHI_p *= -1.0;
				// THETA_m
				Z_m = l3 + l2;
				TVector3 vBoson_Z_m = Z_m.Vect();
				TLorentzVector Lepm_in_Z_m = l3;
				TVector3 boost_Z_m = -(Z_m.BoostVector());
				Lepm_in_Z_m.Boost(boost_Z_m);
				TVector3 vLepm_Z_m = Lepm_in_Z_m.Vect();
				THETA_m = vBoson_Z_m.Angle(vLepm_Z_m);
				//THETA_m = Lepm_in_Z_m.Theta();
				//PHI_m = Lepm_in_Z_m.Phi();
				TVector3 v3_n_m = v3_zaxis_lab.Cross(vBoson_Z_m);
				double angle_zaxis_to_Z_m = v3_zaxis_lab.Angle(vBoson_Z_m);
				v3_n_m *= (1.0/TMath::Sin(angle_zaxis_to_Z_m));
				TVector3 v3_r_m = (1.0/TMath::Sin(angle_zaxis_to_Z_m))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_m))*vBoson_Z_m);
				TVector3 unit_vBoson_Z_m = vBoson_Z_m.Unit();
				TVector3 ProjectTovBoson_vLepm_Z_m = unit_vBoson_Z_m;
				ProjectTovBoson_vLepm_Z_m *= (unit_vBoson_Z_m.Dot(vLepm_Z_m));
				TVector3 PerpendicularComponent_vLepm_Z_m = ( vLepm_Z_m - ProjectTovBoson_vLepm_Z_m );
				PHI_m = v3_n_m.Angle(PerpendicularComponent_vLepm_Z_m);
				if( (v3_r_m.Dot(PerpendicularComponent_vLepm_Z_m) < 0) ) PHI_m *= -1.0;
			}
			//if(deltaM1234<deltaM1432){ 
			if(mother1_l1 == mother1_l2 && mother1_l3 == mother1_l4){
				deltaM=deltaM1234;
				Mll1=(l1+l2).M();
				Ptll1=(l1+l2).Pt();
				Etall1=(l1+l2).Eta();
				deltaR_ll1=l1.DeltaR(l2);
				Mll2=(l3+l4).M();
				Ptll2=(l3+l4).Pt();
				Etall2=(l3+l4).Eta();
				deltaR_ll2=l3.DeltaR(l4);
				//--------------
				ptLepp_Boson1 = l1.Pt();
				etaLepp_Boson1 = l1.Eta();
				phiLepp_Boson1 = l1.Phi();
				massLepp_Boson1 = ElMass;
				ptLepm_Boson1 = l2.Pt();
				etaLepm_Boson1 = l2.Eta();
				phiLepm_Boson1 = l2.Phi();
				massLepm_Boson1 = ElMass;
				ptLepp_Boson2 = l3.Pt();
				etaLepp_Boson2 = l3.Eta();
				phiLepp_Boson2 = l3.Phi();
				massLepp_Boson2 = ElMass;
				ptLepm_Boson2 = l4.Pt();
				etaLepm_Boson2 = l4.Eta();
				phiLepm_Boson2 = l4.Phi();
				massLepm_Boson2 = ElMass;
				// THETA_p
				Z_p = l1 + l2;
				TVector3 vBoson_Z_p = Z_p.Vect();
				Cos_Scatter_Angle = v3_zaxis_lab.Dot(vBoson_Z_p.Unit());
				TLorentzVector Lepp_in_Z_p = l1;
				TVector3 boost_Z_p = -(Z_p.BoostVector());
				Lepp_in_Z_p.Boost(boost_Z_p);
				TVector3 vLepp_Z_p = Lepp_in_Z_p.Vect();
				THETA_p = vBoson_Z_p.Angle(vLepp_Z_p);
				//THETA_p = Lepp_in_Z_p.Theta();
				//PHI_p = Lepp_in_Z_p.Phi();
				TVector3 v3_n_p = v3_zaxis_lab.Cross(vBoson_Z_p);
				double angle_zaxis_to_Z_p = v3_zaxis_lab.Angle(vBoson_Z_p);
				v3_n_p *= (1.0/TMath::Sin(angle_zaxis_to_Z_p));
				TVector3 v3_r_p = (1.0/TMath::Sin(angle_zaxis_to_Z_p))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_p))*vBoson_Z_p);
				TVector3 unit_vBoson_Z_p = vBoson_Z_p.Unit();
				TVector3 ProjectTovBoson_vLepp_Z_p = unit_vBoson_Z_p;
				ProjectTovBoson_vLepp_Z_p *= (unit_vBoson_Z_p.Dot(vLepp_Z_p));
				TVector3 PerpendicularComponent_vLepp_Z_p = ( vLepp_Z_p - ProjectTovBoson_vLepp_Z_p );
				PHI_p = v3_n_p.Angle(PerpendicularComponent_vLepp_Z_p);
				if( (v3_r_p.Dot(PerpendicularComponent_vLepp_Z_p) < 0) ) PHI_p *= -1.0;
				// THETA_m
				Z_m = l3 + l4;
				TVector3 vBoson_Z_m = Z_m.Vect();
				TLorentzVector Lepm_in_Z_m = l3;
				TVector3 boost_Z_m = -(Z_m.BoostVector());
				Lepm_in_Z_m.Boost(boost_Z_m);
				TVector3 vLepm_Z_m = Lepm_in_Z_m.Vect();
				THETA_m = vBoson_Z_m.Angle(vLepm_Z_m);
				//THETA_m = Lepm_in_Z_m.Theta();
				//PHI_m = Lepm_in_Z_m.Phi();
				TVector3 v3_n_m = v3_zaxis_lab.Cross(vBoson_Z_m);
				double angle_zaxis_to_Z_m = v3_zaxis_lab.Angle(vBoson_Z_m);
				v3_n_m *= (1.0/TMath::Sin(angle_zaxis_to_Z_m));
				TVector3 v3_r_m = (1.0/TMath::Sin(angle_zaxis_to_Z_m))*(v3_zaxis_lab - (TMath::Cos(angle_zaxis_to_Z_m))*vBoson_Z_m);
				TVector3 unit_vBoson_Z_m = vBoson_Z_m.Unit();
				TVector3 ProjectTovBoson_vLepm_Z_m = unit_vBoson_Z_m;
				ProjectTovBoson_vLepm_Z_m *= (unit_vBoson_Z_m.Dot(vLepm_Z_m));
				TVector3 PerpendicularComponent_vLepm_Z_m = ( vLepm_Z_m - ProjectTovBoson_vLepm_Z_m );
				PHI_m = v3_n_m.Angle(PerpendicularComponent_vLepm_Z_m);
				if( (v3_r_m.Dot(PerpendicularComponent_vLepm_Z_m) < 0) ) PHI_m *= -1.0;
			}
			ptEl_ptlead = l1.Pt();
			etaEl_ptlead = (l1.Eta());
			phiEl_ptlead = l1.Phi();
			ptEl_pttail = l2.Pt();
			etaEl_pttail = (l2.Eta());
			phiEl_pttail = l2.Phi();
			if(l2.Pt() > l1.Pt()) {
				ptEl_ptlead = l2.Pt();
				etaEl_ptlead = (l2.Eta());
				phiEl_ptlead = l2.Phi();
				ptEl_pttail = l1.Pt();
				etaEl_pttail = (l1.Eta());
				phiEl_pttail = l1.Phi();
			}
			if(l3.Pt() > ptEl_ptlead) {
				ptEl_ptlead = l3.Pt();
				etaEl_ptlead = (l3.Eta());
				phiEl_ptlead = l3.Phi();
			}
			if(l3.Pt() < ptEl_pttail) {
				ptEl_pttail = l3.Pt();
				etaEl_pttail = (l3.Eta());
				phiEl_pttail = l3.Phi();
			}
			if(l4.Pt() > ptEl_ptlead) {
				ptEl_ptlead = l4.Pt();
				etaEl_ptlead = (l4.Eta());
				phiEl_ptlead = l4.Phi();
			}
			if(l4.Pt() < ptEl_pttail) {
				ptEl_pttail = l4.Pt();
				etaEl_pttail = (l4.Eta());
				phiEl_pttail = l4.Phi();
			}

			{
				ptLep_ptlead  = ptEl_ptlead ; 
				etaLep_ptlead = etaEl_ptlead;
				phiLep_ptlead = phiEl_ptlead;
				ptLep_ptlead  = ptEl_pttail ; 
				etaLep_ptlead = etaEl_pttail;
				phiLep_ptlead = phiEl_pttail;
			}
		}



		Pt4l=(l1+l2+l3+l4).Pt();
		P4l=(l1+l2+l3+l4).P();
		M4l=(l1+l2+l3+l4).M();
		E4l=(l1+l2+l3+l4).E();
		Mrecoil=sqrt(pow(1000.0-E4l,2)-pow(P4l,2));
		deltaR12=l1.DeltaR(l2);
		deltaR34=l3.DeltaR(l4);
		deltaR13=l1.DeltaR(l3);
		deltaR24=l2.DeltaR(l4);

		P1_plus_Boson1=(sqrt(2))*(TMath::Sin(THETA_p))*(5*TMath::Cos(THETA_p)+1)*(TMath::Cos(PHI_p));
		P2_plus_Boson1=(sqrt(2))*(TMath::Sin(THETA_p))*(5*TMath::Cos(THETA_p)+1)*(TMath::Sin(PHI_p));
		P3_plus_Boson1=0.25*(5+4*TMath::Cos(THETA_p)+15*TMath::Cos(2*THETA_p));
		P4_plus_Boson1=5*(TMath::Sin(THETA_p)*TMath::Sin(THETA_p))*(TMath::Cos(2*PHI_p));
		P5_plus_Boson1=5*(TMath::Sin(THETA_p)*TMath::Sin(THETA_p))*(TMath::Sin(2*PHI_p));
		P6_plus_Boson1=(sqrt(2))*(TMath::Sin(THETA_p))*((-5)*TMath::Cos(THETA_p)+1)*(TMath::Cos(PHI_p));
		P7_plus_Boson1=(sqrt(2))*(TMath::Sin(THETA_p))*((-5)*TMath::Cos(THETA_p)+1)*(TMath::Sin(PHI_p));
		P8_plus_Boson1=(1.0/4.0/sqrt(3))*((-5)+12*TMath::Cos(THETA_p)-15*TMath::Cos(2*THETA_p));

		P1_plus_Boson2=(sqrt(2))*(TMath::Sin(THETA_m))*(5*TMath::Cos(THETA_m)-1)*(TMath::Cos(PHI_m));
		P2_plus_Boson2=(sqrt(2))*(TMath::Sin(THETA_m))*(5*TMath::Cos(THETA_m)-1)*(TMath::Sin(PHI_m));
		P3_plus_Boson2=0.25*(5-4*TMath::Cos(THETA_m)+15*TMath::Cos(2*THETA_m));
		P4_plus_Boson2=5*(TMath::Sin(THETA_m)*TMath::Sin(THETA_m))*(TMath::Cos(2*PHI_m));
		P5_plus_Boson2=5*(TMath::Sin(THETA_m)*TMath::Sin(THETA_m))*(TMath::Sin(2*PHI_m));
		P6_plus_Boson2=(sqrt(2))*(TMath::Sin(THETA_m))*((-5)*TMath::Cos(THETA_m)-1)*(TMath::Cos(PHI_m));
		P7_plus_Boson2=(sqrt(2))*(TMath::Sin(THETA_m))*((-5)*TMath::Cos(THETA_m)-1)*(TMath::Sin(PHI_m));
		P8_plus_Boson2=(1.0/4.0/sqrt(3))*((-5)-12*TMath::Cos(THETA_m)-15*TMath::Cos(2*THETA_m));

		P1_tilde_Boson1=0;
		P2_tilde_Boson1=0;
		P3_tilde_Boson1=0;
		P4_tilde_Boson1=0;
		P5_tilde_Boson1=0;
		P6_tilde_Boson1=0;
		P7_tilde_Boson1=0;
		P8_tilde_Boson1=0;

		P1_tilde_Boson2=0;
		P2_tilde_Boson2=0;
		P3_tilde_Boson2=0;
		P4_tilde_Boson2=0;
		P5_tilde_Boson2=0;
		P6_tilde_Boson2=0;
		P7_tilde_Boson2=0;
		P8_tilde_Boson2=0;

		double* P_p[8] = {&P1_plus_Boson1, &P2_plus_Boson1, &P3_plus_Boson1, &P4_plus_Boson1, &P5_plus_Boson1, &P6_plus_Boson1, &P7_plus_Boson1, &P8_plus_Boson1}; 
		double* P_tilde_Boson1[8] = {&P1_tilde_Boson1, &P2_tilde_Boson1, &P3_tilde_Boson1, &P4_tilde_Boson1, &P5_tilde_Boson1, &P6_tilde_Boson1, &P7_tilde_Boson1, &P8_tilde_Boson1}; 

		for(int nindex=1; nindex<=8 ; nindex++) {
			for(int mindex=1; mindex<=8 ; mindex++) {
				(*(P_tilde_Boson1[nindex-1])) += ( (a_matrix[nindex-1][mindex-1]) * (*(P_p[mindex-1])) );
			}
		}

		double* P_m[8] = {&P1_plus_Boson2, &P2_plus_Boson2, &P3_plus_Boson2, &P4_plus_Boson2, &P5_plus_Boson2, &P6_plus_Boson2, &P7_plus_Boson2, &P8_plus_Boson2}; 
		double* P_tilde_Boson2[8] = {&P1_tilde_Boson2, &P2_tilde_Boson2, &P3_tilde_Boson2, &P4_tilde_Boson2, &P5_tilde_Boson2, &P6_tilde_Boson2, &P7_tilde_Boson2, &P8_tilde_Boson2}; 

		for(int nindex=1; nindex<=8 ; nindex++) {
			for(int mindex=1; mindex<=8 ; mindex++) {
				(*(P_tilde_Boson2[nindex-1])) += ( (a_matrix[nindex-1][mindex-1]) * (*(P_m[mindex-1])) );
			}
		}

		//////////////////////////////////////////////////////////////////////////////
		tree->Fill();
	}
	bool event_exist = ( (tree->GetEntries() > 0) ? true : false );
	if(event_exist) cout << filename << " contains " << tree->GetEntries() << " entries " << endl;
	file.Write();
	file.Close();

	return event_exist;
}



std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> OptimizeCuts(std::vector<TString> filesin, std::vector<int> sample_types, string& optcut, double ndivisions=20) {

	gStyle->SetOptStat(0000);

	int nfiles = filesin.size();
	int signfiles=0;
	for(int fi=0; fi<nfiles; fi++){
		if(sample_types.at(fi)==0 || sample_types.at(fi)==2) signfiles++;
	}
	int bkgnfiles=0;
	for(int fi=0; fi<nfiles; fi++){
		if(sample_types.at(fi)==1) bkgnfiles++;
	}

	TFile fin1(filesin.at(0), "READ");
	TTree* tin1 = (TTree*)fin1.Get("tree");
	TObjArray* objarr = tin1->GetListOfBranches();
	TObjArray* objarr1 = (TObjArray*)objarr->Clone();
	int objnum = objarr1->GetEntries();
	delete tin1;
	fin1.Close();

	std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> output;
	std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> output_sorted;
	std::pair<double, double> cut_limits;
	std::pair<double, double> n_sig_and_n_bkg;
	double significance=0;
	double n_sig_for_significance=0;
	double n_bkg_for_significance=0;
	string var_name="";
	double left_cut=0;
	double right_cut=0;
	double xmin=0;
	double xmax=100;
	double xrange=0;
	double n_sig=0;
	double n_bkg=0;
	double integral_left=0;
	double integral_right=100;
	double width_step=0;
	for(int k=0; k<objnum; k++){
		std::vector<TH1D*> opt_temps;
		xmin = 0;
		xmax = 100;
		xrange = xmax - xmin;
		width_step = xrange/ndivisions;
		for(int fid=0; fid<nfiles; fid++){ 
			if( !((TString(filesin.at(fid))).Contains("mumuTozzzTo")) ) continue;
			if( sample_types.at(fid) != 0 ) continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 0);
			chin.SetBranchStatus(((objarr1->At(k)))->GetName(), 1);
			chin.SetBranchStatus("SF", 1);
			chin.SetBranchStatus("event_weight", 1);
			if( (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) ) {
				chin.Draw(Form("%s>>opt_temp%i", (objarr1->At(k))->GetName(), k), Form("(%s > -99)*(event_weight)",((objarr1->At(k)))->GetName()), "goff");
				//double xmean = ((TH1F*)(gDirectory->Get(Form("th_temp%i", k))))->GetMean();
				//double xstaddev = ((TH1F*)(gDirectory->Get(Form("th_temp%i", k))))->GetStdDev();
				//xmin = chin.GetMinimum(((objarr1->At(k)))->GetName());
				//xmax = chin.GetMaximum(((objarr1->At(k)))->GetName());
				if( ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k)))) != 0 ) {
					xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge(1);
					xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetNcells() );
					xrange = xmax - xmin;
					width_step = xrange/ndivisions;
					xmin -= 2.*width_step;
					xmax += 2.*width_step;
				}
				else{
					xmin = 0;
					xmax = 100;
					xrange = xmax - xmin;
					width_step = xrange/ndivisions;
					xmin -= 2.0*width_step;
					xmax += 2.0*width_step;
				}
			}
			else if( (string(((objarr1->At(k)))->GetName()) == string("lep_fst")) ) {
				xmin = 0;
				xmax = 100;
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
		}
		for(int fid=0; fid<nfiles; fid++){ 
			if( sample_types.at(fid) > 1 ) continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 0);
			chin.SetBranchStatus(((objarr1->At(k)))->GetName(), 1);
			chin.SetBranchStatus("SF", 1);
			chin.SetBranchStatus("event_weight", 1);
			opt_temps.push_back(new TH1D(Form("opt_file%i_%s",fid,((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			if(string(((objarr1->At(k)))->GetName()) != string("lep_fst")) {
				if( ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k)))) != 0 ) {
					chin.Draw(Form("%s>>%s", ((objarr1->At(k)))->GetName(), Form("opt_file%i_%s",fid,((objarr1->At(k)))->GetName())), Form("10000000*SF*(%s > -99)*(event_weight)",((objarr1->At(k)))->GetName()), "goff");
				}
				else {
					((TH1F*)(gDirectory->Get(Form("opt_file%i_%s",fid,((objarr1->At(k)))->GetName()))))->Fill(1E7);
				}
			}
			else { 
				chin.Draw(Form("%s>>%s", ((objarr1->At(k)))->GetName(), Form("opt_file%i_%s",fid,((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)"), "goff");
			}

		}
		significance = 0;
		n_sig_for_significance = 0;
		n_bkg_for_significance = 0;
		var_name = string( ((objarr1->At(k)))->GetName() );
		//for(int nwidth=1; nwidth<=ndivisions; nwidth++) 
		//for(int nstart=1; (nstart+nwidth-1)<=(ndivisions+4); nstart++) 

		for(int nstart=3; nstart<=(ndivisions+2); nstart++) {
			for(int nwidth=1; (nstart+nwidth-1)<=(ndivisions+4); nwidth++) {
				if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("met")) || ((TString(var_name)).Contains("deltaR_ll1")) || ((TString(var_name)).Contains("delta_ll2")) ) && !( ((TString(var_name)).Contains("met_byhand")) || ((TString(var_name)).Contains("meta")) ) ) {
					nwidth=(ndivisions+4+1-nstart);
				}
				n_sig = 0;
				n_bkg = 0;
				for(int fid=0; fid<nfiles; fid++) { 
					if( sample_types.at(fid) == 0 ) { 
						n_sig += ( (opt_temps.at(fid))->Integral(nstart,(nstart+nwidth-1)) );
					}
					if( sample_types.at(fid) == 1 ) { 
						n_bkg += ( (opt_temps.at(fid))->Integral(nstart,(nstart+nwidth-1)) );
					}
				}
				if( (n_sig>0 && n_bkg>0) && ((n_sig/sqrt(n_bkg)) > significance) ) { 
					significance = (n_sig/sqrt(n_bkg));
					n_sig_for_significance = n_sig;
					n_bkg_for_significance = n_bkg;
					integral_left = (xmin+((nstart-1)*width_step));
					integral_right = (xmin+((nstart-1)*width_step)+(nwidth*width_step));
				}
				if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("met")) || ((TString(var_name)).Contains("deltaR_ll1")) || ((TString(var_name)).Contains("deltaR_ll2")) ) && !( ((TString(var_name)).Contains("met_byhand")) || ((TString(var_name)).Contains("meta")) ) ) {
					break;
				}
			}
		}
		cut_limits.first  = integral_left;
		cut_limits.second = integral_right;
		n_sig_and_n_bkg.first  = n_sig_for_significance;
		n_sig_and_n_bkg.second = n_bkg_for_significance;
		(output[var_name]).first.first  = cut_limits;
		(output[var_name]).first.second = n_sig_and_n_bkg;
		(output[var_name]).second = significance;
	}

	//=======================================================================
	string cut_series = "( ";
	//-----------------------------------------------------------------------
	unsigned int n_var = output.size();
	std::vector<int> index(n_var);
	std::vector<double> significances;
	std::vector<string> var_names;
	std::vector<string> opt_cut_var_names;
	for (const auto& [key, value] : output) {
		significances.push_back(double(value.second));
		var_names.push_back(string(key));
	}
	TMath::SortItr(significances.begin(), significances.end(), index.begin(), true );

	cout<<"****************************************************"<<endl;
	for(int i=0; i<n_var; i++) {
		unsigned int j = index[i];
		(output_sorted[string(var_names.at(j))]).first.first.first    = (output[var_names.at(j)]).first.first.first  ;
		(output_sorted[string(var_names.at(j))]).first.first.second   = (output[var_names.at(j)]).first.first.second ;
		(output_sorted[string(var_names.at(j))]).first.second.first   = (output[var_names.at(j)]).first.second.first  ;
		(output_sorted[string(var_names.at(j))]).first.second.second  = (output[var_names.at(j)]).first.second.second ;
		(output_sorted[string(var_names.at(j))]).second = (output[var_names.at(j)]).second;
		if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("met")) || ((TString(var_name)).Contains("deltaR_ll1")) || ((TString(var_name)).Contains("delta_ll2")) ) && !( ((TString(var_name)).Contains("met_byhand")) || ((TString(var_name)).Contains("meta")) ) ) {
			if( (output_sorted[string(var_names.at(j))]).first.first.first   < xmin+2*width_step ) (output_sorted[string(var_names.at(j))]).first.first.first  = -1E7 ;
			if( (output_sorted[string(var_names.at(j))]).first.first.second  > xmax-2*width_step ) (output_sorted[string(var_names.at(j))]).first.first.second = 1E7 ;
		}
		if( ( ((TString(var_names.at(j))).Contains("M4l")) || ((TString(var_names.at(j))).Contains("Pt4l")) || ((TString(var_names.at(j))).Contains("met")) || ((TString(var_names.at(j))).Contains("deltaR_ll1")) || ((TString(var_names.at(j))).Contains("deltaR_ll2")) ) && !( ((TString(var_names.at(j))).Contains("met_byhand")) || ((TString(var_names.at(j))).Contains("meta")) ) ) {
			cut_series += string(Form("(%s >= %f && %s <= %f)", var_names.at(j).c_str(), (output[var_names.at(j)]).first.first.first, var_names.at(j).c_str(), (output[var_names.at(j)]).first.first.second));
			if( !(cut_series.empty()) ) cut_series += " && ";
			std::cout << '{' << var_names.at(j) << "}: cuts = [" << (output[var_names.at(j)]).first.first.first << ", " << (output[var_names.at(j)]).first.first.second << "] ; significance = " << (output[var_names.at(j)]).second << endl;
			opt_cut_var_names.push_back(var_names.at(j).c_str());
		}
	}
	//-----------------------------------------------------------------------
	//cut_series += " (1==1) ";
	//for (const auto& [key, value] : output) {
	//	std::cout << '{' << key << "}: cuts = [" << value.first.first.first << ", " << value.first.first.second << "] ; significance = " << value.second << endl;
	//	if( ((TString(key)).Contains("deltaM1")) ) continue;
	//	//if( (value.second > 0.01) && (!((TString(key)).Contains("num")) && !((TString(key)).Contains("SF")) && !((TString(key)).Contains("event")) && !((TString(key)).Contains("lepton_fst")) && !((TString(key)).Contains("deltaM1"))) ) 
	//	if( ((TString(key)).Contains("M4l")) || ((TString(key)).Contains("Pt4l")) || ((TString(key)).Contains("deltaM")) || ((TString(key)).Contains("met")) || ((TString(key)).Contains("Mrecoil")) || ((TString(key)).Contains("deltaR_ll1")) || ((TString(key)).Contains("deltaR_ll2")) || ((TString(key)).Contains("Mll1")) || ((TString(key)).Contains("Mll2")) || ((TString(key)).Contains("Ptll1")) || ((TString(key)).Contains("Ptll2")) ) 
	//	{
	//		cut_series += string(Form(" && (%s >= %f && %s <= %f)", key.c_str(), value.first.first.first, key.c_str(), value.first.first.second));
	//	}
	//}
	//-----------------------------------------------------------------------
	cut_series += " ( (Mrecoil>0 && Mrecoil<200) && (Mll1>80 && Mll1<100) && (Mll2>80 && Mll2<100) && (deltaM<20) ) ";
	cut_series += " )";
	//cout<<"****************************************************"<<endl;
	//cout<<"cut_series is :"<<endl<<endl;
	//cout<<cut_series<<endl<<endl;
	//cout<<"****************************************************"<<endl;

	//=======================================================================

	std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> series_output;
	std::pair<double, double> series_cut_limits;
	std::pair<double, double> series_n_sig_and_n_bkg;
	//string opt_cut_series = " (1==1) ";
	string opt_cut_series = " ( (Mrecoil>0 && Mrecoil<200) && (Mll1>80 && Mll1<100) && (Mll2>80 && Mll2<100) && (deltaM<20) ) ";
	for(int k=0; k<opt_cut_var_names.size(); k++){
		std::vector<TH1D*> hist_for_cuts;
		xmin = 0;
		xmax = 100;
		xrange = xmax - xmin;
		width_step = xrange/ndivisions;
		for(int fid=0; fid<nfiles; fid++){ 
			if( sample_types.at(fid) != 0 ) continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 1);
			chin.SetBranchStatus("SF", 1);
			chin.SetBranchStatus("event_weight", 1);
			chin.Draw(Form("%s>>hist_for_cut%i", (opt_cut_var_names.at(k)).c_str(), k), Form("10000000*SF*(%s > -99)*(event_weight)", (opt_cut_var_names.at(k)).c_str()), "goff");
			//double xmean = ((TH1F*)(gDirectory->Get(Form("th_temp%i", k))))->GetMean();
			//double xstaddev = ((TH1F*)(gDirectory->Get(Form("th_temp%i", k))))->GetStdDev();
			//xmin = chin.GetMinimum((opt_cut_var_names.at(k)).c_str());
			//xmax = chin.GetMaximum((opt_cut_var_names.at(k)).c_str());
			if( ((TH1F*)(gDirectory->Get(Form("hist_for_cut%i", k)))) != 0 ) {
				xmin = ((TH1F*)(gDirectory->Get(Form("hist_for_cut%i", k))))->GetBinLowEdge(1);
				xmax = ((TH1F*)(gDirectory->Get(Form("hist_for_cut%i", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("hist_for_cut%i", k))))->GetNcells() );
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.*width_step;
				xmax += 2.*width_step;
			}
			else{
				xmin = 0;
				xmax = 100;
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
		}
		for(int fid=0; fid<nfiles; fid++){ 
			if( sample_types.at(fid) > 1 ) continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 1);
			chin.SetBranchStatus("SF", 1);
			chin.SetBranchStatus("event_weight", 1);
			hist_for_cuts.push_back(new TH1D(Form("opt_hist_for_cuts%i_%s",fid,(opt_cut_var_names.at(k)).c_str()), Form("%s",(opt_cut_var_names.at(k)).c_str()), ndivisions+4, xmin, xmax));
			if( ((TH1F*)(gDirectory->Get(Form("hist_for_cut%i", k)))) != 0 ) {
				chin.Draw(Form("%s>>%s", (opt_cut_var_names.at(k)).c_str(), Form("opt_hist_for_cuts%i_%s",fid,(opt_cut_var_names.at(k)).c_str())), Form("10000000*SF*(%s > -99)*(event_weight)*(%s)",(opt_cut_var_names.at(k)).c_str(),opt_cut_series.c_str()), "goff");
			}
			else {
				((TH1F*)(gDirectory->Get(Form("opt_hist_for_cuts%i_%s",fid,(opt_cut_var_names.at(k)).c_str()))))->Fill(1E7);
			}
		}
		significance = 0;
		n_sig_for_significance = 0;
		n_bkg_for_significance = 0;
		var_name = string( opt_cut_var_names.at(k) );
		for(int nstart=3; nstart<=(ndivisions+2); nstart++) {
			for(int nwidth=1; (nstart+nwidth-1)<=(ndivisions+4); nwidth++) {
				if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("met")) || ((TString(var_name)).Contains("deltaR_ll1")) || ((TString(var_name)).Contains("delta_ll2")) ) && !( ((TString(var_name)).Contains("met_byhand")) || ((TString(var_name)).Contains("meta")) ) ) {
					nwidth=(ndivisions+4+1-nstart);
				}
				n_sig = 0;
				n_bkg = 0;
				for(int fid=0; fid<nfiles; fid++) { 
					if( sample_types.at(fid) == 0 ) { 
						n_sig += ( (hist_for_cuts.at(fid))->Integral(nstart,(nstart+nwidth-1)) );
					}
					if( sample_types.at(fid) == 1 ) { 
						n_bkg += ( (hist_for_cuts.at(fid))->Integral(nstart,(nstart+nwidth-1)) );
					}
				}
				if( (n_sig>0 && n_bkg>0) && ((n_sig/sqrt(n_bkg)) > significance) ) { 
					significance = (n_sig/sqrt(n_bkg));
					n_sig_for_significance = n_sig;
					n_bkg_for_significance = n_bkg;
					integral_left = (xmin+((nstart-1)*width_step));
					integral_right = (xmin+((nstart-1)*width_step)+(nwidth*width_step));
				}
				if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("met")) || ((TString(var_name)).Contains("deltaR_ll1")) || ((TString(var_name)).Contains("delta_ll2")) ) && !( ((TString(var_name)).Contains("met_byhand")) || ((TString(var_name)).Contains("meta")) ) ) {
					break;
				}
			}
		}
		if( integral_left < xmin+2*width_step )  integral_left = -1E7 ;
		if( integral_right > xmax-2*width_step ) integral_right = 1E7 ;
		series_cut_limits.first  = integral_left;
		series_cut_limits.second = integral_right;
		series_n_sig_and_n_bkg.first  = n_sig_for_significance;
		series_n_sig_and_n_bkg.second = n_bkg_for_significance;
		(series_output[var_name]).first.first  = cut_limits;
		(series_output[var_name]).first.second = n_sig_and_n_bkg;
		(series_output[var_name]).second = significance;
		opt_cut_series += string(Form(" && (%s >= %f && %s <= %f)", var_name.c_str(), integral_left, var_name.c_str(), integral_right));
	}




	cout<<"****************************************************"<<endl;
	cout<<"optimized cut_series is :"<<endl<<endl;
	cout<<opt_cut_series<<endl<<endl;
	cout<<"****************************************************"<<endl;

	//optcut = cut_series;
	optcut = opt_cut_series;
	return output;
}



double ComputeSignificance(std::vector<TString> filesin, std::vector<int> sample_types, std::vector<int> sig_types, std::vector<int> bkg_types, string cut="") {

	gStyle->SetOptStat(0000);

	int nfiles = filesin.size();
	int signfiles=0;
	for(int fi=0; fi<nfiles; fi++){
		if(sample_types.at(fi)==0 || sample_types.at(fi)==2) signfiles++;
	}
	int bkgnfiles=0;
	for(int fi=0; fi<nfiles; fi++){
		if(sample_types.at(fi)==0) bkgnfiles++;
	}

	bool is_sig=false;
	bool is_bkg=false;
	double significance=0;
	double n_sig=0;
	double n_bkg=0;
	double total_n_bkg=0;
	std::map<string, double> samples_yield_after_cut;
	cout<<cut.c_str()<<endl;
	for(int fid=0; fid<nfiles; fid++) { 
		is_sig=false;
		is_bkg=false;
		for(int k=0; k<sig_types.size(); k++) {
			if(sample_types.at(fid)==sig_types.at(k)) {
				is_sig=true;
			}
		}
		for(int k=0; k<bkg_types.size(); k++) {
			if(sample_types.at(fid)==bkg_types.at(k)) {
				is_bkg=true;
			}
		}
		if( !is_sig && !is_bkg ) continue;
		TChain chin("tree");
		chin.Add(filesin.at(fid));
		chin.SetBranchStatus("*", 1);
		cout<<filesin.at(fid)<<" = "<<chin.GetEntries(cut.c_str())<<endl;
		string process((filesin.at(fid)).Data());
		//process.erase(process.find("_delphes_preselected.root",0), string("_delphes_preselected.root").size());
		process.erase(process.find(".root",0), string(".root").size());
		//if(sample_types.at(fi)==0)  
		//if((TString(process)).Contains("mumuTozzzTo")) 
		if(is_sig)
		{
			chin.Draw(Form("%s>>event_th_file%i", "event", fid), Form("10000000*SF*(event_weight)*(%s)", cut.c_str()), "goff");
			if( ((TH1F*)(gDirectory->Get(Form("event_th_file%i", fid)))) != NULL ) {
				n_sig += ((TH1F*)(gDirectory->Get(Form("event_th_file%i", fid))))->Integral();
			}
			cout<<"n_sig = "<<n_sig<<endl;
			samples_yield_after_cut[process] = n_sig;
		}
		//if(sample_types.at(fi)==1)  
		//if( !((TString(process)).Contains("mumuTozzzTo")) ) 
		if(is_bkg)
		{
			chin.Draw(Form("%s>>event_th_file%i", "event", fid), Form("10000000*SF*(event_weight)*(%s)", cut.c_str()), "goff");
			if( ((TH1F*)(gDirectory->Get(Form("event_th_file%i", fid)))) != NULL ) {
				n_bkg = ((TH1F*)(gDirectory->Get(Form("event_th_file%i", fid))))->Integral();
				samples_yield_after_cut[process] = n_bkg;
				total_n_bkg += n_bkg;
				cout<<"total_n_bkg = "<<total_n_bkg<<endl;
			}
		}
	}
	cout<<"Summary of yields: ==============================================="<<endl;
	int n_sample = samples_yield_after_cut.size();
	std::vector<int> eveindex(n_sample);
	std::vector<double> events_yield;
	std::vector<string> sample_name;
	for (const auto& [key, value] : samples_yield_after_cut) {
		events_yield.push_back(double(value));
		sample_name.push_back(string(key));
	}
	TMath::SortItr(events_yield.begin(), events_yield.end(), eveindex.begin(), true );

	//----------------------- Print yields summary with sorted order
	for(int k=0; k<events_yield.size(); k++){
		if((TString(sample_name.at(eveindex.at(k)))).Contains("mumuTozzzTo")) {
			cout<<sample_name.at(eveindex.at(k))<<"+++ : signal yields after cut = "<<events_yield.at(eveindex.at(k))<<endl;
		}
		if( !((TString(sample_name.at(eveindex.at(k)))).Contains("mumuTozzzTo")) ) {
			cout<<sample_name.at(eveindex.at(k))<<"--- : background yields after cut = "<<events_yield.at(eveindex.at(k))<<endl;
		}
	}
	cout<<"=================================================================="<<endl;
	if( (n_sig>0 && total_n_bkg>0) ) { 
		significance = (n_sig/sqrt(total_n_bkg));
	}

	return significance;
}



void PlotEvents(std::vector<TString> filesin, std::vector<bool> if_event_exist, std::vector<int> sample_types, string plot_prefix, string cut, int ndivisions, string para_name, double RWstart, double RWstep, int RWnstep) {

	gStyle->SetOptStat(0000);
	for(int ifi=0; ifi<if_event_exist.size(); ifi++) {
		if(filesin.empty()) return;
		if( !(if_event_exist.at(ifi)) ) {
			filesin.erase(filesin.begin()+ifi);
			sample_types.erase(sample_types.begin()+ifi);
		}
	}
	if(filesin.empty()) return;
	int nfiles = filesin.size();
	TFile fin1(filesin.at(0), "READ");
	TTree* tin1 = (TTree*)fin1.Get("tree");
	TObjArray* objarr = tin1->GetListOfBranches();
	TObjArray* objarr1 = (TObjArray*)objarr->Clone();
	int objnum = objarr1->GetEntries();
	delete tin1;
	fin1.Close();

	TCanvas* c_hstacks = new TCanvas("c_hstacks","c_hstacks",2500,1900);
	double xmin=0;
	double xmax=100;
	double xrange=0;
	double width_step=0;
	double signalhistmaxy = 0;

	TH1F* sum_event_weight_temps;
	TH1F* ph_temps[8];
	TH1F* mh_temps[8];
	TH2F* pmh2d_temps[8][8];
	double fa_temps[8];
	double ga_temps[8];
	double hab_temps[8][8];

	TH2F* MvvCosSA_ph_temps[8]; //P~ of boson1
	TH2F* MvvCosSA_mh_temps[8]; //P~ of boson2
	TH2F* MvvCosSA_pmh2d_temps[8][8]; ///P~ of boson1 * P~ of boson2
	TH2F* MvvCosSA_sum_event_weight_temps = new TH2F("MvvCosSA_sum_event_weight_temps", "Event weight ; M_{VV}(M_{4l}) ; cos(#Theta)", 4, 950, 1050, 20, -1.0, 1.0);
	TH2F* MvvCosSA_L_2_temps = new TH2F("MvvCosSA_L_2_temps", "C_{2} ; M_{VV}(M_{4l}) ; cos(#Theta)", 4, 950, 1050, 20, -1.0, 1.0);
	TH2F* MvvCosSA_I_3_temps = new TH2F("MvvCosSA_I_3_temps", "I_{3} ; M_{VV}(M_{4l}) ; cos(#Theta)", 4, 950, 1050, 20, -1.0, 1.0);

	for(int index=1; index<=8; index++) {
		fa_temps[index-1] = 0;
		ga_temps[index-1] = 0;
		MvvCosSA_ph_temps[index-1] = new TH2F(Form("MvvCosSA_ph_temps_%i",index), Form("P%i_tilde_Boson1 ; M_{VV}(M_{4l}) ; cos(#Theta)",index), 4, 950, 1050, 20, -1.0, 1.0);
		MvvCosSA_mh_temps[index-1] = new TH2F(Form("MvvCosSA_mh_temps_%i",index), Form("P%i_tilde_Boson2 ; M_{VV}(M_{4l}) ; cos(#Theta)",index), 4, 950, 1050, 20, -1.0, 1.0);
		for(int index2=1; index2<=8; index2++) {
			hab_temps[index-1][index2-1] = 0;
			MvvCosSA_pmh2d_temps[index-1][index2-1] = new TH2F(Form("MvvCosSA_pmh2d_temps_%i_%i",index,index2), Form("(P%i_tilde_Boson1 #times P%i_tilde_Boson2) ; M_{VV}(M_{4l}) ; cos(#Theta)",index,index2), 4, 950, 1050, 20, -1.0, 1.0);
		}
	}

	for(int fid=0; fid<nfiles; fid++){ 
		TChain chin("tree");
		chin.Add(filesin.at(fid));
		chin.SetBranchStatus("*", 1);
		chin.Draw(Form("event_weight>>+sum_event_weight_temps"), Form("10000000*SF*(event_weight)"), "goff");
		chin.Draw(Form("Cos_Scatter_Angle:M4l>>+MvvCosSA_sum_event_weight_temps"), Form("10000000*SF*(event_weight)"), "goff");
		sum_event_weight_temps = (TH1F*)(gDirectory->Get(Form("sum_event_weight_temps")));
		for(int index=1; index<=8; index++) {
			//chin.Draw(Form("%s>>+ph_temps_%i", Form("P%i_plus_Boson1",index), index), Form("10000000*SF*(event_weight)*(P%i_plus_Boson1)",index), "goff");
			//chin.Draw(Form("%s>>+mh_temps_%i", Form("P%i_plus_Boson2",index), index), Form("10000000*SF*(event_weight)*(P%i_plus_Boson2)",index), "goff");
			chin.Draw(Form("%s>>+ph_temps_%i", Form("P%i_tilde_Boson1",index), index), Form("10000000*SF*(event_weight)*(P%i_tilde_Boson1)*0.5",index), "goff");
			chin.Draw(Form("%s>>+mh_temps_%i", Form("P%i_tilde_Boson2",index), index), Form("10000000*SF*(event_weight)*(P%i_tilde_Boson2)*0.5",index), "goff");
			chin.Draw(Form("Cos_Scatter_Angle:M4l>>+MvvCosSA_ph_temps_%i", index), Form("10000000*SF*(event_weight)*(P%i_tilde_Boson1)*0.5",index), "goff");
			chin.Draw(Form("Cos_Scatter_Angle:M4l>>+MvvCosSA_mh_temps_%i", index), Form("10000000*SF*(event_weight)*(P%i_tilde_Boson2)*0.5",index), "goff");
			for(int index2=1; index2<=8; index2++) {
				//chin.Draw(Form("%s:%s>>+pmh2d_temps_%i_%i", Form("P%i_plus_Boson2",index2), Form("P%i_plus_Boson1",index), index, index2), Form("10000000*SF*(event_weight)*(P%i_plus_Boson1)*(P%i_plus_Boson2)",index,index2), "goff"); // %s:%s is act as y:x, which means P%i_plus_Boson1 in X axis and P%i_plus_Boson2 in Yaxis
				chin.Draw(Form("%s:%s>>+pmh2d_temps_%i_%i", Form("P%i_tilde_Boson2",index2), Form("P%i_tilde_Boson1",index), index, index2), Form("10000000*SF*(event_weight)*(P%i_tilde_Boson1)*(P%i_tilde_Boson2)*0.25",index,index2), "goff");
				chin.Draw(Form("Cos_Scatter_Angle:M4l>>+MvvCosSA_pmh2d_temps_%i_%i",index,index2), Form("10000000*SF*(event_weight)*(P%i_tilde_Boson1)*(P%i_tilde_Boson2)*0.25",index,index2), "goff");
			}
		}
	}

	double sum_event_weight = sum_event_weight_temps->Integral();
	cout<<"sum_event_weight = "<<sum_event_weight<<endl;
	for(int index=1; index<=8; index++) {
		ph_temps[index-1] = ((TH1F*)(gDirectory->Get(Form("ph_temps_%i", index))));
		mh_temps[index-1] = ((TH1F*)(gDirectory->Get(Form("mh_temps_%i", index))));
		fa_temps[index-1] = ((ph_temps[index-1])->Integral())/sum_event_weight;
		ga_temps[index-1] = ((mh_temps[index-1])->Integral())/sum_event_weight;
		(MvvCosSA_ph_temps[index-1])->Divide(MvvCosSA_sum_event_weight_temps);
		(MvvCosSA_mh_temps[index-1])->Divide(MvvCosSA_sum_event_weight_temps);
		for(int index2=1; index2<=8; index2++) {
			pmh2d_temps[index-1][index2-1] = ((TH2F*)(gDirectory->Get(Form("pmh2d_temps_%i_%i", index, index2))));
			hab_temps[index-1][index2-1] = ((pmh2d_temps[index-1][index2-1])->Integral())/sum_event_weight;
			(MvvCosSA_pmh2d_temps[index-1][index2-1])->Divide(MvvCosSA_sum_event_weight_temps);
		}
	}


	double SUM_SQ_fa = 0;
	double SUM_SQ_ga = 0;
	double SUM_SQ_hab = 0;
	for(int ibin=1; ibin<=(MvvCosSA_sum_event_weight_temps->GetNcells()-2); ibin++){
		SUM_SQ_fa = 0;
		SUM_SQ_ga = 0;
		SUM_SQ_hab = 0;
		for(int index=1; index<=8; index++) {
			SUM_SQ_fa += ( (MvvCosSA_ph_temps[index-1]->GetBinContent(ibin))*(MvvCosSA_ph_temps[index-1]->GetBinContent(ibin)) );
			SUM_SQ_ga += ( (MvvCosSA_mh_temps[index-1]->GetBinContent(ibin))*(MvvCosSA_mh_temps[index-1]->GetBinContent(ibin)) );
			for(int index2=1; index2<=8; index2++) {
				SUM_SQ_hab += ((MvvCosSA_pmh2d_temps[index-1][index2-1]->GetBinContent(ibin))*(MvvCosSA_pmh2d_temps[index-1][index2-1]->GetBinContent(ibin)));
			}
		}
		MvvCosSA_L_2_temps->SetBinContent(ibin, (2 * TMath::Max( ( (-2.0/9.0) - (12.0*SUM_SQ_fa) + (6.0*SUM_SQ_ga) + (4.0*SUM_SQ_hab) ) , ( (-2.0/9.0) - (12.0*SUM_SQ_ga) + (6.0*SUM_SQ_fa) + (4.0*SUM_SQ_hab) ) )));
		MvvCosSA_I_3_temps->SetBinContent(ibin, (4.0*((MvvCosSA_pmh2d_temps[4-1][4-1]->GetBinContent(ibin)) + (MvvCosSA_pmh2d_temps[5-1][5-1]->GetBinContent(ibin))) - (4.0/sqrt(3))*((MvvCosSA_pmh2d_temps[6-1][1-1]->GetBinContent(ibin)) + (MvvCosSA_pmh2d_temps[6-1][6-1]->GetBinContent(ibin)) + (MvvCosSA_pmh2d_temps[7-1][2-1]->GetBinContent(ibin)) + (MvvCosSA_pmh2d_temps[7-1][7-1]->GetBinContent(ibin)) + (MvvCosSA_pmh2d_temps[1-1][1-1]->GetBinContent(ibin)) + (MvvCosSA_pmh2d_temps[1-1][6-1]->GetBinContent(ibin)) + (MvvCosSA_pmh2d_temps[2-1][2-1]->GetBinContent(ibin)) + (MvvCosSA_pmh2d_temps[2-1][7-1]->GetBinContent(ibin)))));
	}

	TFile MvvCosSA_hists("MvvCosSA_hists.root", "RECREATE");
	MvvCosSA_hists.cd();
	for(int index=1; index<=8; index++) {
		MvvCosSA_ph_temps[index-1]->Write();
		MvvCosSA_mh_temps[index-1]->Write();
		for(int index2=1; index2<=8; index2++) {
			MvvCosSA_pmh2d_temps[index-1][index2-1]->Write();
		}
	}
	MvvCosSA_sum_event_weight_temps->Write();
	MvvCosSA_L_2_temps->Write();
	MvvCosSA_I_3_temps->Write();
	MvvCosSA_hists.Close();

	for(int xbin=1; xbin<=(MvvCosSA_sum_event_weight_temps->GetXaxis()->GetNbins()); xbin++){
		for(int ybin=1; ybin<=(MvvCosSA_sum_event_weight_temps->GetYaxis()->GetNbins()); ybin++){
			ofstream text_ouput_MvvCosSA_ph_temps(Form("./text_ouput_MvvCosSA_fa_xbin%i_ybin%i.txt",xbin,ybin), ios::ate);
			ofstream text_ouput_MvvCosSA_mh_temps(Form("./text_ouput_MvvCosSA_ga_xbin%i_ybin%i.txt",xbin,ybin), ios::ate);
			ofstream text_ouput_MvvCosSA_pmh2d_temps(Form("./text_ouput_MvvCosSA_hab_xbin%i_ybin%i.txt",xbin,ybin), ios::ate);
			text_ouput_MvvCosSA_ph_temps<<"Mvv"<<"\t"<<"cos(Theta)"<<"\t"<<endl;
			text_ouput_MvvCosSA_ph_temps<<(MvvCosSA_ph_temps[1-1]->GetXaxis()->GetBinCenter(xbin))<<"\t"<<(MvvCosSA_ph_temps[1-1]->GetYaxis()->GetBinCenter(ybin))<<"\t"<<endl<<"=========================================="<<endl<<endl;
			text_ouput_MvvCosSA_ph_temps<<"\t"<<"1         "<<"\t"<<"2         "<<"\t"<<"3         "<<"\t"<<"4         "<<"\t"<<"5         "<<"\t"<<"6         "<<"\t"<<"7         "<<"\t"<<"8         "<<"\t"<<endl;
			text_ouput_MvvCosSA_mh_temps<<"Mvv"<<"\t"<<"cos(Theta)"<<"\t"<<endl;
			text_ouput_MvvCosSA_mh_temps<<(MvvCosSA_mh_temps[1-1]->GetXaxis()->GetBinCenter(xbin))<<"\t"<<(MvvCosSA_mh_temps[1-1]->GetYaxis()->GetBinCenter(ybin))<<"\t"<<endl<<"=========================================="<<endl<<endl;
			text_ouput_MvvCosSA_mh_temps<<"\t"<<"1         "<<"\t"<<"2         "<<"\t"<<"3         "<<"\t"<<"4         "<<"\t"<<"5         "<<"\t"<<"6         "<<"\t"<<"7         "<<"\t"<<"8         "<<"\t"<<endl;
			text_ouput_MvvCosSA_pmh2d_temps<<"Mvv"<<"\t"<<"cos(Theta)"<<"\t"<<endl;
			text_ouput_MvvCosSA_pmh2d_temps<<(MvvCosSA_pmh2d_temps[1-1][1-1]->GetXaxis()->GetBinCenter(xbin))<<"\t"<<(MvvCosSA_pmh2d_temps[1-1][1-1]->GetYaxis()->GetBinCenter(ybin))<<"\t"<<endl<<"=========================================="<<endl<<endl;
			text_ouput_MvvCosSA_pmh2d_temps<<"\t"<<"1         "<<"\t"<<"2         "<<"\t"<<"3         "<<"\t"<<"4         "<<"\t"<<"5         "<<"\t"<<"6         "<<"\t"<<"7         "<<"\t"<<"8         "<<"\t"<<endl;

			for(int index=1; index<=8; index++) {
				text_ouput_MvvCosSA_ph_temps<<"\t"<<(MvvCosSA_ph_temps[index-1])->GetBinContent(xbin,ybin);
				text_ouput_MvvCosSA_mh_temps<<"\t"<<(MvvCosSA_mh_temps[index-1])->GetBinContent(xbin,ybin);
				text_ouput_MvvCosSA_pmh2d_temps<<index;
				for(int index2=1; index2<=8; index2++) {
					text_ouput_MvvCosSA_pmh2d_temps<<"\t"<<(MvvCosSA_pmh2d_temps[index-1][index2-1])->GetBinContent(xbin,ybin);
				}
				text_ouput_MvvCosSA_pmh2d_temps<<"\t"<<endl;
			}
			text_ouput_MvvCosSA_ph_temps<<"\t"<<endl;
			text_ouput_MvvCosSA_mh_temps<<"\t"<<endl;

			text_ouput_MvvCosSA_ph_temps.close();
			text_ouput_MvvCosSA_mh_temps.close();
			text_ouput_MvvCosSA_pmh2d_temps.close();
		}
		//c_hstacks->cd();
		//(ph_temps[index-1])->SetTitle(Form(" ; P^{+}_{%i}; Events", index));
		//(ph_temps[index-1])->Draw("HIST");
		//c_hstacks->SaveAs(Form("./plots/P+_%i.png", index));
		//cout<<"ga_"<<index<<" = "<<ga_temps[index-1]<<endl;
		//cout<<"Clear Pad----------------------------------"<<endl;
		//c_hstacks->Clear();
		//c_hstacks->cd();
		//(mh_temps[index-1])->SetTitle(Form(" ; P^{-}_{%i}; Events", index));
		//(mh_temps[index-1])->Draw("HIST");
		//c_hstacks->SaveAs(Form("./plots/P-_%i.png", index));
	}

	for(int index=1; index<=8; index++) {
		TCanvas tc_MvvCosSA_ph_temps(Form("tc_%s",MvvCosSA_ph_temps[index-1]->GetName()),"");
		MvvCosSA_ph_temps[index-1]->Draw("COLZ");
		tc_MvvCosSA_ph_temps.SaveAs(Form("tc_%s.png",MvvCosSA_ph_temps[index-1]->GetName()));
		TCanvas tc_MvvCosSA_mh_temps(Form("tc_%s",MvvCosSA_mh_temps[index-1]->GetName()),"");
		MvvCosSA_mh_temps[index-1]->Draw("COLZ");
		tc_MvvCosSA_mh_temps.SaveAs(Form("tc_%s.png",MvvCosSA_mh_temps[index-1]->GetName()));
		for(int index2=1; index2<=8; index2++) {
			TCanvas tc_MvvCosSA_pmh2d_temps(Form("tc_%s",MvvCosSA_pmh2d_temps[index-1][index2-1]->GetName()),"");
			MvvCosSA_pmh2d_temps[index-1][index2-1]->Draw("COLZ");
			tc_MvvCosSA_pmh2d_temps.SaveAs(Form("tc_%s.png",MvvCosSA_pmh2d_temps[index-1][index2-1]->GetName()));
			//cout<<"hab_"<<index<<index2<<" = "<<hab_temps[index-1][index2-1]<<endl;
			//cout<<"Clear Pad----------------------------------"<<endl;
			//c_hstacks->Clear();
			//c_hstacks->cd();
			//(pmh2d_temps[((index-1)*8+(index2-1))])->SetTitle(Form(" ; P^{+}_{%i} ; P^{-}_{%i}", index, index2));
			//(pmh2d_temps[((index-1)*8+(index2-1))])->Draw("COLZ");
			//c_hstacks->SaveAs(Form("./plots/P_+%i_-%i.png", index, index2));
		}
	}
	TCanvas tc_MvvCosSA_L_2_temps(Form("tc_%s",MvvCosSA_L_2_temps->GetName()),"");
	MvvCosSA_L_2_temps->Draw("COLZ");
	tc_MvvCosSA_L_2_temps.SaveAs(Form("tc_%s.png",MvvCosSA_L_2_temps->GetName()));
	TCanvas tc_MvvCosSA_I_3_temps(Form("tc_%s",MvvCosSA_I_3_temps->GetName()),"");
	MvvCosSA_I_3_temps->Draw("COLZ");
	tc_MvvCosSA_I_3_temps.SaveAs(Form("tc_%s.png",MvvCosSA_I_3_temps->GetName()));


	for(int index=1; index<=8; index++) {
		SUM_SQ_fa += (fa_temps[index-1] * fa_temps[index-1]);
		SUM_SQ_ga += (ga_temps[index-1] * ga_temps[index-1]);
		for(int index2=1; index2<=8; index2++) {
			SUM_SQ_hab += (hab_temps[index-1][index2-1] * hab_temps[index-1][index2-1]);
		}
	}
	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Output evaluations of Quantum entanglement witness observables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	double L_2 = 2 * TMath::Max( ( (-2.0/9.0) - (12.0*SUM_SQ_fa) + (6.0*SUM_SQ_ga) + (4.0*SUM_SQ_hab) ) , ( (-2.0/9.0) - (12.0*SUM_SQ_ga) + (6.0*SUM_SQ_fa) + (4.0*SUM_SQ_hab) ) );
	L_2 = TMath::Max( L_2 , 0.0 );
	double I_3 = 4.0*(hab_temps[4-1][4-1] + hab_temps[5-1][5-1]) - (4.0/sqrt(3))*(hab_temps[6-1][1-1] + hab_temps[6-1][6-1] + hab_temps[7-1][2-1] + hab_temps[7-1][7-1] + hab_temps[1-1][1-1] + hab_temps[1-1][6-1] + hab_temps[2-1][2-1] + hab_temps[2-1][7-1]);
	cout<<"  Total average L_2 = "<<L_2<<endl;
	cout<<"  Total average I_3 = "<<I_3<<endl;
	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	for(int index=1; index<=8; index++) {
		cout<<"Total average fa_"<<index<<" = "<<fa_temps[index-1]<<endl;
		if(index>1) {
			cout<<"Clear Pad----------------------------------"<<endl;
			c_hstacks->Clear();
		}
		c_hstacks->cd();
		(ph_temps[index-1])->SetTitle(Form(" ; P^{+}_{%i}; Events", index));
		(ph_temps[index-1])->Draw("HIST");
		c_hstacks->SaveAs(Form("./plots/P+_%i.png", index));
		cout<<"Total average ga_"<<index<<" = "<<ga_temps[index-1]<<endl;
		cout<<"Clear Pad----------------------------------"<<endl;
		c_hstacks->Clear();
		c_hstacks->cd();
		(mh_temps[index-1])->SetTitle(Form(" ; P^{-}_{%i}; Events", index));
		(mh_temps[index-1])->Draw("HIST");
		c_hstacks->SaveAs(Form("./plots/P-_%i.png", index));
	}
	for(int index=1; index<=8; index++) {
		for(int index2=1; index2<=8; index2++) {
			cout<<"Total average hab_"<<index<<index2<<" = "<<hab_temps[index-1][index2-1]<<endl;
			cout<<"Clear Pad----------------------------------"<<endl;
			c_hstacks->Clear();
			c_hstacks->cd();
			(pmh2d_temps[index-1][index2-1])->SetTitle(Form(" ; P^{+}_{%i} ; P^{-}_{%i}", index, index2));
			(pmh2d_temps[index-1][index2-1])->Draw("COLZ");
			c_hstacks->SaveAs(Form("./plots/P_+%i_-%i.png", index, index2));
		}
	}




	// -------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------
	// +++++++ Draw the plots to show correaltion between the P functions and leptons' dynamics +++++++
	TH1F* temp_h1f_PvsD = 0;
	TH2F* temp_h2f_PvsD = 0;// temp_h2f_PvsD = ((TH2F*)(gDirectory->Get(Form("temp_h2f_PvsD"))));
	TH3F* temp_h3f_PvsD = 0;
	string LepDynamicsNmae = "";
	for(int k=0; k<objnum; k++){
		LepDynamicsNmae = (((objarr1->At(k)))->GetName());
		if((TString(LepDynamicsNmae)).Contains("J")) continue;
		if( !(TString(LepDynamicsNmae)).Contains("pt") && !(TString(LepDynamicsNmae)).Contains("eta") && !(TString(LepDynamicsNmae)).Contains("phi") ) continue;
		temp_h1f_PvsD = 0;
		temp_h2f_PvsD = 0;// temp_h2f_PvsD = ((TH2F*)(gDirectory->Get(Form("temp_h2f_PvsD"))));
		temp_h3f_PvsD = 0;
		for(int fid=0; fid<nfiles; fid++){ 
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 1);
			for(int index=1; index<=8; index++) {
				//chin.Draw(Form("%s>>ph_temps_%i", Form("P%i_plus_Boson1",index), index), Form("10000000*SF*(event_weight)*(P%i_plus_Boson1)",index), "goff");
				//chin.Draw(Form("%s>>mh_temps_%i", Form("P%i_plus_Boson2",index), index), Form("10000000*SF*(event_weight)*(P%i_plus_Boson2)",index), "goff");
				//chin.Draw(Form("%s>>ph_temps_%i", Form("P%i_plus_Boson1",index), index), Form("10000000*SF*(event_weight)*(P%i_tilde_Boson1)*0.5",index), "goff");
				//chin.Draw(Form("%s>>mh_temps_%i", Form("P%i_plus_Boson2",index), index), Form("10000000*SF*(event_weight)*(P%i_tilde_Boson2)*0.5",index), "goff");
				// ------------------------------------------------ 1D
				c_hstacks->Clear();
				c_hstacks->cd();
				chin.Draw(Form("%s>>temp_h1f_PvsD", Form("P%i_plus_Boson1",index)), "", "goff");
				temp_h1f_PvsD = ((TH1F*)(gDirectory->Get(Form("temp_h1f_PvsD"))));
				temp_h1f_PvsD->SetTitle(Form(" ; P^{+}_{%i}; Events", index));
				temp_h1f_PvsD->Draw("HIST");
				c_hstacks->SaveAs(Form("./plots/P+%i.png", index));
				temp_h1f_PvsD->Reset("ICES M");
				// ------------------------------------------------ 1D
				c_hstacks->Clear();
				c_hstacks->cd();
				chin.Draw(Form("%s>>temp_h1f_PvsD", Form("P%i_plus_Boson2",index)), "", "goff");
				temp_h1f_PvsD = ((TH1F*)(gDirectory->Get(Form("temp_h1f_PvsD"))));
				temp_h1f_PvsD->SetTitle(Form(" ; P^{-}_{%i}; Events", index));
				temp_h1f_PvsD->Draw("HIST");
				c_hstacks->SaveAs(Form("./plots/P-%i.png", index));
				temp_h1f_PvsD->Reset("ICES M");
				// ------------------------------------------------ 1D
				c_hstacks->Clear();
				c_hstacks->cd();
				chin.Draw(Form("%s>>temp_h1f_PvsD", Form("P%i_tilde_Boson1",index)), "", "goff");
				temp_h1f_PvsD = ((TH1F*)(gDirectory->Get(Form("temp_h1f_PvsD"))));
				temp_h1f_PvsD->SetTitle(Form(" ; P^{~(Boson1)}_{%i}; Events", index));
				temp_h1f_PvsD->Draw("HIST");
				c_hstacks->SaveAs(Form("./plots/P~Boson1%i.png", index));
				temp_h1f_PvsD->Reset("ICES M");
				// ------------------------------------------------ 1D
				c_hstacks->Clear();
				c_hstacks->cd();
				chin.Draw(Form("%s>>temp_h1f_PvsD", Form("P%i_tilde_Boson2",index)), "", "goff");
				temp_h1f_PvsD = ((TH1F*)(gDirectory->Get(Form("temp_h1f_PvsD"))));
				temp_h1f_PvsD->SetTitle(Form(" ; P^{~(Boson2)}_{%i}; Events", index));
				temp_h1f_PvsD->Draw("HIST");
				c_hstacks->SaveAs(Form("./plots/P~Boson2%i.png", index));
				temp_h1f_PvsD->Reset("ICES M");
				// ------------------------------------------------ 2D
				c_hstacks->Clear();
				c_hstacks->cd();
				chin.Draw(Form("%s:%s>>temp_h2f_PvsD", LepDynamicsNmae.c_str(), Form("P%i_plus_Boson1",index)), "", "goff");
				temp_h2f_PvsD = ((TH2F*)(gDirectory->Get(Form("temp_h2f_PvsD"))));
				temp_h2f_PvsD->SetTitle(Form(" ; P^{+}_{%i}; %s; Events", index, LepDynamicsNmae.c_str()));
				temp_h2f_PvsD->Draw("COLZ");
				c_hstacks->SaveAs(Form("./plots/P+%i_VS_%s.png", index, LepDynamicsNmae.c_str()));
				temp_h2f_PvsD->Reset("ICES M");
				// ------------------------------------------------ 2D
				c_hstacks->Clear();
				c_hstacks->cd();
				chin.Draw(Form("%s:%s>>temp_h2f_PvsD", LepDynamicsNmae.c_str(), Form("P%i_plus_Boson2",index)), "", "goff");
				temp_h2f_PvsD = ((TH2F*)(gDirectory->Get(Form("temp_h2f_PvsD"))));
				temp_h2f_PvsD->SetTitle(Form(" ; P^{-}_{%i}; %s; Events", index, LepDynamicsNmae.c_str()));
				temp_h2f_PvsD->Draw("COLZ");
				c_hstacks->SaveAs(Form("./plots/P-%i_VS_%s.png", index, LepDynamicsNmae.c_str()));
				temp_h2f_PvsD->Reset("ICES M");
				// ------------------------------------------------ 2D
				c_hstacks->Clear();
				c_hstacks->cd();
				chin.Draw(Form("%s:%s>>temp_h2f_PvsD", LepDynamicsNmae.c_str(), Form("P%i_tilde_Boson1",index)), "", "goff");
				temp_h2f_PvsD = ((TH2F*)(gDirectory->Get(Form("temp_h2f_PvsD"))));
				temp_h2f_PvsD->SetTitle(Form(" ; P^{~(Boson1)}_{%i}; %s; Events", index, LepDynamicsNmae.c_str()));
				temp_h2f_PvsD->Draw("COLZ");
				c_hstacks->SaveAs(Form("./plots/P~Boson1%i_VS_%s.png", index, LepDynamicsNmae.c_str()));
				temp_h2f_PvsD->Reset("ICES M");
				// ------------------------------------------------ 2D
				c_hstacks->Clear();
				c_hstacks->cd();
				chin.Draw(Form("%s:%s>>temp_h2f_PvsD", LepDynamicsNmae.c_str(), Form("P%i_tilde_Boson2",index)), "", "goff");
				temp_h2f_PvsD = ((TH2F*)(gDirectory->Get(Form("temp_h2f_PvsD"))));
				temp_h2f_PvsD->SetTitle(Form(" ; P^{~(Boson2)}_{%i}; %s; Events", index, LepDynamicsNmae.c_str()));
				temp_h2f_PvsD->Draw("COLZ");
				c_hstacks->SaveAs(Form("./plots/P~Boson2%i_VS_%s.png", index, LepDynamicsNmae.c_str()));
				temp_h2f_PvsD->Reset("ICES M");
				for(int index2=1; index2<=8; index2++) {
					//chin.Draw(Form("%s:%s>>pmh2d_temps_%i_%i", Form("P%i_plus_Boson2",index2), Form("P%i_plus_Boson1",index), index, index2), Form("10000000*SF*(event_weight)*(P%i_plus_Boson1)*(P%i_plus_Boson2)",index,index2), "goff"); // %s:%s is act as y:x, which means P%i_plus_Boson1 in X axis and P%i_plus_Boson2 in Yaxis
					//chin.Draw(Form("%s:%s>>pmh2d_temps_%i_%i", Form("P%i_tilde_Boson1",index2), Form("P%i_tilde_Boson1",index), index, index2), Form("10000000*SF*(event_weight)*(P%i_tilde_Boson1)*(P%i_tilde_Boson1)*0.25",index,index2), "goff"); // %s:%s is act as y:x, which means P%i_plus_Boson1 in X axis and P%i_plus_Boson2 in Yaxis
					// ------------------------------------------------ 3D
					c_hstacks->Clear();
					c_hstacks->cd();
					chin.Draw(Form("%s:%s:%s>>temp_h3f_PvsD", LepDynamicsNmae.c_str(), Form("P%i_tilde_Boson2",index2), Form("P%i_tilde_Boson1",index)), "", "goff"); // %s:%s is act as y:x, which means P%i_plus_Boson1 in X axis and P%i_plus_Boson2 in Yaxis
					temp_h3f_PvsD = ((TH3F*)(gDirectory->Get(Form("temp_h3f_PvsD"))));
					temp_h3f_PvsD->SetTitle(Form(" ; P^{~Boson2}_{%i}; P^{~Boson1}_{%i}; %s", index, index2, LepDynamicsNmae.c_str()));
					temp_h3f_PvsD->Draw("COLZ");
					c_hstacks->SaveAs(Form("./plots/P~Boson1%i_VS_P~Boson2%i_VS_%s.png", index, index2, LepDynamicsNmae.c_str()));
					temp_h3f_PvsD->Reset("ICES M");
					// ------------------------------------------------ 3D
					c_hstacks->Clear();
					c_hstacks->cd();
					chin.Draw(Form("%s:%s:%s>>temp_h3f_PvsD", LepDynamicsNmae.c_str(), Form("P%i_plus_Boson2",index2), Form("P%i_plus_Boson1",index)), "", "goff"); // %s:%s is act as y:x, which means P%i_plus_Boson1 in X axis and P%i_plus_Boson2 in Yaxis
					temp_h3f_PvsD = ((TH3F*)(gDirectory->Get(Form("temp_h3f_PvsD"))));
					temp_h3f_PvsD->SetTitle(Form(" ; P^{+}_{%i}; P^{-}_{%i}; %s", index, index2, LepDynamicsNmae.c_str()));
					temp_h3f_PvsD->Draw("COLZ");
					c_hstacks->SaveAs(Form("./plots/P+%i_VS_P-%i_VS_%s.png", index, index2, LepDynamicsNmae.c_str()));
					temp_h3f_PvsD->Reset("ICES M");
				}
			}
		}
	}
	// Using GetSumOfWeights()
//---	cout<<"##################################### Using GetSumOfWeights() #####################################"<<endl;
//---	sum_event_weight = sum_event_weight_temps->GetSumOfWeights();
//---	cout<<"sum_event_weight = "<<sum_event_weight<<endl;
//---	for(int index=1; index<=8; index++) {
//---		fa_temps[index-1] = ((ph_temps[index-1])->GetSumOfWeights())/sum_event_weight;
//---		ga_temps[index-1] = ((mh_temps[index-1])->GetSumOfWeights())/sum_event_weight;
//---		for(int index2=1; index2<=8; index2++) {
//---			hab_temps[index-1][index2-1] = ((pmh2d_temps[index-1][index2-1])->GetSumOfWeights())/sum_event_weight;
//---		}
//---	}
//---	cout<<"===================================== QE elements ====================================="<<endl;
//---	for(int index=1; index<=8; index++) {
//---		cout<<"fa_"<<index<<" = "<<fa_temps[index-1]<<endl;
//---		cout<<"ga_"<<index<<" = "<<ga_temps[index-1]<<endl;
//---	}
//---	for(int index=1; index<=8; index++) {
//---		for(int index2=1; index2<=8; index2++) {
//---			cout<<"hab_"<<index<<index2<<" = "<<hab_temps[index-1][index2-1]<<endl;
//---		}
//---	}
//---	SUM_SQ_fa = 0;
//---	SUM_SQ_ga = 0;
//---	SUM_SQ_hab = 0;
//---	for(int index=1; index<=8; index++) {
//---		SUM_SQ_fa += (fa_temps[index-1] * fa_temps[index-1]);
//---		SUM_SQ_ga += (ga_temps[index-1] * ga_temps[index-1]);
//---		for(int index2=1; index2<=8; index2++) {
//---			SUM_SQ_hab += (hab_temps[index-1][index2-1] * hab_temps[index-1][index2-1]);
//---		}
//---	}
//---	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Output evaluations of Quantum entanglement witness observables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
//---	L_2 = 2 * TMath::Max( ( (-2.0/9.0) - (12.0*SUM_SQ_fa) + (6.0*SUM_SQ_ga) + (4.0*SUM_SQ_hab) ) , ( (-2.0/9.0) - (12.0*SUM_SQ_ga) + (6.0*SUM_SQ_fa) + (4.0*SUM_SQ_hab) ) );
//---	L_2 = TMath::Max( L_2 , 0.0 );
//---	I_3 = 4.0*(hab_temps[4-1][4-1] + hab_temps[5-1][5-1]) - (4.0/sqrt(3))*(hab_temps[6-1][1-1] + hab_temps[6-1][6-1] + hab_temps[7-1][2-1] + hab_temps[7-1][7-1] + hab_temps[1-1][1-1] + hab_temps[1-1][6-1] + hab_temps[2-1][2-1] + hab_temps[2-1][7-1]);
//---	cout<<"   L_2 = "<<L_2<<endl;
//---	cout<<"   I_3 = "<<I_3<<endl;
//---	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	// -------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------





///---	cout<<"================================== Drawing variables distributions =================================="<<endl;
///---	for(int k=0; k<objnum; k++){
///---		if(k>0) {
///---			cout<<"Clear Pad----------------------------------"<<endl;
///---			c_hstacks->Clear();
///---		}
///---		c_hstacks->cd();
///---		c_hstacks->SetRightMargin(0.25);
///---		TLegend lg(0.75,0.70,0.99,0.99);
///---		//lg.SetNColumns(2);
///---		cout<<"Draw "<<k<<" Hists--------------------------"<<endl;
///---		cout<<"Draw Hists---------------------------------- Branch Name: "<< ((objarr1->At(k)))->GetName() <<endl;
///---
///---		std::vector<TH1D*> smh_temps;
///---		xmin = 0;
///---		xmax = 100;
///---		xrange = xmax - xmin;
///---		width_step = xrange/ndivisions;
///---		signalhistmaxy = 0;
///---		for(int fid=0; fid<nfiles; fid++){ 
///---			if(sample_types.at(fid) == 1)continue;
///---			TChain chin("tree");
///---			chin.Add(filesin.at(fid));
///---			chin.SetBranchStatus("*", 1);
///---			chin.SetBranchStatus("SF", 1);
///---			if( (string(((objarr1->At(k)))->GetName()) != string("event_weight")) && (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) ) {
///---				chin.Draw(Form("%s>>opt_temp%i", (objarr1->At(k))->GetName(), k), Form("10000000*SF*(event_weight)*((%s > -99) && %s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
///---				xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge(1);
///---				xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetNcells() );
///---				signalhistmaxy += ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetMaximum();
///---				xrange = xmax - xmin;
///---				width_step = xrange/ndivisions;
///---				xmin -= 2.0*width_step;
///---				xmax += 2.0*width_step;
///---			}
///---			else if( (string(((objarr1->At(k)))->GetName()) == string("event_weight")) ) {
///---				chin.Draw(Form("%s>>opt_temp%i", (objarr1->At(k))->GetName(), k), Form("10000000*SF*(%s)",cut.c_str()), "goff");
///---				xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge(1);
///---				xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetNcells() );
///---				signalhistmaxy += ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetMaximum();
///---				xrange = xmax - xmin;
///---				width_step = xrange/ndivisions;
///---				xmin -= 2.0*width_step;
///---				xmax += 2.0*width_step;
///---			}
///---			else if( (string(((objarr1->At(k)))->GetName()) == string("lep_fst")) ) {
///---				xmin = 0;
///---				xmax = 3;
///---				signalhistmaxy += (chin.GetEntries()/2);
///---				xrange = xmax - xmin;
///---				width_step = xrange/ndivisions;
///---			}
///---		}
///---		cout<<"======================================================================== Decide histos ranges"<<endl;
///---
///---		//-------------------------------------- histograms for SM samples
///---		{
///---			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
///---			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
///---			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
///---			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
///---			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
///---			cout<<"histo vectors prepared"<<endl;
///---		}
///---
///---
///---		//============================================================================== Fill histograms
///---		for(int fid=0; fid<nfiles; fid++){ 
///---			TChain Chain("tree");
///---			Chain.Add(filesin.at(fid));
///---			Chain.SetBranchStatus("*", 1);
///---			Chain.SetBranchStatus("SF", 1);
///---			//------------------------------------------------- SM signal and backgrounds
///---			if( !((filesin.at(fid)).Contains("reweightscan")) && !((filesin.at(fid)).Contains("benchmark")) ) {
///---				if(string(((objarr1->At(k)))->GetName()) == string("lep_fst")) {
///---					if((filesin.at(fid)).Contains("QE_signal")) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
///---					}
///---					else if((filesin.at(fid)).Contains("mumuTomumu")) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
///---					}
///---					else if((filesin.at(fid)).Contains("mumuTovmvm")) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
///---					}
///---					else if((filesin.at(fid)).Contains("ttTo")) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
///---					}
///---					else if((filesin.at(fid)).Contains("mumuTow") || (filesin.at(fid)).Contains("mumuToz") || (filesin.at(fid)).Contains("mumuToh") ) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
///---					}
///---				}
///---				if(string(((objarr1->At(k)))->GetName()) != string("lep_fst")) {
///---					if((filesin.at(fid)).Contains("QE_signal")) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
///---					}
///---					else if((filesin.at(fid)).Contains("mumuTomumu")) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
///---					}
///---					else if((filesin.at(fid)).Contains("mumuTovmvm")) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
///---					}
///---					else if((filesin.at(fid)).Contains("ttTo")) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
///---					}
///---					else if((filesin.at(fid)).Contains("mumuTow") || (filesin.at(fid)).Contains("mumuToz") || (filesin.at(fid)).Contains("mumuToh") ) {
///---						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
///---					}
///---				}
///---			}
///---		}
///---
///---		cout<<"======================================================================== Preapare hist plots"<<endl;
///---
///---
///---		c_hstacks->SetTopMargin(3.0);
///---		c_hstacks->cd();
///---
///---		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
///---			if(((smh_temps.at(hindex))->Integral() > 1) && string((smh_temps.at(hindex))->GetName())==string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
///---				string process((smh_temps.at(hindex))->GetName());
///---				process.erase(process.find("th_",0), string("th_").size());
///---				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
///---				cout<<"-------------------------------------------- Draw SM signal["<<process<<"] plots"<<endl;
///---				cout<<"Integral = "<<(smh_temps.at(hindex))->Integral()<<endl;
///---				cout<<"Mean = "<<(smh_temps.at(hindex))->GetMean()<<endl;
///---				cout<<"RMS = "<<(smh_temps.at(hindex))->GetRMS()<<endl;
///---				(smh_temps.at(hindex))->SetTitle(Form(" ; %s; ",((objarr1->At(k)))->GetName()));
///---				(smh_temps.at(hindex))->SetTitle(Form(" ; %s; ",((objarr1->At(k)))->GetName()));
///---				(smh_temps.at(hindex))->SetTitleOffset(1.2);
///---				(smh_temps.at(hindex))->SetLineWidth(2*3);
///---				(smh_temps.at(hindex))->SetLineStyle(1);
///---				(smh_temps.at(hindex))->SetMarkerStyle(20);
///---				(smh_temps.at(hindex))->SetMarkerSize(0.9*3);
///---				(smh_temps.at(hindex))->SetLineColor(kBlack);
///---				(smh_temps.at(hindex))->SetMarkerColor(kBlack);
///---				//(smh_temps.at(hindex))->SetFillColorAlpha(kBlack, 10);
///---				//(smh_temps.at(hindex))->SetFillStyle(3001);
///---				if(hindex==0) (smh_temps.at(hindex))->Draw("HIST");
///---				else (smh_temps.at(hindex))->Draw("HIST SAME");
///---				(smh_temps.at(hindex))->SetMaximum( signalhistmaxy*10.0 );
///---				lg.AddEntry(smh_temps.at(hindex), process.c_str(), "L");
///---			}
///---		}
///---		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
///---			if(((smh_temps.at(hindex))->Integral() > 1) && string((smh_temps.at(hindex))->GetName())!=string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
///---				string process((smh_temps.at(hindex))->GetName());
///---				process.erase(process.find("th_",0), string("th_").size());
///---				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
///---				cout<<"-------------------------------------------- Draw SM background["<<process<<"] plots"<<endl;
///---				cout<<"Integral = "<<(smh_temps.at(hindex))->Integral()<<endl;
///---				cout<<"Mean = "<<(smh_temps.at(hindex))->GetMean()<<endl;
///---				cout<<"RMS = "<<(smh_temps.at(hindex))->GetRMS()<<endl;
///---				(smh_temps.at(hindex))->SetTitle(Form(" ; %s; ",((objarr1->At(k)))->GetName()));
///---				(smh_temps.at(hindex))->SetTitle(Form(" ; %s; ",((objarr1->At(k)))->GetName()));
///---				(smh_temps.at(hindex))->SetTitleOffset(1.2);
///---				(smh_temps.at(hindex))->SetLineWidth(2*3);
///---				(smh_temps.at(hindex))->SetLineStyle(hindex);
///---				(smh_temps.at(hindex))->SetMarkerStyle(20+hindex);
///---				(smh_temps.at(hindex))->SetMarkerSize(0.9*3);
///---				(smh_temps.at(hindex))->SetLineColor(1+hindex);
///---				(smh_temps.at(hindex))->SetMarkerColor(1+hindex);
///---				//(smh_temps.at(hindex))->SetFillColorAlpha(1+hindex, 10);
///---				//(smh_temps.at(hindex))->SetFillStyle(3001);
///---				if(hindex==0) (smh_temps.at(hindex))->Draw("HIST");
///---				else (smh_temps.at(hindex))->Draw("HIST SAME");
///---				(smh_temps.at(hindex))->SetMaximum( signalhistmaxy*10.0 );
///---				lg.AddEntry(smh_temps.at(hindex), process.c_str(), "LF");
///---			}
///---		}
///---		lg.Draw("SAME");
///---
///---		cout<<"======================================================================== Draw histos"<<endl;
///---		//c_hstacks->SetLogy();
///---		c_hstacks->Update();
///---		cout<<"Save drawing----------------------------------"<<endl;
///---		if( access(Form("./plots/"),0) ) mkdir(Form("./plots/"),0755);
///---		if( access(Form("./plots/%s/", para_name.c_str()),0) ) mkdir(Form("./plots/%s/", para_name.c_str()),0755);
///---		c_hstacks->SaveAs(Form("./plots/%s/%s_%s_comparison_%s.png", para_name.c_str(), para_name.c_str(), plot_prefix.c_str(), ((objarr1->At(k)))->GetName()));
///---		c_hstacks->SaveAs(Form("./plots/%s/%s_%s_comparison_%s.pdf", para_name.c_str(), para_name.c_str(), plot_prefix.c_str(), ((objarr1->At(k)))->GetName()));
///---
///---
///---		cout<<"======================================================================== Prepare histos stacks"<<endl;
///---		cout<<"Clear pad for hstack----------------------------------"<<endl;
///---		c_hstacks->Clear();
///---		c_hstacks->cd();
///---		c_hstacks->SetFrameFillColor(19);
///---		c_hstacks->SetTopMargin(3.0);
///---		c_hstacks->SetRightMargin(0.25);
///---		THStack *hs_sm = new THStack("hs_sm",Form(" ; %s; ",((objarr1->At(k)))->GetName())); 
///---		TLegend lg_sm(0.75,0.80,0.99,0.99,"SM");
///---		// Setup SM signal legend
///---		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
///---			if(((smh_temps.at(hindex))->Integral() > 1) && ((smh_temps.at(hindex))->Integral() > 1) && string((smh_temps.at(hindex))->GetName())==string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
///---				(smh_temps.at(hindex))->SetLineWidth(2*3);
///---				(smh_temps.at(hindex))->SetLineStyle(1);
///---				(smh_temps.at(hindex))->SetMarkerStyle(0);
///---				(smh_temps.at(hindex))->SetMarkerSize(0);
///---				(smh_temps.at(hindex))->SetLineColor(kBlack);
///---				(smh_temps.at(hindex))->SetMarkerColor(0);
///---				//(smh_temps.at(hindex))->SetFillColorAlpha(0, 10);
///---				//(smh_temps.at(hindex))->SetFillStyle(3001);
///---				////// Only for legen, not for stack				hs_sm->Add((smh_temps.at(hindex)));
///---				string process((smh_temps.at(hindex))->GetName());
///---				process.erase(process.find("th_",0), string("th_").size());
///---				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
///---				lg_sm.AddEntry(smh_temps.at(hindex), process.c_str(), "L");
///---			}
///---		}
///---		// Setup SM background legend
///---		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
///---			if(((smh_temps.at(hindex))->Integral() > 1) && string((smh_temps.at(hindex))->GetName())!=string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
///---				(smh_temps.at(hindex))->SetLineWidth(1*3);
///---				(smh_temps.at(hindex))->SetLineStyle(1);
///---				(smh_temps.at(hindex))->SetMarkerStyle(0);
///---				(smh_temps.at(hindex))->SetMarkerSize(0);
///---				(smh_temps.at(hindex))->SetLineColor(1+hindex);
///---				(smh_temps.at(hindex))->SetMarkerColor(0);
///---				(smh_temps.at(hindex))->SetFillColorAlpha(1+hindex, 10);
///---				(smh_temps.at(hindex))->SetFillStyle(3001);
///---				////// Only for legen, not for stack				hs_sm->Add((smh_temps.at(hindex)));
///---				string process((smh_temps.at(hindex))->GetName());
///---				process.erase(process.find("th_",0), string("th_").size());
///---				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
///---				lg_sm.AddEntry(smh_temps.at(hindex), process.c_str(), "LF");
///---			}
///---		}
///---		// Add SM background into stack
///---		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
///---			if(((smh_temps.at(hindex))->Integral() > 1) && string((smh_temps.at(hindex))->GetName())!=string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
///---				hs_sm->Add((smh_temps.at(hindex)));
///---			}
///---		}
///---		// Add SM signal into stack
///---		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
///---			if(((smh_temps.at(hindex))->Integral() > 1) && string((smh_temps.at(hindex))->GetName())==string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
///---				hs_sm->Add((smh_temps.at(hindex)));
///---			}
///---		}
///---		cout<<"======================================================================== Draw histos stacks"<<endl;
///---		hs_sm->Draw("HIST");
///---		lg_sm.Draw("SAME");
///---		//c_hstacks->SetLogy();
///---		c_hstacks->Update();
///---		cout<<"Save drawing----------------------------------"<<endl;
///---		if( access(Form("./plots/"),0) ) mkdir(Form("./plots/"),0755);
///---		if( access(Form("./plots/%s/", para_name.c_str()),0) ) mkdir(Form("./plots/%s/", para_name.c_str()),0755);
///---		c_hstacks->SaveAs(Form("./plots/%s/%s_%s_stacks_%s.png", para_name.c_str(), para_name.c_str(), plot_prefix.c_str(), ((objarr1->At(k)))->GetName()));
///---		c_hstacks->SaveAs(Form("./plots/%s/%s_%s_stacks_%s.pdf", para_name.c_str(), para_name.c_str(), plot_prefix.c_str(), ((objarr1->At(k)))->GetName()));
///---
///---
///---
///---	}
}



std::pair<std::vector<TH1D*>, std::vector<int>> obtain_histograms(std::vector<TString> filesin, std::vector<bool> if_event_exist, std::vector<int> sample_types, string plot_prefix, string cut, int ndivisions, string para_name, double RWstart, double RWstep, int RWnstep, string Xvar="deltaM") {

	std::vector<TH1D*> output_histograms;
	std::vector<int> output_histograms_types;
	std::pair<std::vector<TH1D*>, std::vector<int>> output;


	assert(!(filesin.empty()));
	gStyle->SetOptStat(0000);
	for(int ifi=0; ifi<if_event_exist.size(); ifi++) {
		if( !(if_event_exist.at(ifi)) ) {
			filesin.erase(filesin.begin()+ifi);
			sample_types.erase(sample_types.begin()+ifi);
		}
	}
	int nfiles = filesin.size();
	TFile fin1(filesin.at(0), "READ");
	TTree* tin1 = (TTree*)fin1.Get("tree");
	TObjArray* objarr = tin1->GetListOfBranches();
	TObjArray* objarr1 = (TObjArray*)objarr->Clone();
	int objnum = objarr1->GetEntries();
	delete tin1;
	fin1.Close();

	TCanvas* c_hstacks = new TCanvas("c_hstacks","c_hstacks",2500,1900);
	double xmin=0;
	double xmax=100;
	double xrange=0;
	double width_step=0;
	double signalhistmaxy = 0;
	for(int k=0; k<objnum; k++){
		if(string(((objarr1->At(k)))->GetName()) != string(Xvar)) continue;
		if(k>0) {
			cout<<"Clear Pad----------------------------------"<<endl;
			c_hstacks->Clear();
		}
		c_hstacks->cd();
		c_hstacks->SetRightMargin(0.25);
		TLegend lg(0.75,0.70,0.99,0.99);
		//lg.SetNColumns(2);
		cout<<"Draw "<<k<<" Hists--------------------------"<<endl;
		cout<<"Draw Hists---------------------------------- Branch Name: "<< ((objarr1->At(k)))->GetName() <<endl;

		std::vector<TH1D*> ewh_temps;
		std::vector<TH1D*> rwh_temps;
		std::vector<TH1D*> bmh_temps;
		std::vector<TH1D*> smh_temps;
		std::vector<TH1D*> rw_sub_sm_h_temps;
		std::map<int, string> benchmarks_map;
		std::map<int, double> benchmarks_values_map;
		xmin = 0;
		xmax = 100;
		xrange = xmax - xmin;
		width_step = xrange/ndivisions;
		signalhistmaxy = 0;
		for(int fid=0; fid<nfiles; fid++){ 
			if(sample_types.at(fid) == 1)continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 1);
			chin.SetBranchStatus("SF", 1);
			if( (string(((objarr1->At(k)))->GetName()) != string("weights")) && (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) ) {
				chin.Draw(Form("%s>>opt_temp%i", (objarr1->At(k))->GetName(), k), Form("10000000*SF*(event_weight)*((%s > -99) && %s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
				xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge(1);
				xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetNcells() );
				signalhistmaxy += ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetMaximum();
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
			else if( (string(((objarr1->At(k)))->GetName()) == string("weights")) && (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) ) {
				chin.Draw(Form("%s>>opt_temp%i", (objarr1->At(k))->GetName(), k), Form("10000000*SF*(event_weight)*((%s > -99) && %s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
				xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge(1);
				xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%i", k))))->GetNcells() );
				signalhistmaxy += (chin.GetEntries());
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
			else if( (string(((objarr1->At(k)))->GetName()) == string("lep_fst")) ) {
				xmin = 0;
				xmax = 100;
				signalhistmaxy += (chin.GetEntries()/2);
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
		}
		cout<<"======================================================================== Decide histos ranges"<<endl;

		//-------------------------------------- histograms for SM samples
		{
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			cout<<"histo vectors prepared"<<endl;
		}


		//============================================================================== Fill histograms
		for(int fid=0; fid<nfiles; fid++){ 
			TChain Chain("tree");
			Chain.Add(filesin.at(fid));
			Chain.SetBranchStatus("*", 1);
			Chain.SetBranchStatus("SF", 1);
			//------------------------------------------------- SM signal and backgrounds
			if( !((filesin.at(fid)).Contains("reweightscan")) && !((filesin.at(fid)).Contains("benchmark")) ) {
				if(string(((objarr1->At(k)))->GetName()) == string("lep_fst")) {
					if((filesin.at(fid)).Contains("mumuTozzzTollllvlvl")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTomumu")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTovmvm")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("ttTo")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTow") || (filesin.at(fid)).Contains("mumuToz") || (filesin.at(fid)).Contains("mumuToh") ) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
				}
				if(string(((objarr1->At(k)))->GetName()) != string("lep_fst")) {
					if((filesin.at(fid)).Contains("mumuTozzzTollllvlvl")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTomumu")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTovmvm")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("ttTo")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTow") || (filesin.at(fid)).Contains("mumuToz") || (filesin.at(fid)).Contains("mumuToh") ) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
				}
			}
		}

		for(int hi=0; hi<smh_temps.size(); hi++){
			output_histograms.push_back(smh_temps[hi]);
			if((TString((output_histograms.back())->GetName())).Contains("signal")) output_histograms_types.push_back(0);
			else output_histograms_types.push_back(1);
		}

		output.first = output_histograms;
		output.second = output_histograms_types;

	}

	cout<<"===== Summary of output histograms obtained ====="<<endl;
	for(int hi=0; hi<(output.first).size(); hi++){
		cout<<"Histogram name: "<<(output.first)[hi]->GetName()<<endl;
		cout<<"Histogram integral: "<<(output.first)[hi]->Integral(1, (output.first)[hi]->GetNcells()-2)<<endl;
	}
	cout<<"===== Summary end ==============================="<<endl;

	return output;
}







void selection_draw(string para_name, double RWstart, double RWstep, int RWnstep, string hardcut="( (Mrecoil>0 && Mrecoil<200) && (Mll1>80 && Mll1<100) && (Mll2>80 && Mll2<100) && (deltaM<20) )", double CLv=0.95, int TStype=3){

	gSystem->Load("/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/ExRootAnalysis/libExRootAnalysis.so");
	gSystem->Load("libPhysics");

	gROOT->SetBatch();


	std::vector<TString> Infiles_names;
	std::vector<TString> Outfiles_names_preselected;
	std::vector<bool>    If_Event_exist;
	std::vector<int>    Sample_Types; // 0 for SM signal, 1 for SM background, 2 for AQGC signal

	//--------------------------------------------------------------------------------------------------------------
	Infiles_names.push_back(TString(Form("QE_mumuToZZTollll_SM_1TeV.root"))); Sample_Types.push_back(0);


	//--------------------------------------------------------------------------------------------------------------
	Outfiles_names_preselected.push_back(TString(Form("output_QE_mumuToZZTollll_SM_1TeV.root")));


	if(Infiles_names.size() != Outfiles_names_preselected.size()) {
		cout<<"Files number does not match! Exit!"<<endl;
		return;
	}

	for(size_t ifile=0; ifile<Infiles_names.size(); ifile++) {
		TChain* chain = new TChain("LHEF");
		chain->Add(Infiles_names.at(ifile));
		ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
		Long64_t numberOfEntries = treeReader->GetEntries();
		If_Event_exist.push_back(AnalyseEvents(treeReader, Outfiles_names_preselected.at(ifile), RWnstep));

	}

	string cut="";
	//std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> cut_information;
	//cut_information = OptimizeCuts(Outfiles_names_preselected, Sample_Types, cut);
	cut="(1==1)";

	//cout<<"Computing significance: ----------------------------"<<endl;
	//std::vector<int> Sig_Types;
	//std::vector<int> Bkg_Types;
	//double significance_aftersel = 0;
	//cout<<"++++++++++++++++++++++++++++++++++++++++"<<endl;
	//Sig_Types.push_back(0);
	//Bkg_Types.push_back(1);
	//cout<<"Signal types: ";
	//for(int k=0; k<Sig_Types.size(); k++) {
	//	cout<<Sig_Types.at(k)<<"; ";
	//}
	//cout<<endl;
	//cout<<"Bkg types: ";
	//for(int k=0; k<Bkg_Types.size(); k++) {
	//	cout<<Bkg_Types.at(k)<<"; ";
	//}
	//significance_aftersel = ComputeSignificance(Outfiles_names_preselected, Sample_Types, Sig_Types, Bkg_Types, cut);

	//cout<<endl;
	//cout<<"Signal significance S/sqrt(B) is "<<significance_aftersel<<endl;
	//cout<<"----------------------------------------"<<endl;
	//cout<<"=================================================================="<<endl;
	//cout<<"Optimized cut-based selection: "<<cut<<endl;
	//cout<<"Significance after selection: "<<significance_aftersel<<endl;
	//cout<<"=================================================================="<<endl;

	PlotEvents(Outfiles_names_preselected, If_Event_exist, Sample_Types, "bfcut", "(1==1)", 20, para_name, RWstart, RWstep,  RWnstep);
	//PlotEvents(Outfiles_names_preselected, If_Event_exist, Sample_Types, "afcut", cut, 20, para_name, RWstart, RWstep,  RWnstep);


}
