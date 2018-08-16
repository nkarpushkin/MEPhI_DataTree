#include <fstream>
#include <iostream>
#include <TNtuple.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <Rtypes.h>
#include <vector>
#include <TBrowser.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMath.h>
#include <TKey.h>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include <TSystemFile.h>
#include <TSystemDirectory.h>


#include "/home/nikolay/MEPhI/DataTree_code_and_files/DataTreeQA/cuts/DataTreeCutsConfig.h"
#include "/home/nikolay/MEPhI/DataTree_code_and_files/DataTreeQA/cuts/DataTreeCuts.h"



void DataTree_chain_reader_2()
{


	TFile *QAConfig = TFile::Open("/home/nikolay/MEPhI/DataTree_code_and_files/DataTreeQA/macro/QAConfigurations.root");
	cuts::DataTreeCutsConfig* cc = (cuts::DataTreeCutsConfig*) QAConfig->Get("pbpb_30agev_data_config_cuts");
	cuts::DataTreeCuts _cuts(cc);



	TString source_path = "/home/nikolay/MEPhI/DataTree_code_and_files/20180625_221329/";
	TObjArray files_names;

	TChain *data_chain = new TChain("DataTree");

	TSystemDirectory dir(source_path, source_path);
 	TList *files = dir.GetListOfFiles();
 	if (files) 
	{
		TSystemFile *file; 
		TString fname; 
		TIter next(files); 
		while ((file=(TSystemFile*)next())) 
		{ 
			fname = file->GetName(); 
			if (!file->IsDirectory() && fname.EndsWith("DataTree.root")) 
			//if ((!file->IsDirectory()) && (fname.EndsWith("DataTree.root")) && (files_names.GetLast() < 40))
			{
				files_names.Add(new TObjString(fname));
				data_chain->AddFile(  (source_path + fname).Data() );
			}
		}
	}
	Int_t total_files = files_names.GetLast()+1;
	printf("Total files: %i\n", total_files);
	
	DataTreeEvent* event = new DataTreeEvent;
	data_chain->SetBranchAddress("DTEvent", &event);


	TH1 *th1_hist_ptr = NULL;
	TH2 *th2_hist_ptr = NULL;
	th1_hist_ptr = new TH1F("Vertex_tracks","Vertex_tracks", 100, -5, 5);

	th1_hist_ptr = new TH1F("PSD_EP","PSD_EP", 500, -4, 4);

	th2_hist_ptr = new TH2F("Qvectors_TPC","Qvectors_TPC", 70, -4, 4, 70, -4, 4);
	th2_hist_ptr = new TH2F("Qvectors_PSD","Qvectors_PSD", 70, -4, 4, 70, -4, 4);

	TProfile *tprof_ptr = NULL;
	tprof_ptr = new TProfile("Delta_psi","Delta_psi", 100, 0, 650);

	for(Int_t entry = 0; entry <= data_chain->GetEntries(); entry++)
	{
		data_chain->GetEntry(entry);

		//cuts			
		if(_cuts.IsGoodEvent(*event))
		{

			//Q-vector

			Double_t Q1x_TPC=0;
			Double_t Q1y_TPC=0;

			Double_t Q2x_TPC=0;
			Double_t Q2y_TPC=0;

			for(Int_t track_iter = 0; track_iter < event->GetNVertexTracks(); track_iter++)
			{
				th1_hist_ptr = ((TH1*)(gDirectory->FindObjectAny( "Vertex_tracks" )));
				th1_hist_ptr->Fill( event->GetVertexTrack(track_iter)->GetPhi() );


				//Q-vector
				Q1x_TPC += event->GetVertexTrack(track_iter)->GetPx();
				Q1y_TPC += event->GetVertexTrack(track_iter)->GetPy();
				//Q2x_TPC += (event->GetVertexTrack(track_iter)->GetPt()) * cos( 2*(event->GetVertexTrack(track_iter)->GetPhi()) );
				//Q2y_TPC += (event->GetVertexTrack(track_iter)->GetPt()) * sin( 2*(event->GetVertexTrack(track_iter)->GetPhi()) );
				Q2x_TPC += (event->GetVertexTrack(track_iter)->GetPt()) * (cos( (event->GetVertexTrack(track_iter)->GetPhi()) ) );
				Q2y_TPC += (event->GetVertexTrack(track_iter)->GetPt()) * (sin( (event->GetVertexTrack(track_iter)->GetPhi()) ) );
			}

			Double_t psi1_TPC_EP = atan2(Q1y_TPC,Q1x_TPC);
			Double_t psi2_TPC_EP = atan2(Q2y_TPC,Q2x_TPC);
			th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "Qvectors_TPC" )));
			th2_hist_ptr->Fill( psi1_TPC_EP,psi2_TPC_EP );

			//Q-vector PSD

			Double_t Q1x_PSD=0;
			Double_t Q1y_PSD=0;

			Double_t Q2x_PSD=0;
			Double_t Q2y_PSD=0;
			
			for(Int_t module_iter = 0; module_iter < event->GetNPSDModules(); module_iter++)
			{
				if((event->GetPSDModule(module_iter)->GetId()) <= 45)
				{
					Q1x_PSD +=  (event->GetPSDModule(module_iter)->GetEnergy()) * (event->GetPSDModule(module_iter)->GetPositionComponent(0)) ;
					Q1y_PSD +=  (event->GetPSDModule(module_iter)->GetEnergy()) * (event->GetPSDModule(module_iter)->GetPositionComponent(1)) ;

					//cout<<(event->GetPSDModule(module_iter)->GetPositionComponent(2))<<" "<<module_iter<<endl;
				}
			}
			


			Double_t psi1_PSD_EP =  atan2(Q1y_PSD,Q1x_PSD);

			//cout<<psi1_PSD_EP<<endl;

	
			th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "Qvectors_PSD" )));
			th2_hist_ptr->Fill( psi1_PSD_EP,psi2_TPC_EP );	

			th1_hist_ptr = ((TH1*)(gDirectory->FindObjectAny( "PSD_EP" )));
			th1_hist_ptr->Fill( psi1_PSD_EP );

			tprof_ptr = ((TProfile*)(gDirectory->FindObjectAny( "Delta_psi" )));
			tprof_ptr->Fill( ( event->GetNVertexTracks() ), (psi1_PSD_EP-psi2_TPC_EP) );
			
		}//cut



	}



	TCanvas *canv_tracks_phi = new TCanvas( "VertexTrack_Phi", "VertexTrack_Phi" );
	th1_hist_ptr = ((TH1*)(gDirectory->FindObjectAny( "Vertex_tracks" )));
	th1_hist_ptr->Draw();

	TCanvas *canv_tracks_qvect = new TCanvas( "Qvect_Phi", "Qvect_Phi" );
	th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "Qvectors_TPC" )));
	th2_hist_ptr->Draw("colz");

	TCanvas *canv_tracks_qvect2 = new TCanvas( "Qvect_Phi2", "Qvect_Phi2" );
	th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "Qvectors_PSD" )));
	th2_hist_ptr->Draw("colz");

	TCanvas *canv_tracks_qvect24 = new TCanvas( "Qvect_Phi24", "Qvect_Phi24" );
	th1_hist_ptr = ((TH1*)(gDirectory->FindObjectAny( "PSD_EP" )));
	th1_hist_ptr->Draw();

	TCanvas *canv_tracks_qvect25 = new TCanvas( "Qvect_Phi25", "Qvect_Phi25" );
	tprof_ptr = ((TProfile*)(gDirectory->FindObjectAny( "Delta_psi" )));
	tprof_ptr->Draw();


}
