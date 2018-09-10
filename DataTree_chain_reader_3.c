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



void DataTree_chain_reader_6()
{

	Bool_t IsQxy_TPCmult_PSDen_comparison = 1;
	Bool_t IsQQ_after_recentering_comparison = 1;
	Bool_t IsRecentering_needed = 1;


	TFile *QAConfig = TFile::Open("/home/nikolay/MEPhI/DataTree_code_and_files/DataTreeQA/macro/QAConfigurations.root");
	cuts::DataTreeCutsConfig* cc = (cuts::DataTreeCutsConfig*) QAConfig->Get("pbpb_30agev_data_config_cuts");
	cuts::DataTreeCuts _cuts(cc);

	TString source_path = "/home/nikolay/MEPhI/DataTree_code_and_files/20180625_221329/";//30agev with magn
	//TString source_path = "/home/nikolay/MEPhI/DataTree_code_and_files/20180824_153409/";//30agev no magn
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

	th1_hist_ptr = new TH1F("PSD_EP","PSD_EP", 100, -4, 4);

	th2_hist_ptr = new TH2F("Qvectors_TPC","Qvectors_TPC", 70, -4, 4, 70, -4, 4);
	th2_hist_ptr = new TH2F("Qvectors_PSD","Qvectors_PSD", 70, -4, 4, 70, -4, 4);

	TProfile *tprof_ptr = NULL;
	tprof_ptr = new TProfile("Delta_psi","Delta_psi", 30, 0, 650);

	if(IsQxy_TPCmult_PSDen_comparison)
	{
		th2_hist_ptr = new TH2F("QxTPC_Mult","QxTPC_Mult", 70, -1, 1, 70, 0, 300);
		th2_hist_ptr = new TH2F("QyTPC_Mult","QyTPC_Mult", 70, -1, 1, 70, 0, 300);
		th2_hist_ptr = new TH2F("QxTPC_En","QxTPC_En", 70, -1, 1, 70, 1000, 8000);
		th2_hist_ptr = new TH2F("QyTPC_En","QyTPC_En", 70, -1, 1, 70, 1000, 8000);
		th2_hist_ptr = new TH2F("QxPSD_Mult","QxPSD_Mult", 70, -10, 10, 70, 0, 300);
		th2_hist_ptr = new TH2F("QyPSD_Mult","QyPSD_Mult", 70, -10, 10, 70, 0, 300);
		th2_hist_ptr = new TH2F("QxPSD_En","QxPSD_En", 70, -10, 10, 70, 1000, 8000);
		th2_hist_ptr = new TH2F("QyPSD_En","QyPSD_En", 70, -10, 10, 70, 1000, 8000);
	}

	if(IsQQ_after_recentering_comparison)
	{
		tprof_ptr = new TProfile("QxTPC_QxPSD","QxTPC_QxPSD", 20, 0, 250);
		tprof_ptr = new TProfile("QyTPC_QyPSD","QyTPC_QyPSD", 20, 0, 250);
		tprof_ptr = new TProfile("QxTPC_QyPSD","QxTPC_QyPSD", 20, 0, 250);
		tprof_ptr = new TProfile("QyTPC_QxPSD","QyTPC_QxPSD", 20, 0, 250);
	}

/////////////////Average q vector calculation for further recentering

	Double_t Av_Q1x_TPC=0;
	Double_t Av_Q1y_TPC=0;

	Double_t Av_Q2x_TPC=0;
	Double_t Av_Q2y_TPC=0;
	Double_t Temp_Av_Q2x_TPC=0;
	Double_t Temp_Av_Q2y_TPC=0;

	Double_t Av_Q1x_PSD=0;
	Double_t Av_Q1y_PSD=0;
	Double_t Temp_Av_Q1x_PSD=0;
	Double_t Temp_Av_Q1y_PSD=0;

	Double_t Av_Q2x_PSD=0;
	Double_t Av_Q2y_PSD=0;

	Double_t Total_PSD_Energy;
	Int_t TPC_iter_counter=0;
	Int_t PSD_iter_counter=0;
	Int_t Event_iter_counter=0;
	Int_t good_track_counter=0;

	for(Int_t entry = 0; entry <= data_chain->GetEntries(); entry++)
	{
		data_chain->GetEntry(entry);
		

		//cuts			
		if(_cuts.IsGoodEvent(*event))
		{
			Event_iter_counter++;
			Av_Q2x_TPC=0;
			Av_Q2y_TPC=0;
			good_track_counter = 0;
			for(Int_t track_iter = 0; track_iter < event->GetNVertexTracks(); track_iter++)
			{
				if( _cuts.IsGoodTrack( *(event->GetVertexTrack(track_iter)) ) )
				{
					Av_Q2x_TPC += cos( (event->GetVertexTrack(track_iter)->GetPhi()) );
					Av_Q2y_TPC += sin( (event->GetVertexTrack(track_iter)->GetPhi()) );

					good_track_counter++;
				}

			}
			Av_Q2x_TPC /= good_track_counter;
			Av_Q2y_TPC /= good_track_counter;

			Temp_Av_Q2x_TPC += Av_Q2x_TPC;
			Temp_Av_Q2y_TPC += Av_Q2y_TPC;
	
			Total_PSD_Energy=0;
			Av_Q1x_PSD=0;
			Av_Q1y_PSD=0;
			for(Int_t module_iter = 0; module_iter < event->GetNPSDModules(); module_iter++)
			{
				if((event->GetPSDModule(module_iter)->GetId()) <= 45)
				{
					Av_Q1x_PSD +=  (event->GetPSDModule(module_iter)->GetEnergy()) * (event->GetPSDModule(module_iter)->GetPositionComponent(0)) ;
					Av_Q1y_PSD +=  (event->GetPSDModule(module_iter)->GetEnergy()) * (event->GetPSDModule(module_iter)->GetPositionComponent(1)) ;
					Total_PSD_Energy += event->GetPSDModule(module_iter)->GetEnergy();

				}
			}
			
			Av_Q1x_PSD /= Total_PSD_Energy;
			Av_Q1y_PSD /= Total_PSD_Energy;
			
			Temp_Av_Q1x_PSD += Av_Q1x_PSD;
			Temp_Av_Q1y_PSD += Av_Q1y_PSD;

			if(IsQxy_TPCmult_PSDen_comparison)
			{
				th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxTPC_Mult" )));
				th2_hist_ptr->Fill( Av_Q2x_TPC, good_track_counter );
				th2_hist_ptr->SetTitle("QxTPC_vs_Multiplicity");
				th2_hist_ptr->GetXaxis()->SetTitle("QxTPC");
				th2_hist_ptr->GetYaxis()->SetTitle("Multiplicity");

				th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyTPC_Mult" )));
				th2_hist_ptr->Fill( Av_Q2y_TPC, good_track_counter );
				th2_hist_ptr->SetTitle("QyTPC_vs_Multiplicity");
				th2_hist_ptr->GetXaxis()->SetTitle("QyTPC");
				th2_hist_ptr->GetYaxis()->SetTitle("Multiplicity");

				th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxTPC_En" )));
				th2_hist_ptr->Fill( Av_Q2x_TPC, Total_PSD_Energy );
				th2_hist_ptr->SetTitle("QxTPC_vs_PSDEnergy");
				th2_hist_ptr->GetXaxis()->SetTitle("QxTPC");
				th2_hist_ptr->GetYaxis()->SetTitle("PSDEnergy");

				th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyTPC_En" )));
				th2_hist_ptr->Fill( Av_Q2y_TPC, Total_PSD_Energy );
				th2_hist_ptr->SetTitle("QyTPC_vs_PSDEnergy");
				th2_hist_ptr->GetXaxis()->SetTitle("QyTPC");
				th2_hist_ptr->GetYaxis()->SetTitle("PSDEnergy");

				th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxPSD_Mult" )));
				th2_hist_ptr->Fill( Av_Q1x_PSD, good_track_counter );
				th2_hist_ptr->SetTitle("QxPSD_vs_Multiplicity");
				th2_hist_ptr->GetXaxis()->SetTitle("QxPSD");
				th2_hist_ptr->GetYaxis()->SetTitle("Multiplicity");

				th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyPSD_Mult" )));
				th2_hist_ptr->Fill( Av_Q1y_PSD, good_track_counter );
				th2_hist_ptr->SetTitle("QyPSD_vs_Multiplicity");
				th2_hist_ptr->GetXaxis()->SetTitle("QyPSD");
				th2_hist_ptr->GetYaxis()->SetTitle("Multiplicity");

				th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxPSD_En" )));
				th2_hist_ptr->Fill( Av_Q1x_PSD, Total_PSD_Energy );
				th2_hist_ptr->SetTitle("QxPSD_vs_PSDEnergy");
				th2_hist_ptr->GetXaxis()->SetTitle("QxPSD");
				th2_hist_ptr->GetYaxis()->SetTitle("PSDEnergy");

				th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyPSD_En" )));
				th2_hist_ptr->Fill( Av_Q1y_PSD, Total_PSD_Energy );
				th2_hist_ptr->SetTitle("QyPSD_vs_PSDEnergy");
				th2_hist_ptr->GetXaxis()->SetTitle("QyPSD");
				th2_hist_ptr->GetYaxis()->SetTitle("PSDEnergy");
			}

//cout<<Av_Q1x_PSD<<" "<<Av_Q1y_PSD<<" "<<Total_PSD_Energy<<" "<<Av_Q2x_TPC<<" "<<Av_Q2y_TPC<<" "<<event->GetNVertexTracks()<<endl;

		}
	}

	Av_Q2x_TPC = Temp_Av_Q2x_TPC/Event_iter_counter ;
	Av_Q2y_TPC = Temp_Av_Q2y_TPC/Event_iter_counter ;
	
	Av_Q1x_PSD = Temp_Av_Q1x_PSD/Event_iter_counter ;
	Av_Q1y_PSD = Temp_Av_Q1y_PSD/Event_iter_counter ;


	if(IsQxy_TPCmult_PSDen_comparison)
	{
		TCanvas *canv_av_q = new TCanvas( "Av_q", "Av_q" );
		canv_av_q->DivideSquare(8);
		canv_av_q->cd(1);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxTPC_Mult" )));
		th2_hist_ptr->Draw();
		canv_av_q->cd(2);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyTPC_Mult" )));
		th2_hist_ptr->Draw();
		canv_av_q->cd(3);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxTPC_En" )));
		th2_hist_ptr->Draw();
		canv_av_q->cd(4);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyTPC_En" )));
		th2_hist_ptr->Draw();
		canv_av_q->cd(5);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxPSD_Mult" )));
		th2_hist_ptr->Draw();
		canv_av_q->cd(6);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyPSD_Mult" )));
		th2_hist_ptr->Draw();
		canv_av_q->cd(7);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxPSD_En" )));
		th2_hist_ptr->Draw();
		canv_av_q->cd(8);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyPSD_En" )));
		th2_hist_ptr->Draw();
	}


//cout<<Av_Q1x_PSD<<" "<<Av_Q1y_PSD<<" "<<atan2(Av_Q1y_PSD, Av_Q1x_PSD)<<" "<<Av_Q2x_TPC<<" "<<Av_Q2y_TPC<<" "<<atan2(Av_Q2y_TPC, Av_Q2x_TPC)<<endl;

//////////////////////////////////////////////////////////////////////////////////////////////

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
			good_track_counter=0;
			for(Int_t track_iter = 0; track_iter < event->GetNVertexTracks(); track_iter++)
			{
				if( _cuts.IsGoodTrack( *(event->GetVertexTrack(track_iter)) ) )
				{
					th1_hist_ptr = ((TH1*)(gDirectory->FindObjectAny( "Vertex_tracks" )));
					th1_hist_ptr->Fill( event->GetVertexTrack(track_iter)->GetPhi() );


					//Q-vector
					Q1x_TPC += event->GetVertexTrack(track_iter)->GetPx();
					Q1y_TPC += event->GetVertexTrack(track_iter)->GetPy();
					//Q2x_TPC += (event->GetVertexTrack(track_iter)->GetPt()) * cos( 2*(event->GetVertexTrack(track_iter)->GetPhi()) );
					//Q2y_TPC += (event->GetVertexTrack(track_iter)->GetPt()) * sin( 2*(event->GetVertexTrack(track_iter)->GetPhi()) );
					Q2x_TPC += cos( (event->GetVertexTrack(track_iter)->GetPhi()) );
					Q2y_TPC += sin( (event->GetVertexTrack(track_iter)->GetPhi()) );

					good_track_counter++;
				}

			}
			Q2x_TPC /= good_track_counter;
			Q2y_TPC /= good_track_counter;

			(IsRecentering_needed) ? (Q2x_TPC -= Av_Q2x_TPC):(Q2x_TPC = Q2x_TPC) ;
			(IsRecentering_needed) ? (Q2y_TPC -= Av_Q2y_TPC):(Q2y_TPC = Q2y_TPC) ;

			Double_t psi1_TPC_EP = atan2(Q1y_TPC,Q1x_TPC);
			Double_t psi2_TPC_EP = atan2(Q2y_TPC,Q2x_TPC);
			th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "Qvectors_TPC" )));
			th2_hist_ptr->Fill( psi1_TPC_EP,psi2_TPC_EP );

			//Q-vector PSD

			Double_t Q1x_PSD=0;
			Double_t Q1y_PSD=0;

			Double_t Q2x_PSD=0;
			Double_t Q2y_PSD=0;
			
			Total_PSD_Energy=0;
			PSD_iter_counter=0;
			for(Int_t module_iter = 0; module_iter < event->GetNPSDModules(); module_iter++)
			{
				if((event->GetPSDModule(module_iter)->GetId()) <= 45)
				{
					Q1x_PSD +=  (event->GetPSDModule(module_iter)->GetEnergy()) * (event->GetPSDModule(module_iter)->GetPositionComponent(0)) ;
					Q1y_PSD +=  (event->GetPSDModule(module_iter)->GetEnergy()) * (event->GetPSDModule(module_iter)->GetPositionComponent(1)) ;
					Total_PSD_Energy += (event->GetPSDModule(module_iter)->GetEnergy());
					//cout<<(event->GetPSDModule(module_iter)->GetPositionComponent(2))<<" "<<module_iter<<endl;
				}
			}
			
			Q1x_PSD /= Total_PSD_Energy;
			Q1y_PSD /= Total_PSD_Energy;

			(IsRecentering_needed) ? (Q1x_PSD -= Av_Q1x_PSD):(Q1x_PSD = Q1x_PSD) ;
			(IsRecentering_needed) ? (Q1y_PSD -= Av_Q1y_PSD):(Q1y_PSD = Q1y_PSD) ;

			Double_t psi1_PSD_EP =  atan2(Q1y_PSD,Q1x_PSD);

			//cout<<psi1_PSD_EP<<endl;

	
			th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "Qvectors_PSD" )));
			th2_hist_ptr->Fill( psi1_PSD_EP,psi2_TPC_EP );	

			th1_hist_ptr = ((TH1*)(gDirectory->FindObjectAny( "PSD_EP" )));
			th1_hist_ptr->Fill( psi1_PSD_EP );

			tprof_ptr = ((TProfile*)(gDirectory->FindObjectAny( "Delta_psi" )));
			tprof_ptr->Fill( ( event->GetNVertexTracks() ), (psi1_PSD_EP-psi2_TPC_EP) );

			if(IsQQ_after_recentering_comparison)
			{

				tprof_ptr = ((TProfile*)(gDirectory->FindObjectAny( "QxTPC_QxPSD" )));
				tprof_ptr->Fill( good_track_counter, Q2x_TPC*Q1x_PSD );
				tprof_ptr->SetTitle("QxTPC*QxPSD_vs_multiplicity");
				tprof_ptr->GetXaxis()->SetTitle("multiplicity");
				tprof_ptr->GetYaxis()->SetTitle("QxTPC*QxPSD");

				tprof_ptr = ((TProfile*)(gDirectory->FindObjectAny( "QyTPC_QyPSD" )));
				tprof_ptr->Fill( good_track_counter, Q2y_TPC*Q1y_PSD );
				tprof_ptr->SetTitle("QyTPC*QyPSD_vs_multiplicity");
				tprof_ptr->GetXaxis()->SetTitle("multiplicity");
				tprof_ptr->GetYaxis()->SetTitle("QyTPC*QyPSD");

				tprof_ptr = ((TProfile*)(gDirectory->FindObjectAny( "QxTPC_QyPSD" )));
				tprof_ptr->Fill( good_track_counter, Q2x_TPC*Q1y_PSD );
				tprof_ptr->SetTitle("QxTPC*QyPSD_vs_multiplicity");
				tprof_ptr->GetXaxis()->SetTitle("multiplicity");
				tprof_ptr->GetYaxis()->SetTitle("QxTPC*QyPSD");

				tprof_ptr = ((TProfile*)(gDirectory->FindObjectAny( "QyTPC_QxPSD" )));
				tprof_ptr->Fill( good_track_counter, Q2y_TPC*Q1x_PSD );
				tprof_ptr->SetTitle("QyTPC*QxPSD_vs_multiplicity");
				tprof_ptr->GetXaxis()->SetTitle("multiplicity");
				tprof_ptr->GetYaxis()->SetTitle("QyTPC*QxPSD");	
			
			}
		
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


	if(IsQQ_after_recentering_comparison)
	{
		TCanvas *canv_qq = new TCanvas( "qq", "qq" );
		canv_qq->DivideSquare(4);
		canv_qq->cd(1);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxTPC_QxPSD" )));
		th2_hist_ptr->Draw();
		canv_qq->cd(2);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyTPC_QyPSD" )));
		th2_hist_ptr->Draw();
		canv_qq->cd(3);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QxTPC_QyPSD" )));
		th2_hist_ptr->Draw();
		canv_qq->cd(4);
		th2_hist_ptr = ((TH2*)(gDirectory->FindObjectAny( "QyTPC_QxPSD" )));
		th2_hist_ptr->Draw();
	}

}
