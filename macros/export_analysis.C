#include "../analysis_functions.C"



int main (int argc, char *argv[])
{
	int run_number = -1;
	int event_number = -1;

	std::string input_filename = "";

	for (int iarg=1; iarg<argc; ++iarg)
	{
      		std::string arg (argv[iarg]);
		if (arg[0] == '-')
		{
			if (arg=="-i" || arg=="--input")
			input_filename = std::string(argv[++iarg]);

			else if (arg=="-r" || arg=="--run")
			run_number = atoi(argv[++iarg]);

			else if (arg=="-e" || arg=="--event")
			event_number = atoi(argv[++iarg]);
		}
	}

////open root file TTree
	string run_num=to_string(run_number);
	const char *file = Form("/sps/nemo/scratch/spratt/analysis/vertical_position_reconstruction/run_output_analysis_files/%d/export_file%d.root", run_number,run_number);
	open_file(file);


///////////////////////////////////////
/////////////Histograms////////////////
///////////////////////////////////////

TH1D* counter_2 = new TH1D("counter2","counter2",100,-1,1);
TH1D* sum_2 = new TH1D("sum2","sum2",100,-1,1);
TH1D* output_2 = new TH1D("angle_vs_height","angle_vs_height",100,-1,1);

TH1D* t6_vs_z_counter = new TH1D("t6_vs_z_counter","t6_vs_z_counter",20,-1,1);
TH1D* t6_vs_z_sum = new TH1D("t6_vs_z_sum","t6_vs_z_sum",20,-1,1);
TH1D* t6_vs_z_output = new TH1D("t6_vs_z_output","t6_vs_z_output",20,-1,1);

TH1D* t6_vs_z_counter2 = new TH1D("t6_vs_z_counter2","t6_vs_z_counter2",200,-1,1);
TH1D* t6_vs_z_sum2 = new TH1D("t6_vs_z_sum2","t6_vs_z_sum2",200,-1,1);
TH1D* t6_vs_z_output2 = new TH1D("t6_vs_z_output2","t6_vs_z_output2",200,-1,1);

TH1D* t6_mode_hist_graph = new TH1D("t6_mode_hist_graph","t6_mode_hist_graph",100,-1,1);

TH1D* counter = new TH1D("counter","counter",100,-1,1);
TH1D* sum = new TH1D("sum","sum",100,-1,1);
TH1D* output = new TH1D("output","output",100,-1,1);


TH1D* counter_crate_1 = new TH1D("counter_crate_1","counter_crate_1",30,-1,1);
TH1D* sum_crate_1 = new TH1D("sum_crate_1","sum_crate_1",30,-1,1);
TH1D* output_crate_1 = new TH1D("output_crate_1","output_crate_1",30,-1,1);



TH1D* counter_angle = new TH1D("counter_angle","counter_angle",100,-1,1);
TH1D* sum_angle = new TH1D("sum_angle","sum_angle",100,-1,1);
TH1D* output_angle = new TH1D("output_angle","output_angle",100,-1,1);

TH2D* prop_vs_z = new TH2D("prop_vs_z","prop_vs_z",100,-1,1,100,40000,80000);

TH2D* t5_vs_t6_total = new TH2D("t5_vs_t6_total","t5_vs_t6_total",100,0,60000,100,0,60000);

TH2D* t6_vs_z_2d = new TH2D("t6_vs_z_2d","t6_vs_z_2d",70,-1,1,400,0,60000);

TH2D* prop_vs_z_crate_1 = new TH2D("prop_vs_z_crate_1","prop_vs_z_crate_1",70,-1,1,400,40000,80000);

TGraphErrors* t6_dec_vs_z = new TGraphErrors(200);
TGraph* t6_dec_vs_z_grad = new TGraph(0);
TGraphErrors* t6_width_graph = new TGraphErrors(70);
TGraphErrors* prop_vs_z_graph_width = new TGraphErrors(70);


TH1D** angle_vs_height_counter = make_array_of_hist("angle_vs_height_counter",100,-1,1);
TH1D** angle_vs_height_sum = make_array_of_hist("angle_vs_height_sum",100,-1,1);
TH1D** angle_vs_height_output = make_array_of_hist("angle_vs_height_output",100,-1,1);


TH1D** prop_vs_z_low_angle_counter = make_array_of_hist("prop_vs_z_low_angle_counter",25,-1,1);
TH1D** prop_vs_z_low_angle_sum = make_array_of_hist("prop_vs_z_low_angle_sum",25,-1,1);
TH1D** prop_vs_z_low_angle_output = make_array_of_hist("prop_vs_z_low_angle_output",25,-1,1);
TH2D** prop_vs_z_low_angle_2d = make_array_of_hist("prop_vs_z_low_angle_2d",25,-1,1,1000,30000,80000);

TH2D** t5_vs_t6 = make_array_of_hist("t5_vs_t6",50,0,60000,50,0,60000);

double previous_TS_for_function_modified[2034][3]={0};

TH1D** prop_time_hist = make_array_of_hist("prop_time_hist",1000,30000,80000);
TH1D** prop_time_hist_selection = make_array_of_hist("prop_time_hist_selection",1000,30000,80000);

make_dir(run_num+"/export_analysis");

for (int i=0;i<get_entries;i++){
    event_tree->GetEntry(i); 

double previous_ts_hit_modified_R0 = previous_TS_for_function_modified[cell_num][0];
double previous_ts_hit_modified_R5 = previous_TS_for_function_modified[cell_num][1];
double previous_ts_hit_modified_R6 = previous_TS_for_function_modified[cell_num][2];



if(t5>0 && t6>0 && t5<1000000 && t6<1000000){
     double z = (t5-t6)/(t5+t6);
     double prop_time = t5+t6;
     double deadtime = t0-previous_ts_hit_modified_R0;
     prop_time_hist[cell_num]->Fill(prop_time,1);

     angle_vs_height_counter[cell_num]->Fill(z,1);
     angle_vs_height_sum[cell_num]->Fill(z,theta);


if (deadtime>100000000){
if (t_drift<1000){

if (check_cell_on_edge(cell_num)==false  && is_cell_bad(cell_num)==false && is_next_to_bad_cell(cell_num)==false){
if (what_is_crate_number(cell_num)==1){

counter_crate_1->Fill(z,1);
sum_crate_1->Fill(z,t5+t6);

t6_vs_z_counter->Fill(z,1);
t6_vs_z_sum->Fill(z,t6);
t6_vs_z_2d->Fill(z,t6,1);

t6_vs_z_counter2->Fill(z,1);
t6_vs_z_sum2->Fill(z,t6);


prop_vs_z_crate_1->Fill(z,t5+t6,1);

}
}}
}



t5_vs_t6_total->Fill(t5,t6,1);

t5_vs_t6[cell_num]->Fill(t5,t6,1);

 //if (check_cell_on_edge(cell_num)==false  && is_cell_bad(cell_num)==false && is_next_to_bad_cell(cell_num)==false){

//}


counter_2->Fill(z,1);
sum_2->Fill(z,theta);


 //    cout<<deadtime<<endl;
 //    cout<<t_drift<<endl;
 //    cout<<"----"<<endl;

prop_vs_z->Fill(z,t5+t6);

if (t_drift<400){
if (deadtime>100000000){
if (theta>-0.2 && theta <0.2){
prop_time_hist_selection[cell_num]->Fill(prop_time,1);
}}}


counter->Fill(z,1);
sum->Fill(z,t5+t6);

if (theta>-0.2 && theta <0.2){

counter_angle->Fill(z,1);
sum_angle->Fill(z,t5+t6);

}






prop_vs_z_low_angle_counter[cell_num]->Fill(z,1);
prop_vs_z_low_angle_sum[cell_num]->Fill(z,prop_time);
prop_vs_z_low_angle_2d[cell_num]->Fill(z,prop_time,1);

//}
}

previous_TS_for_function_modified[cell_num][0] = t0;
previous_TS_for_function_modified[cell_num][1] = t5; 
previous_TS_for_function_modified[cell_num][2] = t6;

}


cout<<"averaging"<<endl;
average_over_hist(prop_vs_z_low_angle_sum,prop_vs_z_low_angle_counter,prop_vs_z_low_angle_output);
average_over_hist(sum_angle,counter_angle,output_angle);
average_over_hist(angle_vs_height_sum,angle_vs_height_counter,angle_vs_height_output);
average_over_hist(sum,counter,output);
average_over_hist(sum_2,counter_2,output_2);
average_over_hist(t6_vs_z_sum,t6_vs_z_counter,t6_vs_z_output);
average_over_hist(t6_vs_z_sum2,t6_vs_z_counter2,t6_vs_z_output2);
average_over_hist(sum_crate_1,counter_crate_1,output_crate_1);

double bin_content=0;
for (int l=0;l<output_crate_1->GetNbinsX();l++){
bin_content=0;
bin_content = counter_crate_1->GetBinContent(l);
output_crate_1->SetBinError(l,0.5*((output_crate_1->GetBinContent(l))/(sqrt(bin_content))));
}

////////fitting






////////t6 vs z fitting and analysis !!!


TH1D** y_projection = new TH1D*[t6_vs_z_2d->GetNbinsX()];
TF1** gaus_fits = new TF1*[t6_vs_z_2d->GetNbinsX()];
//TH1D** y_projections_scaled = new TH1D*[hist_2d->GetNbinsX()];

//TF1** fit_gaussian_functions_array = new TF1*[hist_2d->GetNbinsX()];
//TF1** fit_gaussian_functions_array_with_shift = new TF1*[hist_2d->GetNbinsX()];
for (int i=0; i<t6_vs_z_2d->GetNbinsX();i++){
y_projection[i] = t6_vs_z_2d->ProjectionY(("t6_width"+to_string(i)).c_str(),i,i+1);

  for(int k=0; k < y_projection[i]->GetNbinsX();  k++) {
      y_projection[i]->SetBinError(k,TMath::Sqrt(y_projection[i]->GetBinContent(k)));
   }

t6_mode_hist_graph->SetBinContent(i,y_projection[i]->GetXaxis()->GetBinCenter(y_projection[i]->GetMaximumBin()));

gaus_fits[i] = new TF1(("gaus_fits"+to_string(i)).c_str(),"gaus",y_projection[i]->GetXaxis()->GetBinCenter(y_projection[i]->GetMaximumBin())-2000,y_projection[i]->GetXaxis()->GetBinCenter(y_projection[i]->GetMaximumBin())+400);
y_projection[i]->Fit(gaus_fits[i],"RQ");

if (gaus_fits[i]->GetParError(2)<1000){

//t6_width_graph->SetPointError(i,0,gaus_fits[i]->GetParError(2));

t6_width_graph->SetPoint(i,i,gaus_fits[i]->GetParameter(2));
t6_width_graph->SetPointError(i,0,gaus_fits[i]->GetParError(2));
}}



////end of t6 fitting and analysis !!




////prop_vs_z_fitting_analysis





TH1D** y_projectionpz = new TH1D*[prop_vs_z_crate_1->GetNbinsX()];
TF1** gaus_fitspz = new TF1*[prop_vs_z_crate_1->GetNbinsX()];
//TH1D** y_projections_scaled = new TH1D*[hist_2d->GetNbinsX()];

//TF1** fit_gaussian_functions_array = new TF1*[hist_2d->GetNbinsX()];
//TF1** fit_gaussian_functions_array_with_shift = new TF1*[hist_2d->GetNbinsX()];
for (int i=0; i<prop_vs_z_crate_1->GetNbinsX();i++){
y_projectionpz[i] = prop_vs_z_crate_1->ProjectionY(("prop_vs_z_crate_1"+to_string(i)).c_str(),i,i+1);

  for(int k=0; k < y_projectionpz[i]->GetNbinsX();  k++) {
      y_projectionpz[i]->SetBinError(k,TMath::Sqrt(y_projectionpz[i]->GetBinContent(k)));
   }

gaus_fitspz[i] = new TF1(("gaus_fits_prop_vs_z"+to_string(i)).c_str(),"gaus",y_projectionpz[i]->GetXaxis()->GetBinCenter(y_projectionpz[i]->GetMaximumBin())-2000,y_projectionpz[i]->GetXaxis()->GetBinCenter(y_projectionpz[i]->GetMaximumBin())+400);
y_projectionpz[i]->Fit(gaus_fitspz[i],"RQ");

if (gaus_fitspz[i]->GetParError(2)<1000){

//t6_width_graph->SetPointError(i,0,gaus_fits[i]->GetParError(2));

prop_vs_z_graph_width->SetPoint(i,i,gaus_fitspz[i]->GetParameter(2));
prop_vs_z_graph_width->SetPointError(i,0,gaus_fitspz[i]->GetParError(2));
}}
/////////// end of fitting analysis















TF1* straight_line_fit = new TF1("straight_line_fit","([0]*x)+[1]",0.7, 0.85);
t6_vs_z_output2->Fit(straight_line_fit,"R");
//t6_mode_hist_graph->Fit(straight_line_fit,"R");
for (int i=1;i<20;i++){

double a =straight_line_fit->GetParameter(0);
double b =straight_line_fit->GetParameter(1);
double x = t6_vs_z_output->GetBinCenter(i);

t6_dec_vs_z_grad->SetPoint(t6_dec_vs_z_grad->GetN(),x,(t6_vs_z_output->GetBinContent(i)-t6_vs_z_output->GetBinContent(i-1))/(t6_vs_z_output->GetBinCenter(i)-t6_vs_z_output->GetBinCenter(i-1)));
//t6_mode_hist_graph
t6_dec_vs_z->SetPoint(i,i,t6_vs_z_output->GetBinContent(i)-((a*x)+b));
t6_dec_vs_z->SetPointError(i,0,(t6_vs_z_output->GetBinContent(i)/(sqrt(t6_vs_z_counter->GetBinContent(i))))*0.5);

}



  TF1 *fitFunction_iteration_eps_vs_z = new TF1("fitFunction_iteration_eps_vs_z", get_tp_given_eps ,-0.6,0.6,2);
  fit_iteration_eps_vs_z(output_crate_1,fitFunction_iteration_eps_vs_z);




/*
for (int j=1;j<8;j++){
      double counter_sum =0;
for (int i=(j*10);i<((j+1)*10);i++){
double a =straight_line_fit->GetParameter(0);
double b =straight_line_fit->GetParameter(1);
double x = t6_vs_z_output->GetBinCenter(i);
counter_sum += t6_vs_z_output->GetBinContent(i)-((a*x)+b);
}
t6_dec_vs_z->SetPoint(t6_dec_vs_z->GetN(),j,counter_sum*0.1);
}
*/



cout<<"printing"<<endl;


save_all_cell_hists(t5_vs_t6,"export_analysis/t5_vs_t6UDD",run_num);
save_all_cell_hists(angle_vs_height_output,"export_analysis/angle_vs_height_outputUDD",run_num);
save_all_cell_hists(prop_time_hist,"export_analysis/prop_time_hist_cellUDD",run_num);
save_all_cell_hists(prop_time_hist_selection,"export_analysis/prop_time_hist_selectionUDD",run_num);

save_all_cell_hists(prop_vs_z_low_angle_output,"export_analysis/prop_vs_z_average_low_angle_cellUDD",run_num);
save_all_cell_hists(prop_vs_z_low_angle_counter,"export_analysis/prop_vs_z_counter_low_angle_cellUDD",run_num);
save_all_cell_hists(prop_vs_z_low_angle_sum,"export_analysis/prop_vs_z_sum_low_angle_cellUDD",run_num);
save_all_cell_hists(prop_vs_z_low_angle_2d,"export_analysis/prop_vs_z_2d_low_angle_cellUDD",run_num);


/////////////////////////Load new root file///////////////////////////////////////
string folder = name_folder(run_num);
string title_root_file = "/export_analysisUDD.root";
TFile *rootfile = new TFile((folder+title_root_file).c_str(), "RECREATE");//saving all histograms
rootfile->cd();
//////////////////////////////////////////////////////////////////////////////////

prop_vs_z_graph_width->Write("prop_vs_z_graph_width");
t6_mode_hist_graph->Write();
output_crate_1->Write();

t6_width_graph->Write("t6_width");

for (int i=0;i<69;i++){
y_projection[i]->Write();
}

for (int i=0;i<69;i++){
y_projectionpz[i]->Write();
}
fitFunction_iteration_eps_vs_z->Write();

counter_crate_1->Write();

t6_vs_z_2d->Write();

t5_vs_t6_total->Write();

t6_dec_vs_z_grad->Write("Grad");
t6_vs_z_counter->Write();
t6_dec_vs_z->Write("t6 vs straight line difference");
t6_vs_z_output->Write();
t6_vs_z_output2->Write();




output_2->Write();
output_angle->Write();
output->Write();

prop_vs_z->Write();

///////////////////////////Close Files////////////////////////////////////////////
rootfile->Close();
//////////////////////////////////////////////////////////////////////////////////



}
