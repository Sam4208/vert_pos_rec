#include "analysis_functions.h"

TFile *run_file = nullptr;
TTree *event_tree = nullptr;

int run_number = -1;
int run_start = -1;

long long int current_entry = -1;

    int           event;
    int           calo_hit;
    int           tracker_hit;
    int           cell_num;
    int           om_num;
    double          theta;
    double           phi;
    double           h;
    double           r;
    double           tc;
    double           t0;
    double           t1;
    double           t2;
    double           t3;
    double           t4;
    double           t5;
    double           t6;
    double 		   t_drift;
    double           get_entries;


bool open_file(const char *filename) //<
{
  gStyle->SetOptTitle(1);

  if (run_file != nullptr)
    {
     run_file->Close();
     event_tree = nullptr;
     run_number = -1;
     run_start = -1;
    }

  run_file = TFile::Open(filename, "READ");
  if (run_file == nullptr) return false;

  gROOT->cd();

  event_tree = (TTree*)(run_file->Get("Data"));

  if (event_tree != nullptr)
    {
    event_tree->SetBranchAddress("event", &event);
    event_tree->SetBranchAddress("calo_hit", &calo_hit);
    event_tree->SetBranchAddress("tracker_hit", &tracker_hit);
    event_tree->SetBranchAddress("cell_num", &cell_num);
    event_tree->SetBranchAddress("om_num", &om_num);
    event_tree->SetBranchAddress("theta", &theta);
    event_tree->SetBranchAddress("phi", &phi);
    event_tree->SetBranchAddress("h", &h);
    event_tree->SetBranchAddress("r", &r);
    event_tree->SetBranchAddress("tc", &tc);
    event_tree->SetBranchAddress("t0", &t0);
    event_tree->SetBranchAddress("t1", &t1);
    event_tree->SetBranchAddress("t2", &t2);
    event_tree->SetBranchAddress("t3", &t3);
    event_tree->SetBranchAddress("t4", &t4);
    event_tree->SetBranchAddress("t5", &t5);
    event_tree->SetBranchAddress("t6", &t6);
    event_tree->SetBranchAddress("t_drift", &t_drift);

    cout << "Opened File " << filename << ", " << event_tree->GetEntries() << " events available." << endl << endl;

     get_entries = (event_tree->GetEntries() - 1);
 //   get_entries = 10000;


    }


/*

  else printf("+++ no 'event_tree' found in %s\n", filename);

  // autodetect run number
  std::string run_filename (filename);
  std::string run_prefix ("run-");
  std::string run_suffix ("_");
  size_t pos_prefix = run_filename.find(run_prefix) + run_prefix.size();
  size_t pos_suffix = run_filename.find(run_suffix, pos_prefix);
  run_number = atoi(run_filename.substr(pos_prefix, pos_suffix).c_str());
  printf("+++ run_number = %d\n", run_number);

  // look for run start time

  const char *cbd_base_path = "/sps/nemo/snemo/snemo_data/raw_data/CBD";

  std::vector<std::string> log_paths;
  log_paths.push_back(Form("%s/run-%d/snemo_trigger_run-%d.log", cbd_base_path, run_number, run_number));

  for (int crate=6; crate>=0; crate--)
    log_paths.push_back(Form("%s/run-%d/snemo_crate-%d_run-%d.log", cbd_base_path, run_number, crate, run_number));

  for (const std::string & log_path : log_paths)
    {
      std::ifstream log_file (log_path);
      if (!log_file.is_open()) continue;

      std::string log_line;

      while (getline(log_file, log_line))
	{
	  size_t unixtime_index = log_line.find("run.run_unixtime_ms=");

	  if (unixtime_index == std::string::npos)
	    continue;

	  int unixtime = std::stoi(log_line.substr(20));

	  if (unixtime > run_start)
	    run_start = unixtime;
	}
    }

  TDatime run_start_datime (run_start);
  printf("+++ run_start = %s\n", run_start_datime.AsSQLString());

  */

  return true;
}


void combine_runs(int run_num, TH2D** prop_vs_z_2D_hist_array,TH1D** prop_vs_z_1D_hist_sum_array,TH1D** prop_vs_z_1D_hist_counter_array){

cout<<"Combining Run "<<run_num<<endl;

string folder = name_folder(to_string(run_num).c_str());

string prop_time_vs_z_sum = (folder+"/export_analysis/prop_vs_z_sum_low_angle_cell.root").c_str();
string prop_time_vs_z_counter = (folder+"/export_analysis/prop_vs_z_counter_low_angle_cell.root").c_str();
string prop_time_vs_z_2d = (folder+"/export_analysis/prop_vs_z_2d_low_angle_cell.root").c_str();
string prop_time_vs_z_average = (folder+"/export_analysis/prop_vs_z_average_low_angle_cell.root").c_str();

TFile *prop_time_vs_z_sum_file = new TFile(prop_time_vs_z_sum.c_str(), "READ");
TFile *prop_time_vs_z_counter_file = new TFile(prop_time_vs_z_counter.c_str(), "READ");
TFile *prop_time_vs_z_2d_file = new TFile(prop_time_vs_z_2d.c_str(), "READ");
TFile *prop_time_vs_z_average_file = new TFile(prop_time_vs_z_average.c_str(), "READ");

int no_cells=2034;

for (int i=0;i<no_cells;i++){

TH1D* hist = (TH1D*)prop_time_vs_z_average_file->Get(("prop_vs_z_low_angle_output"+to_string(i)).c_str());

double sum=0;
for(int j=20;j<30;j++){
sum += hist->GetBinContent(j);}

double normalise_value = sum/10;

if (normalise_value>0){
TH1D* hist2 = (TH1D*)prop_time_vs_z_sum_file->Get(("prop_vs_z_low_angle_sum"+to_string(i)).c_str());
hist2->Scale(1/normalise_value);

prop_vs_z_2D_hist_array[i]->Add((TH2D*)prop_time_vs_z_2d_file->Get(("prop_vs_z_low_angle_2d"+to_string(i)).c_str()));
prop_vs_z_1D_hist_sum_array[i]->Add(hist2);
prop_vs_z_1D_hist_counter_array[i]->Add((TH1D*)prop_time_vs_z_counter_file->Get(("prop_vs_z_low_angle_counter"+to_string(i)).c_str()));}
}
}






//name a folder
string name_folder(string run_num){
  
    string directory = "/sps/nemo/scratch/spratt/analysis/vertical_position_reconstruction/run_output_analysis_files/";
    std::string folder = directory+run_num;

    return folder;
}











////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////Check Type of Cell///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool check_cell_on_edge(int i){

bool output = false;

int side_;
int column_;
int layer_;

int cell_num;

for (int side=0;side<2;side++){
    for (int column=0;column<113;column++){
        for (int layer=0;layer<9;layer++){
        cell_num =  side*9*113 + column*9 + layer;

//cout<<cell_num<<endl;

if (cell_num == i){
side_=side;
column_=column;
layer_=layer;
}


}}}

if (layer_ == 0 || layer_== 8 ){
output=true;
}
if (column_==0 || column_==112 ){
output=true;
}

return output;

}




bool is_cell_new(double side, double column, double layer){

      std::vector<double> bad_cell_column={9,106,99,53,82,102,99,26};
      std::vector<double> bad_cell_layer={7,2,7,1,1,2,8,3};
      std::vector<double> bad_cell_side={1,0,0,0,1,0,0,0};
      bool output=false;
   
      for (int i=0; i< bad_cell_layer.size(); i++){
         
      if (side==bad_cell_side.at(i)  &&  layer==bad_cell_layer.at(i)   &&  column  == bad_cell_column.at(i)){

         output=true;}}
       
      return output;

}




bool is_cell_other(double side, double column, double layer){

      std::vector<double> bad_cell_column={56,56,56,56};
      std::vector<double> bad_cell_layer={4,5,4,5};
      std::vector<double> bad_cell_side={0,0,1,1};
      bool output=false;
   
      for (int i=0; i< bad_cell_layer.size(); i++){
         
      if (side==bad_cell_side.at(i)  &&  layer==bad_cell_layer.at(i)   &&  column  == bad_cell_column.at(i)){

         output=true;}}
       
      return output;

}


bool is_cell_noisy(double side, double column, double layer){

      std::vector<double> bad_cell_column={1,2,4};
      std::vector<double> bad_cell_layer={1,1,3};
      std::vector<double> bad_cell_side={0,0,0};
      bool output=false;

      for (int i=0; i< bad_cell_layer.size(); i++){
      if (side==bad_cell_side[i]  &&  layer==bad_cell_layer[i]   &&  column  == bad_cell_column[i]  ){
         output=true;}}

      return output;
}



bool is_cell_HV_trip(double side, double column, double layer){

      std::vector<double> bad_cell_column={3,9,21,84,32,70,75,73,80,100};
      std::vector<double> bad_cell_layer={0,0,0,0,0,0,4,2,2,2};
      std::vector<double> bad_cell_side={0,0,0,0,1,0,0,1,1,1};
      bool output=false;

      for (int i=0; i< bad_cell_layer.size(); i++){
      if (side==bad_cell_side[i]  &&  layer==bad_cell_layer[i]   &&  column  == bad_cell_column[i]  ){
         output=true;}
         

         
         }


      return output;
}



bool is_cell_internal_cabling(double side, double column, double layer){

      std::vector<double> bad_cell_column={86,87,63,57,56,47,79,84,91,99,107,110};
      std::vector<double> bad_cell_layer={8,0,8,0,3,8,6,6,2,8,8,0};
      std::vector<double> bad_cell_side={0,0,0,0,0,1,1,1,1,1,1,1};
      bool output=false;

      for (int i=0; i< bad_cell_layer.size(); i++){
      if (side==bad_cell_side[i]  &&  layer==bad_cell_layer[i]   &&  column  == bad_cell_column[i]  ){
         output=true;}}

      return output;
}



bool is_cell_broken_anode(double side, double column, double layer){

      std::vector<double> bad_cell_column={11,77,89,101};
      std::vector<double> bad_cell_layer={0,8,1,3};
      std::vector<double> bad_cell_side={0,0,0,0};
      bool output=false;

      for (int i=0; i< bad_cell_layer.size(); i++){
      if (side==bad_cell_side[i]  &&  layer==bad_cell_layer[i]   &&  column  == bad_cell_column[i]  ){
         output=true;}}

      return output;
}




bool is_cell_bad(double side, double column, double layer){
      bool output = false;

      bool broken_anode = is_cell_broken_anode(side, column, layer);
      bool internal_cab = is_cell_internal_cabling(side, column, layer);
      bool HV_trip      = is_cell_HV_trip(side, column, layer);
      bool noisy        = is_cell_noisy(side, column, layer);
      bool other        = is_cell_other(side, column, layer);
      bool new_         = is_cell_new(side, column, layer);

      if (broken_anode == true  ||  internal_cab  == true    ||    HV_trip == true     ||   noisy == true    ||      other == true ||  new_ == true){
         output = true;
      }

      
      return output;
}




bool is_cell_bad(int i){
      bool output = false;
      int side_;
      int column_;
      int layer_;
      int cell_num;

      for (int side=0;side<2;side++){
         for (int column=0;column<113;column++){
            for (int layer=0;layer<9;layer++){
            cell_num =  side*9*113 + column*9 + layer;


      if (cell_num == i){
      side_=side;
      column_=column;
      layer_=layer;
      }

}}}
      int side=side_;
      int column=column_;
      int layer=layer_;

      bool broken_anode = is_cell_broken_anode(side, column, layer);
      bool internal_cab = is_cell_internal_cabling(side, column, layer);
      bool HV_trip      = is_cell_HV_trip(side, column, layer);
      bool noisy        = is_cell_noisy(side, column, layer);
      bool other        = is_cell_other(side, column, layer);
      bool new_         = is_cell_new(side, column, layer);

      if (broken_anode == true  ||  internal_cab  == true    ||    HV_trip == true     ||   noisy == true    ||      other == true ||  new_ == true){
         output = true;
      }

      return output;
}




bool is_next_to_bad_cell(int i){

 bool output = false;
      int side_;
      int column_;
      int layer_;
      int cell_num;

      for (int side=0;side<2;side++){
         for (int column=0;column<113;column++){
            for (int layer=0;layer<9;layer++){
            cell_num =  side*9*113 + column*9 + layer;


      if (cell_num == i){
      side_=side;
      column_=column;
      layer_=layer;
      }

}}}
      int side=side_;
      int column=column_;
      int layer=layer_;


         if (is_cell_bad(side, column+1, layer) == true ||  is_cell_bad(side, column-1, layer)  ==true  ||   is_cell_bad(side, column, layer+1) ==true   ||      is_cell_bad(side, column, layer-1)==true ){

               output = true;

         }
   return output;

}




bool is_next_to_bad_cell(double side, double column, double layer){

   bool output = false;
         if (is_cell_bad(side, column+1, layer) == true ||  is_cell_bad(side, column-1, layer)  ==true  ||   is_cell_bad(side, column, layer+1) ==true   ||      is_cell_bad(side, column, layer-1)==true ){

               output = true;

         }
   return output;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////End/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






void make_dir(string run_num){

    string folder = name_folder(run_num);
    if(mkdir((folder).c_str(), 0777) == -1);
}

void make_root_direc(string name,TFile *file){

  TDirectory *direc = file->mkdir(name.c_str());
    direc->cd();

}


TFile* make_root_file(string name,string run_num){

  string folder = name_folder(run_num);
  TFile *rootfile = new TFile((folder+"/"+name+".root").c_str(), "RECREATE");//saving all histograms
  rootfile->cd();

  return rootfile;

}





void save_all_cell_hists(TH1D** output_hists,string name,string run_num){

TFile* root_file =  make_root_file(name,run_num);

int num_cells = 2034;
for (int i=0;i<num_cells;i++){
output_hists[i]->Write();
}
root_file->Close();

cout<<"File "<<name<<" saved"<<endl;
}


void save_all_cell_hists(TH2D** output_hists,string name,string run_num){

TFile* root_file =  make_root_file(name,run_num);

int num_cells = 2034;
for (int i=0;i<num_cells;i++){
output_hists[i]->Write();
}
root_file->Close();
cout<<"File "<<name<<" saved"<<endl;
}  



void save_all_cell_hists(TCanvas** canvas,string name,string run_num){

TFile* root_file =  make_root_file(name,run_num);

int num_cells = 2034;
for (int i=0;i<num_cells;i++){
canvas[i]->Write();
}
root_file->Close();
cout<<"File "<<name<<" saved"<<endl;
}  




/*

void save_snd(std::string run_num,sndisplay::tracker *snd_prop, string name, int lower_limit, int upper_limit){
    string folder =  "/sps/nemo/scratch/spratt/analysis/vertical_position_reconstruction/run_output_analysis_files/";
    string title ="/snd_";
    string file_type =".png";
    snd_prop->setrange(lower_limit,upper_limit);
    snd_prop->draw();//saving all sndisplay plits
    snd_prop->update(false);
    snd_prop->canvas->SaveAs((folder+run_num+title+name+file_type).c_str());
}
*/



void average_over_full_2d_hists(TH2D* h_full, TH1D* output){

   
   int hist_bin_no = h_full->GetNbinsX();


   TH1D** y_projections = new TH1D*[hist_bin_no];

	for(int i=0;i<hist_bin_no;i++){

	y_projections[i]  = h_full->ProjectionY(("velocity_1D_hist_"+to_string(i+1)).c_str(),i,i+1);
	double set_value = y_projections[i]->GetMean();
	output->SetBinContent(i,set_value);
	}}



void average_over_hist(TH1D* h_full,TH1D* h_counter,TH1D* h_output){
   
   double number_bins = h_full->GetNbinsX();

	for(int i=0;i<number_bins;i++){

	double total = h_full->GetBinContent(i);
	double no = h_counter->GetBinContent(i);
	if (no>0){
	double set_value = total/no;
	h_output->SetBinContent(i,set_value);
	}}}







void average_over_hist(TH1D** h_full,TH1D** h_counter,TH1D** h_output){
   
   double number_bins = h_full[0]->GetNbinsX();

for (int j=0;j<2034;j++){

	for(int i=0;i<number_bins;i++){

	double total = h_full[j]->GetBinContent(i);
	double no = h_counter[j]->GetBinContent(i);
	if (no>0){
	double set_value = total/no;
	h_output[j]->SetBinContent(i,set_value);
	}}}}





double calc_average_of_hist(TH1D *hist,int bin_lower, int bin_higher){
  double sum=0;
  double no_bins =  hist->GetNbinsX();
  for (int i=bin_lower;i<bin_higher;i++){
  sum = sum + hist->GetBinContent(i);
  }
  double average = sum/(bin_higher-bin_lower);
  
return average;
}



void set_hist_limits(TH1D* hist, double lower, double upper){

hist->SetMaximum(8000);
hist->SetMinimum(4000);


}



/////////


TH1D** make_array_of_hist(string name_of_hist, double bins,double lower_bound,double upper_bound){
    
      int no_cells=2034;
      TH1D** h = new TH1D*[no_cells];
      for (int i=0;i<no_cells;i++){
            string cell_num = to_string(i);
            string name = (name_of_hist+cell_num);
            h[i] = new TH1D(name.c_str(),name.c_str(), bins, lower_bound, upper_bound);
                        }

      return h;
}




TH2D** make_array_of_hist(string name_of_hist, double bins,double lower_bound,double upper_bound,double bins2,double lower_bound2,double upper_bound2){
      
      int no_cells=2034;
      TH2D** h = new TH2D*[no_cells];
      for (int i=0;i<no_cells;i++){
            string cell_num = to_string(i);
            string name = (name_of_hist+cell_num);
            h[i] = new TH2D(name.c_str(),name.c_str(), bins, lower_bound, upper_bound ,bins2, lower_bound2, upper_bound2);
                        }

      return h;
}




TGraph** make_array_of_tgraph(string name_of_hist){
      
      int no_cells=2034;
      TGraph** h = new TGraph*[no_cells];
      for (int i=0;i<no_cells;i++){
            string cell_num = to_string(i);
            string name = (name_of_hist+cell_num);
            h[i] = new TGraph(0);
                        }

      return h;
}




void add_to_hist(TH1D* hist_average, TH1D* hist_counter, TH1D* hist_full_output,TH1D* hist_full_counter){

  int no_bins = hist_average->GetNbinsX();
  double average;
  double count;
  double input;
  for (int i=0;i<no_bins;i++){
       average = hist_average->GetBinContent(i);
       count = hist_counter->GetBinContent(i);
       input = average*count;
       hist_full_output->AddBinContent(i,input);
       hist_full_counter->AddBinContent(i,count);
  }}


double gaussian_fit_function(double* x, double* par){
  double X=x[0];
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = (pow(X,1) - m) / s;
  double y = h * std::exp(-0.5 * a * a);
  return y;

}


TCanvas** make_array_of_canvas(string name_of_hist){
    
      int no_cells=2034;
      TCanvas** h = new TCanvas*[no_cells];
      for (int i=0;i<no_cells;i++){
            string cell_num = to_string(i);
            string name = (name_of_hist+cell_num);
            h[i] = new TCanvas(name.c_str(),name.c_str());
                        }

      return h;
}





TH2D* rebin_2d_hist(TH2D* old,int rebin_factor_x,int rebin_factor_y){

 //the original histogram //create a new TH2 with your bin arrays spec 

  TAxis *xaxis = old->GetXaxis(); 
  TAxis *yaxis = old->GetYaxis(); 
  double xmin = xaxis->GetXmin();
  double xmax = xaxis->GetXmax();
  double ymin = yaxis->GetXmin();
  double ymax = yaxis->GetXmax();
  double ybins = yaxis->GetNbins();
  double xbins = xaxis->GetNbins();

  int new_x_bins = round(xbins/rebin_factor_x);
  int new_y_bins = round(ybins/rebin_factor_y);

  string name = old->GetTitle();

  TH2D *h = new TH2D((name+"rebin").c_str(),(name+"rebin").c_str(),new_x_bins,xmin,xmax,new_y_bins,ymin,ymax);

  for (int j=1; j<=ybins;j++) { 
      for (int i=1; i<=xbins;i++) { 
            h->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),old->GetBinContent(i,j)); } } 

   return h;

}






TH1D** generate_full_distribution(TH1D** average_hists,TH1D** counter_hists){
int no_cells = 2034;
TH1D** output_hist_array = new TH1D*[4];

TH1D* bad = new TH1D("bad","bad", 100,0,1);
TH1D* bad_counter = new TH1D("bad_counter","bad_counter", 100,0,1);
TH1D* bad_output = new TH1D("bad_output","bad_output", 100,0,1);

TH1D* edge = new TH1D("edge","edge", 100,0,1);
TH1D* edge_counter = new TH1D("edge_counter","edge_counter", 100,0,1);
TH1D* edge_output = new TH1D("edge_output","edge_output", 100,0,1);

TH1D* inner = new TH1D("inner","inner", 100,0,1);
TH1D* inner_counter = new TH1D("inner_counter","inner_counter", 100,0,1);
TH1D* inner_output = new TH1D("inner_output","inner_output", 100,0,1);

TH1D* near = new TH1D("near","near", 100,0,1);
TH1D* near_counter = new TH1D("near_counter","near_counter", 100,0,1);
TH1D* near_output = new TH1D("near_output","near_output", 100,0,1);

    for (int i=0; i<no_cells;i++){
        
        if (is_cell_bad(i)==true){
         add_to_hist(average_hists[i],counter_hists[i],bad,bad_counter);
        }

        if (is_next_to_bad_cell(i)==true && is_cell_bad(i)==false){
         add_to_hist(average_hists[i],counter_hists[i],near,near_counter);
        }

        if (check_cell_on_edge(i)==false  && is_cell_bad(i)==false && is_next_to_bad_cell(i)==false){
         add_to_hist(average_hists[i],counter_hists[i],inner,inner_counter);
        }

        if (check_cell_on_edge(i)==true  && is_cell_bad(i)==false && is_next_to_bad_cell(i)==false){
         add_to_hist(average_hists[i],counter_hists[i],edge,edge_counter);
        }}

       // set_hist_limits(inner_output,4000,8000);
       // set_hist_limits(edge_output,4000,8000);
       // set_hist_limits(near_output,4000,8000);
       // set_hist_limits(bad_output,4000,8000);

        average_over_hist(inner,inner_counter,inner_output);
        average_over_hist(edge,edge_counter,edge_output);
        average_over_hist(near,near_counter,near_output);
        average_over_hist(bad,bad_counter,bad_counter);

        output_hist_array[0]=inner_output;
        output_hist_array[1]=edge_output;
        output_hist_array[2]=near_output;
        output_hist_array[3]=bad_output;
        return output_hist_array;
        }





double straight_line_function(double* x, double* par){

   double X=x[0];
  double a=par[0];
  double b=par[1];

  double y = (a*X)+b;

  return y;

}



int cell_num_calc(double side, double column, double layer){

    int cell_num =  side*9*113 + column*9 + layer;
    return cell_num;

}

/*

void make_sndisplay_canvas(){

sndisplay::demonstrator *sndemonstrator = new sndisplay::demonstrator ("demonstrator_test");

}
*/



int calo_track_corresponder(int calo_column, int track_layer){
  if (calo_column == 0 && track_layer >= 0 && track_layer <= 6) return 1;
  if (calo_column == 1 && track_layer >= 2 && track_layer <= 11) return 1;
  if (calo_column == 2 && track_layer >= 8 && track_layer <= 17) return 1;
  if (calo_column == 3 && track_layer >= 14 && track_layer <= 24) return 1;
  if (calo_column == 4 && track_layer >= 19 && track_layer <= 29) return 1;
  if (calo_column == 5 && track_layer >= 25 && track_layer <= 34) return 1;
  if (calo_column == 6 && track_layer >= 31 && track_layer <= 40) return 1;
  if (calo_column == 7 && track_layer >= 37 && track_layer <= 46) return 1;
  if (calo_column == 8 && track_layer >= 43 && track_layer <= 52) return 1;
  if (calo_column == 9 && track_layer >= 49 && track_layer <= 58) return 1;
  if (calo_column == 10 && track_layer >= 55 && track_layer <= 64) return 1;
  if (calo_column == 11 && track_layer >= 61 && track_layer <= 70) return 1;
  if (calo_column == 12 && track_layer >= 66 && track_layer <= 75) return 1;
  if (calo_column == 13 && track_layer >= 72 && track_layer <= 81) return 1;
  if (calo_column == 14 && track_layer >= 78 && track_layer <= 87) return 1;
  if (calo_column == 15 && track_layer >= 84 && track_layer <= 93) return 1;
  if (calo_column == 16 && track_layer >= 90 && track_layer <= 99) return 1;
  if (calo_column == 17 && track_layer >= 96 && track_layer <= 105) return 1;
  if (calo_column == 18 && track_layer >= 101 && track_layer <= 110) return 1;
  if (calo_column == 19 && track_layer >= 107 && track_layer <= 112) return 1;
  else return 0;
}



void fit_gaus_function(TH1D *hist1,TF1 *fitFunction, double p0_lower, double p0_higher, double p1_lower, double p1_higher, double p2_lower, double p2_higher){

    fitFunction->SetParLimits(0,p0_lower,p0_higher);//mean
    fitFunction->SetParLimits(1,p1_lower,p1_higher);//sd
    fitFunction->SetParLimits(2,p2_lower,p2_lower);//hight

    hist1->Fit(fitFunction,"QR0"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    double par_2 = fitFunction->GetParameter(2);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "mean: " << par_0 << std::endl;
    std::cout << "sd: " << par_1 << std::endl;
    std::cout << "height_constant: " << par_2 << std::endl;

}




int what_is_crate_number(int cell_num){
double output=-1;

if (cell_num <342){
output=1;
}
if (cell_num >1017 && cell_num<1359){
output=1;
}

return output;
}




double get_tp_given_eps(Double_t *x,Double_t *par) {
    double L = 1;
    double eps = ((x[0]*0.5));

    double vel = par[0];
    double acc = par[1];

    double best_d;


    for (int i=0;i<=10000;i++){
    double d = i*0.0001;
    double eps_test_value = calc_eps(d,vel,acc);
    if ((eps_test_value-eps)>0){
      best_d = d;
     // cout<<d<<endl;
      break;
    }
    }


    double tp = calc_tp(best_d,vel,acc);

    return tp;
   }



void fit_iteration_eps_vs_z(TH1D* graph1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,0, 0.2);
    fitFunction->SetParLimits(1,0, 0.01);

    graph1->Fit(fitFunction, "R"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "vel: " << par_0 << std::endl;
    std::cout << "acc: " << par_1 << std::endl;

}
   






double calc_tp(double d,double vel,double acc){

double L=1;

double tp = (L/vel) + 2*(acc/pow(vel,3))*(pow(d,2)+pow((L-d),2));

return tp;

}


double calc_eps(double d,double vel,double acc){

double L=1;


double term1 = ((2*d)-L)*pow(vel,2);
double term2 = ((2*acc)*(pow(d,2)-pow((L-d),2)));
double term3 = (L*pow(vel,2));
double term4 = ((2*acc)*(pow(d,2)+pow((L-d),2)));

double eps = (term1+term2)/(term3+term4);

//cout<<term1<<"    "<<term2<<"    "<<term3<<"    "<<term4<<endl;

return eps;

}



