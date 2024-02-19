#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <vector>
#include <TMultiGraph.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include<TApplication.h>
#include <THStack.h>
#include <stdio.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TParameter.h>
#include "TImageDump.h"
#include "TImage.h"
#include "TMath.h"
#include "TPoint.h"
#include "TColor.h"
#include "TVirtualPad.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TText.h"
#include "RStipples.h"
#include "TList.h"
//#include <bits/stdc++.h> 
#include <iostream>
#include <chrono>
#include <ctime>   
#include "Math/ProbFunc.h"
#include <TGraph.h>
#include <random>
//#include "/sps/nemo/scratch/spratt/analysis/analysis_library.h"

using namespace std;


//////////////////////////////////////

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



double calc_tp(double d,double vel,double acc){

double L=1;

double tp = (L/vel) + 2*(acc/pow(vel,3))*(pow(d,2)+pow((L-d),2));

return tp;

}




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





string name_folder(string run_num){
  
    string directory = "/sps/nemo/scratch/spratt/analysis/timestamp_investigation/run_output_analysis_files/";
    std::string folder = directory+run_num;

    return folder;
}





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


void save_all_cell_hists(TGraph** output_hists,string name,string run_num){

TFile* root_file =  make_root_file(name,run_num);

int num_cells = 2034;
for (int i=0;i<num_cells;i++){
output_hists[i]->Write();
}
root_file->Close();

cout<<"File "<<name<<" saved"<<endl;
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




///////////////////////////////////







void average_over_hist(TH1D* h_full,TH1D* h_counter,TH1D* h_output){
   
   double number_bins = h_full->GetNbinsX();

	for(int i=0;i<number_bins;i++){

	double total = h_full->GetBinContent(i);
	double no = h_counter->GetBinContent(i);
	if (no>0){
	double set_value = total/no;
	h_output->SetBinContent(i,set_value);
	}}}





double get_tp_given_eps(Double_t *x,Double_t *par) {
    double L = 1;
    double eps = x[0]-0.5;

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




void fit_iteration_eps_vs_z(TGraph *graph1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,0, 0.1);
    fitFunction->SetParLimits(1,0, pow(10,-8));

    graph1->Fit(fitFunction, "R"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "vel: " << par_0 << std::endl;
    std::cout << "acc: " << par_1 << std::endl;

}



void fit_straight_line(TH1D *hist1,TF1 *fitFunction){

    

   // double standard_difference = pow(4,-1)-pow(5,-1);  
   // fitFunction->SetParLimits(0,pow(5,-1),pow(4,-1));//mean
   // fitFunction->SetParLimits(1,0,100);//mean
   // fitFunction->SetParLimits(3,0,0.1);//mean
    //fitFunction->SetParLimits(1,9.92*pow(10,-6), 9.920000001*pow(10,-6));//sd
    //fitFunction->SetParLimits(2,0, 1000000);//hight
    //fitFunction->SetParLimits(3,-10, 10);//phi

   // fitFunction->SetParameter(0,pow(4.3,-1));
   // fitFunction->SetParameter(1,standard_difference);
   // fitFunction->SetParameter(2,0.001);
   // fitFunction->SetParameter(3,0.01);

    hist1->Fit(fitFunction,"RQ"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);

   // cout<<"Iteration Fit"<<endl;

   // std::cout << "Fit Parameters:" << std::endl;
   // std::cout << "gradient: " << par_0 << std::endl;
   // std::cout << "intercept: " << par_1 << std::endl;

}




void fit_gaussian_function_with_shift(TH1D *hist1,TF1 *fitFunction){

    

    double standard_difference = pow(4,-1)-pow(5,-1);  
    fitFunction->SetParLimits(0,pow(5,-1),pow(4,-1));//mean
    fitFunction->SetParLimits(1,0,100);//mean
    fitFunction->SetParLimits(3,0,0.1);//mean
    //fitFunction->SetParLimits(1,9.92*pow(10,-6), 9.920000001*pow(10,-6));//sd
    //fitFunction->SetParLimits(2,0, 1000000);//hight
    //fitFunction->SetParLimits(3,-10, 10);//phi

    fitFunction->SetParameter(0,pow(4.3,-1));
    fitFunction->SetParameter(1,standard_difference);
    fitFunction->SetParameter(2,0.001);
    fitFunction->SetParameter(3,0.01);

    hist1->Fit(fitFunction,"Q","0"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    double par_2 = fitFunction->GetParameter(2);
    double par_3 = fitFunction->GetParameter(3);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "mean: " << par_0 << std::endl;
    std::cout << "sd: " << par_1 << std::endl;
    std::cout << "height_constant: " << par_2 << std::endl;
    std::cout << "shift : " << par_3 << std::endl;

}









void fit_gaussian_function(TH1D *hist1,TF1 *fitFunction){

    double height_estimate = hist1->GetMaximum();
    cout<<height_estimate<<endl;

    double standard_difference = -pow(4.1,1)+pow(4.5,1);

    //fitFunction->SetParLimits(0,pow(4.0,1),pow(4.5,1));//mean
   // fitFunction->SetParLimits(1,0, 10000);//sd
    //fitFunction->SetParLimits(2,height_estimate*0.5,height_estimate*2);//height

    fitFunction->SetParameter(0,pow(4.2,1));
    fitFunction->SetParameter(1,standard_difference*10);
    fitFunction->SetParameter(2,height_estimate);


    hist1->Fit(fitFunction,"Q","0"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    double par_2 = fitFunction->GetParameter(2);
 
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "mean: " << par_0 << std::endl;
    std::cout << "sd: " << par_1 << std::endl;
    std::cout << "height: " << par_2 << std::endl;


}



void fit_simple_gaussian_function(TH1D *hist1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,2.4,6);//mean
    fitFunction->SetParLimits(1,50, 1000);//sd
    fitFunction->SetParLimits(2,0, 1000000);//hight

    hist1->Fit(fitFunction,"Q","0"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    double par_2 = fitFunction->GetParameter(2);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "mean: " << par_0 << std::endl;
    std::cout << "sd: " << par_1 << std::endl;
    std::cout << "height_constant: " << par_2 << std::endl;

}



void fit_double_gaussian_function(TH1D *hist1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,4000,6000);
    fitFunction->SetParLimits(1,10, 1000);
    fitFunction->SetParLimits(2,10, 10000000);

    fitFunction->SetParLimits(3,4000,6000);
    fitFunction->SetParLimits(4,10, 1000);
    fitFunction->SetParLimits(5,10, 10000000);

    hist1->Fit(fitFunction,"Q"); // "R" for range-based fit

}



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




double pick_from_norm_distribution(double v_central, double width){
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::normal_distribution<double> distribution(v_central,width);
  double vel = distribution(generator);
  return vel;
}



double pick_random_height(){
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution<int> distribution(150,850);
  double height = distribution(generator);
  height = height*0.001;
  return height;
}



void fill_model_2d_hist(TH2D* model_prop_vs_z_hist,double acc_dis_width){

double v_central = 2.17*pow(10,-4);
double v_width = 1*pow(10,-5);
double height = pick_random_height();
double vel    = pick_from_norm_distribution(v_central,v_width);
double k    = pick_from_norm_distribution(5.10710,5.10710*acc_dis_width);
double acc = k*pow(10,-6)*vel;
double eps = calc_eps(height,vel,acc);
double tp =calc_tp(height,vel,acc);
model_prop_vs_z_hist->Fill(eps,tp,1);

//cout<<"  height: "<<height<<"  vel: "<<vel<<"  eps: "<<eps<<"  tp: "<<tp<<endl;

}


void fill_model_1d_hist(TH1D* model_prop_1D_hist,double height,double acc_dis_width,double acc_factor){

double v_central = 2.17*pow(10,-4);
double acc_cent = 5.10710*acc_factor*pow(10,-6)*v_central;
double v_width = 2.5*pow(10,-6);
height = height*0.05;
double vel    = pick_from_norm_distribution(v_central,v_width);
double acc    = pick_from_norm_distribution(acc_cent,acc_cent*acc_dis_width);
double eps = calc_eps(height,vel,acc);
double tp =calc_tp(height,vel,acc);
model_prop_1D_hist->Fill(tp,1);

}



double gaussian_fit_function_for_plot(double X, double m, double s, double h, double shift){
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = ((pow((X-shift),-1)) - pow(m,-1)) *s;
  double y = h*(inv_sqrt_2pi *s) * std::exp(-0.5f * a * a);
  return y;
}

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

double exponential_function(double* x, double* par){
  double X=x[0];
  double A=par[0];
  double B=par[1]; 
  double C=par[2]; 

  double y = (A*X*X)+(B*X)+C;

  return y;

}




double gaussian_fit_function_new(double* x, double* par){
  double X=x[0];
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double phi=par[3];
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = ((((1+phi)/(X)) - (m))/s);
  double y = h*(inv_sqrt_2pi) * std::exp(-0.5f * a * a);
  return y;

}



double simple_gaussian_fit_function(double* x, double* par){
  double X=x[0];
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = ((X) - (m)) / s;
  double y = h*(inv_sqrt_2pi) * std::exp(-0.5f * a * a);
  return y;

}



double gaussian_fit_function_with_shift(double* x, double* par){
  double X=x[0];
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double shift=par[3];
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = ((pow((X-shift),-1)) - (m)) / s;
  double y = h*(inv_sqrt_2pi) * std::exp(-0.5f * a * a);
  return y;

}


double double_gaussian_fit_function(double* x, double* par){
  double X=x[0];
  double m1=par[0];
  double s1=par[1];
  double h1=par[2];

  double m2=par[3];
  double s2=par[4];


  double h2=par[5];

  double inv_sqrt_2pi = 0.3989422804014327;
  double a1 = (X - m1) / s1;
  double a2 = (X - m2) / s2;
  double y1 = h1*(inv_sqrt_2pi / s1) * std::exp(-0.5f * a1 * a1);
  double y2 = h2*(inv_sqrt_2pi / s2) * std::exp(-0.5f * a2 * a2);

  double y_output = y1+y2;

  return y_output;

}

//adds 2d hists for each cell from input file specified in code to the input histogram
void combine_runs_straight_line_fit(int run_num, TH1D** prop_vs_z_1D_hist_array,string file, string name_of_hists){


cout<<"Combining Run "<<run_num<<endl;

string folder = name_folder(to_string(run_num).c_str());
string general_inves = (folder+"/general_inves.root").c_str();

string prop_time_vs_z_2d = (folder+"/individual_cell_data/"+file).c_str();
TFile *prop_time_vs_z_2d_file = new TFile(prop_time_vs_z_2d.c_str(), "READ");

int no_cells=2034;
TH1D** heatmap_hists = new TH1D*[no_cells];
for (int i=0;i<no_cells;i++){
heatmap_hists[i] =  (TH1D*)prop_time_vs_z_2d_file->Get((name_of_hists+to_string(i)).c_str());
prop_vs_z_1D_hist_array[i]->Add(heatmap_hists[i]);

//cout<<prop_vs_z_1D_hist_array[i]->GetMean()<<endl;


}
}


//adds 2d hists for each cell from input file specified in code to the input histogram
void combine_runs(int run_num, TH2D** prop_vs_z_2D_hist_array){

cout<<"Combining Run "<<run_num<<endl;

string folder = name_folder(to_string(run_num).c_str());

string prop_time_vs_z_2d = (folder+"/individual_cell_data/bismuth_loc.root").c_str();
TFile *prop_time_vs_z_2d_file = new TFile(prop_time_vs_z_2d.c_str(), "READ");

int no_cells=2034;
TH2D** heatmap_hists = new TH2D*[no_cells];
for (int i=0;i<no_cells;i++){
heatmap_hists[i] =  (TH2D*)prop_time_vs_z_2d_file->Get(("bismuth_loc"+to_string(i)).c_str());
prop_vs_z_2D_hist_array[i]->Add(heatmap_hists[i]);
}
}




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




TH2D** make_array_of_hist_2d(string name_of_hist, double bins,double lower_bound,double upper_bound,double bins2,double lower_bound2,double upper_bound2){
      
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






void comapare_runs_dummy_code(){


TH2D** heatmap1 = make_array_of_hist_2d("prop_vs_z_2D_hist_array1",100,-1,1,1000,3000,7000);
TH1D** output1 =  make_array_of_hist("prop_vs_z_1D_hist_output_array1",100,-1,1);
TH1D** counter1 = make_array_of_hist("prop_vs_z_1D_hist_counter_array1",100,-1,1);

TH2D** heatmap2 = make_array_of_hist_2d("prop_vs_z_2D_hist_array2",100,-1,1,1000,3000,7000);
TH1D** output2 =  make_array_of_hist("prop_vs_z_1D_hist_output_array2",100,-1,1);
TH1D** counter2 = make_array_of_hist("prop_vs_z_1D_hist_counter_array2",100,-1,1);

TH2D** heatmap3 = make_array_of_hist_2d("prop_vs_z_2D_hist_array3",100,-1,1,1000,3000,7000);
TH1D** output3 =  make_array_of_hist("prop_vs_z_1D_hist_output_array3",100,-1,1);
TH1D** counter3 = make_array_of_hist("prop_vs_z_1D_hist_counter_array3",100,-1,1);

TH2D** heatmap4 = make_array_of_hist_2d("prop_vs_z_2D_hist_array4",100,-1,1,1000,3000,7000);
TH1D** output4 =  make_array_of_hist("prop_vs_z_1D_hist_output_array4",100,-1,1);
TH1D** counter4 = make_array_of_hist("prop_vs_z_1D_hist_counter_array4",100,-1,1);

//combine_runs(1051,heatmap1,output1,counter1);
//combine_runs(1054,heatmap2,output2,counter2);
//combine_runs(1058,heatmap3,output3,counter3);
//combine_runs(1062,heatmap4,output4,counter4);




TH1D** y_projection_1 = make_array_of_hist("y_projection_1",1000,3000,7000);
TH1D** y_projection_2 = make_array_of_hist("y_projection_2",1000,3000,7000);
TH1D** y_projection_3 = make_array_of_hist("y_projection_3",1000,3000,7000);
TH1D** y_projection_4 = make_array_of_hist("y_projection_4",1000,3000,7000);

for (int i=0;i<2034;i++){
for (int j=0;j<1000;j++){
      y_projection_1[i]->Add(heatmap1[i]->ProjectionY(("vel1"+to_string(j+1)).c_str(),j,j+1));
      y_projection_2[i]->Add(heatmap2[i]->ProjectionY(("vel2"+to_string(j+1)).c_str(),j,j+1));
      y_projection_3[i]->Add(heatmap3[i]->ProjectionY(("vel3"+to_string(j+1)).c_str(),j,j+1));
      y_projection_4[i]->Add(heatmap4[i]->ProjectionY(("vel4"+to_string(j+1)).c_str(),j,j+1));
}
y_projection_1[i]->Rebin(4);
y_projection_1[i]->Scale(1/(y_projection_1[i]->Integral()));

y_projection_2[i]->Rebin(4);
y_projection_2[i]->Scale(1/(y_projection_2[i]->Integral()));

y_projection_3[i]->Rebin(4);
y_projection_3[i]->Scale(1/(y_projection_3[i]->Integral()));

y_projection_4[i]->Rebin(4);
y_projection_4[i]->Scale(1/(y_projection_4[i]->Integral()));


string run_num = "10";
string folder = name_folder(run_num);
string title_root_file = "/fits.root";
TFile *rootfile = new TFile((folder+title_root_file).c_str(), "RECREATE");//saving all histograms
rootfile->cd();

TCanvas** compare_distrs = new TCanvas*[2034];
for (int i=0;i<2034;i++){
compare_distrs[i] = new TCanvas(("compare_distrs"+to_string(i)).c_str(),("compare_distrs"+to_string(i)).c_str());

y_projection_1[i]->SetLineColor(kBlue);
y_projection_2[i]->SetLineColor(kRed);
y_projection_3[i]->SetLineColor(kOrange);
y_projection_4[i]->SetLineColor(kGreen);

y_projection_1[i]->Draw();
y_projection_2[i]->Draw("Same");
y_projection_3[i]->Draw("Same");
y_projection_4[i]->Draw("Same");

compare_distrs[i]->Update();
compare_distrs[i]->Write();
}

rootfile->Close();


}}







TF1** fit_straight_line_to_average_plot(TH1D** hist_array){


TF1** fit_staright_line_array = new TF1*[2034];
for (int i =0;i<2034;i++){

      string title = hist_array[i]->GetName();

    //  cout<<hist_array[i]->GetMean()<<endl;

fit_staright_line_array[i] = new TF1(("straight_line_fit"+title+to_string(i+1)).c_str(), straight_line_function,0, 2000, 2);
fit_straight_line(hist_array[i],fit_staright_line_array[i]);

}
return fit_staright_line_array;
}






TGraph** plot_grad_against_height(TF1** fit_staright_line_array2,TF1**fit_staright_line_array3,TF1**fit_staright_line_array4,TF1**fit_staright_line_array5){


double grad_2;
double grad_3;
double grad_4;
double grad_5;


TGraph** graph_array = make_array_of_tgraph("grad_vs_height");

for (int i=0;i<2034;i++){

      TAxis *axis = graph_array[i]->GetXaxis();
   axis->SetLimits(0,1.1);    

grad_2 = fit_staright_line_array2[i]->GetParameter(0);
grad_3 = fit_staright_line_array3[i]->GetParameter(0);
grad_4 = fit_staright_line_array4[i]->GetParameter(0);
grad_5 = fit_staright_line_array5[i]->GetParameter(0);


cout<<grad_2<<"  "<<grad_3<<"  "<<grad_4<<"  "<<grad_5<<endl;

graph_array[i]->SetPoint(graph_array[i]->GetN(),0,grad_2);
graph_array[i]->SetPoint(graph_array[i]->GetN(),1,grad_3);
graph_array[i]->SetPoint(graph_array[i]->GetN(),2,grad_4);
graph_array[i]->SetPoint(graph_array[i]->GetN(),3,grad_5);



}
return graph_array;
}



















/////////////////////////////////////Main Function//////////////////////////////
int main() {



int no_cells=2034;

////////////////////////////////////Name Files//////////////////////////////////
string run_num = "1048";
string folder = name_folder(run_num);
/////////////////////////////////////////////////////////////////////////////////


///////////////////////////Load File one////////////////////////////////////
TH2D** heatmap = make_array_of_hist_2d("bismuth_tracks_2d_",1000,-5000,5000,1000,-1,1);

combine_runs(1048,heatmap);
combine_runs(1055,heatmap);
combine_runs(1059,heatmap);
combine_runs(1063,heatmap);

//loop
//////

TH2D* hist_2d = heatmap[217];
//int rebin_factor = 2;
//int no_bins_1d_hist = round((1000/rebin_factor));
//hist_2d->Rebin2D(5,rebin_factor);

//TH1D** y_projections = new TH1D*[hist_2d->GetNbinsX()];
//TH1D** y_projections_scaled = new TH1D*[hist_2d->GetNbinsX()];

//TF1** fit_gaussian_functions_array = new TF1*[hist_2d->GetNbinsX()];
//TF1** fit_gaussian_functions_array_with_shift = new TF1*[hist_2d->GetNbinsX()];

TH1D* y_projection = hist_2d->ProjectionY(("2d_bismuth"+to_string(100)).c_str(),370,380);

TF1** gaus_fits = new TF1*[5];
double rangel[5]={-0.7,-0.4,-0.1,0.2,0.5};
double rangeu[5]={-0.5,-0.2,0.1,0.4,0.7};

for (int i=0;i<5;i++){
gaus_fits[i] = new TF1(("gaus_fits"+to_string(i)).c_str(),"gaus",rangel[i],rangeu[i]);
y_projection->Fit(gaus_fits[i],"RQN");
}

double mean[5]={0};

for (int i =0;i<5;i++){
mean[i] = gaus_fits[i]->GetParameter(1);

cout<<mean[i]<<endl;
}




//y_projections_scaled[i] = new TH1D(("PDF_tp_z_value_"+to_string(i+1)).c_str(),("PDF_tp_z_value_"+to_string(i+1)).c_str(),no_bins_1d_hist,3,7);
//y_projections[i]->Scale(1/(y_projections[i]->Integral()));

//for (int j=0;j<no_bins_1d_hist;j++){
//y_projections_scaled[i]->SetBinContent(j,y_projections[i]->GetBinContent(j));
//}



//TF1** fit_staright_line_array = fit_straight_line_to_average_plot(hist_array_1d);
//double grad = fit_staright_line_array[i]->GetParameter(0);
//double intercept = fit_staright_line_array[i]->GetParameter(1);

//if (grad>0 && grad<1 && intercept> 10000){

//grad_vs_prop->SetPoint(grad_vs_prop->GetN(),intercept,grad);
//grad_vs_intercept_counter->Fill(intercept,1);
//grad_vs_intercept_full->Fill(intercept,grad);
//}

//// 
/*
for (int i=0;i<grad_vs_intercept_counter->GetNbinsX();i++){

      if (grad_vs_intercept_counter->GetBinContent(i)>0){
double output_value = grad_vs_intercept_full->GetBinContent(i)/grad_vs_intercept_counter->GetBinContent(i);
grad_vs_intercept_output->SetPoint(grad_vs_intercept_output->GetN(),grad_vs_intercept_counter->GetBinCenter(i),output_value);}

}
*/


/////////////////////////Load new root file///////////////////////////////////////
string title_root_file = "/fits_bismuth.root";
TFile *rootfile = new TFile((folder+title_root_file).c_str(), "RECREATE");//saving all histograms
rootfile->cd();
//////////////////////////////////////////////////////////////////////////////////

TCanvas* can1 = new TCanvas("1");
can1->cd();

y_projection->Draw();
for (int i=0;i<5;i++){
gaus_fits[i]->Draw("Same");
}

can1->Write();


///////////////////////////Close Files////////////////////////////////////////////
rootfile->Close();
//////////////////////////////////////////////////////////////////////////////////

  return 0;
}















