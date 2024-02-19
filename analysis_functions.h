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
#include <iostream>
#include <chrono>
#include <ctime>   
#include "Math/ProbFunc.h"
#include <TGraph.h>
#include <random>
//#include "/sps/nemo/scratch/spratt/sndisplay/sndisplay.cc"
using namespace std;


bool open_file(const char *filename);

void combine_runs(int run_num, TH2D** prop_vs_z_2D_hist_array,TH1D** prop_vs_z_1D_hist_sum_array,TH1D** prop_vs_z_1D_hist_counter_array);

string name_folder(string run_num);

///////////Check cell type

bool check_cell_on_edge(int i);

bool is_cell_new(double side, double column, double layer);

bool is_cell_other(double side, double column, double layer);

bool is_cell_noisy(double side, double column, double layer);

bool is_cell_HV_trip(double side, double column, double layer);

bool is_cell_internal_cabling(double side, double column, double layer);

bool is_cell_broken_anode(double side, double column, double layer);

bool is_cell_bad(int i);

bool is_cell_bad(double side, double column, double layer);

bool is_next_to_bad_cell(int i);

bool is_next_to_bad_cell(double side, double column, double layer);


////////////////////



//////

void make_dir(string run_num);

void make_root_direc(string name,TFile *file);

TFile* make_root_file(string name,string run_num);

void save_all_cell_hists(TH1D** output_hists,string name,string run_num);

void save_all_cell_hists(TH2D** output_hists,string name,string run_num);




void average_over_full_2d_hists(TH2D* h_full, TH1D* output);

void average_over_hist(TH1D* h_full,TH1D* h_counter,TH1D* h_output);

void average_over_hist(TH1D** h_full,TH1D** h_counter,TH1D** h_output);

double calc_average_of_hist(TH1D *hist,int bin_lower, int bin_higher);

void set_hist_limits(TH1D* hist, double lower, double upper);

TH1D** make_array_of_hist(string name_of_hist, double bins,double lower_bound,double upper_bound);

TH2D** make_array_of_hist(string name_of_hist, double bins,double lower_bound,double upper_bound,double bins2,double lower_bound2,double upper_bound2);


TGraph** make_array_of_tgraph(string name_of_hist);

void make_sndisplay_canvas();

void save_all_cell_hists(TCanvas** canvas,string name,string run_num);

void add_to_hist(TH1D* hist_average, TH1D* hist_counter, TH1D* hist_full_output,TH1D* hist_full_counter);

TCanvas** make_array_of_canvas(string name_of_hist);

TH2D* rebin_2d_hist(TH2D* old,int rebin_factor_x,int rebin_factor_y);

TH1D** generate_full_distribution(TH1D** average_hists,TH1D** counter_hists);

double straight_line_function(double* x, double* par);

int cell_num_calc(double side, double column, double layer);

//void make_sndisplay_canvas();

//void save_snd(std::string run_num,sndisplay::tracker *snd_prop, string name, int lower_limit, int upper_limit);

int calo_track_corresponder(int calo_column, int track_layer);

double gaussian_fit_function(double* x, double* par);

void fit_gaus_function(TH1D *hist1,TF1 *fitFunction, double p0_lower, double p0_higher, double p1_lower, double p1_higher, double p2_lower, double p2_higher);


int what_is_crate_number(int cell_num);

double get_tp_given_eps(Double_t *x,Double_t *par);

double calc_tp(double d,double vel,double acc);

double calc_eps(double d,double vel,double acc);

void fit_iteration_eps_vs_z(TH1D* graph1,TF1 *fitFunction);

