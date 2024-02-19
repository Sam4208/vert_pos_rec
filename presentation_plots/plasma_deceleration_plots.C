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
#include <vector>
//#include "/sps/nemo/scratch/spratt/analysis/analysis_library.h"

using namespace std;



double calc_d1_exp(double t1){

double k = -1*pow(10,2);
double v0 = 6.024*pow(10,-5);
double a0 = 5*pow(10,-4)*pow(v0,2);

double term1 = a0*pow(t1,3)*pow(6*k,-1);
double term2 = a0*pow(t1,2)*0.5;
double term3 = -v0*t1;

double d1 = -term1-term2-term3;

return d1;

}


double calc_d2_exp(double t2, double t1, double d1){

double k = -1*pow(10,2);
double v0 = 6.024*pow(10,-5);
double a0 = 5*pow(10,-4)*pow(v0,2)*(1-(t1*pow(k,-1)));

double term1 = a0*(pow(t2,3)-pow(t1,3))*pow(12*k,-1);
double term2 = 0.25*a0*(pow(t2,2)-pow(t1,2));
double term3 = -(v0-(a0*t1)-(a0*pow(t1,2)*pow(2*k,-1)))*(t2-t1);

double d2 = d1 - term1-term2-term3;

return d2;

}





double find_closest(vector<double> d_array, vector<double> t_array){

int element =0;

    double min_x,min_y,max_x,max_y;
    auto temp = std::minmax_element(d_array.begin(),d_array.end());
    min_x = *(temp.first);
    max_x = *(temp.second);


double min = min_x;
//cout<<d_array.size()<<endl;
for (int i=0;i<d_array.size();i++){
//cout<<d_array.at(i)<<endl;
if (d_array.at(i)==min){
element =i;

}

}

double time = t_array.at(element);
//cout<<time<<"---"<<endl;

return time;
}


void calc_t1_t2(double d1, double d2, double& t1, double& t2){

vector<double> t1_test_array={10000};
vector<double> t2_test_array={10000};
vector<double> d1_test_array={10000};
vector<double> d2_test_array={10000};
//cout<<"--------"<<endl;

for (int i=1;i<60000;i++){

double t1_test = i*1;

double d1_test = calc_d1_exp(t1_test);

t1_test_array.push_back(t1_test);
d1_test_array.push_back(abs(d1_test-d1));
//cout<<(d1_test)<<"    "<<d1<<endl;
}
t1 = find_closest(d1_test_array,t1_test_array);


for (int i=1;i<60000;i++){

double t2_test = i*1;

double d2_test = calc_d2_exp(t2_test,t1,d1);

t2_test_array.push_back(t2_test);
d2_test_array.push_back(abs(d2_test-d2));

}
t2 = find_closest(d2_test_array,t2_test_array);
//cout<<t1<<"   "<<t2<<endl;

}






void calc_t_expo_acc(double d,double& t1, double& t2){


if(d<1.5){
double d1=d;
double d2=3-d;
calc_t1_t2(d1,d2,t1,t2);

}



if(d>1.5){
double d1=3-d;
double d2=d;
calc_t1_t2(d1,d2,t2,t1);

}

//cout<<t1<<"   "<<t2<<endl;

}






















double calc_t_dec(double d){

double v=6.024*pow(10,-5);
double a=0.005367003367*pow(v,2);
   
   double term1 = pow(v,2)-(2*(a*d));
   double term2 = v-sqrt(term1);
   double time = term2/a;
 //  double     time = (abs(d)/v) + ((2*a)/pow(v,3))*pow(d,2);

     return time;

}


double calc_t_con(double d){

double v= 6.024*pow(10,-5);;
   
   double     time = (abs(d)/v);

     return time;

}



double calc_t_nemo_3(double d){

double v=6.024*pow(10,-5);
double acc=0.005367003367*pow(v,2);
double time=0;
if (d>1.5){
double d2=3-d;
double t1 = calc_t_dec(d2);
 
double a = 0.25*acc;
double b = (((acc*0.5)*t1)-v);
double c = (t1*(v-((0.75)*acc*t1)))+d-d2;

double term1 = -b-sqrt(pow(b,2)-(4*a*c));

double t2 = term1*(pow(2*a,-1));
time = t2;
}

if (d<1.5){
double d1=d;
double t1 = t1 = calc_t_dec(d1);
time = t1;
}


//cout<<time<<endl;
return time;
}






int main(){



TH1D* t6_dec_vs_z = new TH1D("histdec","histdec",300,-1.5,1.5);
TGraph* t5_dec_vs_z = new TGraph(0);
TH1D* t6_nemo_vs_z = new TH1D("histnemo","histnemo",300,-1.5,1.5);
TGraph* t5_nemo_vs_z = new TGraph(0);

TGraph* prop_exp_vs_z_exp = new TGraph(0);
TGraph* z_exp_vs_z_true = new TGraph(0);

TGraph* t5_vs_t6_con = new TGraph(0);
TGraph* t5_vs_t6_nemo_3 = new TGraph(0);
TGraph* t5_vs_t6_dec = new TGraph(0);
TGraph* prop_a2_vs_h_a2 = new TGraph(0);
TGraph* prop_meas_vs_z_meas = new TGraph(0);
TGraph* prop_cons_vs_z_cons = new TGraph(0);
TGraph* prop_meas_vs_z_true = new TGraph(0);
TGraph* z_measured_min_z_true_vs_z_true = new TGraph(0);
TGraph* z_measured_min_z_true_vs_z_true_nemo_3 = new TGraph(0);

for (int i=0;i<300;i++){

double d=i*0.01;



double t5_exp=0;
double t6_exp=0;
calc_t_expo_acc(d,t5_exp,t6_exp);

//cout<<t5_exp<<"   "<<t6_exp<<endl;





double t5_nemo_3 = calc_t_nemo_3(d);
double t6_nemo_3 = calc_t_nemo_3(3-d);
double tp_nemo_3 = t5_nemo_3+t6_nemo_3;




double t5_dec = calc_t_dec(d);
double t6_dec = calc_t_dec(3-d);


double t5_con = calc_t_con(d);
double t6_con = calc_t_con(3-d);

double z_exp = 1.5*(t5_exp-t6_exp)/(t5_exp+t6_exp);

cout<<"t5 "<<t5_exp<<"      t6 "<<t6_exp<<"     h "<<z_exp<<endl;

double z_nemo_3 = 1.5*(t5_nemo_3-t6_nemo_3)/(t5_nemo_3+t6_nemo_3);
double z_dec = 1.5*(t5_dec-t6_dec)/(t5_dec+t6_dec);
double z_con = 1.5*(t5_con-t6_con)/(t5_con+t6_con);
double z_true = 1.5*(d-(3-d))/(3);

double tp_dec = t5_dec+t6_dec;
double tp_con = t5_con+t6_con;
double tp_exp = t5_exp+t6_exp;


if (z_exp>-10 && z_exp<10){
prop_exp_vs_z_exp->SetPoint(prop_exp_vs_z_exp->GetN(),z_exp,tp_exp);
z_exp_vs_z_true->SetPoint(z_exp_vs_z_true->GetN(),z_true,z_exp-z_true);
}




prop_meas_vs_z_meas->SetPoint(prop_meas_vs_z_meas->GetN(),z_dec,tp_dec);
prop_cons_vs_z_cons->SetPoint(prop_cons_vs_z_cons->GetN(),z_con,tp_con);
prop_meas_vs_z_true->SetPoint(prop_meas_vs_z_true->GetN(),z_true,tp_dec);
z_measured_min_z_true_vs_z_true->SetPoint(z_measured_min_z_true_vs_z_true->GetN(),z_true,z_dec-z_true);

if (z_nemo_3>-10 && z_nemo_3<10){
t5_vs_t6_dec->SetPoint(t5_vs_t6_dec->GetN(),t5_dec-t5_con,t6_dec-t6_con);
t5_vs_t6_nemo_3->SetPoint(t5_vs_t6_nemo_3->GetN(),t5_nemo_3-t5_con,t6_nemo_3-t6_con);
t5_vs_t6_con->SetPoint(t5_vs_t6_con->GetN(),t5_con,t6_con);


t6_dec_vs_z->SetBinContent(i,t6_dec);
t5_dec_vs_z->SetPoint(t5_dec_vs_z->GetN(),z_con,t5_dec-t5_con);
t6_nemo_vs_z->SetBinContent(i,t6_nemo_3);
t5_nemo_vs_z->SetPoint(t5_nemo_vs_z->GetN(),z_con,t5_nemo_3-t5_con);


z_measured_min_z_true_vs_z_true_nemo_3->SetPoint(z_measured_min_z_true_vs_z_true_nemo_3->GetN(),z_true,z_nemo_3-z_true);
prop_a2_vs_h_a2->SetPoint(prop_a2_vs_h_a2->GetN(),z_nemo_3,tp_nemo_3);}
}


/////fits!!!!!

TGraph* t6_nemo_vs_z_attempt = new TGraph(0);
TGraph* t6_dec_vs_z_attempt = new TGraph(0);

TF1* straight_line_fit = new TF1("straight_line_fit","([0]*x)+[1]",0.5*1.5, 0.8*1.5);
t6_nemo_vs_z->Fit(straight_line_fit,"R");
for (int i=30;i<t6_nemo_vs_z->GetNbinsX()-30;i++){
double a =straight_line_fit->GetParameter(0);
double b =straight_line_fit->GetParameter(1);
double x = t6_nemo_vs_z->GetBinCenter(i);
if (t6_nemo_vs_z->GetBinContent(i)>0){
t6_nemo_vs_z_attempt->SetPoint(t6_nemo_vs_z_attempt->GetN(),x,t6_nemo_vs_z->GetBinContent(i)-((a*x)+b));}
}

TF1* straight_line_fit2 = new TF1("straight_line_fit","([0]*x)+[1]",0.5*1.5, 0.8*1.5);
t6_dec_vs_z->Fit(straight_line_fit2,"R");
for (int i=30;i<t6_nemo_vs_z->GetNbinsX()-30;i++){
double a =straight_line_fit2->GetParameter(0);
double b =straight_line_fit2->GetParameter(1);
double x = t6_nemo_vs_z->GetBinCenter(i);
if (t6_dec_vs_z->GetBinContent(i)>0){
t6_dec_vs_z_attempt->SetPoint(t6_dec_vs_z_attempt->GetN(),x,t6_dec_vs_z->GetBinContent(i)-((a*x)+b));}


}

TFile *rootfile = new TFile("/sps/nemo/scratch/spratt/analysis/timestamp_investigation/presentation_plots/plasma_deceleration_plots.root", "RECREATE");//saving all histograms

rootfile->cd();
//////////////////////////////////////////////////////////////////////////////////

t6_dec_vs_z_attempt->Write("t6_dec_attempt");
t6_nemo_vs_z_attempt->Write("t6_nemo_attempt");

t6_nemo_vs_z->Write("t6_vs_z_nemo");
t5_nemo_vs_z->Write("t5_vs_z_nemo");

t6_dec_vs_z->Write("t6_vs_z_dec");
t5_dec_vs_z->Write("t5_vs_z_dec");

prop_exp_vs_z_exp->Write("exp");
z_exp_vs_z_true->Write("exp");

prop_meas_vs_z_meas->Write();
prop_cons_vs_z_cons->Write();
prop_meas_vs_z_true->Write();
z_measured_min_z_true_vs_z_true->Write("z_vs_z_true");
z_measured_min_z_true_vs_z_true_nemo_3->Write("z_vs_z_nemo_3");
prop_a2_vs_h_a2->Write("prop_nemo_3");



TCanvas* canvas0 = new TCanvas("t6_vs_z_attempt");
  t6_dec_vs_z_attempt->Draw();
	t6_nemo_vs_z_attempt->Draw("Same");
	canvas0->Update();
	canvas0->Draw();
	canvas0->Write();

TCanvas* canvas1 = new TCanvas("t5_vs_z");
  t6_nemo_vs_z->Draw();
	t6_dec_vs_z->Draw("Same");
	canvas1->Update();
	canvas1->Draw();
	canvas1->Write();


TCanvas* canvas3 = new TCanvas("z_vs_z");
  z_measured_min_z_true_vs_z_true->Draw();
	z_measured_min_z_true_vs_z_true_nemo_3->Draw("Same");
	canvas3->Update();
	canvas3->Draw();
	canvas3->Write();

TCanvas* canvas2 = new TCanvas("t_vs_t");
  t5_vs_t6_nemo_3->SetLineColor(kOrange);
  t5_vs_t6_dec->SetLineColor(kRed);
  t5_vs_t6_dec->Draw();
	t5_vs_t6_nemo_3->Draw("Same");
 // t5_vs_t6_con->Draw("Same");
	canvas2->Update();
	canvas2->Draw();
	canvas2->Write();


///////////////////////////Close Files////////////////////////////////////////////
rootfile->Close();

    return 1;
}

