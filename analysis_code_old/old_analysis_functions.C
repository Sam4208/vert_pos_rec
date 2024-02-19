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


///////////////////////////////////



double get_tp_given_eps(Double_t *x,Double_t *par) {
    double L = 1;
    double eps = (x[0]*0.5)+0.5;

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





void plot_runs(int run_num, TCanvas** canvas){

cout<<"Combining Run "<<run_num<<endl;

string folder = name_folder(to_string(run_num).c_str());

string prop_time_vs_z_sum = (folder+"/export_analysis/prop_time_hist_cell.root").c_str();
TFile *prop_time_vs_z_sum_file = new TFile(prop_time_vs_z_sum.c_str(), "READ");

int no_cells=2034;

for (int i=0;i<no_cells;i++){
canvas[i]->cd();

TH2* hist = ((TH2D*)prop_time_vs_z_sum_file->Get(("prop_time_hist"+to_string(i)).c_str()))->Rebin(10);

 //  int maxBin = hist->GetMaximumBin();
   

    // Get the bin center of the maximum Y value
    //double maxX = hist->GetBinCenter(maxBin);
    //double maxY = hist->GetBinContent(maxBin);

    hist->Scale(1/(hist->GetBinContent(hist->GetMaximumBin())));
 
    hist->Draw("Same,L");

    // Display the canvas


    //canvas[i]->Update();

    // Clean up memory
  //  hist->Delete();
   // text->Delete();


}
}





