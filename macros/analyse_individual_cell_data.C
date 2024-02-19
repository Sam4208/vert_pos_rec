#include "../analysis_functions.C"

//////////////////////////////////////



/////////////////////////////////////Main Function//////////////////////////////
int main() {



int no_cells=2034;

////////////////////////////////////Name Files//////////////////////////////////
string run_num = "combined";
string folder = name_folder(run_num);


cout<<"Start"<<endl;
/////////////////////////////////////////////////////////////////////////////////


///////////////////////////Load File one////////////////////////////////////

TH1D** hists_sum =  make_array_of_hist("sum",50,-1,+1);
TH1D** hists_counter =  make_array_of_hist("counter",50,-1,+1);
TH1D** hists_output =  make_array_of_hist("output",50,-1,+1);
TH2D** hists_2d =  make_array_of_hist("2d",50,-1,+1,1000,30000,80000);


//10hr non-source runs
combine_runs(1046,hists_2d,hists_sum,hists_counter);
combine_runs(1051,hists_2d,hists_sum,hists_counter);
combine_runs(1054,hists_2d,hists_sum,hists_counter);
combine_runs(1062,hists_2d,hists_sum,hists_counter);
combine_runs(1058,hists_2d,hists_sum,hists_counter);
combine_runs(1087,hists_2d,hists_sum,hists_counter);
combine_runs(1093,hists_2d,hists_sum,hists_counter);
combine_runs(1100,hists_2d,hists_sum,hists_counter);
combine_runs(1108,hists_2d,hists_sum,hists_counter);
combine_runs(1115,hists_2d,hists_sum,hists_counter);
combine_runs(1122,hists_2d,hists_sum,hists_counter);
combine_runs(1129,hists_2d,hists_sum,hists_counter);
combine_runs(1184,hists_2d,hists_sum,hists_counter);
combine_runs(1186,hists_2d,hists_sum,hists_counter);
combine_runs(1191,hists_2d,hists_sum,hists_counter);
combine_runs(1197,hists_2d,hists_sum,hists_counter);
combine_runs(1202,hists_2d,hists_sum,hists_counter);
combine_runs(1207,hists_2d,hists_sum,hists_counter);

/*Bismuth runs
combine_runs(1205,hists_2d,hists_sum,hists_counter);
//combine_runs(1048,hists_2d,hists_sum,hists_counter);
//combine_runs(1055,hists_2d,hists_sum,hists_counter);
//combine_runs(1059,hists_2d,hists_sum,hists_counter);
//combine_runs(1063,hists_2d,hists_sum,hists_counter);
combine_runs(1097,hists_2d,hists_sum,hists_counter);
combine_runs(1105,hists_2d,hists_sum,hists_counter);
combine_runs(1112,hists_2d,hists_sum,hists_counter);
combine_runs(1119,hists_2d,hists_sum,hists_counter);
combine_runs(1126,hists_2d,hists_sum,hists_counter);
combine_runs(1133,hists_2d,hists_sum,hists_counter);
combine_runs(1149,hists_2d,hists_sum,hists_counter);
combine_runs(1156,hists_2d,hists_sum,hists_counter);
combine_runs(1170,hists_2d,hists_sum,hists_counter);
combine_runs(1175,hists_2d,hists_sum,hists_counter);
combine_runs(1180,hists_2d,hists_sum,hists_counter);

combine_runs(1190,hists_2d,hists_sum,hists_counter);
combine_runs(1195,hists_2d,hists_sum,hists_counter);
combine_runs(1200,hists_2d,hists_sum,hists_counter);
*/




average_over_hist(hists_sum, hists_counter,hists_output);

////////////////////////////////////////////////////////////////////////////////

save_all_cell_hists(hists_sum,"sum_prop_vs_z_cell",run_num);
save_all_cell_hists(hists_counter,"counter_prop_vs_z_cell",run_num);
save_all_cell_hists(hists_output,"average_prop_vs_z_cell",run_num);

///////////////////////////Close Files////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////
//}
  return 0;
}















