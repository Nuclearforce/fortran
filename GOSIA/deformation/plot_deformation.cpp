#include <TROOT.h>
#include <iostream>
#include <vector>
#include <cstdio>
#include <stdio.h>

void plot_deformation()
{
    TCanvas *c1 = new TCanvas("c1","",900,900);
    
	/*double Q2[8]={1.6681,2.2462,1.9842,1.2854,0.4853,2.0956,1.502,0.3032};
	double Q2_upper[8]={0.0267,0.0833,0.0997,1.5014,0.3847,0.1375,0.0807,1.4936};
	double Q2_lower[8]={0.0265,0.0815,0.0616,0.4317,0.2732,0.0919,0.0737,0.0765};
	double x_Q2[8]={0,2,4,6,8,10,12,14};*/
	double Q2[8]={1.6681,2.2462,1.9842,2.0956,1.502};
	double Q2_upper[8]={0.0267,0.0833,0.0997,0.1375,0.0807};
	double Q2_lower[8]={0.0265,0.0815,0.0616,0.0919,0.0737};
	double x_Q2[8]={0,2,4,6,8};
    TGraphAsymmErrors *Q2_plot =  new TGraphAsymmErrors(5,x_Q2,Q2,0,0,Q2_upper,Q2_lower);
    Q2_plot->SetTitle("");
    Q2_plot->SetMarkerStyle(3);
    Q2_plot->SetLineColor(1);
    Q2_plot->SetMarkerColor(1);
    Q2_plot->GetXaxis()->SetTitle("State");
    Q2_plot->GetXaxis()->CenterTitle();
    Q2_plot->GetXaxis()->SetTitleSize(0.05);
    Q2_plot->GetXaxis()->SetTitleOffset(0.9);
    Q2_plot->GetXaxis()->CenterTitle();
    Q2_plot->GetYaxis()->SetTitle("<Q^{2}>");
    Q2_plot->GetYaxis()->CenterTitle();
    Q2_plot->GetYaxis()->SetTitleSize(0.05);
    Q2_plot->GetYaxis()->SetTitleOffset(0.9);
    TAxis *xaxis = Q2_plot->GetXaxis();
   	xaxis->SetLabelSize(0);
    Q2_plot->Draw("AP");
    
    TLatex *t0 = new TLatex(0.05,1.30,"0^{+}_{1}");
    t0->SetTextAlign(22);
    t0->SetTextFont(43);
    t0->SetTextSize(24);
    t0->SetTextAngle(0);
    t0->Draw();
    TLatex *t = new TLatex(2.05,1.30,"2^{+}_{1}");
    t->SetTextAlign(22);
    t->SetTextFont(43);
    t->SetTextSize(24);
    t->SetTextAngle(0);
    t->Draw();
    TLatex *t2 = new TLatex(4.05,1.30,"4^{+}_{1}");
    t2->SetTextAlign(22);
    t2->SetTextFont(43);
    t2->SetTextSize(24);
    t2->SetTextAngle(0);
    t2->Draw();
    TLatex *t3 = new TLatex(6.05,1.30,"2^{+}_{2}");
    t3->SetTextAlign(22);
    t3->SetTextFont(43);
    t3->SetTextSize(24);
    t3->SetTextAngle(0);
    t3->Draw();
    TLatex *t4 = new TLatex(8.05,1.30,"4^{+}_{2}");
    t4->SetTextAlign(22);
    t4->SetTextFont(43);
    t4->SetTextSize(24);
    t4->SetTextAngle(0);
    t4->Draw();
    
    /*double Q2_Wu[7]={1.46,1.6,1.87,1.6,1.25,2.14,1.51};
	double Q2_upper_Wu[7]={0.12,0.07,0.69,0.42,0.28,0.82,1.46};
	double Q2_lower_Wu[7]={0.1,0.06,0.25,0.25,0.21,0.57,0.31};
	double x_Q2_Wu[7]={0.2,2.2,4.2,6.2,8.2,10.2,12.2};*/
	double Q2_Wu[7]={1.46,1.6,1.87,2.14,1.51};
	double Q2_upper_Wu[7]={0.12,0.07,0.69,0.82,1.46};
	double Q2_lower_Wu[7]={0.1,0.06,0.25,0.57,0.31};
	double x_Q2_Wu[7]={0.2,2.2,4.2,6.2,8.2};
    TGraphAsymmErrors *Q2_plot_Wu =  new TGraphAsymmErrors(5,x_Q2_Wu,Q2_Wu,0,0,Q2_upper_Wu,Q2_lower_Wu);
    Q2_plot_Wu->SetMarkerStyle(4);
    Q2_plot_Wu->SetLineColor(2);
    Q2_plot_Wu->SetMarkerColor(2);
    Q2_plot_Wu->Draw("P same");
    
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(Q2_plot,"Present work","p");
    legend->AddEntry(Q2_plot_Wu,"Wu","p");
    legend->Draw();
    
    TCanvas *c2 = new TCanvas("c2","",900,900);
    
    double cos3d[8]={-0.5607,-0.4684,-0.4517/*,-0.1227,-0.1323*/};
	double cos3d_upper[8]={0.0546,0.0825,0.2002/*,0.3109,0.2214*/};
	double cos3d_lower[8]={0.0548,0.0742,0.1766/*,0.2769,0.2063*/};
	double x_cos3d[8]={0,2,4};
	TGraphAsymmErrors *cos3d_plot =  new TGraphAsymmErrors(3,x_cos3d,cos3d,0,0,cos3d_upper,cos3d_lower);
    cos3d_plot->SetMarkerStyle(3);
    cos3d_plot->SetLineColor(1);
    cos3d_plot->SetMarkerColor(1);
    cos3d_plot->SetTitle("");
    cos3d_plot->GetXaxis()->SetTitle("State");
    cos3d_plot->GetXaxis()->CenterTitle();
    cos3d_plot->GetXaxis()->SetTitleSize(0.05);
    cos3d_plot->GetXaxis()->SetTitleOffset(0.9);
    cos3d_plot->GetXaxis()->CenterTitle();
    cos3d_plot->GetYaxis()->SetTitle("<Q^{3}cos(3#delta )>");
    cos3d_plot->GetYaxis()->CenterTitle();
    cos3d_plot->GetYaxis()->SetTitleSize(0.05);
    cos3d_plot->GetYaxis()->SetTitleOffset(0.9);
   	cos3d_plot->GetXaxis()->SetLabelSize(0);
    cos3d_plot->Draw("AP");
    
    TLatex *t00 = new TLatex(0.05,-0.71,"0^{+}_{1}");
    t00->SetTextAlign(22);
    t00->SetTextFont(43);
    t00->SetTextSize(24);
    t00->SetTextAngle(0);
    t00->Draw();
    TLatex *t1 = new TLatex(2.05,-0.71,"2^{+}_{1}");
    t1->SetTextAlign(22);
    t1->SetTextFont(43);
    t1->SetTextSize(24);
    t1->SetTextAngle(0);
    t1->Draw();
    TLatex *t21 = new TLatex(4.05,-0.71,"4^{+}_{1}");
    t21->SetTextAlign(22);
    t21->SetTextFont(43);
    t21->SetTextSize(24);
    t21->SetTextAngle(0);
    t21->Draw();
    /*    TLatex *t31 = new TLatex(6.05,-0.76,"2^{+}_{2}");
    t31->SetTextAlign(22);
    t31->SetTextFont(43);
    t31->SetTextSize(24);
    t31->SetTextAngle(0);
    t31->Draw();
    TLatex *t41 = new TLatex(8.05,-0.76,"4^{+}_{2}");
    t41->SetTextAlign(22);
    t41->SetTextFont(43);
    t41->SetTextSize(24);
    t41->SetTextAngle(0);
    t41->Draw();*/
    
    
    double cos3d_Wu[7]={-0.48,-0.47,-0.43,-0.21,-0.53};
	double cos3d_upper_Wu[7]={0.03,0.05,0.17,0.187,0.3};
	double cos3d_lower_Wu[7]={0.06,0.09,0.12,0.17,0.16};
	double x_cos3d_Wu[7]={0.2,2.2,4.2,6.2,8.2};
    TGraphAsymmErrors *cos3d_plot_Wu =  new TGraphAsymmErrors(3,x_cos3d_Wu,cos3d_Wu,0,0,cos3d_upper_Wu,cos3d_lower_Wu);
    cos3d_plot_Wu->SetMarkerStyle(4);
    cos3d_plot_Wu->SetLineColor(2);
    cos3d_plot_Wu->SetMarkerColor(2);
    cos3d_plot_Wu->Draw("P same");
   
    auto legend_2 = new TLegend(0.1,0.7,0.48,0.9);
    legend_2->AddEntry(cos3d_plot,"Present work","p");
    legend_2->AddEntry(cos3d_plot_Wu,"Wu","p");
    legend_2->Draw();
    
    TCanvas *c3 = new TCanvas("c3","",900,900);
    double beta[8]={0.143722,0.1667843,0.156755,0.16109,0.136384};
	double beta_error_low[8]={0.0014372800,0.00316784076,0.00215675,0.0056103,0.00413639};
	double beta_error_up[8]={0.0014372800,0.00316784076,0.00415675,0.0046103,0.00313639};
	double beta_wu[8]={0.1344,0.13899,0.152,0.162,0.1367};
	double beta_error_up_wu[8]={0.005,0.003,0.03,0.023,0.015};
	double beta_error_low_wu[8]={0.005,0.003,0.01,0.028,0.06};
	double x_beta[8]={0,2,4,6,8};
	double x_beta_wu[8]={0.2,2.2,4.2,6.2,8.2};
	TGraphAsymmErrors *beta_plot =  new TGraphAsymmErrors(5,x_beta,beta,0,0,beta_error_up,beta_error_low);
	TGraphAsymmErrors *beta_plot_wu =  new TGraphAsymmErrors(5,x_beta_wu,beta_wu,0,0,beta_error_up_wu,beta_error_low_wu);
    beta_plot->SetMarkerStyle(3);
    beta_plot->SetLineColor(1);
    beta_plot->SetMarkerColor(1);
    beta_plot->SetTitle("");
    beta_plot->GetXaxis()->SetTitle("State");
    beta_plot->GetXaxis()->CenterTitle();
    beta_plot->GetXaxis()->SetTitleSize(0.05);
    beta_plot->GetXaxis()->SetTitleOffset(0.9);
    beta_plot->GetXaxis()->CenterTitle();
    beta_plot->GetYaxis()->SetTitle("#beta");
    beta_plot->GetYaxis()->CenterTitle();
    beta_plot->GetYaxis()->SetTitleSize(0.05);
    beta_plot->GetYaxis()->SetTitleOffset(0.9);
   	beta_plot->GetXaxis()->SetLabelSize(0);
    beta_plot->Draw("AP");
    beta_plot_wu->SetMarkerStyle(4);
    beta_plot_wu->SetLineColor(2);
    beta_plot_wu->SetMarkerColor(2);
    beta_plot_wu->Draw("P same");

    auto legend3 = new TLegend(0.1,0.7,0.48,0.9);
    legend3->AddEntry(beta_plot,"Present work","p");
    legend3->AddEntry(beta_plot_wu,"Wu","p");
    legend3->Draw();
    
    TLatex *t000 = new TLatex(0.05,0.1280,"0^{+}_{1}");
    t000->SetTextAlign(22);
    t000->SetTextFont(43);
    t000->SetTextSize(24);
    t000->SetTextAngle(0);
    t000->Draw();
    TLatex *t10 = new TLatex(2.05,0.1280,"2^{+}_{1}");
    t10->SetTextAlign(22);
    t10->SetTextFont(43);
    t10->SetTextSize(24);
    t10->SetTextAngle(0);
    t10->Draw();
    TLatex *t210 = new TLatex(4.05,0.1280,"4^{+}_{1}");
    t210->SetTextAlign(22);
    t210->SetTextFont(43);
    t210->SetTextSize(24);
    t210->SetTextAngle(0);
    t210->Draw();
    TLatex *t310 = new TLatex(6.05,0.1280,"2^{+}_{2}");
    t310->SetTextAlign(22);
    t310->SetTextFont(43);
    t310->SetTextSize(24);
    t310->SetTextAngle(0);
    t310->Draw();
    TLatex *t410 = new TLatex(8.05,0.1280,"4^{+}_{2}");
    t410->SetTextAlign(22);
    t410->SetTextFont(43);
    t410->SetTextSize(24);
    t410->SetTextAngle(0);
    t410->Draw();
    
   
    TCanvas *c4 = new TCanvas("c4","",900,900);
    double gamma[8]={41.3680737465977,39.3101623681799,38.9509355346866};
	double gamma_error[8]={1.3,1.65,4};
	double gamma_wu[8]={39.561, 39.344, 38.5};
	double gamma_error_up_wu[8]={0.6, 1.1, 3.5};
	double gamma_error_low_wu[8]={1.3, 2.0, 2.6};

	double x_gamma[8]={0,2,4};
	double x_gamma_wu[8]={0.2,2.2,4.2};
	TGraphErrors *gamma_plot =  new TGraphErrors(3,x_gamma,gamma,0,gamma_error);
	//TGraphAsymmErrors *cos3d_plot_Wu =  new TGraphAsymmErrors(5,x_cos3d_Wu,cos3d_Wu,0,0,cos3d_upper_Wu,cos3d_lower_Wu);
	TGraphAsymmErrors *gamma_plot_wu =  new TGraphAsymmErrors(3,x_gamma_wu,gamma_wu,0,0,gamma_error_up_wu,gamma_error_low_wu);
    gamma_plot->SetMarkerStyle(3);
    gamma_plot->SetLineColor(1);
    gamma_plot->SetMarkerColor(1);
    gamma_plot->SetTitle("");
    gamma_plot->GetXaxis()->SetTitle("State");
    gamma_plot->GetXaxis()->CenterTitle();
    gamma_plot->GetXaxis()->SetTitleSize(0.05);
    gamma_plot->GetXaxis()->SetTitleOffset(0.9);
    gamma_plot->GetXaxis()->CenterTitle();
    gamma_plot->GetYaxis()->SetTitle("#gamma");
    gamma_plot->GetYaxis()->CenterTitle();
    gamma_plot->GetYaxis()->SetTitleSize(0.05);
    gamma_plot->GetYaxis()->SetTitleOffset(0.9);
   	gamma_plot->GetXaxis()->SetLabelSize(0);
    
    gamma_plot->Draw("AP");
    gamma_plot_wu->SetMarkerStyle(4);
    gamma_plot_wu->SetLineColor(2);
    gamma_plot_wu->SetMarkerColor(2);
    gamma_plot_wu->Draw("same P");

    auto legend4 = new TLegend(0.1,0.7,0.48,0.9);
    legend4->AddEntry(gamma_plot,"Present work","p");
    legend4->AddEntry(gamma_plot_wu,"Wu","p");
    legend4->Draw();
    
    TLatex *t0000 = new TLatex(0.05,33.8,"0^{+}_{1}");
    t0000->SetTextAlign(22);
    t0000->SetTextFont(43);
    t0000->SetTextSize(24);
    t0000->SetTextAngle(0);
    t0000->Draw();
    TLatex *t100 = new TLatex(2.05,33.8,"2^{+}_{1}");
    t100->SetTextAlign(22);
    t100->SetTextFont(43);
    t100->SetTextSize(24);
    t100->SetTextAngle(0);
    t100->Draw();
    TLatex *t2100 = new TLatex(4.05,33.8,"4^{+}_{1}");
    t2100->SetTextAlign(22);
    t2100->SetTextFont(43);
    t2100->SetTextSize(24);
    t2100->SetTextAngle(0);
    t2100->Draw();
    
}

