/* Written by Iris Ponce
This code is written to analyze the frequency of the different flipped bits.
This will will give an idea of which bit flips the most. The result will be
a root file, a TH2D with the different logarithms
*/


//some standard C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <chrono>
#include <math.h>
#include <cmath>
#include <fstream>
#include "TSpectrum.h"
//some ROOT includes
//#include "-lSpectrum"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2S.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TLine.h"
#include "TStyle.h"
#include "TExec.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMath.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"



//convenient for us! let's not bother with art and std namespaces!
using namespace art;
using namespace std;

using namespace std::chrono;

int main(int argc, char** argv)
{
    vector<string> filenames { argv[1] };
    TFile f_output("flipbitFrequency19021_output.root","RECREATE");
    InputTag wire_tag { "sndaq", "", "SupernovaAssembler" };
    InputTag wire_tag_d { "sndeco", "", "DataRecoSN"};
    // InputTag wire_tag_d { "sndeco", "", "CalDataSN" }; // after deconvolution
    InputTag hit_tag_d { "gaushit", "", "DataRecoSN"};
    //   int similar = 0;
    size_t _maxEvts = 1000;
    size_t evCtr = 0;
    double rois[8256] = {0};
    vector<double> differences, flips;
    ofstream flipsJumps;
   // TH1D
    flipsJumps.open ("/uboone/app/users/iponce/myGallery/cpp/flipbitFrequency19021.txt");
    TH1D flips_histU("Flips_Uplane", " ; Difference from non flipped signals", 4096,0,4096 );
    TH1D flips_histV("Flips_Vplane", " ; Difference from non flipped signals", 4096,0,4096 );
    TH1D flips_histY("Flips_Yplane", " ; Difference from non flipped signals", 4096,0,4096 );
    TH2D flips_Channels("Flips_channels", " ;Wire;ADC",8256,0,8256,4096,0,4096);
    //double signals[8256] = {0};
    
    for(gallery::Event ev(filenames); !ev.atEnd(); ev.next())
    {
        if(evCtr >= _maxEvts)
        {
            break;
        }
        auto t_begin = high_resolution_clock::now();
        // int run = ev.eventAuxiliary().run();
        int event = ev.eventAuxiliary().event();
        
        cout << "Processing "
        //     << "Run " << run << ", "
        << "Event " << event << endl;
        
        //this is getting the handle and accessing the RecobWires from the root file
        auto const& wire_handle_d = ev.getValidHandle< vector<recob::Wire> >(wire_tag_d);
        auto const& wire_vec_d(*wire_handle_d);
        
        auto const& wire_handle = ev.getValidHandle< vector<recob::Wire> >(wire_tag);
        auto const& wire_vec(*wire_handle);
        //   double rois[8256] = {0};
        cout <<  "predecon:"<<  wire_vec.size() << " decon: " << wire_vec_d.size() <<endl;
        for (unsigned int index = 0; index < wire_vec.size();index++)//this is checking through all the
            //different wires
        {
            auto zsROIs_pre = wire_vec[index].SignalROI();
            int channel_pre = wire_vec[index].Channel();
            //SignalROI gives a sparse vector with different ranges. This is going through all these ranges
            for (auto iROI_pre = zsROIs_pre.begin_range(); iROI_pre != zsROIs_pre.end_range(); ++iROI_pre)
            {
                auto ROI_pre = *iROI_pre;
                const size_t firstTick_pre = ROI_pre.begin_index();
                const size_t endTick_pre = ROI_pre.end_index();
                //int similar = 0;
                double difftoint = 0;
                double preSample = 0;
                double postSample = 0;
		int flips = 0;
                //this is looping through the different waveforms and checking for flipped bits
                for (size_t iTick = ROI_pre.begin_index(); iTick < ROI_pre.end_index()-1; ++iTick)
                {
                    int flips =0;
                    //this checks whether iTick is at the begining of the range
                    //takes the presample to be 0
                    if (iTick - ROI_pre.begin_index() == 0)
                    {
                        preSample = 0;
                        postSample = ((ROI_pre[iTick+1] + ROI_pre[iTick+2])/2);
                        difftoint = ROI_pre[iTick] - postSample;
                        if (difftoint > 32)
                        {
                            flips += 1;
                            if (flips>0)
                            {
                                differences.push_back(difftoint);
                                flipsJumps << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " " << iTick << " " << preSample << " " << postSample << " "  <<ROI_pre[iTick] <<" "<< difftoint << "\n";
                                flips_Channels.Fill(channel_pre,difftoint);
                                if( channel_pre < 2400)
                                {
                                    flips_histU.Fill(difftoint,1);
                                }
                                else if(channel_pre >= 2400 && channel_pre < 4800)
                                {
                                    flips_histV.Fill(difftoint,1);
                                }
                                else
                                {
                                    flips_histY.Fill(difftoint,1);
		//		    cout << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " "  << difftoint << " " << preSample << " " << postSample << endl;
                                }
                            }
                        }
                        
                    }//ENDS LOOP FOR FIRST SAMPLE
                    else if (iTick - ROI_pre.begin_index() == 1)
                    {
                        preSample = ROI_pre[iTick-1];
                        postSample = ((ROI_pre[iTick+1] + ROI_pre[iTick+2])/2);
                        if (preSample > postSample)
                        {
                            difftoint = ROI_pre[iTick] - postSample;
                            
                            if (difftoint > 32)
                            {
                                flips += 1;
                                if (flips>0)
                                {
                                    differences.push_back(difftoint);
                                    flipsJumps << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre <<" " << iTick<< " " << preSample << " " << postSample << " "<< ROI_pre[iTick] << " " << difftoint << "\n";
                                    flips_Channels.Fill(channel_pre,difftoint);
                                    if( channel_pre < 2400)
                                    {
                                        flips_histU.Fill(difftoint,1);
                                    }
                                    else if(channel_pre >= 2400 && channel_pre < 4800)
                                    {
                                        flips_histV.Fill(difftoint,1);
                                    }
                                    else
                                    {
                       //                  cout << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " "  << difftoint << " " << preSample << " " << postSample << endl;                    
		                        flips_histY.Fill(difftoint,1);
                                    }
                                }
                            }
                        }
                        
                        else
                        {
                            difftoint = ROI_pre[iTick] - preSample;
                            if ( difftoint > 32)
                            {
                                flips += 1;
                                if (flips>0)
                                {
                                    flipsJumps << event << " " <<channel_pre<< " " << firstTick_pre<< " " << endTick_pre << " " << iTick << " " << preSample << " " << postSample << " "<<  ROI_pre[iTick]<<" "<<difftoint << "\n";
                                    differences.push_back(difftoint);
                                    flips_Channels.Fill(channel_pre,difftoint);
                                    if( channel_pre < 2400)
                                    {
                                        flips_histU.Fill(difftoint,1);
                                    }
                                    else if(channel_pre >= 2400 && channel_pre < 4800)
                                    {
                                        flips_histV.Fill(difftoint,1);
                                    }
                                    else
                                    {
                                        flips_histY.Fill(difftoint,1);
                     //                   cout << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " "  << difftoint << " " << preSample << " " << postSample << endl;
				     }
                                }
                                
                                
                            }
                        }
                    }//ENDS LOOPS FOR SECOND SAMPLE
                    
                    //checks if you're at the end of the ROI
                    else if ((ROI_pre.end_index()-1) - iTick == 0)
                    {
                        postSample = 0;
                        preSample =((ROI_pre[iTick-1] + ROI_pre[iTick-2])/2);
                        difftoint = ROI_pre[iTick] - preSample; //assume presample is corect
                        if (difftoint > 32)
                        {
                            flips += 1;
                            if (flips>0)
                            {
                                flipsJumps <<event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " " << iTick<< " "<< preSample << " " << postSample << " " << ROI_pre[iTick] << " " << difftoint << "\n";
                                differences.push_back(difftoint);
                                flips_Channels.Fill(channel_pre,difftoint);
                                if( channel_pre < 2400)
                                {
                                    flips_histU.Fill(difftoint,1);
                                }
                                else if(channel_pre >= 2400 && channel_pre < 4800)
                                {
                                    flips_histV.Fill(difftoint,1);
                                }
                                else
                                {
                                    flips_histY.Fill(difftoint,1);
	             //               cout << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " "  << difftoint << " " << preSample << " " << postSample << endl;
                                }
                            }
                        }
                    }//ENDS LOOP OVER LAST SAMPLE
                    else if ((ROI_pre.end_index()-1) - iTick == 1)
                    {
                        postSample = ROI_pre[iTick-1];
                        preSample = ((ROI_pre[iTick-1] + ROI_pre[iTick-2])/2);
                        if (preSample > postSample)
                        {
                            difftoint = ROI_pre[iTick] - postSample;
                            if (difftoint > 32)
                            {
                                flips += 1;
                                if (flips>0)
                                {
                                    flipsJumps << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre <<" " << iTick<<  " "<< preSample << " " << postSample << " " << ROI_pre[iTick] << " " << difftoint << "\n";
                                    differences.push_back(difftoint);
                                    flips_Channels.Fill(channel_pre,difftoint);
                                    if( channel_pre < 2400)
                                    {
                                        flips_histU.Fill(difftoint,1);
                                    }
                                    else if(channel_pre >= 2400 && channel_pre < 4800)
                                    {
                                        flips_histV.Fill(difftoint,1);
                                    }
                                    else
                                    {
                       //                 cout << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " "  << difftoint << " " << preSample << " " << postSample << endl;
                                        flips_histY.Fill(difftoint,1);
                                    }
                                }
                            }
                        }
                        else
                        {
                            difftoint = ROI_pre[iTick] - preSample;
                            if (difftoint > 32)
                            {
                                flips += 1;
                                if ( flips>0)
                                {
                                    flipsJumps << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre <<" " << iTick << " "<< preSample << " " << postSample << " "<< ROI_pre[iTick] <<" " <<difftoint << "\n";
                                    differences.push_back(difftoint);
                                    flips_Channels.Fill(channel_pre,difftoint);
                                    if( channel_pre < 2400)
                                    {
                                        flips_histU.Fill(difftoint,1);
                                    }
                                    else if(channel_pre >= 2400 && channel_pre < 4800)
                                    {
                                        flips_histV.Fill(difftoint,1);
                                    }
                                    else
                                    {
                     //                     cout << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " "  << difftoint << " " << preSample << " " << postSample << endl;
                                        flips_histY.Fill(difftoint,1);
                                    }
                                }
                            }
                        }
                    }//ENDS LOOP OVER SECOND TO LAST SAMPLE
                    else
                    {
                        preSample = ((ROI_pre[iTick-2]+ROI_pre[iTick-1])/2);
                        postSample = ((ROI_pre[iTick+2]+ROI_pre[iTick+1])/2);
                        if (preSample > postSample)
                        {
                            difftoint = ROI_pre[iTick] - postSample;
                            if (difftoint > 32)
                            {
                                flips += 1;
                                if (flips>0)
                                {
                                    flipsJumps <<event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " " << iTick << " " << preSample << " " << postSample << " " <<ROI_pre[iTick] << " " << difftoint << "\n";
                                    flips_Channels.Fill(channel_pre,difftoint);
                                    differences.push_back(difftoint);
                                    if( channel_pre < 2400)
                                    {
                                        flips_histU.Fill(difftoint,1);
                                    }
                                    else if(channel_pre >= 2400 && channel_pre < 4800)
                                    {
                                        flips_histV.Fill(difftoint,1);
                                    }
                                    else
                                    { 
                         //           cout << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " "  << difftoint << " " << preSample << " " << postSample << endl;
                                        flips_histY.Fill(difftoint,1);
                                    }
                                }
                            }
                        }
                        else
                        {
                            difftoint = ROI_pre[iTick] - preSample;
                            if (difftoint > 32)
                            {
                                flips +=1;
                                if (flips>0)
                                {
                                    flipsJumps<<event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre<< " " << iTick<< " " << preSample << " " << postSample << " " << ROI_pre[iTick] << " " <<  difftoint << "\n";
                                    differences.push_back(difftoint);
                                    flips_Channels.Fill(channel_pre,difftoint);
                                    if( channel_pre < 2400)
                                    {
                                        flips_histU.Fill(difftoint,1);
                                    }
                                    else if(channel_pre >= 2400 && channel_pre < 4800)
                                    {
                                        flips_histV.Fill(difftoint,1);
                                    }
                                    else
                                    { 
                                    // cout << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " "  << difftoint << " " << preSample << " " << postSample << endl;
                                        flips_histY.Fill(difftoint,1);
                                    }
                                }
                            }
                        }
                        
                    }//ENDS LOOP OVER ANY OTHER SAMPLE IN ROI
                    
                }//ENDS THE LOOP FOR FLIPS
                //checks whether flips were found
                if (flips > 0)
                {
                    rois[channel_pre] = rois[channel_pre]+1;
                    if (event == 2 && channel_pre == 72 && firstTick_pre == 1712 && endTick_pre == 1732)
                    {
                        cout << event << " " <<channel_pre << " " << firstTick_pre << " " << endTick_pre << " " << flips << " diff " << difftoint << endl;
                    }
                    
                }
                
            }//ENDS LOOP OVER ROIS
        }//ENDS LOOP OVER WIRES
        
        //auto const& hit_handle = ev.getValidHandle< vector<recob::Hit> >(hit_tag_d);
        //auto const& hit_vec(*hit_handle);
        
        //this is going through the different wires with the deconvolved data
        for (unsigned int i=0; i<wire_vec_d.size();i++)
        {
            //this is getting the ROI Signals for every wire.
            auto zsROIs_d = wire_vec_d[i].SignalROI();
            int channel_d = wire_vec_d[i].Channel();
            
            auto zsROIs = wire_vec[i].SignalROI();
            //int channel = wire_vec[i].Channel();
            
            //going through the entire ROI Signals excluding the last tick to integrate it
            for (auto iROI_d = zsROIs_d.begin_range(); iROI_d != zsROIs_d.end_range(); ++iROI_d)
            {
                
                auto ROI_d = *iROI_d;
                const size_t firstTick = ROI_d.begin_index();
                const size_t endTick = ROI_d.end_index();
                double integral;
                vector<double> stdDev;
                double sDev;
                TH1D horig("roi_original", "roi_original;Tick;ADC", endTick - firstTick, firstTick, endTick); // new hist of the waveform
                for (size_t iTick = ROI_d.begin_index(); iTick < ROI_d.end_index(); iTick++ )
                {    //not including last sample
                    horig.Fill((int)iTick,ROI_d[iTick]);
                    stdDev.push_back(ROI_d[iTick]);
                }
                integral = horig.Integral();
                sDev = TMath::StdDev(stdDev.begin(),stdDev.end());
                //writing to the txt files for the ROI integrals
 //               ROIfile << event << " " <<channel_d << " " << firstTick << " " << endTick << " " << integral << " " <<  sDev << "\n";
                
                
                //sometimes we want specific waveforms saved change the values as needed.
                if (event ==2 && channel_d == 207 && (int)firstTick ==1678 && (int)endTick == 1732)
                    //if ( (integral < 96 && integral > 95) && channel_d >= 5464 && channel_d <5465 )
                {
                    // cout << "int:" << integral << " startTick: " << firstTick << " endTick:"  <<endTick <<" channel:" << channel_d <<endl;
                   cout << sDev << endl; 
                    // this is checking for a specific RAW waveforms with the parameters given above (could be eraed)
                    for (unsigned int j = 0; j < wire_vec.size(); j++)
                    {
                        auto zsROIs = wire_vec[j].SignalROI();
                        int channel = wire_vec[j].Channel();
                        for (auto iROI = zsROIs.begin_range();iROI!=zsROIs.end_range();++iROI)
                        {
                            auto ROI = *iROI;
                            const size_t oneTick = ROI.begin_index();
                            const size_t lastTick = ROI.end_index();
                            if (channel == channel_d && oneTick == firstTick && endTick == lastTick)
                            {
                                // cout << "decon " << wire_vec_d.size() << " non decon " << wire_vec.size() << endl
                                cout << "YAY " << lastTick << "and " << oneTick<< endl;
                                //      double d_integral;
                                TH1D h_orig("roi_preD","Waveform Pre-Deconvolution;Tick;ADC",lastTick - oneTick,oneTick,lastTick);
                                for (size_t tick = ROI.begin_index(); tick < ROI.end_index(); tick++)
                                {
                                    h_orig.Fill((int)tick,ROI[tick]);
                                }
                                cout << "int:" << integral << " startTick: " << firstTick << " endTick:"  <<endTick <<" channel:" << channel_d <<endl;
                                f_output.cd();
                                h_orig.Write();
                            }
                            
                        }//once the RAW waveform is found, it's saved in the ROOT file
                    }
                    //cout << "int:" << integral << " startTick: " << firstTick << " endTick:"  <<endTick <<" channel:" << channel_d <<endl;
                    f_output.cd();
                    horig.Write();
                    //f_output.Close();
                }
                else if((int)firstTick == 1894 && (int)endTick == 1921 && channel_d == 8058)
                    //else if ( (integral > 15 && integral < 20)&& channel_d > 8050 && channel_d < 8090)
                {
                    //  f_output.cd();
                    //  horig.Write();
                    //  f_output.Close();
                    cout << "noise: " << "int:" << integral << " startTick: " << firstTick << " endTick:"  <<endTick <<" channel:" << channel_d <<endl;
                }
                //now we search for the corresponding Hit for the ROI signal, if applicable.
                
                /*
                 int hit = 0;
                 for (unsigned int j=0; j<hit_vec.size();j++)
                 {
                 if ((double)channel_d == (double)hit_vec[j].Channel() && ((double)hit_vec[j].StartTick() >= (double)firstTick && (double)hit_vec[j].StartTick() <= (double)endTick))
                 {
                 similar++;
                 hit += (double)hit_vec[j].Integral();
                 }
                 
                 else
                 {
                 hitsVsROI.Fill(integral,0);
                 
                 }
                 }
                 similar++;
                 hitsVsROI.Fill(integral,hit);
                 */
                
            }//integrates
            
            
        }//goes through different wires
        
        
        //    cout << hit_vec.size() << " hits" << endl;
        cout << wire_vec_d.size() << " wires" << endl;
        auto t_end = high_resolution_clock::now();
        duration<double,std::milli> time_total_ms(t_end-t_begin);
        //cout << "\tEvent took " << time_total_ms.count() << " ms to process." << endl;
        evCtr++;
    }//goes through all events
    double max = *max_element(differences.begin(),differences.end());
    cout << "Max:" << max << endl;
    double power2 = 0;
    //double flips[differences.size()] = {0};
    TH2D flipfrequency("flip_freq","frequency of signal differences; Power of 2; Frequency",5000,0,5000,1,1,2);

    for (unsigned int i=0; i < differences.size(); i++)
    {
        power2 = log2 (differences[i]);
        flipfrequency.Fill(power2,1);
        flips.push_back(differences[i]);
        	//cout << differences[i] << endl;
    }	
    TCanvas *c1 = new TCanvas;
    
    c1->cd();
    flips_histU.Draw();
    c1->Update();
    for (int i = 0; i < 13; i++)
    {
	TLine *l=new TLine(pow(2,i),c1->GetUymin(),pow(2,i),c1->GetUymax());
	l->Draw("same");
    }
    c1->cd();
    
    TH1D* d_U = new TH1D("background_U","",flips_histU.GetNbinsX(),0,flips_histU.GetNbinsX());
    TSpectrum *s_U = new TSpectrum();
    //int bins = flips_histU.GetNbinsX();
    double source_U[4096] = {0};
    for ( int i =0; i < flips_histU.GetNbinsX(); i++)
    {
	source_U[i] = flips_histU.GetBinContent(i+1);
    }

    s_U->Background(source_U,flips_histU.GetNbinsX(),8,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kFALSE,TSpectrum::kBackSmoothing3,kFALSE);
//","kBackOrder2","kFALSE","kBackSmoothing3","kFALSE");
    for ( int i =0; i<flips_histU.GetNbinsX();i++)    {
	d_U->SetBinContent(i+1,source_U[i]);
    }
    d_U->Draw("SAME L");
    
    c1->Print(".png");
    
    TCanvas *c2 = new TCanvas;
    c2->cd();
    flips_histV.Draw();
    c2->Update();
    for (int i=0; i<13; i++)
    {
	c2->cd();
	TLine *V = new TLine(pow(2,i),c2->GetUymin(),pow(2,i),c2->GetUymax());
        V->Draw("same");
 
    }
    c2->cd();
    
    TH1D* d_V = new TH1D("background_V","",flips_histV.GetNbinsX(),0,flips_histV.GetNbinsX());
    TSpectrum *s_V = new TSpectrum();
    //int bins = flips_histU.GetNbinsX();
    double source_V[4096] = {0};
    for ( int i =0; i < flips_histV.GetNbinsX(); i++)
    {
        source_V[i] = flips_histV.GetBinContent(i+1);
    }
    
    s_V->Background(source_V,flips_histU.GetNbinsX(),8,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kFALSE,TSpectrum::kBackSmoothing3,kFALSE);
    //","kBackOrder2","kFALSE","kBackSmoothing3","kFALSE");
    for ( int i =0; i<flips_histV.GetNbinsX();i++)    {
        d_V->SetBinContent(i+1,source_V[i]);
    }
    d_V->Draw("SAME L");
    c2->Print(".png");
    
    TCanvas *c3 = new TCanvas; 
    c3->cd();
    flips_histY.Draw();
    c3->Update();
    for (int i=0; i<13; i++)
    {
	c3->cd();
        TLine *Y = new TLine(pow(2,i),c3->GetUymin(),pow(2,i),c3->GetUymax());
        Y->Draw("same");
    }
    c3->cd();
    
    TH1D* d_Y = new TH1D("background_Y","",flips_histY.GetNbinsX(),0,flips_histY.GetNbinsX());
    TSpectrum *s_Y = new TSpectrum();
    //int bins = flips_histU.GetNbinsX();
    double source_Y[4096] = {0};
    for ( int i =0; i < flips_histY.GetNbinsX(); i++)
    {
        source_Y[i] = flips_histY.GetBinContent(i+1);
    }
    
    s_Y->Background(source_Y,flips_histU.GetNbinsX(),8,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kFALSE,TSpectrum::kBackSmoothing3,kFALSE);
    //","kBackOrder2","kFALSE","kBackSmoothing3","kFALSE");
    for ( int i =0; i<flips_histY.GetNbinsX();i++)    {
        d_Y->SetBinContent(i+1,source_Y[i]);
    }
    d_Y->Draw("SAME L");
    c3->Print(".png");
    
    
   TCanvas *c4 = new TCanvas;
   c4->cd();
   flips_Channels.Draw("colz");
   c4->Update();
   for (int i = 0; i <13; i++)
   {
	c4->cd();
        TLine *W = new TLine(0,pow(2,i),8256,pow(2,i));
	W->Draw("same");
   }
   c4->cd();
   c4->Print(".png");

    //TSpectrum *s = new TSpectrum();
    //s->Background(flips_histU, flips_histU.GetNbins(),"kBackIncreasingWindow","kBackOrder2","kFALSE","kBackSmoothing3","kFALSE");
    f_output.cd();
    d_U->Write();
    d_V->Write();
    d_Y->Write();
    flipfrequency.Write();
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    flips_histU.Write();
    flips_histV.Write();
    flips_histY.Write();
    flips_Channels.Write();
    f_output.Close();

    flipsJumps.close();
    f_output.Close();
    delete c1;
    delete c2;
    delete c3;
    //delete l;
}
//end









