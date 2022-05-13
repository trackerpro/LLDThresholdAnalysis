import os
import cx_Oracle
import argparse
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join, basename
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
import ROOT
import math
ROOT.PyConfig.DisableRootLogon = True
ROOT.PyConfig.IgnoreCommandLineOptions = False

parser = argparse.ArgumentParser()
parser.add_argument('--partition', help='Specify partition: TEC, TIB, TID or TOB',  type=str, choices=['TECM','TECP','TIB','TOB'])
parser.add_argument('--referenceTemperature', help = 'reference temperature for rescaling',  type=int, default=0, nargs='?')
parser.add_argument('--runNumberMin', help = 'Min run number for good gain scans',  type=int, default=0, nargs='?')
parser.add_argument('--runNumberMax', help = 'Max run number for good gain scans',  type=int, default=999999, nargs='?')
parser.add_argument('--runNumberToExclude', help='bad run numbers that need to be manually excluded', nargs="+", default=[348601,348580,315860,312747,313887,348581,348602,335412,327725,312744,348579,348599,335404,335403]);
parser.add_argument('--outputTrackerMaps', help= 'output folder for tracker maps', type=str, default='outputTrackerMaps', nargs='?')
parser.add_argument('--inputDirectory', help = 'Name of the input direcotry with numpy files',  type=str, default="OPTOSCANS", nargs='?')
parser.add_argument('--dcuChannel', help = 'Express which DCU channel data has to be considered', type=str, default='thybrid', choices=['tsilicon','thybrid','tdcu'])
parser.add_argument('--thresholdGain', help = 'Express which LLD gain to be considered for threshold', type=str, default='Threshold0',  choices=['SelectedThreshold','Threshold0','Threshold1','Threshold2','Threshold3'])
parser.add_argument('--verbose', action='store_true', default=False, help='verbosity flag');
parser.add_argument('--dumpTrackerMap', action='store_true', default=False, help='dump tracker map for each run');
parser.add_argument('-b',"--batch",dest="batch",default=False,action="store_true")
args = parser.parse_args()

x_years_run_tk  = [287069,308166,327740]
x_years_run_pp  = [(271036,284078),(294927,307082),(314472,327564)]
x_years_lumi    = [33.74,75.32,125.125,192.98]


def setTDRStyle ():
  
  ROOT.gStyle.SetCanvasBorderMode(0);
  ROOT.gStyle.SetCanvasColor(0);
  ROOT.gStyle.SetCanvasDefH(600);
  ROOT.gStyle.SetCanvasDefW(600);
  ROOT.gStyle.SetCanvasDefX(0);
  ROOT.gStyle.SetCanvasDefY(0);

  ROOT.gStyle.SetPadBorderMode(0);
  ROOT.gStyle.SetPadColor(0); 
  ROOT.gStyle.SetPadGridX(0);
  ROOT.gStyle.SetPadGridY(0);
  ROOT.gStyle.SetGridColor(0);
  ROOT.gStyle.SetGridStyle(3);
  ROOT.gStyle.SetGridWidth(1);

  ROOT.gStyle.SetFrameBorderMode(0);
  ROOT.gStyle.SetFrameBorderSize(1);
  ROOT.gStyle.SetFrameFillColor(0);
  ROOT.gStyle.SetFrameFillStyle(0);
  ROOT.gStyle.SetFrameLineColor(1);
  ROOT.gStyle.SetFrameLineStyle(1);
  ROOT.gStyle.SetFrameLineWidth(1);
  ROOT.gStyle.SetHistLineColor(1);
  ROOT.gStyle.SetHistLineStyle(0);
  ROOT.gStyle.SetHistLineWidth(1);

  ROOT.gStyle.SetEndErrorSize(2);
  ROOT.gStyle.SetFuncColor(2);
  ROOT.gStyle.SetFuncStyle(1);
  ROOT.gStyle.SetFuncWidth(1);
  ROOT.gStyle.SetOptDate(0);
  
  ROOT.gStyle.SetOptFile(0);
  ROOT.gStyle.SetOptStat(0);
  ROOT.gStyle.SetStatColor(0); 
  ROOT.gStyle.SetStatFont(42);
  ROOT.gStyle.SetStatFontSize(0.04);
  ROOT.gStyle.SetStatTextColor(1);
  ROOT.gStyle.SetStatFormat("6.4g");
  ROOT.gStyle.SetStatBorderSize(1);
  ROOT.gStyle.SetStatH(0.1);
  ROOT.gStyle.SetStatW(0.15);

  ROOT.gStyle.SetPadTopMargin(0.07);
  ROOT.gStyle.SetPadBottomMargin(0.13);
  ROOT.gStyle.SetPadLeftMargin(0.12);
  ROOT.gStyle.SetPadRightMargin(0.05);

  ROOT.gStyle.SetOptTitle(0);
  ROOT.gStyle.SetTitleFont(42);
  ROOT.gStyle.SetTitleColor(1);
  ROOT.gStyle.SetTitleTextColor(1);
  ROOT.gStyle.SetTitleFillColor(10);
  ROOT.gStyle.SetTitleFontSize(0.05);

  ROOT.gStyle.SetTitleColor(1, "XYZ");
  ROOT.gStyle.SetTitleFont(42, "XYZ");
  ROOT.gStyle.SetTitleSize(0.05, "XYZ");
  ROOT.gStyle.SetTitleXOffset(0.9);
  ROOT.gStyle.SetTitleYOffset(1.05);
 
  ROOT.gStyle.SetLabelColor(1, "XYZ");
  ROOT.gStyle.SetLabelFont(42, "XYZ");
  ROOT.gStyle.SetLabelOffset(0.007, "XYZ");
  ROOT.gStyle.SetLabelSize(0.04, "XYZ");

  ROOT.gStyle.SetAxisColor(1, "XYZ");
  ROOT.gStyle.SetStripDecimals(1); 
  ROOT.gStyle.SetTickLength(0.025, "XYZ");
  ROOT.gStyle.SetNdivisions(510, "XYZ");
  ROOT.gStyle.SetPadTickX(1); 
  ROOT.gStyle.SetPadTickY(1);

  ROOT.gStyle.SetOptLogx(0);
  ROOT.gStyle.SetOptLogy(0);
  ROOT.gStyle.SetOptLogz(0);

  ROOT.gStyle.SetPaperSize(20.,20.);
  ROOT.gStyle.SetPaintTextFormat(".2f");


def plotGraph(graph_list,x_min,x_max,y_min,y_max,years,name,years_time=[]):

  color = [ROOT.kBlack,ROOT.kRed,ROOT.kBlue,ROOT.kGreen+1,ROOT.kOrange,ROOT.kViolet+1,ROOT.kAzure+10,ROOT.kMagenta,ROOT.kGray+1];

  canvas = ROOT.TCanvas("canvas","canvas",1000,600);
  canvas.cd();
  x_min = x_min-0.001*x_min;
  x_max = x_max+0.001*x_max;
  y_min = y_min*1.15 if y_min < 0 else y_min*0.85;
  y_max = y_max*0.85 if y_max < 0 else y_max*1.15;
  frame = ROOT.TH1F("frame","frame",1,x_min,x_max);
  frame.GetYaxis().SetRangeUser(y_min,y_max);
  if "run" in name:
    frame.GetXaxis().SetTitle("CMS Run Number");
  elif "lumi" in name:
    frame.GetXaxis().SetTitle("Integrated luminosity (fb^-1)");
  frame.GetYaxis().SetTitle("LLD Threshold Change [mA]");
  frame.Draw("*");
    
  boxes = [];
  for entries in years_time:
    g = ROOT.TBox(entries[0],y_min,entries[1],y_max);
    g.SetLineWidth(0);
    g.SetFillStyle(1001);
    g.SetFillColorAlpha(ROOT.kGray+1,0.3)
    g.Draw("same");
    boxes.append(g);

  lines = [];  
  for entry in years:
    line = ROOT.TLine();
    line.SetLineColor(ROOT.kBlack);
    line.SetLineWidth(2);
    line.SetLineStyle(7);
    line.SetX1(entry);
    line.SetX2(entry);
    line.SetY1(y_min);
    line.SetY2(y_max);
    lines.append(line);
    line.Draw("Lsame");

  latex = ROOT.TLatex();
  latex.SetNDC();
  latex.SetTextAlign(11);
  latex.SetTextFont(62);
  latex.SetTextSize(0.04);
  if "lumi" in name:
    latex.DrawLatex(0.22,0.85,"2016");
    latex.DrawLatex(0.50,0.85,"2017");
    latex.DrawLatex(0.78,0.85,"2018");
  else:
    latex.DrawLatex(0.18,0.85,"2015+2016");
    latex.DrawLatex(0.45,0.85,"2017");
    latex.DrawLatex(0.65,0.85,"2018");
    latex.DrawLatex(0.82,0.85,"LS2");

  latex.SetTextSize(0.05);
  latex.SetTextFont(62);
  latex.SetTextAlign(11);
  latex.DrawLatex(0.170,0.95,"CMS");
  latex.SetTextSize(0.05);
  latex.SetTextFont(52);
  latex.SetTextAlign(11);
  latex.DrawLatex(0.24,0.95,"Preliminary");

  latex.SetTextSize(0.05);
  latex.SetTextFont(62);
  latex.SetTextAlign(11);
  if "TIB" in name:
    latex.DrawLatex(0.22,0.72,"TIB");
  elif "TOB" in name:
    latex.DrawLatex(0.22,0.72,"TOB");
  elif "TECP" in name:
    latex.DrawLatex(0.22,0.72,"TEC+");
  elif "TECM" in name:
    latex.DrawLatex(0.22,0.72,"TEC-");

  base_label = "Ring"
  if "TIB" in name or "TOB" in name:
    base_label = "Layer"

  legend = ROOT.TLegend(0.15,0.45,0.4,0.7);
  legend.SetFillStyle(0);
  legend.SetFillColor(0);
  legend.SetBorderSize(0);
  for i,entry in enumerate(graph_list):
    entry.SetMarkerStyle(20);
    entry.SetMarkerSize(1);
    entry.SetLineColor(color[i]);
    entry.SetMarkerColor(color[i]);
    entry.SetLineWidth(1)
    entry.Draw("PLsame");
    legend.AddEntry(entry,base_label+" "+str(i+1),"PL");
  legend.Draw("same");

  #canvas.SaveAs(name+".png","png");
  canvas.SaveAs(name+".pdf","pdf");

def error_median(x):
  if len(x) > 1:
    return 1.253*x.std()/math.sqrt(len(x));
  else:
    return 0;

partName = ''
if    args.partition == 'TECM':
    partName = 'TM_09-JUN-2009_1'
elif  args.partition == 'TIB':
    partName = 'TI_27-JAN-2010_2'
elif  args.partition == 'TECP':
    partName = 'TP_09-JUN-2009_1'
elif  args.partition == 'TOB':
    partName = 'TO_30-JUN-2009_1'


WorkingDirectory = os.getcwd()
outputTrackerMapPath = join(WorkingDirectory,args.outputTrackerMaps)
os.system('mkdir -p '+outputTrackerMapPath);

### Constant parameters for all DCUs
### DCU calibration: https://gitlab.cern.ch/giassi/cmstrackerdaq/-/blob/trackerDAQ_5_0_0/FecSoftwareV3_0/generic/include/TkDcuConversionFactors.h
### DCU calibration: https://gitlab.cern.ch/giassi/cmstrackerdaq/-/blob/trackerDAQ_5_0_0/FecSoftwareV3_0/generic/src/common/TkDcuConversionFactors.cc
BETA = 3280.
T25  = 298.1 
T0   = 273.1
TC   = 55. # characteristic temperature
RTH0 = 10000. 

conn_str = os.path.expandvars("$CONFDB")
conn     = cx_Oracle.connect(conn_str)
e        = conn.cursor()    
e.execute("select distinct DCUHARDID,ADCGAIN0,ADCOFFSET0,I20,I10,TSGAIN,TSOFFSET from viewdcuconversion where PARTITIONNAME='%s' and dcuType='FEH'"%(partName));
q_dcu_conversion = e.fetchall();
df_const = pd.DataFrame(q_dcu_conversion,columns=['dcuHardId','adcGain0','adcOffset0','i20','i10','tsGain','tsOffset'])

        
### Selection of input files
mypath = WorkingDirectory + '/' + args.inputDirectory
string_det = ""
if args.partition == 'TIB':  string_det = "_TI"
elif args.partition == 'TOB':  string_det = "_TO_"
elif args.partition == 'TECP': string_det = "_TP_"
elif args.partition == 'TEC':  string_det = "_TM_"

#### Make list of temperature and threshold files
files_threshold   = [f for f in listdir(mypath) if isfile(join(mypath, f)) and 'Temp' not in f and string_det in f ]
files_temperature = [f for f in listdir(mypath) if isfile(join(mypath, f)) and 'Temp' in f and string_det in f ] 
print('Number of input threshold files = ',len(files_threshold));
print('Number of input temperature files = ',len(files_temperature));

ret_df_vs_run = pd.DataFrame()


for file_ in files_threshold:

    run,part = basename(file_).split('_')[1:3]
    name     = basename(file_).split('.')[0]
    
    ### run selection
    if int(run) < args.runNumberMin: continue;
    if int(run) > args.runNumberMax: continue;
    if int(run) in args.runNumberToExclude: continue;
    ### open the file
    print ('file name ',name,' run number = ',run,' partition = ',part);
    thresholds_array = np.load(join(mypath,file_),allow_pickle=True)
    
    file_temp = [f for f in files_temperature if name in f]    
    if len(file_temp)==0 : 
        print ('No Temperature file found for run ',run)
        continue;

    if len(file_temp)>1: 
        raise Exception('There are two temperature files! This shoudn t happen!')
    else: 
        temperature_array = np.load(join(mypath,file_temp[0]),allow_pickle=True)
        
    ### runs for which we have both threshold and temperature arrays / measurements
    df_thr  = pd.DataFrame(thresholds_array)
    df_thr  = df_thr.reset_index(drop=True);
    df_thr  = df_thr[df_thr['IsValid'] == 1]
    df_thr  = df_thr[(df_thr['Threshold0'] > 0) & (df_thr['Threshold0'] < 65535.0)]
    df_thr  = df_thr[(df_thr['Threshold1'] > 0) & (df_thr['Threshold1'] < 65535.0)]
    df_thr  = df_thr[(df_thr['Threshold2'] > 0) & (df_thr['Threshold2'] < 65535.0)]
    df_thr  = df_thr[(df_thr['Threshold3'] > 0) & (df_thr['Threshold3'] < 65535.0)]    

    ### exclude empty dataframes
    if df_thr.shape[0]==0:
        print ("Empty threshold dataframe for run --> ",run);
        continue;

    def condition(row):
      if row['Gain'] == 0: return row['Threshold0'];
      elif row['Gain'] == 1: return row['Threshold1'];
      elif row['Gain'] == 2: return row['Threshold2'];
      else : return row['Threshold3'];

    df_thr['SelectedThreshold'] = df_thr.apply(condition,axis=1);
    df_thr  = df_thr.dropna(subset=[args.thresholdGain])
    
    df_temp = pd.DataFrame(temperature_array)
    df_temp = df_temp.dropna(subset=[args.dcuChannel])
    df_temp = df_temp[((df_temp.tsilicon > 0) & (df_temp.tsilicon < 4096))];
    df_temp = df_temp[((df_temp.thybrid > 0) & (df_temp.thybrid < 4096))];
    df_temp = df_temp[((df_temp.tdcu > 0) & (df_temp.tdcu < 4096))];
    df_temp = df_temp.reset_index(drop=True);
    
    if df_temp.shape[0]==0:
        print ("Empty temperature dataframe for run --> ",run);
        continue;
        
    ### Grouping the various fibers belonging to the same module and take the average
    df_thr = df_thr.groupby('dcuHardId',as_index=False);
    df_thr_m = df_thr[['dcuHardId','SelectedThreshold','Threshold0','Threshold1','Threshold2','Threshold3','Detector','Layer','Cl','Detid','Mod','luminosity']].mean();
    df_thr_e = df_thr[['dcuHardId','SelectedThreshold','Threshold0','Threshold1','Threshold2','Threshold3']].agg([error_median]);
    df_thr_m = df_thr_m.reset_index(drop=True) 
    df_thr_e = df_thr_e.reset_index(drop=True) 
    df_thr   = df_thr_m;
    df_thr["SelectedThresholdErr"] = df_thr_e.SelectedThreshold;
    df_thr["Threshold0Err"] = df_thr_e.Threshold0;
    df_thr["Threshold1Err"] = df_thr_e.Threshold1;
    df_thr["Threshold2Err"] = df_thr_e.Threshold2;
    df_thr["Threshold3Err"] = df_thr_e.Threshold3;
    df_thr["Ring"] = df_thr.Mod*0.1;
    df_thr["Ring"] = df_thr["Ring"].astype('int');    
    df_thr   = df_thr[(df_thr["Threshold0Err"]>0) & (df_thr["Threshold1Err"]>0) & (df_thr["Threshold2Err"]>0) & (df_thr["Threshold3Err"]>0)];

    if df_thr.shape[0]==0:
        print ("Empty threshold dataframe for run --> ",run);
        continue;

    ### if there are multiple temperature measurements for a given device in the time range considered, take the median value
    df_temp = df_temp.groupby('dcuHardId');
    df_temp_m = df_temp[['dcuHardId','tsilicon','thybrid','tdcu']].median();
    df_temp_e = df_temp[['dcuHardId','tsilicon','thybrid','tdcu']].agg([error_median]);
    df_temp_m = df_temp_m.reset_index(drop=True) 
    df_temp_e = df_temp_e.reset_index(drop=True) 
    df_temp = df_temp_m;
    df_temp['tsiliconErr'] = df_temp_e.tsilicon;
    df_temp['thybridErr']  = df_temp_e.thybrid;
    df_temp['tdcuErr']     = df_temp_e.tdcu;
    ### Merge temperature and threshold data frames
    df_all  = pd.merge(df_thr,df_temp, on='dcuHardId', how='outer') 
    df_all  = df_all.dropna(subset=[args.thresholdGain,args.dcuChannel])    
    df_all  = df_all.reset_index(drop=True) 
    df_all  = pd.merge(df_all,df_const,on='dcuHardId', how='outer')
    df_all  = df_all.dropna(subset=['tsGain'])
    df_all  = df_all.dropna(subset=[args.thresholdGain,args.dcuChannel])    
    df_all  = df_all.reset_index(drop=True)    

    if df_all.shape[0]==0:
        print ("Empty overall dataframe for run --> ",run);
        continue;
    
    #tranform temperature from ADC count to celsius
    df_all['Threshold'] = df_all[args.thresholdGain];
    df_all['ThresholdErr1sup'] = df_all[args.thresholdGain]+df_all[args.thresholdGain+"Err"];
    df_all['ThresholdErr1sdw'] = df_all[args.thresholdGain]-df_all[args.thresholdGain+"Err"];
    if args.dcuChannel == "tsilicon":        
        if args.partition == "TOB" or args.partition == "TECP" or args.partition == "TECM":                
            df_all['val'] = np.log(2.*(abs(df_all[args.dcuChannel]-df_all.adcOffset0))/(df_all.i20*df_all.adcGain0*RTH0));
            df_all['valErr1sup'] = np.log(2.*(abs(df_all[args.dcuChannel]+df_all[args.dcuChannel+"Err"]-df_all.adcOffset0))/(df_all.i20*df_all.adcGain0*RTH0));
            df_all['valErr1sdw'] = np.log(2.*(abs(df_all[args.dcuChannel]-df_all[args.dcuChannel+"Err"]-df_all.adcOffset0))/(df_all.i20*df_all.adcGain0*RTH0));
        else:
            df_all['val'] = np.log(1.*(abs(df_all[args.dcuChannel]-df_all.adcOffset0))/(df_all.i20*df_all.adcGain0*RTH0));
            df_all['valErr1sup'] = np.log(1.*(abs(df_all[args.dcuChannel]+df_all[args.dcuChannel+"Err"]-df_all.adcOffset0))/(df_all.i20*df_all.adcGain0*RTH0));
            df_all['valErr1sdw'] = np.log(1.*(abs(df_all[args.dcuChannel]-df_all[args.dcuChannel+"Err"]-df_all.adcOffset0))/(df_all.i20*df_all.adcGain0*RTH0));
        df_all['Temperature'] = 1./(1./T25+df_all.val/BETA)-T0;
        df_all['TemperatureErr1sup'] = 1./(1./T25+df_all.valErr1sup/BETA)-T0;
        df_all['TemperatureErr1sdw'] = 1./(1./T25+df_all.valErr1sdw/BETA)-T0;

    elif  args.dcuChannel == "thybrid":
        df_all['val'] = np.log(abs((df_all[args.dcuChannel]-df_all.adcOffset0))/(df_all.adcGain0*df_all.i10*RTH0));    
        df_all['valErr1sup'] = np.log((abs(df_all[args.dcuChannel]+df_all[args.dcuChannel+"Err"]-df_all.adcOffset0))/(df_all.adcGain0*df_all.i10*RTH0));        
        df_all['valErr1sdw'] = np.log((abs(df_all[args.dcuChannel]-df_all[args.dcuChannel+"Err"]-df_all.adcOffset0))/(df_all.adcGain0*df_all.i10*RTH0));        
        df_all['Temperature'] = 1./(1./T25+df_all.val/BETA)-T0;
        df_all['TemperatureErr1sup'] = 1./(1./T25+df_all.valErr1sup/BETA)-T0;
        df_all['TemperatureErr1sdw'] = 1./(1./T25+df_all.valErr1sdw/BETA)-T0;

    elif  args.dcuChannel == "tdcu":
        df_all['val'] = (abs(df_all[args.dcuChannel]-df_all.tsOffset))/df_all.tsGain;            
        df_all['valErr1sup'] = (abs(df_all[args.dcuChannel]+df_all[args.dcuChannel+"Err"]-df_all.tsOffset))/df_all.tsGain;            
        df_all['valErr1sdw'] = (abs(df_all[args.dcuChannel]-df_all[args.dcuChannel+"Err"]-df_all.tsOffset))/df_all.tsGain;            
        df_all['Temperature'] = df_all.val+T25-T0;
        df_all['TemperatureErr1sup'] = df_all.valErr1sup+T25-T0;
        df_all['TemperatureErr1sdw'] = df_all.valErr1sdw+T25-T0;        
    else:
        print("Temperature not found --> skip");
        continue;

    df_all = df_all.dropna(subset=['Temperature','TemperatureErr1sup','TemperatureErr1sdw'])
    df_all = df_all.reset_index(drop=True)
    df_all['ThresholdRescaled'] = df_all.Threshold*np.exp(((args.referenceTemperature+T0)-(df_all.Temperature+T0))/TC)
    df_all['ThresholdRescaledErr1sup'] = df_all.ThresholdErr1sup*np.exp(((args.referenceTemperature+T0)-(df_all.TemperatureErr1sup+T0))/TC)
    df_all['ThresholdRescaledErr1sdw'] = df_all.ThresholdErr1sdw*np.exp(((args.referenceTemperature+T0)-(df_all.TemperatureErr1sdw+T0))/TC)    

    #### Plot tracker map for threshold and Temperature
    if args.dumpTrackerMap:
        df_forMap       = df_all.loc[:,('Detid','Threshold')]
        df_forMap.Detid = df_forMap.Detid.astype(int)        
        df_forMap.to_csv(outputTrackerMapPath+"/trackerMap.txt",sep=' ',index=False,header=False)
        os.system(('print_TrackerMap '+outputTrackerMapPath+'/trackerMap.txt \"Threshold Current\" '+outputTrackerMapPath+"/threshold_map_%s_%s.png 1250 False False 2 6")%(run,args.partition));

        df_forMap = df_all.loc[:,('Detid','Temperature')]
        df_forMap.Detid = df_forMap.Detid.astype(int) 
        df_forMap.to_csv(outputTrackerMapPath+"/trackerMap.txt", sep=' ', index=False, header=False); 
        os.system(('print_TrackerMap '+outputTrackerMapPath+'/trackerMap.txt \"Temperature\" '+outputTrackerMapPath+"/temperature_map_%s_%s.png 1250 False False -15 25")%(run,args.partition));
        
        df_forMap       = df_all.loc[:,('Detid','ThresholdRescaled')]
        df_forMap.Detid = df_forMap.Detid.astype(int)        
        df_forMap.to_csv(outputTrackerMapPath+"/trackerMap.txt",sep=' ',index=False,header=False)
        os.system(('print_TrackerMap '+outputTrackerMapPath+'/trackerMap.txt \"Threshold Current\" '+outputTrackerMapPath+"/threshold_map_%s_%s_rescaled.png 1250 False False 2 6")%(run,args.partition));
        os.system("rm "+outputTrackerMapPath+"/*.coor");

    ### Regroup in TIB layer 5,6,7,8 into a single one
    if args.partition == 'TIB':
        df_all.loc[df_all.Layer==6,"Layer"] = 5
        df_all.loc[df_all.Layer==7,"Layer"] = 5
        df_all.loc[df_all.Layer==8,"Layer"] = 5

    #### Prepare the data frame for the final plot
    if args.partition == 'TOB' or args.partition == 'TIB':
        grouped = df_all.groupby('Layer')
    elif args.partition == 'TECP' or args.partition == 'TECM':
        grouped = df_all.groupby('Ring')

    median_vs_run  = grouped[['luminosity','Threshold','ThresholdRescaled']].median().reset_index();    
    median_error_vs_run = grouped[['ThresholdErr1sup','ThresholdErr1sdw','ThresholdRescaledErr1sup','ThresholdRescaledErr1sdw']].sem().reset_index();
    median_vs_run["ThresholdErr1sup"] = median_error_vs_run.ThresholdErr1sup*1.253;
    median_vs_run["ThresholdErr1sdw"] = median_error_vs_run.ThresholdErr1sdw*1.253;
    median_vs_run["ThresholdRescaledErr1sup"] = median_error_vs_run.ThresholdRescaledErr1sup*1.253;
    median_vs_run["ThresholdRescaledErr1sdw"] = median_error_vs_run.ThresholdRescaledErr1sdw*1.253;
    median_vs_run["run"] = int(run);
    ret_df_vs_run = pd.concat([ret_df_vs_run, median_vs_run])

if args.partition == 'TIB' or args.partition == 'TOB':
  ret_df_vs_lumi = ret_df_vs_run.groupby(['Layer','luminosity'])[['ThresholdRescaled','Threshold']].mean().reset_index();
  ret_df_vs_lumi_err = ret_df_vs_run.groupby(['Layer','luminosity'])[['ThresholdRescaledErr1sup','ThresholdRescaledErr1sdw','ThresholdErr1sup','ThresholdErr1sdw']].sem().reset_index();
elif args.partition == 'TECP' or args.partition == 'TECM':
  ret_df_vs_lumi = ret_df_vs_run.groupby(['Ring','luminosity'])[['ThresholdRescaled','Threshold']].mean().reset_index();
  ret_df_vs_lumi_err = ret_df_vs_run.groupby(['Ring','luminosity'])[['ThresholdRescaledErr1sup','ThresholdRescaledErr1sdw','ThresholdErr1sup','ThresholdErr1sdw']].sem().reset_index();


ret_df_vs_lumi["ThresholdRescaledErr1sup"] = ret_df_vs_lumi_err.ThresholdRescaledErr1sup;
ret_df_vs_lumi["ThresholdRescaledErr1sdw"] = ret_df_vs_lumi_err.ThresholdRescaledErr1sdw;
ret_df_vs_lumi["ThresholdErr1sup"] = ret_df_vs_lumi_err.ThresholdErr1sup;
ret_df_vs_lumi["ThresholdErr1sdw"] = ret_df_vs_lumi_err.ThresholdErr1sdw;


## rescale to zero threshold
first_run  = ret_df_vs_run.run.min() 
first_lumi = ret_df_vs_lumi.luminosity.min() 

if args.partition == 'TIB' or args.partition == 'TOB':        
    layers = set(ret_df_vs_run.Layer);
elif args.partition == 'TECP' or args.partition == 'TECM':
    layers = set(ret_df_vs_run.Ring);


for layer in layers:
    if args.partition == 'TIB' or args.partition == 'TOB':        

        zero_vs_run           = ret_df_vs_run[(ret_df_vs_run.run == first_run) & (ret_df_vs_run.Layer == layer)].Threshold.min()
        zero_rescaled_vs_run  = ret_df_vs_run[(ret_df_vs_run.run == first_run) & (ret_df_vs_run.Layer == layer)].ThresholdRescaled.min() 
        zero_vs_lumi          = ret_df_vs_lumi[(ret_df_vs_lumi.luminosity == first_lumi) & (ret_df_vs_lumi.Layer == layer)].Threshold.min() 
        zero_rescaled_vs_lumi = ret_df_vs_lumi[(ret_df_vs_lumi.luminosity == first_lumi) & (ret_df_vs_lumi.Layer == layer)].ThresholdRescaled.min() 

        ret_df_vs_run.loc[ret_df_vs_run.Layer == layer, 'ThresholdRescaled'] -= zero_rescaled_vs_run
        ret_df_vs_run.loc[ret_df_vs_run.Layer == layer, 'ThresholdRescaledErr1sup'] -= zero_rescaled_vs_run
        ret_df_vs_run.loc[ret_df_vs_run.Layer == layer, 'ThresholdRescaledErr1sdw'] -= zero_rescaled_vs_run

        ret_df_vs_run.loc[ret_df_vs_run.Layer == layer, 'Threshold'] -= zero_vs_run
        ret_df_vs_run.loc[ret_df_vs_run.Layer == layer, 'ThresholdErr1sup'] -= zero_vs_run
        ret_df_vs_run.loc[ret_df_vs_run.Layer == layer, 'ThresholdErr1sdw'] -= zero_vs_run
 
        ret_df_vs_lumi.loc[ret_df_vs_lumi.Layer == layer, 'ThresholdRescaled'] -= zero_rescaled_vs_lumi
        ret_df_vs_lumi.loc[ret_df_vs_lumi.Layer == layer, 'ThresholdRescaledErr1sup'] -= zero_rescaled_vs_lumi
        ret_df_vs_lumi.loc[ret_df_vs_lumi.Layer == layer, 'ThresholdRescaledErr1sdw'] -= zero_rescaled_vs_lumi

        ret_df_vs_lumi.loc[ret_df_vs_lumi.Layer == layer, 'Threshold'] -= zero_vs_lumi
        ret_df_vs_lumi.loc[ret_df_vs_lumi.Layer == layer, 'ThresholdErr1sup'] -= zero_vs_lumi
        ret_df_vs_lumi.loc[ret_df_vs_lumi.Layer == layer, 'ThresholdErr1sdw'] -= zero_vs_lumi

    elif args.partition == 'TECP' or args.partition == 'TECM':

        zero_vs_run           = ret_df_vs_run[(ret_df_vs_run.run == first_run) & (ret_df_vs_run.Ring == layer)].Threshold.min() 
        zero_rescaled_vs_run  = ret_df_vs_run[(ret_df_vs_run.run == first_run) & (ret_df_vs_run.Ring == layer)].ThresholdRescaled.min() 
        zero_vs_lumi          = ret_df_vs_lumi[(ret_df_vs_lumi.luminosity == first_lumi) & (ret_df_vs_lumi.Ring == layer)].Threshold.min() 
        zero_rescaled_vs_lumi = ret_df_vs_lumi[(ret_df_vs_lumi.luminosity == first_lumi) & (ret_df_vs_lumi.Ring == layer)].ThresholdRescaled.min() 

        ret_df_vs_run.loc[ret_df_vs_run.Ring == layer, 'ThresholdRescaled'] -= zero_rescaled_vs_run
        ret_df_vs_run.loc[ret_df_vs_run.Ring == layer, 'ThresholdRescaledErr1sup'] -= zero_rescaled_vs_run
        ret_df_vs_run.loc[ret_df_vs_run.Ring == layer, 'ThresholdRescaledErr1sdw'] -= zero_rescaled_vs_run

        ret_df_vs_run.loc[ret_df_vs_run.Ring == layer, 'Threshold'] -= zero_vs_run 
        ret_df_vs_run.loc[ret_df_vs_run.Ring == layer, 'ThresholdErr1sup'] -= zero_vs_run 
        ret_df_vs_run.loc[ret_df_vs_run.Ring == layer, 'ThresholdErr1sdw'] -= zero_vs_run 

        ret_df_vs_lumi.loc[ret_df_vs_lumi.Ring == layer, 'ThresholdRescaled'] -= zero_rescaled_vs_lumi
        ret_df_vs_lumi.loc[ret_df_vs_lumi.Ring == layer, 'ThresholdRescaledErr1sup'] -= zero_rescaled_vs_lumi
        ret_df_vs_lumi.loc[ret_df_vs_lumi.Ring == layer, 'ThresholdRescaledErr1sdw'] -= zero_rescaled_vs_lumi

        ret_df_vs_lumi.loc[ret_df_vs_lumi.Ring == layer, 'Threshold'] -= zero_vs_lumi
        ret_df_vs_lumi.loc[ret_df_vs_lumi.Ring == layer, 'ThresholdErr1sup'] -= zero_vs_lumi
        ret_df_vs_lumi.loc[ret_df_vs_lumi.Ring == layer, 'ThresholdErr1sdw'] -= zero_vs_lumi

### style option for the final plot
print("Produce final plot for "+args.partition);

x_min_run = ret_df_vs_run.run.min() 
x_max_run = ret_df_vs_run.run.max() 
x_min_lumi = ret_df_vs_lumi.luminosity.min() 
x_max_lumi = ret_df_vs_lumi.luminosity.max() 

y_min_run_thr = ret_df_vs_run.Threshold.min() 
y_max_run_thr = ret_df_vs_run.Threshold.max() 
y_min_lumi_thr = ret_df_vs_lumi.Threshold.min() 
y_max_lumi_thr = ret_df_vs_lumi.Threshold.max() 

y_min_run_thr_rescaled = ret_df_vs_run.ThresholdRescaled.min() 
y_max_run_thr_rescaled = ret_df_vs_run.ThresholdRescaled.max() 
y_min_lumi_thr_rescaled = ret_df_vs_lumi.ThresholdRescaled.min() 
y_max_lumi_thr_rescaled = ret_df_vs_lumi.ThresholdRescaled.max() 

graph_thr_vs_run = [];
graph_thr_vs_lumi = [];
graph_thr_rescaled_vs_run = [];
graph_thr_rescaled_vs_lumi = [];

for label, df in ret_df_vs_run.groupby('Layer') if args.partition == 'TIB' or args.partition == 'TOB' else ret_df_vs_run.groupby('Ring'):         
    x_val = df.run.values.tolist();
    y_val = df.Threshold.tolist();
    y_val_err_1sup = df.ThresholdErr1sup.tolist()
    y_val_err_1sdw = df.ThresholdErr1sdw.tolist()
    z_val = df.ThresholdRescaled.tolist();
    z_val_err_1sup = df.ThresholdRescaledErr1sup.tolist();
    z_val_err_1sdw = df.ThresholdRescaledErr1sdw.tolist();
    graph_1 = ROOT.TGraphAsymmErrors();
    graph_2 = ROOT.TGraphAsymmErrors();
    for i in range(0,len(df.index)):
      graph_1.SetPoint(i,x_val[i],y_val[i]);
      if label == 1:
        print(i," ",x_val[i]," ",y_val[i]," ",z_val[i])
      graph_1.SetPointError(i,0,0,y_val_err_1sdw[i],y_val_err_1sup[i]);
      graph_2.SetPoint(i,x_val[i],z_val[i]);
      graph_1.SetPointError(i,0,0,z_val_err_1sdw[i],z_val_err_1sup[i]);
    graph_thr_vs_run.append(graph_1)
    graph_thr_rescaled_vs_run.append(graph_2)

for label, df in ret_df_vs_lumi.groupby('Layer') if args.partition == 'TIB' or args.partition == 'TOB' else ret_df_vs_lumi.groupby('Ring'):         
    x_val = df.luminosity.values.tolist();
    y_val = df.Threshold.tolist();
    y_val_err_1sup = df.ThresholdErr1sup.tolist()
    y_val_err_1sdw = df.ThresholdErr1sdw.tolist()
    z_val = df.ThresholdRescaled.tolist();
    z_val_err_1sup = df.ThresholdRescaledErr1sup.tolist();
    z_val_err_1sdw = df.ThresholdRescaledErr1sdw.tolist();
    graph_1 = ROOT.TGraphAsymmErrors();
    graph_2 = ROOT.TGraphAsymmErrors();
    for i in range(0,len(df.index)):
      graph_1.SetPoint(i,x_val[i],y_val[i]);
      graph_1.SetPointError(i,0,0,y_val_err_1sdw[i],y_val_err_1sup[i]);
      graph_2.SetPoint(i,x_val[i],z_val[i]);
      graph_1.SetPointError(i,0,0,z_val_err_1sdw[i],z_val_err_1sup[i]);
    graph_thr_vs_lumi.append(graph_1)
    graph_thr_rescaled_vs_lumi.append(graph_2)

### plotting part
setTDRStyle();
plotGraph(graph_thr_vs_run,x_min_run,x_max_run,y_min_run_thr,y_max_run_thr,x_years_run_tk,args.thresholdGain+"_"+args.partition+"_vs_run",x_years_run_pp);
plotGraph(graph_thr_vs_lumi,x_min_lumi,x_max_lumi,y_min_lumi_thr,y_max_lumi_thr,x_years_lumi,args.thresholdGain+"_"+args.partition+"_vs_lumi");
plotGraph(graph_thr_rescaled_vs_run,x_min_run,x_max_run,y_min_run_thr_rescaled,y_max_run_thr_rescaled,x_years_run_tk,args.thresholdGain+"_rescaled_"+args.dcuChannel+"_"+args.partition+"_vs_run",x_years_run_pp);
plotGraph(graph_thr_rescaled_vs_lumi,x_min_lumi,x_max_lumi,y_min_lumi_thr_rescaled,y_max_lumi_thr_rescaled,x_years_lumi,args.thresholdGain+"_rescaled_"+args.dcuChannel+"_"+args.partition+"_vs_lumi");

