import pandas as pd
from pdb import set_trace
import numpy as np
import os
import csv
import matplotlib.pyplot as plt
import os,sys,getopt,glob,cx_Oracle,subprocess
import argparse
import os.path
import time,calendar
from datetime import datetime
from datetime import timedelta

# OPTSCAN query from tkCommissioner
def printQuery(analysisId, runNumber): 
    myQuery = ( "select distinct "
                "TKF.DETECTOR Detector,"
                "TKF.SIDE Side,"
                "TKF.LAYER Layer,"
                "TKF.CL Cl,"
                "TKF.CR Cr,"
                "TKF.MOD Mod,"
                "TKF.RACK Rack,"
                "TKF.DETID Detid,"
                "DCU.DCUHARDID dcuHardId,"
                "FEC.CRATESLOT Crate,"
                "FEC.FECSLOT   Fec,"
                "RING.RINGSLOT Ring,"
                "CCU.CCUADDRESS Ccu,"
                "CCU.ARRANGEMENT CcuArrangement,"
                "HYBRID.I2CCHANNEL I2CChannel,"
                "DEVICE.I2CADDRESS I2CAddress,"
                "ROUND((DEVICE.I2CADDRESS-.5)/2)-16 lasChan,"
                " AOS.DEVICEID DeviceId,"
                " AOS.FEDID FedId,"
                " AOS.FEUNIT FeUnit,"
                " AOS.FECHAN FeChan,"
                " AOS.FEDAPV FeApv,"
                " AOS.GAIN Gain,"
                " AOS.BIAS0 Bias0,"
                " AOS.MEASGAIN0 MeasGain0,"
                " AOS.ZEROLIGHT0 ZeroLight0,"
                " AOS.LINKNOISE0 LinkNoise0,"
                " AOS.LIFTOFF0 LiftOff0,"
                " AOS.TICKHEIGHT0 TickHeight0,"
                " AOS.THRESHOLD0 Threshold0,"
                " AOS.BASELINESLOP0 BaselinesSlop0,"
                " AOS.BIAS1 Bias1,"
                " AOS.MEASGAIN1 MeasGain1,"
                " AOS.ZEROLIGHT1 ZeroLight1,"
                " AOS.LINKNOISE1 LinkNoise1,"
                " AOS.LIFTOFF1 LiftOff1,"
                " AOS.TICKHEIGHT1 TickHeight1,"
                " AOS.THRESHOLD1 Threshold1,"
                " AOS.BASELINESLOP1 BaselinesSlop1,"
                " AOS.BIAS2 Bias2,"
                " AOS.MEASGAIN2 MeasGain2,"
                " AOS.ZEROLIGHT2 ZeroLight2,"
                " AOS.LINKNOISE2 LinkNoise2,"
                " AOS.LIFTOFF2 LiftOff2,"
                " AOS.TICKHEIGHT2 TickHeight2,"
                " AOS.THRESHOLD2 Threshold2,"
                " AOS.BASELINESLOP2 BaselinesSlop2,"
                " AOS.BIAS3 Bias3,"
                " AOS.MEASGAIN3 MeasGain3,"
                " AOS.ZEROLIGHT3 ZeroLight3,"
                " AOS.LINKNOISE3 LinkNoise3,"
                " AOS.LIFTOFF3 LiftOff3,"
                " AOS.TICKHEIGHT3 TickHeight3,"
                " AOS.THRESHOLD3 Threshold3,"
                " AOS.BASELINESLOP3 BaselinesSlop3,"
                #" CASE WHEN AOS.GAIN=0 THEN AOS.THRESHOLD0 WHEN AOS.GAIN=1 THEN AOS.THRESHOLD1 WHEN AOS.GAIN=2 THEN AOS.THRESHOLD2 ELSE AOS.THRESHOLD3 END Threshold,"
                " AOS.ISVALID IsValid"
                " from"
                " ANALYSISOPTOSCAN AOS join"
                " ANALYSIS on AOS.ANALYSISID = ANALYSIS.ANALYSISID join"
                " RUN on RUN.RUNNUMBER  = ANALYSIS.RUNNUMBER  join"
                " STATEHISTORY on STATEHISTORY.STATEHISTORYID = RUN.STATEHISTORYID join"
                " DEVICE on AOS.DEVICEID=DEVICE.DEVICEID join"
                " HYBRID on DEVICE.HYBRIDID=HYBRID.HYBRIDID join"
                " CCU on HYBRID.CCUID=CCU.CCUID join"
                " RING on CCU.RINGID=RING.RINGID join"
                " FEC on RING.FECID=FEC.FECID   join"
                " DEVICE b on b.HYBRIDID = HYBRID.HYBRIDID join"
                " DCU on b.DEVICEID = DCU.DEVICEID and DCU.VERSIONMAJORID = STATEHISTORY.FECVERSIONMAJORID and DCU.VERSIONMINORID = STATEHISTORY.FECVERSIONMINORID left outer join"
                " tk_fibers tkf on DCU.DCUHARDID = tkf.dcuid and "
                " mod( AOS.FECHAN,3) = mod(fiber,3) "
                " where "
                " AOS.ANALYSISID=%s and RUN.RUNNUMBER=%s")
    #print(myQuery %(analysisId,runNumber))
    return myQuery %(analysisId,runNumber) 
    

# function used to write in Pandas dataframe a subset of the info downloaded from gain scans analysis (only for a subset of selected runs) 
def writeDF(infos_q3, start_time, delta_time, luminosity):
   
    df = pd.DataFrame(infos_q3)
    df.columns = ['Detector','Side','Layer','Cl','Cr','Mod','Rack','Detid','dcuHardId','Crate','Fec','Ring','Ccu','CcuArrangement','I2CChannel','I2CAddress','lasChan','DeviceId','FedId','FedUnit','FeChan','FeApv','Gain','Bias0','MeasGain0','ZeroLight0','LinkNoise0','LiftOff0','TickHeight0','Threshold0','BaselinesSlop0','Bias1','MeasGain1','ZeroLight1','LinkNoise1','LiftOff1','TickHeight1','Threshold1','BaselinesSlop1','Bias2','MeasGain2','ZeroLight2','LinkNoise2','LiftOff2','TickHeight2','Threshold2','BaselinesSlop2','Bias3','MeasGain3','ZeroLight3','LinkNoise3','LiftOff3','TickHeight3','Threshold3','BaselinesSlop3','IsValid']
    df['start_time'] = start_time
    df['delta_time'] = delta_time
    df['luminosity'] = luminosity
    rr = df.to_records(index=False)
    rr.dtype.names = [str(i) for i in rr.dtype.names]
    return rr

########## Parsing options to riun the code
parser = argparse.ArgumentParser()
parser.add_argument('--partition', help='Specify partition: TEC, TIB, TID or TOB',  type=str, choices=['TECM','TECP','TIB','TOB'])
parser.add_argument('--minTimeGlobalRuns', help  = 'Minimum time extension for global runs --> in seconds',  type=int, default=600, nargs='?')
parser.add_argument('--timeWindowMatching', help = 'Max time range in order to match a gain scan with either a pedestal or a GR',  type=int, default=12*3600, nargs='?')
parser.add_argument('--runNumberMin', help = 'Min run number for good gain scans',  type=int, default=0, nargs='?')
parser.add_argument('--runNumberMax', help = 'Max run number for good gain scans',  type=int, default=999999, nargs='?')
parser.add_argument('--outputDirectory', help = 'Name of the output direcotry with numpy files',  type=str, default="OPTOSCANS", nargs='?')
parser.add_argument('--forceFileRecreation', help = 'Foce re-creation of the files in case they already exsists',  action='store_true', default=False)
parser.add_argument('--verbose', action='store_true', default=False, help='verbosity flag');
parser.add_argument('--lumiByDayFile', help = 'recorded luminosity from https://cern.ch/cmslumi/publicplots/lumiByDay.csv', type=str, default='input_files/lumiByDay.csv',nargs='?')
parser.add_argument('--dcuquerytimewindown', help = 'in case no temperature measurements are found by matching with a run use a standard \pm window in time (sec)', type=int, default=12*3600, nargs='?')
parser.add_argument('--mindcuthresholdwithreadings', help = 'in case there are less dcu readings than this fraction of DCUs abort', type=float, default=0.3, nargs='?')
options = parser.parse_args()

### Partition Name
print("saveThresholdFiles.py execution for partition ",options.partition);
partName = ''
if    options.partition == 'TECM': 
    partName = 'TM_09-JUN-2009_1'
elif  options.partition == 'TIB':  
    partName = 'TI_27-JAN-2010_2'  
elif  options.partition == 'TECP': 
    partName = 'TP_09-JUN-2009_1'
elif  options.partition == 'TOB':  
    partName = 'TO_30-JUN-2009_1'
 
#connect to Oracle
conn_str = os.path.expandvars("$CONFDB")
conn     = cx_Oracle.connect(conn_str)
e        = conn.cursor()

##### Query to download all runs for a given partition with RUNMODE 18 --> GlobalRuns 
e.execute("select runPar.string_value,a.RUNNUMBER,runPar.TIME,a.STARTTIME,a.ENDTIME from RUN a join PARTITION b on a.PARTITIONID=b.PARTITIONID join  CMS_RUNINFO.runsession_parameter runPar ON  runPar.RUNNUMBER=a.RUNNUMBER  where a.RUNMODE = 18 and runPar.NAME='CMS.LVL0:TRIGGER_MODE_AT_START' and b.PARTITIONNAME='%s' ORDER BY a.RUNNUMBER" %partName)
infos_GR = e.fetchall()

time_GR=[]
for result in infos_GR:
    if not ("collisions" in result[0] or "cosmics" in result[0] or "circulating" in result[0]): continue; ### take only interesting global runs
    if result[3] and result[4]: ## check whether start and end time are in the DB table
        time_GR.append((result[0],result[1],result[3],result[4],np.abs(result[4]-result[3]).total_seconds()))
    else:
        print ("For run ",result[1]," either no start or stop time found");

dfGR_time = pd.DataFrame(time_GR,columns=['runType','runNumber','timeGR_start','timeGR_stop','deltaT'])        
print ("Found a total of ",len(time_GR)," Global runs matching the required conditions")

##### Query to download all runs for a given partition with RUNMODE 2 --> Pedestal runs
e.execute("select runPar.string_value,a.RUNNUMBER,runPar.TIME,a.STARTTIME,a.ENDTIME from RUN a join PARTITION b on a.PARTITIONID=b.PARTITIONID join CMS_RUNINFO.runsession_parameter runPar ON  runPar.RUNNUMBER=a.RUNNUMBER  where a.RUNMODE = 2 and b.PARTITIONNAME='%s' ORDER BY a.RUNNUMBER" %partName)
infos_PED = e.fetchall()

time_PED=[]
for result in infos_PED:
    if "http://cmsrc-tracker.cms:18000" in result[0]: ### check the entry in the RUN table written by the levelOne function manager
        if result[3] and result[4]: ## check whether start and end time are in the DB table
            time_PED.append((result[0],result[1],result[3],result[4],np.abs(result[4]-result[3]).total_seconds()))
        else:
            print ("For run ",result[1]," either no start or stop time found");
                

dfPED_time = pd.DataFrame(time_PED,columns=['runType','runNumber','timePED_start','timePED_stop','deltaT'])
print ("Found a total of ",len(time_PED)," pedestal runs matching the required conditions")

##### Query to download the gain scan runs --> RUN MODE 4
e.execute("select runPar.string_value,a.RUNNUMBER,runPar.TIME,a.STARTTIME,a.ENDTIME from RUN a join PARTITION b on a.PARTITIONID=b.PARTITIONID join CMS_RUNINFO.runsession_parameter runPar ON  runPar.RUNNUMBER=a.RUNNUMBER  where a.RUNMODE = 4 and b.PARTITIONNAME='%s' ORDER BY a.RUNNUMBER" %partName)
infos_GAIN = e.fetchall()

runs_GAIN = []
for result in infos_GAIN:
    if "http://cmsrc-tracker.cms:18000" in result[0]:
        if result[4] and result[3]:
            runs_GAIN.append((result[0],result[1],result[3],result[4],np.abs(result[4]-result[3]).total_seconds()))
        else:
            print ("For run ",result[1]," either no start or stop time found");

dfGAIN_time = pd.DataFrame(runs_GAIN,columns=['runType','runNumber','timeGAIN_start','timeGAIN_stop','deltaT'])
print ("Found a total of ",len(runs_GAIN)," gain scans matching the required conditions")

####### Divide global run in collisions from global run not in collision
dfGR_time_coll = dfGR_time.loc[(dfGR_time['runType'].str.contains("collisions")) & (dfGR_time['deltaT']>options.minTimeGlobalRuns) ] 
dfGR_time_coll.reset_index(drop=True, inplace=True)
dfGR_time_noColl = dfGR_time.loc[(~dfGR_time['runType'].str.contains("collisions")) & (dfGR_time['deltaT']>options.minTimeGlobalRuns)]
dfGR_time_noColl.reset_index(drop=True, inplace=True)

print ("Number of global runs of type collision passing selection --> ",dfGR_time_coll.shape[0]);
print ("Number of global runs of type cosmics/circulating passing selection --> ",dfGR_time_noColl.shape[0]);

### read lumi by day file:
integrated_luminosity_recorded = 0.;
luminosity_record = [];
if options.lumiByDayFile:
    with open(options.lumiByDayFile,'r') as csvfile:
        csvreader = csv.reader(csvfile,delimiter=',', quotechar='|');
        for idx,row in enumerate(csvreader):
            if idx == 0: continue;
            integrated_luminosity_recorded = integrated_luminosity_recorded + float(row[2])*0.001*0.001*0.001 ## conversion in fb
            luminosity_record.append((datetime.strptime(row[0]+" 00:00:00", '%Y-%m-%d %H:%M:%S'),integrated_luminosity_recorded));                        
dfLumi_time = pd.DataFrame(luminosity_record,columns=['start_time','luminosity'])


#### Save thresholds and info of good Gain scans matched with a run containing temperature measurements
good_runs=[]
n_pedestals = 0;
n_noncollisions = 0;
n_collisions = 0;
#### match gain scans with either pedestal runs or global runs in order to retrive temperature measurements
for run  in runs_GAIN:
    if run[1] and run[2] and run[3] : ### check that the run number exists        

        ## difference in luminosity
        diff_luminosity = np.abs(dfLumi_time.start_time-run[2]);

        ## time difference
        delta_start_pedestal = np.abs(dfPED_time.timePED_start-run[2]).min();
        delta_stop_pedestal  = np.abs(dfPED_time.timePED_stop-run[2]).min();

        if options.verbose: 
            print ('GAIN scan ',run[1],' --> delta_start_pedestal = ',delta_start_pedestal, " delta_stop_pedestal = ",delta_stop_pedestal)
        if delta_start_pedestal.total_seconds() < delta_stop_pedestal.total_seconds() and delta_start_pedestal.total_seconds() < options.timeWindowMatching: 
            diff_start = np.abs(dfPED_time.timePED_start-run[2]);
            good_runs.append(
                (run[0], ## run-mode
                 run[1], ## run-number
                 run[2], ## start-time
                 'Pedestal', 
                 dfPED_time['runNumber'].iloc[diff_start.argmin()], ## run number of the pedestal run matched i.e. minimizes diff_start
                 dfPED_time['timePED_start'].iloc[diff_start.argmin()], ## star time of the pedestal run which minimizes diff_start
                 dfPED_time['deltaT'].iloc[diff_start.argmin()],
                 dfLumi_time['luminosity'].iloc[diff_luminosity.argmin()] ## integrated luminosity
             )) ## time extension of the pedestal run
            n_pedestals = n_pedestals+1;
            continue;
        elif delta_stop_pedestal.total_seconds() < delta_start_pedestal.total_seconds() and delta_stop_pedestal.total_seconds() < options.timeWindowMatching:
            diff_start = np.abs(dfPED_time.timePED_stop-run[2]);
            good_runs.append(
                (run[0], ## run-mode
                 run[1], ## run-number
                 run[2], ## start-time
                 'Pedestal', 
                 dfPED_time['runNumber'].iloc[diff_start.argmin()], ## run number of the pedestal run matched i.e. minimizes diff_start
                 dfPED_time['timePED_start'].iloc[diff_start.argmin()], ## stop time of the pedestal run which minimizes diff_start
                 dfPED_time['deltaT'].iloc[diff_start.argmin()],
                 dfLumi_time['luminosity'].iloc[diff_luminosity.argmin()] ## integrated luminosity
             )) ## time extension of the pedestal run
            n_pedestals = n_pedestals+1;
            continue;
        else:
            delta_start_noColl = np.abs(dfGR_time_noColl.timeGR_start-run[2]).min()
            delta_stop_noColl  = np.abs(dfGR_time_noColl.timeGR_stop-run[2]).min()
            if options.verbose: 
                print ('GAIN scan ',run[1],' --> delta_start_noColl = ',delta_start_noColl, " delta_stop_noColl = ",delta_stop_noColl)
            if delta_start_noColl.total_seconds() < delta_stop_noColl.total_seconds() and delta_start_noColl.total_seconds() < options.timeWindowMatching: 
                diff_start = np.abs(dfGR_time_noColl.timeGR_start-run[2]);
                good_runs.append((run[0], ## run-mode                                                                                                                                            
                                  run[1], ## run-number                                                                                                                                       
                                  run[2], ## start-time                                                                                                                                          
                                  'noCollision',
                                  dfGR_time_noColl['runNumber'].iloc[diff_start.argmin()], ## run number of the non-collision run matched i.e. minimizes diff_start               
                                  dfGR_time_noColl['timeGR_start'].iloc[diff_start.argmin()], ## star time of the non-collision run which minimizes diff_start                    
                                  dfGR_time_noColl['deltaT'].iloc[diff_start.argmin()],
                                  dfLumi_time['luminosity'].iloc[diff_luminosity.argmin()] ## integrated luminosity
                              )) ## time extension of the non-collision run                                                
                n_noncollisions = n_noncollisions+1;
                continue;
            elif delta_stop_noColl.total_seconds() <  delta_start_noColl.total_seconds() and delta_stop_noColl.total_seconds() < options.timeWindowMatching:
                diff_start = np.abs(dfGR_time_noColl.timeGR_stop-run[2]);
                good_runs.append((run[0], ## run-mode                                                                                                                                            
                                  run[1], ## run-number                                                                                                                                       
                                  run[2], ## start-time                                                                                                                                          
                                  'noCollision',
                                  dfGR_time_noColl['runNumber'].iloc[diff_start.argmin()], ## run number of the non-collision run matched i.e. minimizes diff_start               
                                  dfGR_time_noColl['timeGR_start'].iloc[diff_start.argmin()], ## star time of the non-collision run which minimizes diff_start                    
                                  dfGR_time_noColl['deltaT'].iloc[diff_start.argmin()],
                                  dfLumi_time['luminosity'].iloc[diff_luminosity.argmin()] ## integrated luminosity
                              )) ## time extension of the non-collision run                                                
                n_noncollisions = n_noncollisions+1;
                continue;

            else:                
                delta_start_coll = np.abs(dfGR_time_coll.timeGR_start-run[2]).min()
                delta_stop_coll  = np.abs(dfGR_time_coll.timeGR_stop-run[2]).min()
                if options.verbose: 
                    print ('GAIN scan ',run[1],' --> delta_start_coll = ',delta_start_coll, " delta_stop_coll = ",delta_stop_coll)
                if delta_start_coll.total_seconds() < delta_stop_coll.total_seconds() and delta_start_coll.total_seconds() < options.timeWindowMatching: 
                    diff_start = np.abs(dfGR_time_coll.timeGR_start-run[2]);
                    good_runs.append((run[0], ## run-mode                                                                                                                                
                                      run[1], ## run-number                                                                                                                                   
                                      run[2], ## start-time                                                                                                                               
                                      'collision',
                                      dfGR_time_coll['runNumber'].iloc[diff_start.argmin()], ## run number of the collision run matched i.e. minimizes diff_start               
                                      dfGR_time_coll['timeGR_start'].iloc[diff_start.argmin()], ## star time of the collision run which minimizes diff_start                    
                                      dfGR_time_coll['deltaT'].iloc[diff_start.argmin()],
                                      dfLumi_time['luminosity'].iloc[diff_luminosity.argmin()] ## integrated luminosity
                                  )) ## time extension of the collision run                                                  
                    n_collisions = n_collisions+1;
                    continue;
                elif delta_stop_coll.total_seconds() - delta_start_coll.total_seconds() and delta_start_coll.total_seconds() < options.timeWindowMatching:
                    diff_start = np.abs(dfGR_time_coll.timeGR_stop-run[2]);
                    good_runs.append((run[0], ## run-mode                                                                                                                                
                                      run[1], ## run-number                                                                                                                                   
                                      run[2], ## start-time                                                                                                                               
                                      'collision',
                                      dfGR_time_coll['runNumber'].iloc[diff_start.argmin()], ## run number of the collision run matched i.e. minimizes diff_start               
                                      dfGR_time_coll['timeGR_start'].iloc[diff_start.argmin()], ## star time of the collision run which minimizes diff_start                    
                                      dfGR_time_coll['deltaT'].iloc[diff_start.argmin()],
                                      dfLumi_time['luminosity'].iloc[diff_luminosity.argmin()] ## integrated luminosity
                                  )) ## time extension of the collision run                                                  
                    n_collisions = n_collisions+1;
                    continue;                        
                else:
                    print ("No matching found for run ", run[1]);

print ("##### Summary of matching --> n-runs matched with pedestal = ",n_pedestals," over ",len(runs_GAIN))
print ("##### Summary of matching --> n-runs matched with non collision = ",n_noncollisions," over ",len(runs_GAIN))
print ("##### Summary of matching --> n-runs matched with collision = ",n_collisions," over ",len(runs_GAIN))
         
### loop over runs (pedestal or global) that have been matched to a gain-scan within the aferementioned time window
list_selected_runs = [];
list_selected_long_dcu = [];
list_missing_analyses = [];
list_missing_temperatures = [];
number_good_analyses = 0;
number_good_temperatures = 0;

for run in good_runs:
    
    if run[1] < options.runNumberMin or run[1] > options.runNumberMax: 
        continue;    
    if options.verbose: 
        print ("Good GAIN scan ",run[1]," query for temperature values .. ");
    
    list_selected_runs.append(run[1]);

    #### for each gain scan selected for which we have a possible good temperature measurement, we check if the analysis results are in the DB
    e.execute("select max(analysisid), ANALYSISTYPE, RUNNUMBER, PARTITIONNAME from analysis a join partition b on a.PARTITIONID = b.PARTITIONID where PARTITIONNAME='%s' and RUNNUMBER=%s group by ANALYSISTYPE,RUNNUMBER,PARTITIONNAME" %(partName, run[1])) 
    info_q2 = e.fetchall()

    if len(info_q2)==0:
        list_missing_analyses.append(run[1]);
        continue;

    analysisid    = info_q2[0][0]
    partitionName = info_q2[0][3]
    analysisType  = info_q2[0][1]

    if str(analysisType) != 'OPTOSCAN': continue;

    ##### Create file and output directory
    fileName = analysisType+'_'+str(run[1])+'_'+partitionName+'.npy'
    WorkingDirectory = os.getcwd()

    if not os.path.exists(options.outputDirectory):
        os.mkdir(options.outputDirectory)

    OPTOSCAN_path = WorkingDirectory + "/" + options.outputDirectory;
    docs_dir = os.path.expanduser(OPTOSCAN_path)

    ##### in case a file is already presence and not forced, skip
    if os.path.isfile(docs_dir+'/'+fileName) and not options.forceFileRecreation: continue;

    e.execute(printQuery(analysisid,run[1]))
    infos_q3 = e.fetchall()

    arraytosave = writeDF(infos_q3, int(calendar.timegm(run[5].timetuple())), int(run[6]), run[7])
    np.save(os.path.join(docs_dir,fileName),arraytosave)
    number_good_analyses = number_good_analyses + 1; 
    print (fileName,' done...')

    ##### check if the Temperature file for that given run doesn't exist yet
    fileNameTemp = analysisType+'_'+str(run[1])+'_'+partitionName+'_Temp.npy'
    
    if os.path.isfile(docs_dir+'/'+fileNameTemp) and not options.forceFileRecreation: continue;

    ##### Start and end time of the run used to measure temperature
    time_start =  run[5]-timedelta(seconds=3600)
    time_end   =  run[5]+timedelta(seconds=run[6]+3600);
    unixtime_start = calendar.timegm(time_start.timetuple())  
    unixtime_end = calendar.timegm(time_end.timetuple())

    e.execute("select distinct dcuhardid from cms_trk_TKCC.viewalldcuhardid where DCUTYPE='FEH' and partitionname='%s'"%(partName));
    info_qdcu = e.fetchall();

    #e.execute("with pgroup as ( select distinct dcu_id, detid, cable, detector from  CMS_TRK_DCS_CONF.DCUS join cms_trk_TKCC.tk_fibers on dcu_id=tk_fibers.dcuid) Select distinct dcutimestamp, dcuhardid, channel0, channel4, channel7 from cms_trk_TKCC.dcuchanneldata dd join pgroup on (dd.dcuhardid=pgroup.dcu_id)  where dcutimestamp between %s and %s and channel4 <> 0 and channel0 <> 0 and channel7 <> 0 and detector='%s'" %(unixtime_start,unixtime_end,options.partition))
    e.execute("select distinct dcutimestamp, dcuhardid, channel0, channel4, channel7 from cms_trk_TKCC.viewalldcuvalues where dcutimestamp between %s and %s and channel4 <> 0 and channel0 <> 0 and channel7 <> 0 and partitionname='%s' and dcutype='FEH'"%(unixtime_start,unixtime_end,partName))
    info_qTemp = e.fetchall()

    if not info_qTemp or len(info_qTemp) < len(info_qdcu)*options.mindcuthresholdwithreadings:  
        print("Fall back in making a query for DCU values within a mindcuthresholdwithreadings time window");
        time_start =  run[2]-timedelta(seconds=options.dcuquerytimewindown)
        time_end   =  run[2]+timedelta(seconds=options.dcuquerytimewindown)
        unixtime_start = calendar.timegm(time_start.timetuple())  
        unixtime_end   = calendar.timegm(time_end.timetuple())

        #e.execute("with pgroup as ( select distinct dcu_id, detid, cable, detector from  CMS_TRK_DCS_CONF.DCUS join cms_trk_TKCC.tk_fibers on dcu_id=tk_fibers.dcuid) Select distinct dcutimestamp, dcuhardid, channel0, channel4, channel7 from cms_trk_TKCC.dcuchanneldata dd join pgroup on (dd.dcuhardid=pgroup.dcu_id)  where dcutimestamp between %s and %s and channel4 <> 0 and channel0 <> 0 and channel7 <> 0 and detector='%s'" %(unixtime_start,unixtime_end,options.partition))
        e.execute("select distinct dcutimestamp, dcuhardid, channel0, channel4, channel7 from cms_trk_TKCC.viewalldcuvalues where dcutimestamp between %s and %s and channel4 <> 0 and channel0 <> 0 and channel7 <> 0 and partitionname='%s' and dcutype='FEH'"%(unixtime_start,unixtime_end,partName))
        info_qTemp = e.fetchall()        
        if info_qTemp and len(info_qTemp) >= len(info_qdcu)*options.mindcuthresholdwithreadings:
            list_selected_long_dcu.append(run[1]);

    if not info_qTemp: 
        print ('file temp empty for run ',run[1],' ',docs_dir+'/'+fileNameTemp);
        list_missing_temperatures.append(run[1]);
        continue;
    elif len(info_qTemp) < len(info_qdcu)*options.mindcuthresholdwithreadings:        
        print ('file temp has not enough readings for run ',run[1],' ',docs_dir+'/'+fileNameTemp);
        list_missing_temperatures.append(run[1]);
        continue;
    else :
        #save query output to numpy array
        df_Temp = pd.DataFrame(info_qTemp)
        df_Temp.columns = ['dcuTimeStamp','dcuHardId','tsilicon','thybrid', 'tdcu']
        df_Temp['typeRun']   = run[3]
        df_Temp['run_start'] = run[5]
        df_Temp['run_delta'] = run[6]
        rr_Temp = df_Temp.to_records(index=False)
        rr_Temp.dtype.names = [str(i) for i in rr_Temp.dtype.names]
        np.save(os.path.join(docs_dir,fileNameTemp),rr_Temp)
        print (fileNameTemp, ' done...')
        number_good_temperatures = number_good_temperatures + 1;
        continue

conn.close()
print ('The script has finished running.')
print ('Number of missing analyses for good Gain Scans --> ',len(list_missing_analyses),' out of ',len(list_selected_runs));
print ('Number of missing temperatures for good Gain Scans --> ',len( list_missing_temperatures),' out of ',len(list_selected_runs));
print ('Number of good gain scan files --> ',number_good_analyses);
print ('Number of good temperature files --> ',number_good_temperatures);
print ('Number of temperature files obtained from larger DCU window time -->',len(list_selected_long_dcu));
