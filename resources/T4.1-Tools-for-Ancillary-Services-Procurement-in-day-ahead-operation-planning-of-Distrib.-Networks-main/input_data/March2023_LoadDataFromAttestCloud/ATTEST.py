get_ipython().magic('%reset -f')
get_ipython().magic('%cls')

import pandas as pd
import numpy as np
import os, sys, shutil, re, time
from difflib import SequenceMatcher
from win10toast import ToastNotifier



np.set_printoptions(suppress=True, precision=2)
pd.set_option("display.max_rows"   , 1000,
              "display.max_columns", 50,
              "display.width"      , 2000,
              "display.precision"  , 3)
np.set_printoptions(threshold=np.inf, linewidth=np.inf)

PWD = os.getcwd()
cases_path = r"D:\Work\Produce\TESTIFY\SLA\input_data\March2023_LoadDataFromAttestCloud"
origins    = os.path.dirname(cases_path)

# SubFolders = [i for i in os.listdir(cases_path) if os.path.isdir(os.path.join(cases_path,i))]

# You can't use the dict(key = value, key2 = value2) notation because my keys here contain dots.

cases =  {
          # "es_dx_01_2020.ods"        : ["ES DX 03" ,
          #                              ['ES_Dx_03_2030' , 'ES_Dx_03_2040' , 'ES_Dx_03_2050'],
          #                              [1.020408163     , 1.204081633     ,	1.306122449    ]],

          # "pt_dx_01_2020.ods"         : ["PT DX 01" ,
          #                               ['PT_Dx_01_2030' , 'PT_Dx_01_2040'  , 'PT_Dx_01_2050'],
          #                               [1.168849815     , 1.436098655      , 1.704035874    ]],
          #                            ### ['PT_Dx_01_2030'], [1.01]], # To simulate growth factor (1.01) for one year (2030) only
          #                            ### [ 9, 10, 11, 14, 20, 33, 61, 65, 99, 102] ### Where to add loads

          # "hr_dx_01_2020_brown.ods"  :  ["HR DX 01\\BROWN" ,
          #                               ['HR_Dx_01_2030'  , 'HR_Dx_01_2040' , 'HR_Dx_01_2050'],
          #                               [1.175428234      , 1.305375074     , 1.683402245    ]],

          # "hr_dx_01_2020_red.ods"  :    ["HR DX 01\\RED" ,
          #                               ['HR_Dx_01_2030'  , 'HR_Dx_01_2040' , 'HR_Dx_01_2050'],
          #                               [1.175428234      , 1.305375074     , 1.683402245    ]],

          # "hr_dx_01_2020_green.ods"  :  ["HR DX 01\\GREEN" ,
          #                               ['HR_Dx_01_2030'  , 'HR_Dx_01_2040' , 'HR_Dx_01_2050'],
          #                               [1.175428234      , 1.305375074     , 1.683402245    ]],

          "uk_dx_01_2020.ods"        : ["UK DX 01" ,
                                        ['UK_Dx_01_2030' , 'UK_Dx_01_2040' , 'UK_Dx_01_2050' ],
                                        [1.047405509     , 1.310057655     ,	1.444586803    ]],
                                        ### [1.047405509     , 1.310057655     ,	1.93              ]],
         }

LoadGrowthYears = ['2030','2040','2050']
CustomCustomName = ''
# LoadGrowthPD = pd.DataFrame.from_dict(LoadGrowthDict, columns = LoadGrowthYears, orient='index')

Equal_Split    = True
EV_weight      = 1.0
Load_Weight    = 1.0
Scale_RES      = 1.0
Scale_ESS      = 1.0
Line_Derate    = 1
W_nd_WO_Flex   = False
Bus_Volts      = []#[0.95, 1.05]
#%%
def SurveyAssets(FlexExcel):
     StorageSheet   = FlexExcel.parse('Storage_Addt')
     RESSheet       = FlexExcel.parse('RES_Addt')
     LoadsSheet     = FlexExcel.parse('Loads_Addt')
     PLoadProfile   = FlexExcel.parse('pLoad_Profiles_Addt')
     DR_Loads       = LoadsSheet.query('Status == 1').bus_i.ravel()
     DR_Loads_IDX   = np.argwhere(np.isin(PLoadProfile.bus_i, DR_Loads)).ravel()

     Data2Present   = {'C_RES'     : RESSheet.Pmax.sum(),
                       'LoadGrowth': 1,
                       'C_ESS'     : StorageSheet.energy.sum(),
                       'C_Flex'    : PLoadProfile.loc[DR_Loads_IDX,:].max(axis=1).sum(),
                       }
     return Data2Present

FilterName = lambda x: re.sub('_[dtDT]x_0\d{1}_',':',x)

def  DiscipLines(LinesSheetNew, RegexMatch):
     FromTo = LinesSheetNew.loc[:,['fbus','tbus']].values
     FromTo = np.sort(FromTo,axis=1)
     LinesSheetNew.loc[:,['fbus','tbus']] = FromTo
     LinesSheetNew = LinesSheetNew.sort_values(['fbus','tbus'])
     LinesSheetNew = LinesSheetNew.reset_index(drop=True)
     if   len(RegexMatch)==1:
          RegexMatch = [*map(int,RegexMatch[0])]
          BoolTest   = np.logical_not(np.isin(LinesSheetNew.loc[:,['fbus','tbus']].values, RegexMatch).all(axis=1))
          # since not(A*B) = not(A) + not(B), you can get the same effect with:
          # BoolTest = np.isin(LinesSheetNew.loc[:,['fbus','tbus']].values, RegexMatch, invert=True).any(axis=1)
          LinesSheetNew = LinesSheetNew[BoolTest.tolist()]

     elif len(RegexMatch)>1:
          raise ValueError('More than one match were found for RegEx: \'.+Without_(\d+)_(\d+)\'')

     return LinesSheetNew

def  StatusColumn(DF):
     col_index      = np.argwhere('status' == np.char.lower(DF.columns.values.astype(str))).item()
     col_name       = DF.columns[col_index]
     return col_index, col_name


def  DownSample(DF, n_bundle):
     assert n_bundle%1 == 0, 'n_bundle must be an integer'
     # breakpoint()
     if  np.isin('bus_i', DF.columns):
          restore_bus_i = True
          DF_t                = DF.drop(labels=['bus_i'],axis=1).transpose()
     else:
          restore_bus_i=False
          DF_t                = DF.transpose()

     ColsDict            = DF_t.to_dict(orient='list')
     ColsDict            = {key: np.reshape(value,(-1,n_bundle))   for key, value in ColsDict.items()}
     ColsDict            = {key: np.average(value,axis=1).tolist() for key, value in ColsDict.items()}
     ColsValues          = np.vstack([*ColsDict.values()])
     if ColsValues.shape[1] != 24:
          print('Downsampled profile has length %i != 24' %ColsValues.shape[1])
          # breakpoint()
          []

     columns             = ['t%i' %i for i in np.arange(ColsValues.shape[1])+1]
     # DF_t_downsampled    = pd.DataFrame.from_dict(ColsDict, orient='index', columns = ['t%i' %i for i in range(1,25)])
     DF_t_downsampled    = pd.DataFrame(ColsValues, index = ColsDict.keys(), columns=columns)

     if restore_bus_i==True:
          OutArchive          = pd.concat([DF.loc[:,'bus_i'], DF_t_downsampled],axis=1)
     else:
          OutArchive          = pd.concat([DF_t_downsampled],axis=1)
     OutArchive.index    = DF.index

     # print('Load profile downsampled by combining each %i consecutive instants together.\nNew archive is [%i x %i]' %(n_bundle, *OutArchive.shape))
     return OutArchive

#%%
# EV_PV_ESS_File = r'EV-PV-Storage_Data_for_Simulations_updated_17_MARCH.xlsx' Used up to May 19th
EV_PV_ESS_File = r'EV-PV-Storage_Data_for_Simulations_updated_17_APRIL.xlsx' # Added on May 20th

NewData        = os.path.join(cases_path, EV_PV_ESS_File)
NewDataExcel   = pd.ExcelFile(NewData)
# SheetsNames    = NewDataExcel.sheet_names

Data2Present   = {}
CaseNamesFiles = []
breakpoint()

for  key , values in cases.items(): # key = [*cases.keys()][0]; values = [*cases.values()][0]
     print('='*90,('Doing: %s @<%s>' %(key, time.ctime())).center(80),sep='\n')
     CaseFileName   = os.path.join(origins, key)
     FlexFile       = CaseFileName.replace('.ods','_flex.ods')
     FlexExcel      = pd.ExcelFile(FlexFile, engine="odf")
     FlexSheets     = FlexExcel.sheet_names
     # CaseNamesFiles.append(['BaseCase', CaseFileName])

     Data2Present[FilterName(key).replace('.ods','')] = SurveyAssets(FlexExcel)

     # customname     = re.findall('\w{2}_dx_01_2020_?(\w+)?\.ods',key)
     cuntry, customname     = re.findall('(\w{2})_dx_01_2020(_\w+)?\.ods',key)[0]
     AddtDict       = {}

     for  Addt in ['Loads','RES','Storage']:
          AddtDict[Addt] = FlexExcel.parse(Addt+'_Addt')
          AddtDict[Addt].index+=1
     # print(Loads_Addt[Loads_Addt['Status']==1])
     del Addt

     FolderName, yearSheets, GrowthFactors = values
     LoadProfilesFolder = os.path.join(cases_path, FolderName)
     # LoadProFilesNames = [os.path.join(LoadProfilesFolder,i) for i in os.listdir(LoadProfilesFolder) if i.endswith('.xlsx')]
     LoadProFilesNames = [i for i in os.listdir(LoadProfilesFolder) if (i.endswith('.xlsx') == True) and (i.startswith('~$') == False)]

     # ===================================================================
     customfolder   = 'EditedInputFiles'
     NewPath        = os.path.join(PWD,customfolder,FolderName)
     try:
          os.mkdir(NewPath)
     except FileExistsError as _FileExistsError_:
          print('Path "%s" already exists' %NewPath)

     shutil.copyfile(CaseFileName, os.path.join(NewPath, 'BaseCase.ods'))
     shutil.copyfile(FlexFile    , os.path.join(NewPath, 'BaseCase_flex.ods'))
     CaseNamesFiles.append([FolderName.replace(' ','_') + '_BaseCase', os.path.join(NewPath, 'BaseCase.ods'), time.ctime()])

     # ===================================================================

     if W_nd_WO_Flex==True:
          x,y,z = np.meshgrid(yearSheets, LoadProFilesNames, ['WithFlex','WOFlex'], indexing='ij')
          x=x.ravel(); y=y.ravel(); z=z.ravel()
          xy = [*zip(x,y,z)]
     else:
          x,y = np.meshgrid(yearSheets, LoadProFilesNames, indexing='ij')
          x=x.ravel(); y=y.ravel();
          xy = [*zip(x,y)]

     # FilesNames = np.hstack([np.vstack(x.ravel()),np.vstack(y.ravel())])
     # FilesNames = ['_'.join(i) for i in xy]

     LoadPQ={}
     for  i_profile, loadproFileName in enumerate(LoadProFilesNames): # i_profile = 0; loadproFile = LoadProFiles[i_profile]
          # print(pd.read_excel.__doc__)
          loadprofile = pd.read_excel(os.path.join(LoadProfilesFolder, loadproFileName), index_col=[0,1], header=None) # loadprofile.reorder_levels([1,0]);
          LoadPQ[loadproFileName] = dict(P = loadprofile.loc[(slice(None),'P'),:].droplevel(1),
                                         Q = loadprofile.loc[(slice(None),'Q'),:].droplevel(1))
          del loadprofile

     # https://stackoverflow.com/questions/18715688/find-common-substring-between-two-strings
     # print(SequenceMatcher.__doc__)
     # match = SequenceMatcher(None, LoadProFiles).get_matching_blocks()
     PV_ESS_EV = {}

     # ------------------------------------------------------------------------------------
     for  idx, yearSheet in enumerate(yearSheets): # yearSheet = values[0]
          PV_ESS  = pd.read_excel(NewData, sheet_name = yearSheet, header=0, usecols='C:D',nrows=0)
          PV_cap  = PV_ESS.columns[0] * Scale_RES
          ESS_cap = PV_ESS.columns[1] * Scale_ESS

          EV_data = pd.read_excel(NewData, sheet_name = yearSheet, header=1, usecols='F:I', index_col=0)
          EV_data = EV_data.dropna()

          AddWhere_global  = {}
          if   cuntry in ['pt']: # change this to "FALSE" only for PT and ES cases, where EV and PV are concentrated on one bus.
               DivideBy = 1
               # Look in the EV-PV-Storage file, look at the list of nodes and PV-share on each load, and put the PV/EV/ESS there accordingly
               Additions_Distribution_per_node  = pd.read_excel(NewData, sheet_name = yearSheet, header=1, usecols='A:D')
               AddWherePD                    = Additions_Distribution_per_node[Additions_Distribution_per_node.astype(bool).all(axis=1)]
               # AddWhere_Special              = AddWherePD.Node.values.tolist()
               AddWhere_Special              = AddWherePD.Node.dropna().values.tolist()
               AddWhere_global['RES']        = AddWhere_Special
               AddWhere_global['Storage']    = AddWhere_Special
               AddWhere_global['Loads']      = np.sort([ 9, 10, 11, 14, 20, 33, 61, 65, 99, 102]).tolist() # This is for EV

          else:
               if   '2030' in yearSheet:
                    DivideBy = 2
                    year = '2030'
               elif '2040' in yearSheet:
                    DivideBy = 4
                    year = '2040'
               elif '2050' in yearSheet:
                    DivideBy = 6
                    year = '2050'

          PV_ESS_EV[yearSheet] = dict(PV_cap      = PV_cap,
                                      ESS_cap     = ESS_cap,
                                      EV_data     = EV_data,
                                      LoadGrowth  = GrowthFactors[idx],
                                      DivideBy    = DivideBy)
          del EV_data, ESS_cap, PV_cap, DivideBy, yearSheet

     #%% ------------------------------------------------------------------------------------
     for  itr, permutation in enumerate(xy):
          if   W_nd_WO_Flex==True:
               year, loadprofilename, logical_flex = permutation
          else:
               year, loadprofilename               = permutation
               logical_flex = 'WithFlex'

          # if   any([i in year for i in ['2040','2050']]):
          #       print('Skipping {%s} as instructed' %' - '.join(permutation));
          #       continue #I want to generate only XXX cases.

          print('-'*75,
                ('Doing #%4i/%-4i @<%s>: %s' %(itr+1, len(xy), time.ctime(), '_'.join(permutation).replace('.xlsx',''))).center(60),
                '-'*75,
                sep='\n')
          PermutName     = '_'.join(permutation).replace('.xlsx','')
          loadprofile    = LoadPQ[loadprofilename]
          LoadGrowth     = PV_ESS_EV[year]['LoadGrowth']
          PV_cap         = PV_ESS_EV[year]['PV_cap']
          ESS_cap        = PV_ESS_EV[year]['ESS_cap']
          EV_data        = PV_ESS_EV[year]['EV_data']
          DivideBy       = PV_ESS_EV[year]['DivideBy']

          Data2Present[FilterName(PermutName)] = {'C_RES'      : PV_cap,
                                                  'LoadGrowth' : LoadGrowth,
                                                  'C_ESS'      : ESS_cap,
                                                  'C_Flex'     : EV_data.loc[:,'EV load (MW)'].max(axis=0) ,
                                                 }
          # continue #### i enable this "continue" and skip everything below, when i want to see just Data2Present

          # %% ===================================================================
          #   COPY AND EDIT THE NETWORK DATA (MATPOWER DATA) (e.g. data in mpc.branch, and data in mpc.bus)

          DestinationFile= os.path.join(NewPath, PermutName+customname.upper()+CustomCustomName.upper()+'.ods')
          shutil.copyfile(CaseFileName, DestinationFile)
          # CaseNamesFiles.append('("%s_%s","%s"),' %(year,loadprofilename.replace('.xlsx',''),DestinationFile))
          CaseNamesFiles.append(['_'.join([*permutation, CustomCustomName]).replace('.xlsx',''), DestinationFile, time.ctime()])

          ExcelWriter    = pd.ExcelWriter(DestinationFile, mode='w')
          Write2Excel    = lambda Data, sheet: Data.to_excel(ExcelWriter,
                                                             sheet_name=sheet,
                                                             index=False,
                                                             startrow=0,
                                                             startcol=0,
                                                             # float_format="%0.2f",
                                                             )
          EditedSheets   = ['Lines','Buses']
          CaseExcel      = pd.ExcelFile(CaseFileName, engine="odf")
          # breakpoint()
          # ------------------------------------------------------------------
          BusesSheet     = CaseExcel.parse('Buses')
          BusesSheetNew  = BusesSheet.copy(deep=True)
          if   len(Bus_Volts) > 0:
               BusesSheetNew.loc[:,['Vmin','Vmax']] = Bus_Volts

          Write2Excel(BusesSheetNew, 'Buses')

          # ------------------------------------------------------------------
          LinesSheet     = CaseExcel.parse('Lines')
          LinesSheetNew  = LinesSheet.copy(deep=True)
          if Line_Derate!=1:
                LinesSheetNew.loc[:,['rateA','rateB','rateC']] = LinesSheet.loc[:,['rateA','rateB','rateC']] * Line_Derate
          # breakpoint()
          ######################################################################
          # This "Wihtout_..." concerns creating a line contingency by removing a certain line in a meshed network.
          # we implemented this case for the UK grid only, because it has a single mesh/ring formed by lines 1-7 and 2-7
          RegexMatch = re.findall('.+Without_(\d+)_(\d+)',DestinationFile) # Without lines ... to simulate a line-contingency
          if   True: #len(RegexMatch)>0:
               # Sort Bus_From - Bus_To pairs
               LinesSheetNew = DiscipLines(LinesSheetNew, RegexMatch)
          else:
               RegexMatch = []
          ######################################################################

          Write2Excel(LinesSheetNew, 'Lines')

          # ------------------------------------------------------------------
          for  sheet in np.setdiff1d(CaseExcel.sheet_names, EditedSheets):
               if   np.char.startswith(sheet,'\'file:///'):
                    continue
               TheSheet = CaseExcel.parse(sheet)
               Write2Excel(TheSheet, sheet)

          ExcelWriter.close() # CLOSE THE NETWORK DATA FILE

          #%% ===================================================================
          #   COPY AND EDIT THE FLEXIBILITY DATA
          DestinationFlexFile = DestinationFile.replace('.ods','_flex.ods')
          shutil.copyfile(FlexFile, DestinationFlexFile)
          ExcelWriter = pd.ExcelWriter(DestinationFlexFile, mode='w')
          Write2Excel = lambda Data, sheet: Data.to_excel(ExcelWriter,
                                                          sheet_name  = sheet,
                                                          index       = False,
                                                          startrow    = 0,
                                                          startcol    = 0,
                                                          # float_format="%0.7g",
                                                          )

          AddWhere  = {}
          Shares    = {}
          SlackBus  = FlexExcel.parse('Profile_Slack')
          SlackBus  = SlackBus.dropna()

          Summary_Table = FlexExcel.parse('Buses_Addt').loc[:,['bus_i']].sort_values('bus_i',axis=0).reset_index(drop=True).copy(deep=True)
          Summary_Table.loc[:,['PV (MW)','EES (MWh)','EV (%)']]=0
          t1_t24 = np.char.add('t',[*map(str,range(1,25))])

          #%% ### ExcelWriter.close() ### For debugging
          for  Addt in ['Loads','RES','Storage']:
               # columns_names_lower = np.char.lower(AddtDict[Addt].columns.values.astype(str)) # sometimes he uses Status, and sometimes status
               # col_index      = np.argwhere('status' == columns_names_lower).item()
               # col_name       = AddtDict[Addt].columns[col_index]
               col_index, col_name       = StatusColumn(AddtDict[Addt])
               # DivideBy2      = np.min([AddtDict[Addt].shape[0], DivideBy])

               # ------------------------------------------------------------------
               if   Addt in AddWhere_global.keys():
                    AddWhere[Addt] = AddWhere_global[Addt]
                    Shares[Addt]   = np.ones(len(AddWhere[Addt]))/len(AddWhere[Addt])

               else:
                    if   AddtDict[Addt].shape[0] >= DivideBy:
                         # florin told me to split the EV load over X number of loads / buses.
                         # this is only possible if there are already X many buses > DivideBy
                         # so i do this on buses which already have loads.
                         # buses with no existing loads have a very small chance of getting selected for this.
                         Candidates_List     = AddtDict[Addt].bus_i.values
                         prob                = AddtDict[Addt].iloc[:,col_index].values +0.01
                         Prob2               = prob/sum(prob)

                    else:
                         Buses = FlexExcel.parse('Buses_Addt');
                         Buses.loc[:,'Prob'] = 0.01
                         PreferredBuses = AddtDict[Addt][AddtDict[Addt].iloc[:,col_index]==1].bus_i.values
                         PreferredBuses_idx = np.argwhere(np.in1d(Buses.bus_i, PreferredBuses)).ravel()
                         Buses.loc[PreferredBuses_idx, 'Prob'] = 1
                         SlackBus_idx   = np.argwhere(np.in1d(Buses.bus_i, SlackBus.bus_i.values.item())).ravel()
                         Buses.loc[SlackBus_idx, 'Prob'] = 0
                         prob           = Buses.Prob.values
                         Prob2          = prob/sum(prob)
                         # Candidates_List= np.setdiff1d(Buses.bus_i.values, SlackBus.bus_i.values.item())
                         Candidates_List= Buses.bus_i.values

                    AddWhere[Addt] = np.sort(np.random.choice(Candidates_List, DivideBy, p = Prob2, replace=False)).tolist()
                    share_         = np.random.rand(DivideBy) * (1-Equal_Split) + 0.1
                    Shares[Addt]   = (share_ / sum(share_))
                    # breakpoint() # Record Shares & Add_Where in a sheet, into a separate Excel file. call it EV-PV-Storage_Data. take it and go

               # ------------------------------------------------------------------
               if   Addt in 'Loads':
                    NewAddt        = AddtDict[Addt]
                    col_index, col_name = StatusColumn(NewAddt)
                    NewAddt.loc[:,col_name] = 0
                    FlexIDX        = NewAddt.index[np.flatnonzero(np.isin(NewAddt.bus_i , AddWhere[Addt]))]
                    if   logical_flex == 'WOFlex':
                         col_index, col_name       = StatusColumn(AddtDict[Addt])
                         NewAddt.loc[:,col_name] = 0
                    else:
                         NewAddt.loc[FlexIDX,'Status'] = 1

                    Write2Excel(NewAddt, 'Loads_Addt')

                    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    p_profile           = FlexExcel.parse('pLoad_Profiles_Addt').sort_values('bus_i',axis=0).reset_index(drop=True).copy(deep=True);
                    q_profile           = FlexExcel.parse('qLoad_Profiles_Addt').sort_values('bus_i',axis=0).reset_index(drop=True).copy(deep=True);

                    if True: # Q: why import then edit? why not set p_profile right away? A: to preserve the format of the PD.DataFrame
                         ## assert loadprofile['P'].shape == p_profile.iloc[:,1:].shape, 'The new P load profile and the exisitng P load profile have different shapes: [%i x %i] =/= [%i x %i] ' %(*loadprofile['P'].shape, *p_profile.iloc[:,1:].shape)
                         ## assert loadprofile['Q'].shape == q_profile.iloc[:,1:].shape, 'The new Q load profile and the exisitng Q load profile have different shapes: [%i x %i] =/= [%i x %i] ' %(*loadprofile['P'].shape, *p_profile.iloc[:,1:].shape)
                         CommonIndices,a2,a3 = np.intersect1d(p_profile.bus_i.values,
                                                              loadprofile['P'].index.values,
                                                              return_indices=True)
                         ExcludedBuses       = np.setdiff1d(loadprofile['P'].index, p_profile.bus_i.values)
                         IgnoredLoad         = loadprofile['P'].loc[ExcludedBuses,:].abs().values.sum()
                         if   IgnoredLoad > 0:
                              print('Some load is not incorporated. Press any key to continue');
                              breakpoint()
                              []

                         # assert np.all(p_profile.bus_i.values[a2] == loadprofile['P'].index.values[a3]), 'sum ting wong'
                         assert len(CommonIndices)>0, 'New load profile and old load profile dont have any common buses (or bus names)'
                         if   p_profile.shape[1] > 25: # 25 because 1 column is bus_i, and 24 columns for time.
                              p_profile                = DownSample(p_profile, 4)
                              q_profile                = DownSample(q_profile, 4)

                         if   loadprofile['P'].shape[1] > 24:
                              p_profile.iloc[a2,1:]    = DownSample(loadprofile['P'].iloc[a3],4).values[:,:24]
                              q_profile.iloc[a2,1:]    = DownSample(loadprofile['Q'].iloc[a3],4).values[:,:24]

                         else:
                              p_profile.iloc[a2,1:]    = loadprofile['P'].iloc[a3].values
                              q_profile.iloc[a2,1:]    = loadprofile['Q'].iloc[a3].values

                    Write2Excel(p_profile, 'old_pLoad_Profiles_Addt')
                    # Write2Excel(q_profile, 'old_qLoad_Profiles_Addt')  # To speed up. i never look at these anyway

                    ev_profile          = EV_data.iloc[:,[0]].values
                    ev_profile_split    = np.reshape(Shares[Addt],(-1,1)) @ ev_profile.T # @ is the matrix multiplicaiton operator

                    p_profile_new       = p_profile.copy(deep=True).iloc[:,:ev_profile.size+1]
                    p_profile_new.iloc[:,1:] = p_profile_new.iloc[:,1:] * LoadGrowth
                    Write2Excel(p_profile_new, 'grown_pLoad_Profiles_WO_EV')

                    FlexIDX             = np.argwhere(np.isin(p_profile_new.bus_i , AddWhere[Addt])).ravel()
                    assert len(FlexIDX)>0,'FlexIDX is empty'
                    EV_profile = p_profile_new.iloc[FlexIDX,:ev_profile.size+1].copy(deep=True)
                    EV_profile.iloc[:,1:ev_profile.size+1] = ev_profile_split

                    p_profile_new.iloc[FlexIDX,1:] = (EV_weight   * ev_profile_split) + \
                                                     (Load_Weight * p_profile_new.iloc[FlexIDX,1:ev_profile.size+1].values)

                    p_profile_new.index  = range(p_profile_new.shape[0])
                    # if   ((p_profile_new.iloc[:,1:] / p_profile.iloc[:,1:]).dropna().abs()<=1).any(axis=None):
                    #      print(p_profile_new.iloc[:,1:] / p_profile.iloc[:,1:])
                    #      # assert LoadGrowth > 1, 'LoadGrowth <= 1'
                    #      print('Funny Load Growth = %g' %LoadGrowth)
                    #      breakpoint()
                    # I disabled this breakpoint because in the HR-Brown case, there are negative loads.

                    Write2Excel(p_profile_new, 'pLoad_Profiles_Addt')
                    Write2Excel(EV_profile, 'P_EV_Profiles_Addt')

                    FlexIDX_             = np.argwhere(np.isin(Summary_Table.bus_i , AddWhere[Addt])).ravel()
                    Summary_Table.iloc[FlexIDX_,-1] = Shares[Addt]

                    FlexIDX    = np.argwhere(np.isin(q_profile.bus_i , AddWhere[Addt])).ravel()
                    PF         = q_profile.iloc[:,1:ev_profile.size+1] / p_profile.iloc[:,1:ev_profile.size+1]
                    PF[np.isnan(PF)]    = 0
                    q_profile_new       = q_profile.copy(deep=True).iloc[:,:ev_profile.size+1]
                    q_profile_new.iloc[:,1:] = q_profile_new.iloc[:,1:] * LoadGrowth
                    # Write2Excel(q_profile_new, 'grown_qLoad_Profiles_WO_EV') # To speed up. i never look at these anyway
                    q_profile_new.iloc[FlexIDX,1:] = p_profile_new.iloc[FlexIDX,1:] * PF.iloc[FlexIDX,:]
                    q_profile_new.index  = range(q_profile_new.shape[0])
                    Write2Excel(q_profile_new, 'qLoad_Profiles_Addt')
                    # del p_profile_new, p_profile, q_profile_new, q_profile
                    # ExcelWriter.close(); os.startfile(DestinationFlexFile) # for debugging

               # ------------------------------------------------------------------
               elif Addt in 'RES':
                    # breakpoint()
                    Buses_Addt = FlexExcel.parse('Buses_Addt');
                    # Write2Excel(Buses_Addt, 'old_Buses_Addt')

                    Buses_Addt.loc[:,'hasGEN']=0
                    Buses_Addt.loc[np.isin(Buses_Addt.bus_i, AddWhere[Addt]),'hasGEN'] = 1
                    Buses_Addt.to_excel(ExcelWriter, sheet_name='Buses_Addt',index=False)

                    Gen_cost       = FlexExcel.parse('Gens_cost_Addt');
                    NonSlackBus    = Gen_cost.bus_i!=SlackBus.bus_i.values.item()
                    ModelRow       = Gen_cost[NonSlackBus].iloc[[0],:]
                    RES_cost       = pd.concat([ModelRow]*DivideBy, axis=0)
                    RES_cost.bus_i = AddWhere[Addt]
                    Gen_cost_new   = pd.concat([Gen_cost[Gen_cost.bus_i == SlackBus.bus_i.values.item()],
                                                RES_cost], axis=0)
                    Gen_cost_new.index  = range(Gen_cost_new.shape[0])
                    Write2Excel(Gen_cost_new, 'Gens_cost_Addt')

                    Gen            = FlexExcel.parse('Gens_Addt');
                    # Write2Excel(Gen, 'old_Gens_Addt')
                    NonSlackBus    = Gen.bus_i!=SlackBus.bus_i.values.item()
                    ModelRow       = Gen[NonSlackBus].iloc[[0],:]
                    ModelRow.Dgtype= 1 # 1 is PV, 2 is wind, 0 is not RES but dispatchable / grid
                    RES_new        = pd.concat([ModelRow]*DivideBy, axis=0)
                    RES_new.bus_i  = AddWhere[Addt]
                    RES_new        = pd.concat([Gen[Gen.bus_i == SlackBus.bus_i.values.item()],
                                                RES_new], axis=0)
                    RES_new.index  = range(RES_new.shape[0])
                    Write2Excel(RES_new, 'Gens_Addt')

                    RES_Addt       = FlexExcel.parse('RES_Addt');
                    # Write2Excel(RES_Addt, 'raw_RES_Addt')
                    ModelRow       = RES_Addt.iloc[[0],:]
                    RES_Addt_new   = pd.concat([ModelRow]*DivideBy, axis=0)
                    RES_Addt_new.bus_i = AddWhere[Addt]
                    # breakpoint()
                    RES_Addt_new.Pmax  = Shares[Addt] * PV_cap
                    RES_Addt_new.Qmax  = RES_Addt_new.Pmax*0.5
                    RES_Addt_new.index = range(RES_Addt_new.shape[0])
                    Write2Excel(RES_Addt_new, 'RES_Addt')

                    p_min_profile  = FlexExcel.parse('pGen_Profiles_Min_Addt');
                    # Write2Excel(p_min_profile, 'raw_pGen_Profiles_Min_Addt')
                    NonSlackBus    = p_min_profile.bus_i!=SlackBus.bus_i.values.item()
                    ModelRow       = p_min_profile[NonSlackBus].iloc[[0],:]
                    p_min_new      = pd.concat([ModelRow]*DivideBy, axis=0)
                    p_min_new.bus_i= AddWhere[Addt]
                    p_min_new      = pd.concat([p_min_profile[NonSlackBus == False],
                                                p_min_new], axis=0)
                    p_min_new.index= range(p_min_new.shape[0])
                    Write2Excel(p_min_new, 'pGen_Profiles_Min_Addt')

                    p_max_profile  = FlexExcel.parse('pGen_Profiles_Max_Addt');
                    # Write2Excel(p_max_profile, 'raw_pGen_Profiles_Max_Addt')
                    NonSlackBus    = p_max_profile.bus_i!=SlackBus.bus_i.values.item()
                    ModelRow       = p_max_profile[NonSlackBus].iloc[[0],:]
                    p_max_new      = pd.concat([ModelRow]*DivideBy, axis=0)
                    p_max_new.bus_i= AddWhere[Addt]
                    p_max_new.iloc[:,1:] = np.tile(np.reshape(Shares[Addt] * PV_cap,(-1,1)), (1,24))
                    p_max_new      = pd.concat([p_max_profile[NonSlackBus == False],
                                                p_max_new], axis=0)
                    p_max_new.index= range(p_max_new.shape[0])
                    Write2Excel(p_max_new, 'pGen_Profiles_Max_Addt')
                    FlexIDX    = np.argwhere(np.isin(Summary_Table.bus_i, AddWhere[Addt])).ravel()
                    Summary_Table.iloc[FlexIDX,1] =  Shares[Addt] * PV_cap # column 0 is bus_i


                    q_min_profile  = FlexExcel.parse('qGen_Profiles_Min_Addt');
                    NonSlackBus    = q_min_profile.bus_i!=SlackBus.bus_i.values.item()
                    ModelRow       = q_min_profile[NonSlackBus].iloc[[0],:]
                    q_min_new      = pd.concat([ModelRow]*DivideBy, axis=0)
                    q_min_new.bus_i= AddWhere[Addt]
                    q_min_new.iloc[:,1:] = 0
                    q_min_new      = pd.concat([q_min_profile[NonSlackBus == False],
                                                q_min_new], axis=0)
                    q_min_new.index= range(q_min_new.shape[0])
                    Write2Excel(q_min_new, 'qGen_Profiles_Min_Addt')

                    q_max_profile  = FlexExcel.parse('qGen_Profiles_Max_Addt');
                    NonSlackBus    = q_max_profile.bus_i!=SlackBus.bus_i.values.item()
                    ModelRow       = q_max_profile[NonSlackBus].iloc[[0],:]
                    q_max_new      = pd.concat([ModelRow]*DivideBy, axis=0)
                    q_max_new.bus_i= AddWhere[Addt]
                    q_max_new.iloc[:,1:] = np.tile(np.reshape(Shares[Addt] * PV_cap,(-1,1)), (1,24)) * 0
                    q_max_new      = pd.concat([q_max_profile[NonSlackBus == False],
                                                q_max_new], axis=0)
                    q_max_new.index= range(q_max_new.shape[0])
                    Write2Excel(q_max_new, 'qGen_Profiles_Max_Addt')

               # ------------------------------------------------------------------
               elif Addt in 'Storage':
                    ESS_cost       = FlexExcel.parse('Storage_cost_Addt')
                    ESS_cost_new   = pd.concat([ESS_cost.iloc[[0],:]]*DivideBy,axis=0)
                    ESS_cost_new.bus_i = AddWhere[Addt]
                    ESS_cost_new.index= range(ESS_cost_new.shape[0])
                    Write2Excel(ESS_cost_new, 'Storage_cost_Addt')

                    # FlexExcel.parse('Storage_Addt')
                    ESS            = FlexExcel.parse('Storage_Addt')
                    ModelRow       = ESS.iloc[[0],:]
                    ESS_new        = pd.concat([ModelRow]*DivideBy, axis=0)
                    ESS_new.bus_i  = AddWhere[Addt]
                    ESS_new.eRat   = ESS_cap * Shares[Addt]
                    ESS_new.chRat  = ESS_new.eRat * ESS_new.chRat
                    ESS_new.disRat = ESS_new.eRat * ESS_new.disRat

                    ESS_new.index= range(ESS_new.shape[0])
                    if   logical_flex == 'WOFlex':
                         col_index, col_name = StatusColumn(ESS_new)
                         # col_index      = np.argwhere('status' == np.char.lower(ESS_new.columns.values.astype(str))).item()
                         # col_name       = ESS_new[Addt].columns[col_index]
                         ESS_new.loc[:,col_name] = 0
                    else:
                         []
                         #### ESS_new.loc[:,'Status'] = 1

                    Write2Excel(ESS_new, 'Storage_Addt')
                    FlexIDX    = np.argwhere(np.isin(Summary_Table.bus_i, AddWhere[Addt])).ravel()
                    Summary_Table.iloc[FlexIDX,2] =  Shares[Addt] * ESS_cap # column 0 is bus_i

          # ------------------------------------------------------------------
          EditedSheets = ['Gens_cost_Addt','Gens_Addt','RES_Addt','Buses_Addt',
                          'Storage_cost_Addt','Storage_Addt',
                          'pLoad_Profiles_Addt','qLoad_Profiles_Addt','Loads_Addt',
                          'pGen_Profiles_Min_Addt','pGen_Profiles_Max_Addt',
                          'qGen_Profiles_Min_Addt','qGen_Profiles_Max_Addt']

          # ------------------------------------------------------------------
          if   True: #len(RegexMatch)>0: # let the lines be sorted anyway cuz if the lines are sorted in one file, but not sorted in the other file, i get errors
               Lines_Addt          = FlexExcel.parse('Lines_Addt')
               Lines_Addt          = DiscipLines(Lines_Addt, RegexMatch)
               Write2Excel(Lines_Addt, 'Lines_Addt')
               EditedSheets.append('Lines_Addt')

          Write2Excel(Summary_Table, 'Summary_Table')
          # Summary_Sheet1 = pd.DataFrame.from_dict(AddWhere,orient='index')
          # Summary_Sheet2 = pd.DataFrame.from_dict(Shares,orient='index')
          # Summary_Sheet  = pd.concat([Summary_Sheet1,Summary_Sheet2],axis=0,keys=['AddWhere','Shares'])

          # Load_EV_summary = pd.concat([pd.DataFrame(ev_profile_split, index=AddWhere['Loads'], columns=t1_t24),
          #                              pd.DataFrame(p_profile.iloc[FlexIDX,1:ev_profile.size+1].values, index=AddWhere['Loads'], columns=t1_t24),
          #                              pd.DataFrame(p_profile.iloc[FlexIDX,1:ev_profile.size+1].values * LoadGrowth, index=AddWhere['Loads'], columns=t1_t24),
          #                              ], axis=0, keys=['EV','P_old','P_old_x_LoadGrowth'])


          # Write2Excel(Summary_Sheet, 'Add_Where')

          # sys.exit()
          for  sheet in np.setdiff1d(FlexExcel.sheet_names, EditedSheets):
               if   np.char.startswith(sheet,'\'file:///'):
                    continue

               TheSheet = FlexExcel.parse(sheet)
               TheSheet.to_excel(ExcelWriter, sheet_name=sheet,index=False)

          # Data2Present[key] = SurveyAssets(FlexExcel)

          ExcelWriter.close()
          # os.startfile(DestinationFlexFile) # for debugging
          # ------------------------------------------------------------------
          # breakpoint()
          # sys.exit()
          []

Data2PresentDF = pd.DataFrame.from_dict(Data2Present, orient='index')
print(Data2PresentDF, end='\n'*3)
# print('Cases names and file paths',*CaseNamesFiles,sep='\n')
assert len(CaseNamesFiles) == np.unique(np.char.array(CaseNamesFiles)[:,0]).size, 'Different cases have the same name'
CaseNamesFilesNP         = np.char.array(CaseNamesFiles)
CaseNamesFilesNP         = '"' + CaseNamesFilesNP + '"'
CaseNamesFilesNP[:,1]    = np.char.replace(CaseNamesFilesNP[:,1],'D:\\Work\\Produce\\TESTIFY\\SLA\\','')
CaseNamesFilesNP         = np.char.replace(CaseNamesFilesNP,'\\','/')
CaseNames_strLen         = np.char.str_len(CaseNamesFilesNP).max(axis=0)
CaseNamesFilesNP         = np.char.ljust(CaseNamesFilesNP,CaseNames_strLen)
print(*[('(%s, %s),' %tuple(i[:-1])) for i in CaseNamesFilesNP],sep='\n')

# -------------------------------------------
SummaryPD      = pd.concat([pd.DataFrame(CaseNamesFilesNP,columns=['Case','File','Time']),
                            Data2PresentDF.reset_index(drop=False)],
                           axis=1)
SummaryLogFile = r"D:\Work\Produce\TESTIFY\SLA\input_data\March2023_LoadDataFromAttestCloud\EditedInputFiles\InLog.xlsx"
ExcelWriter    = pd.ExcelFile(SummaryLogFile, engine="openpyxl")
if   len(ExcelWriter.sheet_names) == 0:
     SummaryPD.to_excel(ExcelWriter, sheet_name = 'Sheet1', startrow = 0, index = False, header = True)
else:
     # breakpoint()
     ExistingRows = ExcelWriter.parse('Sheet1')
     pd.concat([SummaryPD,ExistingRows],axis=0).to_excel(ExcelWriter, startrow = 0, index = False, header = True)

ExcelWriter.close()
os.startfile(SummaryLogFile)
# -------------------------------------------

ToastNotifier().show_toast(title='ATTEST case generator', msg = 'Finished %i cases' %len(CaseNamesFiles), duration=15)
