get_ipython().magic('%reset -f')
get_ipython().magic('%cls')

import pandas as pd
import numpy as np
import os, sys, shutil, re, time
# import xlsxwriter as XLS
from difflib import SequenceMatcher
from win10toast import ToastNotifier

np.set_printoptions(suppress  = True,
                    precision = 2)
pd.set_option("display.max_rows"   , 1000,
              "display.max_columns", 50,
              "display.width"      , 2000,
              "display.precision"  , 3)
np.set_printoptions(threshold=np.inf, linewidth=np.inf)

PWD            = os.getcwd()
cases_path     = r"D:\Work\Produce\TESTIFY\SLA\input_data\March2023_LoadDataFromAttestCloud\EditedInputFiles"
EV_PV_ESS_File = r'EV-PV-Storage_Data_for_Simulations_updated_17_APRIL.xlsx' # Added on May 20th
Output_EV_PV_EES_File = r"D:\Work\Produce\TESTIFY\SLA\output_data\March2023_LoadDataFromAttestCloud\EditedInputFiles\Output_EV-PV-Storage_Data.xlsx"

cases =  {
# "Year":[2030, 2040, 2050],
"ES" : [1.020408163     , 1.204081633     , 1.306122449 ],
"PT" : [1.168849815     , 1.436098655     , 1.704035874 ],
"HR" : [1.175428234     , 1.305375074     , 1.683402245 ],
"UK" : [1.047405509     , 1.310057655     , 1.444586803 ],
}

LdGrwthPD = pd.DataFrame.from_dict(cases, orient = 'index', columns = [2030, 2040, 2050])

#%%
def MyWalk(path):
     TheList   = [*os.walk(path)]
     OutFiles  = []
     OutFolders=[]
     for  i, Row in enumerate(TheList): # i=0;  Row=TheList[i]
          # print(i)
          if len(Row[2])>0:
               Row2 = np.hstack([*map(np.vstack, np.meshgrid(Row[0], Row[2]))])
               Row3 = [os.path.join(*i) for i in Row2]
               OutFiles.extend(Row3)

          if  len(Row[1])>0:
               Row4 = np.hstack([*map(np.vstack, np.meshgrid(Row[0], Row[1]))])
               Row5 = [os.path.join(*i) for i in Row4]
               OutFolders.extend(Row5)

     return OutFiles, OutFolders

#%%
EV_PV_ESS_excel     = pd.ExcelFile(EV_PV_ESS_File)
EV_PV_ESS_sheets    = EV_PV_ESS_excel.sheet_names
permutations        = np.char.split(EV_PV_ESS_sheets,'_')

#%%
if False:
     FilesList = !dir /b /s
     # if you get any syntax errors from this line, it means you have a problem somewhere else,
     # but spyder hates this line and will point fingers at it right away.
     # to debug, disable this line and re-run the code to see where the real syntax error is
     ExcelWriter = pd.ExcelWriter(Output_EV_PV_EES_File, mode='w', engine = 'xlsxwriter')

elif os.path.exists(Output_EV_PV_EES_File):
     CopyFile       = Output_EV_PV_EES_File.replace('.xlsx',' - Copy.xlsx')
     assert os.path.isfile(CopyFile),'File does not exist'
     ContentsSheet  = pd.read_excel(CopyFile, sheet_name='Contents')
     FilesList      = ContentsSheet.File.values
     ExcelWriter    = pd.ExcelWriter(Output_EV_PV_EES_File, mode='w', engine = 'xlsxwriter')

else:
     raise Exception('Can\'t fetch files')

FlexFiles = [i for i in FilesList if i.endswith('_flex.ods')]
# Sheet2Read= ['pLoad_Profiles_Addt', 'RES_Addt',' Storage_Addt']
# print(*enumerate(FlexFiles),sep='\n'*2)
# FlexFiles = [FlexFiles[14]]
# FlexFiles = FlexFiles[5:]

AllSheets = []
Skipped   = []
Exported  = [];
breakpoint()

for  i, flex_file in enumerate(FlexFiles[::-1]): # i = 14; flex_file = FlexFiles[i]
     print('%s Doing File #%3i/%i @%s:\n"%s"' %('\n'*3 + '='*120 + '\n', i+1, len(FlexFiles), time.ctime(), flex_file))
     # os.startfile(flex_file)

     # EV_PV_ESS.xlsx file ----------------------------------------------------------------------------------------------
     matching_perm  = [(i,j) for i,j in enumerate(permutations) if all([x in flex_file for x in j])]
     if   len(matching_perm)==0:
          print('\nNo EV-PV-EES data sheet matching the name of the file #%i:\n"%s".\nSkipping...' %(i+1, flex_file))
          Skipped.append(i)
          continue
     else:
          sheet = '_'.join(matching_perm[0][1])
          N_sheet_duplicates = np.isin(AllSheets, sheet).sum()
          AllSheets.append(sheet)

          if   N_sheet_duplicates>0:
                sheet += ' (%i)' %(N_sheet_duplicates+1)
                print('Sheet name changed to: "%s"' %sheet)

          if sheet in ExcelWriter.sheets.keys():
               breakpoint()
               print('Duplicate names!')

          ExcelWriter.write_cells([], sheet_name=sheet) # to create the sheet
          if False:
               ExcelWriter.sheets[sheet].write_string(0,0,str(i))
               ExcelWriter.sheets[sheet].write_string(0,1,flex_file)
               #### ExcelWriter.close();  os.startfile(Output_EV_PV_EES_File) ## for debugging

          ev_pv_ess_sheet = EV_PV_ESS_sheets[matching_perm[0][0]]
          print('\nMatching with EV-PV-EES data of "%s"' %ev_pv_ess_sheet.upper())

     if True:
          pv_ess_data    = EV_PV_ESS_excel.parse(sheet_name = ev_pv_ess_sheet,
                                                 header     = None,
                                                 usecols    = 'C:D',
                                                 nrows      = 1)

          node_data      = EV_PV_ESS_excel.parse(sheet_name = ev_pv_ess_sheet,
                                                 header     = 1,
                                                 usecols    = 'A:E',
                                                 )

          ev_data        = EV_PV_ESS_excel.parse(sheet_name = ev_pv_ess_sheet,
                                                 header     = 1,
                                                 usecols    = 'F:I',
                                                 index_col  = 0).dropna()

     else:
          pv_ess_data    = pd.read_excel(EV_PV_ESS_excel,
                                          sheet_name   = ev_pv_ess_sheet,
                                          header       = None,
                                          usecols      = 'C:D',
                                          nrows        = 1)

          node_data      = pd.read_excel(EV_PV_ESS_excel,
                                          sheet_name   = ev_pv_ess_sheet,
                                          header       = 1,
                                          usecols      = 'A:E',
                                          )

          ev_data        = pd.read_excel(EV_PV_ESS_excel,
                                          sheet_name   = ev_pv_ess_sheet,
                                          header       = 1,
                                          usecols      = 'F:I',
                                          index_col    = 0).dropna()

     # flex.ods file ----------------------------------------------------------------------------------------------
     # os.startfile(flex_file)
     if   True:
          # os.startfile(flex_file)
          try:
               flex_excel     = pd.ExcelFile(flex_file, engine="odf")
          except Exception as Opening_Error:
               print("Failed to open file. skipping #%i/%i:\n\"%s\"" %(i+1, len(FlexFiles)+1, flex_file))
               Skipped.append(i)
               ExcelWriter.sheets[sheet].write_string(0,2,'Skipped because failed to open _flex file')
               continue

          flex_sheets    = flex_excel.sheet_names
          flex_sheets    = [i for i in flex_sheets if i.startswith("\'file:///")==False]
          print('Finished loading XLSX file')

          bus_sheet      = flex_excel.parse('Buses_Addt')
          buses          = np.sort(bus_sheet.bus_i.values)
          res_sheet      = flex_excel.parse('RES_Addt')
          storage_sheet  = flex_excel.parse('Storage_Addt')

          pload_sheet    = flex_excel.parse('pLoad_Profiles_Addt',
                                            header=0,
                                            index_col=0)
          loads_sheet    = flex_excel.parse('Loads_Addt',
                                            header=0,
                                            index_col=0)
          loads_sheet    = loads_sheet[loads_sheet.Status==1]
          i_ev           = loads_sheet.index.tolist()
          N_ev           = len(i_ev); #loads_sheet.Status.sum()

          if N_ev==0:
               print("case has no EVs. Skipping")
               Skipped.append(i)
               ExcelWriter.sheets[sheet].write_string(0,2,'Skipped because N_ev==0, WO_flex')
               continue # This is a case WITHOUTFlex. i can't locate where EVs are because i had set Status=0 on load flexibility

          Full_load_buses_with_EV = pload_sheet.loc[i_ev,:]

          ev_split       = ev_data['EV load (MW)'] / N_ev
          if True:
               ev_data.to_excel(ExcelWriter,
                                sheet_name  = sheet,
                                index       = True,
                                startrow    = 1,
                                startcol    = 5,
                                # float_format="%0.7g",
                                )
          else:
               EV_profile     = pd.concat([ev_split]*N_ev, axis=1)
               EV_profile.columns=i_ev
               EV_profile.index = Full_load_buses_with_EV.columns

               Load_WO_EV = Full_load_buses_with_EV - ev_split.values

               vertical_load = pd.concat([EV_profile, Load_WO_EV.T], axis=1, keys= ['EV_profile','P_load without EV'])
               vertical_load.columns.names = ['Load Type','bus_i']
               vertical_load.to_excel(ExcelWriter,
                                      sheet_name  = sheet,
                                      index       = True,
                                      startrow    = 0,
                                      startcol    = 5,
                                      # float_format="%0.7g",
                                      )

     # else:
     #      res_sheet      = pd.read_excel(flex_file, sheet_name='RES_Addt')
     #      storage_sheet  = pd.read_excel(flex_file, sheet_name='Storage_Addt')
     #      pload_sheet    = pd.read_excel(flex_file, sheet_name='pLoad_Profiles_Addt')

     # compare = pd.concat([old_pload_sheet.loc[loads_sheet.index].T, Load_WO_EV.T], axis=1, keys = ['old','calculated'])
     # calc_grwth = compare.calculated / compare.old

     c_res          = res_sheet.Pmax.sum()
     i_res          = res_sheet.bus_i.values.tolist()

     c_ees          = storage_sheet.eRat.sum()
     i_ees          = storage_sheet.bus_i.values.tolist()

     # iloc_ev        = [loads_sheet.index.get_loc(i) for i in i_ev]

     assert abs(c_res - pv_ess_data.iloc[0,0]) < 1e-3, 'RES sizes don\'t match'
     assert abs(c_ees - pv_ess_data.iloc[0,1]) < 1e-3, 'EES sizes don\'t match'

     country, _, _, year = matching_perm[0][1]
     load_growth_factor  = LdGrwthPD.loc[country,int(year)]
     # calc_grwth/load_growth_factor
     SummaryPD = pd.DataFrame(0,columns = ['Node Ratio','PV (MW)','EES (MWh)','EV (number)'], index = buses)
     SummaryPD.loc[i_ev,'Node Ratio']   = 1/N_ev
     SummaryPD.loc[i_ev,'EV (number)']  = 1
     SummaryPD.loc[i_ees,'EES (MWh)']   = c_ees/len(i_ees)
     SummaryPD.loc[i_res,'PV (MW)']     = c_res/len(i_res)
     TotalPD                            = SummaryPD.iloc[:,1:].sum().to_frame(name='In the network').T

     SummaryPD.to_excel(ExcelWriter,
                        sheet_name  = sheet,
                        index       = True,
                        startrow    = 1,
                        startcol    = 0,
                        # float_format="%0.7g",
                        )

     TotalPD.to_excel(ExcelWriter,
                        sheet_name  = sheet,
                        index       = True,
                        header      = False,
                        startrow    = 0,
                        startcol    = 1,
                        # float_format="%0.7g",
                        )

     Exported.append([i, sheet, flex_file])

     # --------------------------------------------------------------------------------------------------------------
     # if   isinstance(old_pload_sheet, None.__class__):
     #      print('Sheet "old_pLoad_Profiles_Addt" could not be found in file #%i\n"%s".\nSkipping...' %(i, flex_file))
     #      continue
     # if   np.isin('old_pLoad_Profiles_Addt', flex_sheets).all():
     #      old_pload_sheet= flex_excel.parse('old_pLoad_Profiles_Addt',
     #                                        header=0,
     #                                        index_col=0)
     #      # pd.concat([Load_WO_EV, old_pload_sheet.loc[i_ev]], axis=0, keys=['Calculated','recorded']).T

     # else:
     #       old_pload_sheet = None


     print('Finished writing sheet "%s" to XLSX file' %sheet)
     print('Total exported sheets = %i' %len(Exported))
     print('Total skipped  sheets = %i' %len(Skipped))
     del sheet
     # ExcelWriter.save() # This closes the file. you dont want to close it yet. next iterations will not write to the file.
     # breakpoint()
     # ExcelWriter.save()

ToCPD = pd.DataFrame(Exported, columns = ['#','Case','File'])
ToCPD.to_excel(ExcelWriter,
                   sheet_name  = 'Contents',
                   index       = False,
                   startrow    = 0,
                   startcol    = 0,
                   # float_format="%0.7g",
                   )

ExcelWriter.close()
os.startfile(Output_EV_PV_EES_File)

# with open(Output_EV_PV_EES_File.replace('.xlsx','.txt'), 'w') as f:
#     f.write('\n'.join(['%i: %s' %(i,j) for i,j in enumerate(FlexFiles[::-1])]))

sys.exit('Finished Successfully without errors')
#%%
'''
ii=iloc_ev[0];
# ii=2;
print('row#%i, bus#%i' %(ii, loads_sheet.index[ii]))
test = pd.concat([pload_sheet.iloc[ii,:],
                  old_pload_sheet.iloc[ii,:] * load_growth_factor,
                  pload_sheet.iloc[ii,:] - (old_pload_sheet.iloc[ii,:] * load_growth_factor),
                  old_pload_sheet.iloc[ii,:],
                  pload_sheet.iloc[ii,:] / old_pload_sheet.iloc[ii,:],
                  ],
                 axis=1,keys=['pload','old_pload_x_grwth','Err','old_pload','pload/old'])
print(test)
'''

# TheSheet.to_excel(ExcelWriter, sheet_name=sheet,index=False)

#%%

Output_File    = pd.ExcelFile(Output_EV_PV_EES_File)
sheet_names    = Output_File.sheet_names
ToC = []

for  sheet in sheet_names: # i = 0; sheet = sheet_names[i]
     the_data = Output_File.parse(sheet_name = sheet,
                                  header     = None,
                                  usecols    = 'A:B',
                                  nrows      = 1).values.ravel().tolist()
     the_data.insert(1,sheet)
     ToC.append(the_data)

ToCPD = pd.DataFrame(ToC, columns = ['#','Case','File'])
ToCPD.to_excel(Output_EV_PV_EES_File.replace('.xlsx','table of contents.xlsx'), sheet_name = 'Contents')
