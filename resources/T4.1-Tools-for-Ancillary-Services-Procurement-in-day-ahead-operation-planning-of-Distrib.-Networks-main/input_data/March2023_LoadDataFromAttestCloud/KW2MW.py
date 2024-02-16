get_ipython().magic('%reset -f')
get_ipython().magic('%cls')

import pandas as pd
import numpy as np
import os, sys, shutil, re

# FolderPath     = r"D:\Work\Produce\TESTIFY\SLA\input_data\March2023_LoadDataFromAttestCloud\PT DX 01\Original in KW"
FolderPath     = r"D:\Work\Produce\TESTIFY\SLA\input_data\March2023_LoadDataFromAttestCloud\ES DX 03\Original in KW"
DataStartsAt   = 2
ExpectedLength = 96

#%%
OutputPath     = os.path.dirname(FolderPath)

FilesNames     = np.array(os.listdir(FolderPath))
FilesNames     = FilesNames[np.char.endswith(FilesNames,'.xlsx')].tolist()

# print(pd.read_csv.__doc__)
# print(pd.DataFrame.to_excel.__doc__)
N_files = len(FilesNames)
for  i, file in enumerate(FilesNames):
     # breakpoint()
     print('Doing %3i/%-3i: %s' %(i, N_files, file))
     FileFullPath = os.path.join(FolderPath,file)
     FileContents = pd.read_excel(FileFullPath, index_col=[0,1], header=None).squeeze("columns")*1e-3
     FileContents = FileContents.sort_index(axis=0, level=0)

     NewFileFullPath = FileFullPath.replace(FolderPath, OutputPath)
     # print(pd.DataFrame.to_excel.__doc__)
     # breakpoint()
     FileContents.to_excel(NewFileFullPath, header=False)
     # os.startfile(NewFileFullPath)
     del FileContents

     # breakpoint()

# load#, bus#
sys.exit('Finished')
