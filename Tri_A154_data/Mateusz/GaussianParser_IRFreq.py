
# coding: utf-8

# In[1]:


import os
from os.path import join
import re
import pandas as pd


# In[2]:


import os

os.getcwd()


# In[3]:


#creating columns for a dataframe
l = ['Molecule', 'Freq. A', 'Freq. B','Freq. C', 'IR Int. A', 'IR Int. B', 'IR Int. C', 'Atom', 'AN', 'X1','Y1','Z1', 
    'X2','Y2','Z2', 'X3','Y3','Z3']
#creating an empty df
FIR_df = pd.DataFrame(columns=l)

for (dirname, dirs, files) in os.walk('.'):      #'.' indicates to start in the current directory and walk downwards
    
    for filename in files:                       #looping through files in a directory
        #getting the path and folder name
        dirpath = os.getcwd()
        #print("current directory is : " + dirpath)
        foldername = os.path.basename(dirpath)
        #print("Directory name is : " + foldername)
        
        if filename.endswith('.log'):            #choosing only log files
            #reading file     
            file = open(os.path.join(dirname, filename), "r")
            #print(3)
            print('')                            #line to use as a process probe
            print(str(dirname))                  #line to use as a process probe

            #a flag for boolean decision; will help to identify sections
            section_flag = False
        
            for line in file:                              #parsing lines in a log file
                #print(4)
                if re.search('^ Frequencies', line):       #processing Frequencies line 
                    #print('Freq located')
                    freq_line = line.strip()                 #strip line
                    freq_line = freq_line.split()            #splitting the line
                    freq_line.insert(2, str(dirname[2:]))
                    #print(freq_line)
                
                elif re.search('^ IR Inten', line):        #processing IR line
                    #print('IR located')
                    ir_line = line.strip()                   #strip line
                    ir_line = ir_line.split()                #splitting the line
                    #print(ir_line)
                
                elif re.search('^  Atom', line):           #locating Atom lines
                    #print('Atom located')
                    section_flag = True                    #changing the flag, ie, this are valid lines
                    
                elif re.search('^\s*.\d\s\s\s\d', line) and section_flag == True:    #selecting lines with vibration of modes
                    #print('Coords located')
                    #print(line)
                    ser = pd.Series(freq_line[2:]+ir_line[3:]+line.split(), index=l)  #create a pd.Series object...
                    FIR_df = FIR_df.append(ser, ignore_index=True)   #adding the object as row to the dataframe
                
                else:
                    section_flag = False                   #failing in any of the other lines, implies invalid line
                    if re.search('^ - Thermochemistry', line):   #this line tells us the corrdinates are over, no need to parse more
                        break
            file.close()                                   #closing the file
            print('File Done!')
            print('FIR_df contains :'+str(len(FIR_df))+' rows')
        else:
            pass


# In[4]:


len(FIR_df)


# In[5]:


FIR_df


# In[6]:


#IR_df = FIR_df.iloc[:130871]


# In[7]:


len(FIR_df)


# In[8]:


FIR_df.to_csv('DF1.csv')


# In[9]:


#FIR_df.head(5)


# In[10]:


#FIR_df.Tri_A154_ga09

