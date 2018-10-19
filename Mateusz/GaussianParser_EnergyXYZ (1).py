
# coding: utf-8

# In[1]:


import os
from os.path import join
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal 
from matplotlib.ticker import NullFormatter
import sys
from PIL import Image
import os


# In[2]:


os.getcwd()


# In[3]:


"%cd .."


# # I. Create Dataframe for Coordinates  & Energy 
# 
# ## Attributes 
#     1) Center 
#     2) Atomic Number
#     3) Atomic Type
#     4-6) 'X','Y', 'Z' (Coordinates)
#     7) E(RPBE1PBE) ---Energy
#     8) Molecule 

# In[4]:


#creating columns for a dataframe
l = ['Center Number', 'Atomic Number', 'Atomic Type', 'X', 'Y', 'Z', 'E(RPBE1PBE)', 'Molecule']
#creating an empty df
coords_df = pd.DataFrame(columns = l)

for (dirname, dirs, files) in os.walk('.'):
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
            #section_flag = False
            
            tables_num = 0
            energies_num = 0
            for line in file:                              #parsing lines in a log file; 1st time
                #print(4)
                #if re.search('^ Frequencies', line):       #processing Frequencies line 
                if re.search('^ *Standard orientation:', line):  ##counting std orient. tables
                    tables_num += 1
                elif re.search('^ SCF Done', line):         #counting SCF lines
                    energies_num += 1
            #print('UUU')
            con3 = 0
            con4 = 0
            #print('tables_num', tables_num)
            
            file.seek(0)        #rewinding the file
            
            for line in file:
               # print(5)
                if re.search('^ *Standard orientation:', line):
                    con3 += 1
                    #print('con3', con3)
                if con3 == tables_num and re.search('^\s\s\s\s\s.\d\s*\d', line):
                    ##print(6)
                    ser = pd.Series(line.split()+['',''], index=l)                  #create a pd.Series object...
                    coords_df = coords_df.append(ser, ignore_index=True)          #
                elif con3 == tables_num and re.search('^ SCF Done', line):
                    #print(7)
                    line = line.split()        
                    #add to the df
                    coords_df['E(RPBE1PBE)'][len(coords_df)-84:] = pd.Series(((line[4]+' ')*84).split())#, index = coords_df.index)
                    coords_df['Molecule'][len(coords_df)-84:] = pd.Series(((str(dirname[2:])+' ')*84).split()) 
                    break    

            file.close()                                   #closing the file
            print('File Done!')
            print('coords_df contains: '+str(len(coords_df))+' rows')
            #print('FIR_df contains: '+str(energies_num)+' SCF line')
        else:
            pass


# In[5]:


#shape oof coords_df
coords_df.head(5)


# # II.  Clean Dataframe for coordinates-Energy & upload Dataframe w/ IR + Frequencies for A, B, C

# In[6]:


#make copy of original coors_df to manipulate in case later we need the original df 
#examine shape of copied df
df_1= coords_df.copy()
print(df_1.shape)
print(type(df_1))
df_1.head(5)


# In[7]:


#load df1.csv which contains the Intensities + Frequencies of Tri_data00
df_3 = pd.read_csv("DF1.csv")


# In[8]:


#isolate just just teh attributes we need: Intensties + frequencies.
#make a copy so we can manipulate the data and preserve the orginal uploaded df
#examine shape of df_3
df_3= df_3.iloc[:,1:8]
df_4= df_3.copy()
print(df_3.shape)
df_3.head(5)


# In[9]:


#set index and drop the duplictae values in the df and reexamine the shape 
df_4 = df_4.set_index("Molecule", drop = True)
df_4.head(5)
df_4.drop_duplicates(subset=['Freq. A','Freq. B','Freq. C','IR Int. A','IR Int. B','IR Int. C'], keep='first', inplace=True)


# In[10]:


#look at shape again 
print (df_4.shape)
df_4.head(5)


# In[11]:


#set the same index (molecule) on df_1 as df_4
df_1 = df_1.set_index("Molecule", drop = True)
print(df_1.shape)
df_1.head(5)


# In[12]:


#set index again after editiing
df_3 = df_3.set_index("Molecule", drop = True)
print(df_3.shape)
df_3.head(5)


# # III. Join Dataframe of Coordinates + Dataframe of IR-Frequencies; called df_combo

# In[13]:


#we join df_3 and df_1 to make a single useful df.  
#we use df_3 instead of df_4 because df_3 has the repeated values (same shape) as the df that is being joined with, df_1
df_combo= df_3.join(df_1)
print(df_combo.shape)
df_combo.head(5)


# In[14]:


df_combo1= df_combo.drop(df_combo.columns[[6, 8]], axis=1)
df_combo1.head(5)


# # IV. Isolate IR- Freq (A, B, C) of Tri_A154_0000 for Gaussian 
# 
#     1) Isolate IR & Freq A
#       a. index
#       b. lists 
#       c. Gaussian function
#     
#     2) Isolate IR & Freq B
#       a. index
#       b. lists 
#       c. Gaussian function

# In[15]:


#isolate Tri_A154_0000
df_A154_00= df_4.loc["Tri_A154_0000", : ]
print(df_A154_00.shape)
df_A154_00.head(5)


# In[16]:


#Isolate respective IR & Franquencies A, B, and C
dfA24IR_A= df_A154_00.loc["Tri_A154_0000", ['IR Int. A']]
print('1) dfA24IRF_A: ' + str((dfA24IR_A.shape)))
dfA24F_A= df_A154_00.loc["Tri_A154_0000", ['Freq. A']]
print('3) dfA24F_A: ' + str((dfA24F_A.shape)))

dfA24IR_B= df_A154_00.loc["Tri_A154_0000", ['IR Int. B']]
print('4) dfA24IRF_B: ' + str((dfA24IR_B.shape)))
dfA24F_B= df_A154_00.loc["Tri_A154_0000", ['Freq. B']]
print('6) dfA24F_B: ' + str((dfA24F_B.shape)))

dfA24IR_C= df_A154_00.loc["Tri_A154_0000", ['IR Int. C']]
print('7) dfA24IRF_C: ' + str((dfA24IR_C.shape)))
dfA24F_C= df_A154_00.loc["Tri_A154_0000", ['Freq. C']]
print('9) dfA24F_C: ' + str((dfA24F_C.shape)))


# ### 1) Int & Freq A

# In[17]:


#use pandas to create respective index IR & Fre A
dfA24IR_A.index = pd.RangeIndex(len(dfA24IR_A.index))
dfA24F_A.index = pd.RangeIndex(len(dfA24F_A.index))
dfA24IR_A.head(5)
print(dfA24IR_A.dtypes)
dfA24F_A.head(5)


# In[18]:


#use pandas to change the respetive df's from float to INT values for intensity A
dfA24IR_Ai= pd.to_numeric(dfA24IR_A['IR Int. A'])
dfA24IR_Ai= dfA24IR_Ai.values.T.tolist()
dfA24IR_Ai


# In[19]:


#use pandas to change the respetive df's from float to INT values for frequency A
dfA24F_Ai= pd.to_numeric(dfA24F_A['Freq. A'])
dfA24F_Ai= dfA24F_Ai.values.T.tolist()
dfA24F_Ai


# In[20]:


#create a list using fucntion zip.  Because IR and frequeny have the same length, each entry has a tuple 
dfA24IRF_A = np.zeros((len(dfA24IR_Ai), 2) )
for n, ir, i  in zip(range(len(dfA24IR_Ai)), dfA24IR_Ai, dfA24F_Ai): 
     dfA24IRF_A[n,0] = ir; dfA24IRF_A[n,1]=i


# In[21]:


#use gaussian function to obatin our array
broaden=2
def gaussian(X,x0,s):
    return np.exp(-0.5*((X-x0)/s)**2)

def lorentz(X,x0,s):
    return 1/(np.pi*s*(1+((X-x0)/s)**2))

IR_A000A = np.zeros((4001,))
X=np.linspace(0, 4000,4001)
Y=np.zeros((4001,))
for l in range(dfA24IRF_A.shape[0]):
    Y = dfA24IRF_A[l,1]*gaussian(X, dfA24IRF_A[l,0], broaden)
    IR_A000A += Y


# In[22]:


IR_A000A


# ### 1) Int & Freq B

# In[23]:


#use pandas to create respective index IR & Freq B
dfA24IR_B.index = pd.RangeIndex(len(dfA24IR_B.index))
dfA24F_B.index = pd.RangeIndex(len(dfA24F_B.index))
dfA24IR_B.head(5)
print(dfA24IR_B.dtypes)
dfA24F_B.head(5)


# In[24]:


#use pandas to change the respetive df's from float to INT values for intensity B
dfA24IR_Bi= pd.to_numeric(dfA24IR_B['IR Int. B'])
dfA24IR_Bi= dfA24IR_Bi.values.T.tolist()


# In[25]:


#use pandas to change the respetive df's from float to INT values for frequency B
dfA24F_Bi= pd.to_numeric(dfA24F_B['Freq. B'])
dfA24F_Bi= dfA24F_Bi.values.T.tolist()


# In[26]:


#create a list for Int and Freq B
dfA24IRF_B = np.zeros((len(dfA24IR_Bi), 2) )
for n, ir, i  in zip(range(len(dfA24IR_Bi)), dfA24IR_Bi, dfA24F_Bi): 
     dfA24IRF_B[n,0] = ir; dfA24IRF_B[n,1]=i


# In[27]:


#use gaussian function on IR & freq B
def gaussian(X,x0,s):
    return np.exp(-0.5*((X-x0)/s)**2)

def lorentz(X,x0,s):
    return 1/(np.pi*s*(1+((X-x0)/s)**2))

IR_A000B = np.zeros((4001,))
X=np.linspace(0, 4000,4001)
Y=np.zeros((4001,))
for l in range(dfA24IRF_B.shape[0]):
    Y = dfA24IRF_B[l,1]*gaussian(X, dfA24IRF_B[l,0], broaden)
    IR_A000B += Y


# In[28]:


IR_A000B


# ### 1) Int & Freq C

# In[29]:


#use pandas to create respective index IR & Fre C
dfA24IR_C.index = pd.RangeIndex(len(dfA24IR_C.index))
dfA24F_C.index = pd.RangeIndex(len(dfA24F_C.index))
dfA24IR_C.head(5)
print(dfA24IR_C.dtypes)
dfA24F_C.head(5)


# In[30]:


#use pandas to change the respetive df's from float to INT values for intensity C
dfA24IR_Ci= pd.to_numeric(dfA24IR_C['IR Int. C'])
dfA24IR_Ci= dfA24IR_Ci.values.T.tolist()


# In[31]:


#use pandas to change the respetive df's from float to INT values for frequency B
dfA24F_Ci= pd.to_numeric(dfA24F_C['Freq. C'])
dfA24F_Ci= dfA24F_Ci.values.T.tolist()


# In[32]:


#create a list for Int and Freq C
dfA24IRF_C = np.zeros((len(dfA24IR_Ci), 2) )
for n, ir, i  in zip(range(len(dfA24IR_Ci)), dfA24IR_Ci, dfA24F_Ci): 
     dfA24IRF_C[n,0] = ir; dfA24IRF_C[n,1]=i


# In[33]:


#use gaussian function on IR & freq B
def gaussian(X,x0,s):
    return np.exp(-0.5*((X-x0)/s)**2)

def lorentz(X,x0,s):
    return 1/(np.pi*s*(1+((X-x0)/s)**2))

IR_A000C = np.zeros((4001,))
X=np.linspace(0, 4000,4001)
Y=np.zeros((4001,))
for l in range(dfA24IRF_C.shape[0]):
    Y = dfA24IRF_C[l,1]*gaussian(X, dfA24IRF_C[l,0], broaden)
    IR_A000C += Y


# In[34]:


IR_A000C


# # V. Gaussian Broadning

# In[35]:


from scipy import interpolate
import sys
from optparse import OptionParser


# In[36]:


#attempt to load gaussain using our paramaters (still working progress)
def error(msg):
   """ write error message and quit
   """
   sys.stderr.write(msg + "\n")
   sys.exit(3)


# In[37]:


def deriv(spec,h):
   """ calculate first derivative of function 'spec'
       using the central finite difference method up to 6th order,
       and for the first 3 and last 3 grid points the
       forward/backward finite difference method up to 2nd order.
       ...as used in f77-program and suggested by Zanazzi-Jona...
   """ 
   der_spec =[[i[0],0] for i in spec]

   length=len(spec)
   for i in range(3,length-3):
      der_spec[i][1]=(-1*spec[i-3][1]+9*spec[i-2][1]-45*spec[i-1][1]+45*spec[i+1][1]-9*spec[i+2][1]+1*spec[i+3][1])/(60*h)
   for i in range(0,3):
      der_spec[i][1]=(-11*spec[i][1]+18*spec[i+1][1]-9*spec[i+2][1]+2*spec[i+3][1])/(6*h)
   for i in range(length-3,length):
      der_spec[i][1]=(11*spec[i][1]-18*spec[i-1][1]+9*spec[i-2][1]-2*spec[i-3][1])/(6*h)

   return der_spec


# In[38]:


def get_range(IR_A000A,IR_A000B,w_incr,shift,start,stop):
   """ determine wavenumber range within the comparison between theoretical
       and experimental spectrum is performed (depends on the shift)
   """
   de1=start+shift-IR_A000B[0][0]
   if (de1 >= 0 ):
      de1=int((start+shift-IR_A000B[0][0])/w_incr+0.00001)
      enstart=de1
      tnstart=int((start-IR_A000A[0][0])/w_incr+0.00001)
   else:
      de1=int((start+shift-IR_A000B[0][0])/w_incr-0.00001)
      enstart=0
      tnstart=int((start-IR_A000A[0][0])/w_incr-de1+0.00001)
   de2=stop+shift-IR_A000B[-1][0]
   if (de2 <= 0 ):
      de2=int((stop+shift-IR_A000B[-1][0])/w_incr-0.00001)
      enstop=len(IR_A000B)+de2
      tnstop=len(IR_A000A)+int((stop-IR_A000A[-1][0])/w_incr-0.00001) 
   else:
      de2=int((stop+shift-IR_A000B[-1][0])/w_incr+0.00001)
      enstop=len(IR_A000B)
      tnstop=len(IR_A000A)+int((stop-IR_A000A[-1][0])/w_incr-de2-0.00001)
   return tnstart, tnstop, enstart, enstop
 


# In[39]:


def integrate(integrand,delta):
   """ integrate using the trapezoid method as Zanazzi-Jona suggested and was used in the f77-program...
   """
   integral = 0.5*(integrand[0][1]+integrand[-1][1])   
   for i in range(1,len(integrand)-1):
      integral += integrand[i][1]
   return integral*delta


# In[40]:


def ypendry(spec,d1_spec,VI):
   """ calculate the Pendry Y-function: y=l^-1/(l^-2+VI^2) with l=I'/I (logarithmic derivative),
       J.B. Pendry, J. Phys. C: Solid St. Phys. 13 (1980) 937-44
   """
   y=[[i[0],0] for i in spec]

   for i in range(len(spec)):
      if (abs(spec[i][1]) <= 1.E-7):
         if (abs(d1_spec[i][1]) <= 1.E-7):
            y[i][1] = 0 
         else:
            y[i][1] = (spec[i][1]/d1_spec[i][1])/((spec[i][1]/d1_spec[i][1])**2+VI**2)
      else:
         y[i][1] = (d1_spec[i][1]/spec[i][1])/(1+(d1_spec[i][1]/spec[i][1])**2*(VI**2))
   return y


# In[41]:



def _main():

   usage = """ %prog [options] r-fac.in
        Reads two spectra and calculates various R-factors -- FS 2011
        Attention: both spectra have to be given on the same, equidistant grid!
        NOTE: in the f77-program R1 is scaled by 0.75 and R2 is scaled by 0.5; this is not done here
        Please provide a file r-fac.in with the following specifications (without the comment lines!!!) 
        (the numbers are just examples, choose them according to your particular case)
        start=1000       # where to start the comparison
        stop=1800        # where to stop the comparison
        w_incr=0.5       # grid interval of the spectra -- should be 1 or smaller! (otherwise integrations/derivatives are not accurate)
        shift_min=-10    # minimal shift of the theoretical spectrum 
        shift_max=+10    # maximal shift of the experimental spectrum
        shift_incr=1     # shift interval
        r=pendry         # which r-factor should be calculated? options: pendry, ZJ, R1, R2 (give a list of the requested r-factors separated by comma)
        VI=10            # approximate half-width of the peaks (needed for pendry r-factor)
        """
   parser = OptionParser(usage=usage)
   parser.add_option("-t","--theory", 
                     default="theory.dat",dest="theory",metavar="filename",
                     help="input file theoretical spectrum [default: %default]")
   parser.add_option("-e","--exp", 
                     default="exp.dat",dest="exp",metavar="filename",
                     help="input file experimental spectrum [default: %default]")

   options, args = parser.parse_args()
   
# read spectra in
   input_theory=open(options.theory)
   tspec=[map(float,line.split()) for line in input_theory.readlines()]
   input_theory.close()
   input_exp=open(options.exp)
   IR_A000B=[map(float,line.split()) for line in input_exp.readlines()]
   input_exp.close()

# read r-fac.in
   if (len(args) != 1):
      error("Please provide r-fac.in and execute: python calc-r-factors.py r-fac.in")
   for line in open(args[0]):
      line=line.split("=")
      if line[0]=="start":
         start=float(line[1])
      elif line[0]=="stop":
         stop=float(line[1])
      elif line[0]=="w_incr":
         w_incr=float(line[1])
      elif line[0]=="shift_min":
         shift_min=float(line[1])
      elif line[0]=="shift_max":
         shift_max=float(line[1])
      elif line[0]=="shift_incr":
         shift_incr=float(line[1])
      elif line[0]=="r":
         r=line[1].split(",") 
         r=[i.strip() for i in r]  # strip off leading and trailing whitespaces
         sys.stdout.write("Requested r-factors: ")
         for i in r:
            sys.stdout.write("%s " % i)
         sys.stdout.write("\nNOTE: in the f77-program R1 is scaled by 0.75 and R2 is scaled by 0.5; this is not done here\n\n")   
      elif line[0]=="VI":
         VI=float(line[1])

# perform some checks of the input data...
   if (int(shift_incr/w_incr+0.00001) == 0):
      error("Error: shift_incr cannot be smaller than w_incr!")
   if (start-IR_A000B[0][0] < 0) or (IR_A000A[-1][0]-stop < 0):
      error("check experimental spectrum!!")
   if (start-IR_A000A[0][0] < 0) or (IR_A000A[-1][0]-stop < 0):
      error("check theoretical spectrum!!")
   if (int((IR_A000B[-1][0]-IR_A000B[0][0])/w_incr+0.0001) != len(IR_A000B)-1 ) or (int((IR_A000A[-1][0]-IR_A000A[0][0])/w_incr+0.0001) != len(IR_A000A)-1 ):
      error("check w_incr!!")

 
# cut out data points that are not needed in order to save time...
   if (IR_A000B[0][0]-(start+shift_min-w_incr*25) < 0):
         IR_A000B=IR_A000A[-1*int((IR_A000A[0][0]-(start+shift_min-w_incr*25))/w_incr-0.00001):][:]
   if (IR_A000B[-1][0]-(stop+shift_max+w_incr*25) > 0):
         IR_A000B=IR_A000B[:-1*(int((IR_A000B[-1][0]-(stop+shift_max+w_incr*25))/w_incr+0.00001)+1)][:] 
   if (IR_A000A[0][0]-(start-w_incr*25) < 0):
         IR_A000A=IR_A000A[-1*int((IR_A000A[0][0]-(start-w_incr*25))/w_incr-0.00001):][:]
   if (IR_A000A[-1][0]-(stop+w_incr*25) > 0):
         IR_A000A=IR_A000A[:-1*(int((IR_A000A[-1][0]-(stop+w_incr*25))/w_incr+0.00001)+1)][:]

   
# set negative intensity values to zero
   for i in range(0,len(IR_A000B)):
      if (IR_A000B[i][1]<0):
         IR_A000B[i][1]=0
   for i in range(0,len(IR_A000A)):
      if (IR_A000A[i][1]<0):
         IR_A000A[i][1]=0
   
# start calculating derivatives...
   d1_espec = deriv(IR_A000B,w_incr)   
   d1_tspec = deriv(IR_A000A,w_incr)
# calculate the second derivatives if the Zanazzi-Jona R-factor is requested   
   if "ZJ" in r:
      d2_tspec = deriv(d1_tspec,w_incr)
      d2_espec = deriv(d1_espec,w_incr)
# calculate Pendry Y-function if Pendry R-factor is requested      
   if "pendry" in r:
      ye = ypendry(espec,d1_espec,VI)
      yt = ypendry(tspec,d1_tspec,VI)
   


   min_pendry = [1.E100,0]
   min_r1     = [1.E100,0]
   min_r2     = [1.E100,0]
   min_zj     = [1.E100,0]
# start with loop over x-axis shifts
   for shift in numpy.arange(shift_min,shift_max+shift_incr,shift_incr):
      # get the interval within the two spectra are compared
      tnstart,tnstop,enstart,enstop = get_range(IR_A000A,IR_A000B,w_incr,shift,start,stop) 
      sys.stdout.write("\nshift: %9.3f, theory-start: %5d, theory-end: %5d, exp-start: %5d, exp-end: %5d\n" % (shift,tspec[tnstart][0],tspec[tnstop-1][0],espec[enstart][0],espec[enstop-1][0]))
      s_espec = numpy.array(IR_A000B[enstart:enstop]) # cut out the interval within which the comparison takes place
      s_tspec = numpy.array(IR_A000A[tnstart:tnstop])
      s_d1_espec = numpy.array(d1_IR_A000B[enstart:enstop])
      s_d1_tspec = numpy.array(d1_IR_A000A[tnstart:tnstop])
      c_scale=integrate(s_espec,w_incr)/integrate(s_tspec,w_incr)
      if "pendry" in r:
         # see J.B. Pendry, J. Phys. C: Solid St. Phys. 13 (1980) 937-44
         s_yt = numpy.array(yt[tnstart:tnstop]) # cut out the interval within which the comparison takes place
         s_ye = numpy.array(ye[enstart:enstop])
         te2 = integrate((s_yt-s_ye)**2,w_incr) # integrate (yt-ye)^2
         t2e2 = integrate(s_yt**2+s_ye**2,w_incr) # integrate yt^2+ye^2
         r_pend = te2/t2e2
         sys.stdout.write("Pendry R-factor : %f, shift: %f\n" % (r_pend,shift))
         if (r_pend < min_pendry[0] ):
            min_pendry=[r_pend,shift]
      if "R1" in r:
         # see  M.A. van Hove, S.Y. Tong, and M.H. Elconin, Surfac Science 64 (1977) 85-95
         r1 = integrate(abs(s_espec-c_scale*s_tspec),w_incr)/integrate(abs(s_espec),w_incr)
         sys.stdout.write("R1 R-factor     : %f, shift: %f\n" % (r1,shift))
         if (r1 < min_r1[0]):
            min_r1=[r1,shift]
      if "R2" in r:
         # see  M.A. van Hove, S.Y. Tong, and M.H. Elconin, Surfac Science 64 (1977) 85-95
         r2 = integrate((s_espec-c_scale*s_tspec)**2,w_incr)/integrate(s_espec**2,w_incr)
         sys.stdout.write("R2 R-factor     : %f, shift: %f\n" % (r2,shift))
         if (r2 < min_r2[0]):
            min_r2=[r2,shift]
      if "ZJ" in r:      
         # E. Zanazzi, F. Jona, Surface Science 62 (1977), 61-88
         s_d2_tspec = numpy.array(d2_tspec[tnstart:tnstop])
         s_d2_espec = numpy.array(d2_espec[enstart:enstop])

         epsilon = 0
         for i in s_d1_espec:
            if abs(i[1]) > epsilon:
               epsilon = abs(i[1])
         
         integrand = abs(c_scale*s_d2_tspec-s_d2_espec)*abs(c_scale*s_d1_tspec-s_d1_espec)/(abs(s_d1_espec)+epsilon)
         # interpolate integrand onto 10 times denser grid, see publication by Zanazzi & Jona
         incr = 0.1*w_incr
         grid_old = numpy.arange(0,len(integrand))*w_incr
         grid_new = numpy.arange(grid_old[0],grid_old[-1]+incr,incr)
         spl = interpolate.splrep(grid_old,integrand.T[1],k=3,s=0)
         integrand_dense = interpolate.splev(grid_new,spl,der=0)
         integrand_dense = numpy.vstack((grid_new,integrand_dense)).T
         # calculate reduced Zanazzi-Jona R-factor r=r/0.027
         r_zj = integrate(integrand_dense,incr)/(0.027*integrate(abs(s_espec),w_incr))
         sys.stdout.write("red. ZJ R-factor: %f, shift %f\n" % (r_zj,shift))
         if (r_zj < min_zj[0]):
            min_zj=[r_zj,shift]


# find minimal r-factor and write it out
   sys.stdout.write("\nMinimal r-factors:\n")
   if "pendry" in r:
      sys.stdout.write("minimal r-factor: Delta = %8.5f, Pendry R-factor = %7.5f \n" % ( min_pendry[1], min_pendry[0]))
   if "R1" in r:
      sys.stdout.write("minimal r-factor: Delta = %8.5f, R1 R-factor = %7.5f \n" % ( min_r1[1], min_r1[0]))
   if "R2" in r:
      sys.stdout.write("minimal r-factor: Delta = %8.5f, R2 R-factor = %7.5f \n" % ( min_r2[1], min_r2[0]))
   if "ZJ" in r:
      sys.stdout.write("minimal r-factor: Delta = %8.5f, ZJ R-factor = %7.5f \n" % ( min_zj[1], min_zj[0]))


if __name__ == "__main__":
      _main()

