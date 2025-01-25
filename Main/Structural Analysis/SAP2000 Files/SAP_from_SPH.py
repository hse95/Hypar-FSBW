# -*- coding: utf-8 -*-
"""
The code creates SAP2000 model for Kinetic Umbrella based on the parameters given to get the dynamic structural response.
It reads the pressure time history in pressure.csv as the input.
@author: Gaoyuan Wu
"""
#%%
# Preparatons: Importing Basic packages needed
import numpy as np  # Numerical packages
import matplotlib.pyplot as plt  # plotting
import pandas as pd # Pandas, transitions between csv/xlsx to python
import os
import sys
import comtypes.client
from scipy.signal import savgol_filter

gvty_water = 9.7706 # KN per m^3, density times GA g 

# Function that transform local para to Cartesian coordinates for given lambda and beta
def lbdToCart(lbd,beta,gamma,theta_star,b1,b2,c1,c2,r): #Note that lbd is an array
    Cart_coord = np.array(np.ones(3)) #Initialization of Cartesian 
    k1 = np.sqrt((b1**2) + (r*(1-lbd[0]))**2)
    k2 = np.sqrt((b2**2) + (r*(1-lbd[0]))**2) #k1 & k2
    theta_prime = np.arcsin(r*(1-lbd[0])/k1) + np.arcsin(r*(1-lbd[0])/k2) # Angle of upper panel
    theta_1 = theta_star - np.arctan(r * (1 - lbd[0]) / b1)
    theta_2 = theta_1 + theta_prime
    #print(lbd)
    #print(beta)
    #print(gamma)
    if gamma == 1: #Left panel
        c = c1
    if gamma == 2: #Right panel
        c = -c2
    if beta == 1:
        Cart_coord[0] = k1 * lbd[1] * np.cos(theta_1) #X
        Cart_coord[1] = c * lbd[0] #Y
        Cart_coord[2] = k1 * lbd[1] * np.sin(theta_1) #Z
    if beta == 2:
        Cart_coord[0] = k2 * lbd[1] * np.cos(theta_2) + k1 * np.cos(theta_1) #X
        Cart_coord[1] = c * lbd[0] #Y
        Cart_coord[2] = k2 * lbd[1] * np.sin(theta_2) + k1 * np.sin(theta_1)#Z
    return Cart_coord,k1,k2,theta_1,theta_2
        
def Macaulay(A): #Note that lbd is an array
    if A >= 0:
        return A
    else:
        return 0
    
def getCoordNForce(n_c1,n_c2,n_b1,n_b2,dw,theta_star,r,c1,c2,b1,b2,Hw,hb,d,Lw):
    gvty_water = 9.7706 # KN per m^3, density times GA g 
    # Before implementing lbdToCart function ,we have to get our lbd, beta & gamma first based on our SAP2000 settings
    # Nodes indexing, we specify bottom left node as NODE 1, and the indexing increases to the right; when reaching the right edge, it starts from the left at the adjacent row above

    nodes_num = (n_c1 + n_c2 + 1) * (n_b1 + n_b2 + 1) #Number of nodes
    nodes_coord = np.array(np.zeros((nodes_num,3))) #Array storing all nodes' coordinates
    lbd = np.array(np.zeros((nodes_num,2))) #Array storing all nodes' lambda (local para)
    gamma = np.array(np.zeros((nodes_num))) #Array storing all nodes' gamma (local para)
    beta = np.array(np.zeros((nodes_num))) #Array storing all nodes' beta (local para)
    k1 = np.array(np.zeros((nodes_num))) #Array storing all nodes' k1
    k2 = np.array(np.zeros((nodes_num))) #Array storing all nodes' k2
    theta_1 = np.array(np.zeros((nodes_num))) #Array storing all nodes' theta_1
    theta_2 = np.array(np.zeros((nodes_num))) #Array storing all nodes' theta_2
    
    x1 = np.array(np.zeros((nodes_num,3))) #Array storing tributary nodes
    x2 = np.array(np.zeros((nodes_num,3))) #Array storing tributary nodes
    x3 = np.array(np.zeros((nodes_num,3))) #Array storing tributary nodes
    x4 = np.array(np.zeros((nodes_num,3))) #Array storing tributary nodes
    
    #Local vector,initilization
    uk = np.array(np.zeros((nodes_num,3)))
    vk = np.array(np.zeros((nodes_num,3)))
    
    #Unit vector:
    nk = np.array(np.zeros((nodes_num,3)))
    
    #Pressure initialization:
    StaticP_val = np.array(np.zeros((nodes_num)))        #Static component, value
    StaticP_vec = np.array(np.zeros((nodes_num,3)))      #Static Component, vector
    StaticP_Area = np.array(np.zeros((nodes_num,3)))     #Static Component times tributary area
    

    WaveP_val = np.array(np.zeros((nodes_num)))        #Wave component, value
    WaveP_vec = np.array(np.zeros((nodes_num,3)))      #Wave Component, vector
    WaveP_Area = np.array(np.zeros((nodes_num,3)))     #Wave Component times tributary area
    
    Tri_A = np.array(np.zeros((nodes_num))) #Tributary Area
    
    # Get local para for nodes and transform to Cartesian
    for i in range (nodes_num):
        if 0 <= i % (n_c1 + n_c2 + 1) and i % (n_c1 + n_c2 + 1) <= n_c1:
            gamma[i] = 1
            lbd[i,0] = 1 - (i % (n_c1 + n_c2 + 1))/n_c1 #lambda 1
        if i % (n_c1 + n_c2 + 1) > n_c1:
            gamma[i] = 2
            lbd[i,0] = ((i % (n_c1 + n_c2 + 1)) - n_c1)/n_c2 #lambda 2
        if i <= ((n_c1 + n_c2 + 1) * (n_b1+1) - 1):
            beta[i] = 1
            lbd[i,1] = (int(i / (n_c1 + n_c2 + 1)))/ (n_b1)  #which row
        if i > ((n_c1 + n_c2 + 1) * (n_b1+1) - 1):
            beta[i] = 2 
            lbd[i,1] = ((int(i / (n_c1 + n_c2 + 1))) - n_b1) / n_b2
        i_coord,k1[i],k2[i],theta_1[i],theta_2[i] = lbdToCart(lbd[i,:],beta[i],gamma[i],theta_star,b1,b2,c1,c2,r)
        nodes_coord[i,0] = round(i_coord[0] , 3) 
        nodes_coord[i,1] = round(i_coord[1] , 3)
        nodes_coord[i,2] = round(i_coord[2] , 3)# Storing Cartesian coordinates 
        
    ####################### Coordinates Finished#############################
    
    ####################### Let's comupte pressure###########################
    
    
    ####################### Hydrostatic Component#########################
        lbd_10 = np.array([[0],[lbd[i,1]]]) #Set lbd1 = 0
        lbd_11 = np.array([[1],[lbd[i,1]]]) #Set lbd1 = 1
        lbd_21 = np.array([[lbd[i,0]],[1]]) #Set lbd2 = 1
        lbd_20 = np.array([[lbd[i,0]],[0]]) #Set lbd2 = 0
        vk[i,:] = (lbdToCart(lbd_21,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0] - (lbdToCart(lbd_20,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0] #Local vector vk
        
        if gamma[i] == 1: #Left panel
            uk[i,:] = (lbdToCart(lbd_10,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0] - (lbdToCart(lbd_11,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0] #Local vector uk
        if gamma[i] == 2: #Right Panel
            uk[i,:] = (lbdToCart(lbd_11,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0] - (lbdToCart(lbd_10,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0] #Local vector uk
   
        cross_p = np.cross(vk[i,:],uk[i,:])
        nk[i,:] = cross_p / np.linalg.norm(cross_p) # Unit vector  
        #Compute hydrostatic pressure on node
        h = np.sqrt(b1**2+r**2) * np.sin(theta_star - np.arctan(r/b1)) #Verte height
        if beta[i] == 1:
            if (dw - k1[i] * lbd[i,1] * np.sin(theta_1[i])) >= 0:
                StaticP_val[i] = gvty_water * (dw - k1[i] * lbd[i,1] * np.sin(theta_1[i]))
            else:
                StaticP_val[i] = 0
        if beta[i] == 2:
            if (dw - (r * lbd[i,0] * np.cos(theta_star) + h + k2[i] * lbd[i,1] * np.sin(theta_2[i]))) >= 0:
                StaticP_val[i] = gvty_water * (dw - (r * lbd[i,0] * np.cos(theta_star) + h + k2[i] * lbd[i,1] * np.sin(theta_2[i])))
            else:
                StaticP_val[i] = 0
                
        StaticP_vec[i,:] = StaticP_val[i] * nk[i,:] # Pressure vector at nodes
        
        # Compute acting area of the pressure
        if gamma[i] == 1:
            delta_lbd1 = 1 / (2*n_c1) #Left Panel diff
        if gamma[i] == 2:
            delta_lbd1 = 1 / (2*n_c2) #Right Panel diff
        if beta[i] == 1:
            delta_lbd2 = 1 / (2*n_b1) #Lower Panel diff
        if beta[i] == 2:
            delta_lbd2 = 1 / (2*n_b2) #Upper Panel diff
        
        if lbd[i,0] == 1 : #Lbd1 = 1 
            if lbd[i,1] == 1:
                tp_lbd = lbd[i,:]
                x1[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0]],[lbd[i,1] - delta_lbd2]])
                x4[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
            if lbd[i,1] == 0:
                tp_lbd = np.array([[lbd[i,0]],[lbd[i,1] + delta_lbd2]])
                x1[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = lbd[i,:]
                x4[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
            else:
                tp_lbd = np.array([[lbd[i,0]],[lbd[i,1] + delta_lbd2]])
                x1[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0]],[lbd[i,1] - delta_lbd2]])
                x4[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
            
        if lbd[i,0]!= 1:
            if lbd[i,1] == 1:
                tp_lbd = np.array([[lbd[i,0] + delta_lbd1],[lbd[i,1]]])
                x1[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0] + delta_lbd1],[lbd[i,1] - delta_lbd2]])
                x4[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
            if lbd[i,1] == 0:
                tp_lbd = np.array([[lbd[i,0] + delta_lbd1],[lbd[i,1] + delta_lbd2]])
                x1[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0] + delta_lbd1],[lbd[i,1]]])
                x4[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
            else:
                tp_lbd = np.array([[lbd[i,0] + delta_lbd1],[lbd[i,1] + delta_lbd2]])
                x1[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0] + delta_lbd1],[lbd[i,1] - delta_lbd2]])
                x4[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
        
        if lbd[i,0] == 0:
            if lbd[i,1] == 1:
                tp_lbd = lbd[i,:]
                x2[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0]],[lbd[i,1] - delta_lbd2]])
                x3[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
            if lbd[i,1] == 0:
                tp_lbd = np.array([[lbd[i,0]],[lbd[i,1] + delta_lbd2]])
                x2[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = lbd[i,:]
                x3[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
            else:
                tp_lbd = np.array([[lbd[i,0]],[lbd[i,1] + delta_lbd2]])
                x2[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0]],[lbd[i,1] - delta_lbd2]])
                x3[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
        if lbd[i,0] != 0:
            if lbd[i,1] == 1:
                tp_lbd = np.array([[lbd[i,0] - delta_lbd1],[lbd[i,1]]])
                x2[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0] - delta_lbd1],[lbd[i,1] - delta_lbd2]])
                x3[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
            if lbd[i,1] == 0:
                tp_lbd = np.array([[lbd[i,0] - delta_lbd1],[lbd[i,1] + delta_lbd2]])
                x2[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0] - delta_lbd1],[lbd[i,1]]])
                x3[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
            else:
                tp_lbd = np.array([[lbd[i,0] - delta_lbd1],[lbd[i,1] + delta_lbd2]])
                x2[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
                tp_lbd = np.array([[lbd[i,0] - delta_lbd1],[lbd[i,1] - delta_lbd2]])
                x3[i] = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[0]
        
        A_k = 0.5 * np.linalg.norm(np.cross((x2[i] - x4[i]),(x1[i] - x3[i]))) # Area of the quad
        Tri_A[i] = A_k
        StaticP_Area[i,:] = A_k * StaticP_vec[i,:]
        
        ########################## Wave Component########################################
        #print("node",i+1)
        #print("gamma=",gamma[i])
        #print("beta=",beta[i])
        #print("lbd=",lbd[i,:])
        #print("k1=",k1[i])
        #print("k2=",k2[i])
        #print("theta_1=",theta_1[i]*180/np.pi)
        #print("theta_2=",theta_2[i]*180/np.pi)
        #print("z=" ,nodes_coord[i,2])
        #print("dw=",dw)
        #print("dw+yita_star",(dw + 1.5 * Hw))
        # The Goda pressure coefficients
        alpha_1 = 0.6 + 0.5 * ((4 * np.pi * d)/(np.sinh(4 * np.pi * d / Lw) * Lw)) ** 2
        alpha_2 = min(((hb - dw)/(3 * hb)) * (Hw / dw)**2, 2 * dw / Hw)
        alpha_3 = 1 - dw/d * (1 - 1/np.cosh(2 * np.pi * d/Lw))

        # The Goda wave pressure distribution
        p1 = (alpha_1 + alpha_2) * gvty_water * Hw
        p3 = alpha_3 * p1
        #print("alpha3 =",alpha_3)
        # Parameters used based on Geometry 
        tp_lbd = np.array([[0], [lbd[i,1]]]) #Set lbd1 == 0
        tp_k1,none,tp_theta1 = (lbdToCart(tp_lbd,beta[i],gamma[i],theta_star,b1,b2,c1,c2,r))[1:4]
        h = tp_k1 * np.sin(tp_theta1) 
        hv = r * lbd[i,0] * np.cos(theta_star) + h
        #print("h=",h)
        #print("hv=",hv)

        # Parameters used based on wave para
        yita_star = 1.5 * Hw
        l1_prime = Macaulay(min(hv,dw + yita_star) - dw)
        l1_prime = l1_prime / np.sin(theta_1[i])
        l2_prime = Macaulay(dw + yita_star - max(hv,dw))
        l2_prime = l2_prime / np.sin(theta_2[i])
        l1 = min(dw,hv) / np.sin(theta_1[i])
        l2 = Macaulay(dw - hv) / np.sin(theta_2[i])
        #print("l1=",l1)
        #print("l2=",l2)
        #print("l1_prime=",l1_prime)
        #print("l2_prime=",l2_prime,"\n")

        ######## Wave Component Values #######
        if beta[i] == 1: 
            tp_wave1 = p3 + (p1 - p3) / (l1 + l2) * lbd[i,1] * k1[i]
            tp_wave2 = (p1 * (l1  + l2 + l1_prime + l2_prime) - p1 * lbd[i,1] * k1[i]) / (l1_prime + l2_prime)
            tp_waveratio = ((l1  + l2 + l1_prime + l2_prime) - lbd[i,1] * k1[i]) / (l1_prime + l2_prime)
            tp_wave = min(tp_wave1,tp_wave2)
            #print("wave1" ,tp_wave1)
            #print("wave2" ,tp_wave2)
            #print("2ratio",tp_waveratio)
            #print("tp_wave",tp_wave)
            WaveP_val[i] = Macaulay(tp_wave)
            #print(WaveP_val[i],'\n')
        if beta[i] == 2:
            tp_wave1 = p3 + (p1 - p3) / (l1 + l2) * (lbd[i,1] * k2[i] + k1[i])
            tp_wave2 = ((p1 * (l1  + l2 + l1_prime + l2_prime) - p1 * (lbd[i,1] * k2[i] + k1[i])) / (l1_prime + l2_prime))
            tp_waveratio = ((l1  + l2 + l1_prime + l2_prime) - (lbd[i,1] * k2[i] + k1[i])) / (l1_prime + l2_prime)
            tp_wave = min(tp_wave1,tp_wave2)
            #print("wave1" ,tp_wave1)
            #print("wave2" ,tp_wave2)
            #print("2ratio",tp_waveratio)
            #print("tp_wave",tp_wave)
            WaveP_val[i] = Macaulay(tp_wave)
            #print(WaveP_val[i],'\n')

            #WaveP_val[i] = Macaulay(min(p3 + (p1 - p3) / (l1 + l2) * (lbd[i,1] * k2[i] + k1[i]), (p1 * (l1  + l2 + l1_prime + l2_prime) - p1 * (lbd[i,1] * k2[i] + k1[i]) / (l1_prime + l2_prime))))
        
        WaveP_vec[i,:] = WaveP_val[i] * nk[i,:] # Wave Pressure vector at nodes
        WaveP_Area[i,:] = A_k * WaveP_vec[i,:]  # Wave Pressure times tributary area


    return nodes_coord, StaticP_Area, StaticP_val, StaticP_vec, WaveP_Area, WaveP_vec, WaveP_val,nk,Tri_A


def getShellConnectivity(n_c1,n_c2,n_b1,n_b2):
    area_num = (n_c1 + n_c2) * (n_b1 + n_b2) #Area numbers
    area_cnct = np.array(np.zeros((area_num,4))) #Array storing connectivity of Shell elements
    
    for i in range(area_num):
        row = int (i/(n_c1 + n_c2)) + 1 #Which row the area is in
        area_cnct[i,0] = int((row-1)*(n_c1 + n_c2 + 1) + 1 + (i - (row -1)*(n_c1 + n_c2))) #Joint 1 of the shell element
        area_cnct[i,1] = int(area_cnct[i,0] + (n_c1 + n_c2 + 1))#Joint 2 of the shell element
        area_cnct[i,2] = int(area_cnct[i,0] + (n_c1 + n_c2 + 2))#Joint 3 of the shell element
        area_cnct[i,3] = int(area_cnct[i,0] + 1) #Joint 4 of the shell element
        
    return area_cnct

''' 
The following use SAP2000's OAPI to generate the model.
The first is system settings.
'''
def StartSAPModel(Project_Name,APIPath= 'D:\OAPI_Example'):
    #set the following flag to True to attach to an existing instance of the program
    
    #otherwise a new instance of the program will be started
    
    AttachToInstance = False
    
     
    
    #set the following flag to True to manually specify the path to SAP2000.exe
    
    #this allows for a connection to a version of SAP2000 other than the latest installation
    
    #otherwise the latest installed version of SAP2000 will be launched
    
    SpecifyPath = False
    
     
    
    #if the above flag is set to True, specify the path to SAP2000 below
    
    ProgramPath = 'C:\Program Files (x86)\Computers and Structures\SAP2000 20\SAP2000.exe'
    
    
    if not os.path.exists(APIPath):
    
            try:
    
                os.makedirs(APIPath)
    
            except OSError:
    
                pass
    
    ModelPath = APIPath + os.sep + Project_Name
    
    
    if AttachToInstance:
    
        #attach to a running instance of SAP2000
    
        try:
    
            #get the active SapObject
    
            mySapObject = comtypes.client.GetActiveObject("CSI.SAP2000.API.SapObject")
    
        except (OSError, comtypes.COMError):
    
            print("No running instance of the program found or failed to attach.")
    
            sys.exit(-1)
    
    else:
    
        #create API helper object
    
        helper = comtypes.client.CreateObject('SAP2000v20.Helper')
    
        helper = helper.QueryInterface(comtypes.gen.SAP2000v20.cHelper)
    
        if SpecifyPath:
    
            try:
    
                #'create an instance of the SAPObject from the specified path
    
                mySapObject = helper.CreateObject(ProgramPath)
    
            except (OSError, comtypes.COMError):
    
                print("Cannot start a new instance of the program from " + ProgramPath)
    
                sys.exit(-1)
    
        else:
    
            try:
    
                #create an instance of the SAPObject from the latest installed SAP2000
    
                mySapObject = helper.CreateObjectProgID("CSI.SAP2000.API.SapObject")
    
            except (OSError, comtypes.COMError):
    
                print("Cannot start a new instance of the program.")
    
                sys.exit(-1)
    
     
    
        #start SAP2000 application
    
        mySapObject.ApplicationStart()
    
    
    '''
    From the codes above, the SAP2000 Program stars
    '''
    #create SapModel object
    
    SapModel = mySapObject.SapModel
    
     
    
    #initialize model
    
    SapModel.InitializeNewModel()
    
     
    
    #create new blank model
    
    ret = SapModel.File.NewBlank()
    
    return SapModel,ModelPath,mySapObject

'''
Units used: KN,m
'''
def unit(SapModel,unit_index = 6):
    ret = SapModel.SetPresentUnits(unit_index)
    

'''
Concrete Shell Section
SAP Unit:kn,m,c

For Hpar:

    Material:
    Self weight = 23.52 kn/m^3
    E = 37900000
    Poisson U = 0.2
    Thermal Expansion = 9.9e-06
    f'c= 65000
    expected compressive strength = 65000
    
    Shell section:
    Thickness = 100 mm -> 0.1m
    Type: thick
'''
def section_concrete_shell(SapModel,mat_name,section_name,gamma = 23.52, E = 4.4*10**7, U = 0.2, T_E = 9.9e-06,t = 0.1):
    
    #material
    material_conc = 2 #In SAP2000, index 2 if for concrete
    ret = SapModel.PropMaterial.SetMaterial(mat_name, material_conc)
    
    #assign isotropic mechanical properties to material
    ret = SapModel.PropMaterial.SetMPIsotropic(mat_name, E, U, T_E)
    
    #define shell section
    thick_shell = 2 #identifiers for thick shell
    ret = SapModel.PropArea.SetShell_1(section_name,thick_shell,0, MatProp=mat_name,MatAng=0,Thickness = t,Bending = t)
    
    return section_name
    
    
'''
Steel Channel Frame
SAP Unit:kn,m,c

For reinforcement at four edges
    Material:
        Self weight = 76.9729
        E = 1.999E+08
        Poisson U = 0.3
        Thermal Expansion = 1.170E-05
    
    Cross section:
        Angle: 125PFC
        t3 = 0.125
        t2 = 0.065
        tf = 7.5*10**-3
        tw = 4.7*10**-3       
'''

def section_steel_channel(SapModel,mat_name,section_name,gamma = 76.9729, E = 1.999E+08, U = 0.3, T_E = 1.170E-05,t3 = 0.125,t2 = 0.065,tf = 7.5*10**-3,tw = 4.7*10**-3):
    
    #material
    material_steel = 1 #In SAP2000, index 1 if for concrete
    ret = SapModel.PropMaterial.SetMaterial(mat_name, material_steel)
    
    #assign isotropic mechanical properties to material
    ret = SapModel.PropMaterial.SetMPIsotropic(mat_name, E, U, T_E)
    
    
    #define frame-vhannel section
    ret = SapModel.PropFrame.SetChannel(section_name,mat_name,t3,t2,tf,tw)
    
    return section_name
    
'''
Link elements for hinge mechanism.
U1->10000000
'''

def link_hinge(SapModel,hinge_name,DOF=[True,0,0,0,0,0],ke=[10000000,0,0,0,0,0],ce=[10000000,0,0,0,0,0]):
    fixed = [0,0,0,0,0,0]
    ret = SapModel.PropLink.SetLinear(hinge_name, DOF, fixed, ke, ce, 0,0)
    ret = SapModel.PropLink.SetSpringData(hinge_name, 1,1)
    return hinge_name

'''
Gap for base support->Side & Mid
'''
def gap_side(SapModel,gap_side_name,DOF=[True,True,True,0,0,0],Nonlinear=[True,0,0,0,0,0],ce=[0,0,0,0,0,0],ke=[10000000,10000000,10000000,0,0,0],k=[10000000,0,0,0,0,0]):
    fixed = [0,0,0,0,0,0]
    dis = [0,0,0,0,0,0]
    ret = SapModel.PropLink.SetGap(gap_side_name, DOF,fixed, Nonlinear, ke,ce, k, dis,0, 0)
    ret = SapModel.PropLink.SetSpringData(gap_side_name, 1,1)

def gap_mid(SapModel,gap_mid_name,DOF=[True,True,0,0,0,0],Nonlinear=[True,0,0,0,0,0],ce=[0,0,0,0,0,0],ke=[10000000,10000000,0,0,0,0],k=[10000000,0,0,0,0,0]):
    fixed = [0,0,0,0,0,0]
    dis = [0,0,0,0,0,0]
    ret = SapModel.PropLink.SetGap(gap_mid_name, DOF,fixed, Nonlinear, ke,ce, k, dis,0, 0)
    ret = SapModel.PropLink.SetSpringData(gap_mid_name, 1,1)
    
def CreateNode(SapModel,Sym_coord):
    node_num = Sym_coord.shape[0] #number of nodes
    for i in range (node_num):
        x = Sym_coord[i,0]
        y = Sym_coord[i,1]
        z = Sym_coord[i,2]
        ret = SapModel.PointObj.AddCartesian(x,y,z)

def CreateArea(SapModel,Area_cnct):
    area_num = Area_cnct.shape[0] #number of elements
    for i in range (area_num):
        v_1 = str(int(Area_cnct[i,0]))
        v_2 = str(int(Area_cnct[i,1]))
        v_3 = str(int(Area_cnct[i,2]))
        v_4 = str(int(Area_cnct[i,3])) 
        pts= [v_1,v_2,v_3,v_4]
        ret = SapModel.AreaObj.AddByPoint(4,pts)

def AssignProp_Area(SapModel,Area_name,Prop_Name):
    ret = SapModel.AreaObj.SetProperty(Area_name,Prop_Name)
'''
th is an array with time and value[time,value]
'''
def TimeHistory(SapModel,th,th_name):
    number_times = th.shape[0]
    time = np.ndarray.tolist(th[:,0])
    value = np.ndarray.tolist(th[:,1])
    ret = SapModel.Func.FuncTH.SetUser(th_name, number_times, time, value)

'''
Read the pressure_Press.csv to create a np.array.
Get rid of the negative values led by numerical instability 
Smooth the curve: from scipy.signal import savgol_filter

'''
def SPH_csv_To_Array(SapModel,file_name,area_cnct):
    area_num = area_cnct.shape[0] #number of areas
    og_dtframe = pd.read_csv(file_name,sep=";")
    og_array = og_dtframe.to_numpy() #to numpy array; first row time
    dt = round(float(og_array[254,1]),4) # time step
    n_array = og_array[254:,1:] # only time and values, get rid of the first 50 time steps
    smth_array = n_array #before smoothing
    for i in range(area_num):
        for j in range(n_array.shape[0]):
            if i == 0:
                smth_array[j,0]=j*dt
            if float(n_array[j,i+1]) <0:
                n_array[j,i+1] = 0
            else:
                n_array[j,i+1] = float(n_array[j,i+1])
        #temp_smth = savgol_filter(n_array[:,i+1],51,7)
        #smth_array[:,i+1] = temp_smth
    
    return dt,n_array

'''
Read the pressure_Press.csv to create a np.array.
Get the maximum value to create a equivalent static case

'''
def SPH_csv_To_Array_Max(SapModel,file_name,area_cnct):
    area_num = area_cnct.shape[0] #number of areas
    og_dtframe = pd.read_csv(file_name,sep=";")
    og_array = og_dtframe.to_numpy() #to numpy array; first row time
    dt = round(float(og_array[254,1]),4) # time step
    n_array = og_array[254:,2:] # only values, get rid of the first 50 time steps
    max_array = np.zeros(area_num) #array that storing the max pressure value
    for i in range(area_num):
        local_max = 0
        for j in range(n_array.shape[0]):
            if float(n_array[j,i]) >= local_max:
                local_max = float(n_array[j,i])
        max_array[i] = local_max #storing the max pressure
        
    
    return max_array

'''
Read the pressure_Press.csv to create a np.array.
Get the maximum value to create a equivalent static case

'''
def SPH_csv_To_Array_T(SapModel,file_name,area_cnct):
    area_num = area_cnct.shape[0] #number of areas
    og_dtframe = pd.read_csv(file_name,sep=";")
    og_array = og_dtframe.to_numpy() #to numpy array; first row time
    n_array = og_array[117,2:] # 117 row
    for i in range(area_num):
        n_array[i] = float(n_array[i])
        
    
    return n_array

#%%
def Create_SAP_Hypar_Sea(Project_Name,APIPath,n_c1,n_c2,n_b1,n_b2,theta_star,r,c1,c2,b1,b2,h_off,col_b,col_h,press_file,num_dt,time_step,Hw,dw,hb,d,Lw,E_c = 3.79*10**7,gamma_shell = 23.52, E_shell = 3.79*10**7, U_shell = 0.2, T_E_shell = 9.9e-06,t_shell = 0.1):
    SapModel,ModelPath,my_sapobj = StartSAPModel(Project_Name,APIPath) #Start SAP2000
    
    #Unit to KN,M
    unit(SapModel,unit_index = 6)
    
    #Create Material & Section
    section_concrete_shell(SapModel,"CONC_SHELL","SHELL",gamma_shell, E_shell, U_shell, T_E_shell,t_shell)
    #section_steel_channel(SapModel,"STEEL_REIN","REINFORCEMENT",gamma = 76.9729, E = 1.999E+08, U = 0.3, T_E = 1.170E-05,t3 = 0.125,t2 = 0.065,tf = 7.5*10**-3,tw = 4.7*10**-3)
    
    #Links Definition
    link_hinge(SapModel,"hinge",DOF=[True,0,0,0,0,0],ke=[10000000,0,0,0,0,0],ce=[10000000,0,0,0,0,0])
    gap_side(SapModel,"gap_side",DOF=[True,True,True,0,0,0],Nonlinear=[True,0,0,0,0,0],ce=[0,0,0,0,0,0],ke=[10000000,10000000,10000000,0,0,0],k=[10000000,0,0,0,0,0])
    gap_mid(SapModel,"gap_mid",DOF=[True,True,0,0,0,0],Nonlinear=[True,0,0,0,0,0],ce=[0,0,0,0,0,0],ke=[10000000,10000000,0,0,0,0],k=[10000000,0,0,0,0,0])
    
    #Hypar Geometry, FEM discretization.
    Sym_coord,StaticP_Area, StaticP_val, StaticP_vec, WaveP_Area, WaveP_vec, WaveP_val,nk,Tri_A = getCoordNForce(n_c1,n_c2,n_b1,n_b2,dw,theta_star,r,c1,c2,b1,b2,Hw,hb,d,Lw) #Get coordinates & nodes
    Area_cnct = getShellConnectivity(n_c1,n_c2,n_b1,n_b2)
    #Create nodes & area
    CreateNode(SapModel,Sym_coord)
    CreateArea(SapModel,Area_cnct)
    area_num = Area_cnct.shape[0]
    node_num = Sym_coord.shape[0]
    #Area assignment 
    for i in range(area_num):
        AssignProp_Area(SapModel, "{}".format(i+1), "SHELL")
    
    '''
    Assign Base Link.
    Nodes at the base: n_c1+n_c2+1
    '''
    #base_nodes_num = n_c1+n_c2+1
    #for i in range(base_nodes_num):
    #    if i == 0 or i == (base_nodes_num-1): #On two sides,lateral constraints
    #        ret = SapModel.LinkObj.AddByPoint("{0}".format(i+1), "", "base_{}".format(i), True,"gap_side") 
    #    else:
    #        ret = SapModel.LinkObj.AddByPoint("{0}".format(i+1), "", "base_{}".format(i), True,"gap_mid") 
            
    '''
    Assign Hinge Links and Column Stiffness
    '''    
    #Hinge coordinates
    vertex_coord = Sym_coord[(n_c1+n_c2+1)*(n_b1+1)-(n_c2+1),:]
    mid_top_bot = 0.5*(Sym_coord[n_c1,:]+Sym_coord[(n_c1+n_c2+1)*(n_b1+n_b2+1)-n_c2-1])
    unit_vecto_hinge= (vertex_coord-mid_top_bot)/np.linalg.norm(vertex_coord-mid_top_bot)
    hinge_coord = vertex_coord + h_off*unit_vecto_hinge
    
    #Create hinge node
    #ret = SapModel.PointObj.AddCartesian(hinge_coord[0],hinge_coord[1],hinge_coord[2])
    
    #Column stiffness, Assign Springs to the hinge do in x and y
    #K_cx = 3*E_c*(col_h*col_b**3)/(12*hinge_coord[2]**3)
    #K_cy = 3*E_c*(col_b*col_h**3)/(12*hinge_coord[2]**3)
    #hinge_index = (n_c1+n_c2+1)*(n_b1+n_b2+1)+1
    #ret = SapModel.PointObj.SetSpring("{0}".format(hinge_index),[K_cx,K_cy,0,0,0,0]) #Set Spring to the Hinge
    #ret = SapModel.PointObj.SetRestraint("{0}".format(hinge_index),[False,False,True,True,False,True])
    
    '''#Hinge ZoneRigid Link
                row_1_cen = (n_c1+n_c2+1)*(n_b1-1)+(n_c1+1)
                row_2_cen = row_1_cen+(n_c1+n_c2+1)
                row_3_cen = row_2_cen+(n_c1+n_c2+1)
                mid_hinge = [row_1_cen,row_2_cen,row_3_cen]
                row_node_num = n_c1+n_c2+1 # number of nodes each row
            
                hinge_link_index = 0
                for i in range(3):
                    if i != 2:
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]-1), "{0}".format(mid_hinge[i]-1+row_node_num), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]-1), "{0}".format(mid_hinge[i]), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]), "{0}".format(mid_hinge[i]+row_node_num), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]), "{0}".format(mid_hinge[i]+1), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]+1), "{0}".format(mid_hinge[i]+1+row_node_num), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]-1), "{0}".format(hinge_index), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]), "{0}".format(hinge_index), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]+1), "{0}".format(hinge_index), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                    else:
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]-1), "{0}".format(mid_hinge[i]), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]), "{0}".format(mid_hinge[i]+1), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]-1), "{0}".format(hinge_index), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]), "{0}".format(hinge_index), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        hinge_link_index +=1
                        ret = SapModel.LinkObj.AddByPoint("{0}".format(mid_hinge[i]+1), "{0}".format(hinge_index), "hinge_{}".format(hinge_link_index), False,"link_hinge") 
                        '''
    #refresh view, update (initialize) zoom
    ret = SapModel.View.RefreshView(0, False)
    
    ##############################################
    #SPH time history loading case + Dead loading#
    ##############################################
    '''
    Define the time-histories, assign the loads to the structure
    '''          
    # Time history of each element
    area_num = Area_cnct.shape[0] #number of elements
    dt,th_array = SPH_csv_To_Array(SapModel,press_file,Area_cnct)
    for i in range(area_num):
        temp_th = np.vstack((th_array[:,0].T,th_array[:,i+1].T)).T
        TimeHistory(SapModel,temp_th,"element_{}".format(i+1))
    

    
    #Add wave load pattern. Assign each load pattern with a unit pressure on the element.   
    for i in range(area_num):
        ret = SapModel.LoadPatterns.Add("SPH_wave_{}".format(i+1),14,0,False)
        ret = SapModel.AreaObj.SetLoadSurfacePressure("{}".format(i+1),"SPH_wave_{}".format(i+1),-1,1) #bottom face,value is 1KN/m^2


    #Modal analysis -80 modes so that 90% modal participation is achieved
    ret = SapModel.LoadCases.ModalEigen.SetCase("Modal_80")
    ret = SapModel.LoadCases.ModalEigen.SetNumberModes("Modal_80", 80, 1)


    ################################################################# "SPH-Wave-with-damping"
    #Set dynamic cases + Dead Load
    ret = SapModel.LoadCases.DirHistNonlinear.SetCase("SPH_WAVE")
    
    load_type = []
    load_name = []
    th_func = []
    l_scale = []
    t_scale = []
    at = []
    Csys = []
    ang = []
    for i in range(area_num+1):
        load_type.append("Load")
        at.append(0)
        ang.append(0)
        Csys.append("Global")
        if i == 0:
            load_name.append("DEAD")
            th_func.append("UNIFTH")
            l_scale.append(1)
            t_scale.append(1)
        else:
            load_name.append("SPH_wave_{}".format(i))
            th_func.append("element_{}".format(i))
            l_scale.append(0.001) #PA TO KPA
            t_scale.append(1)
    #

    #Set simulation time and time step
    ret = SapModel.LoadCases.DirHistNonlinear.SetTimeStep("SPH_WAVE",num_dt,time_step)

    #Set proportional damping
    ret = SapModel.LoadCases.DirHistNonlinear.SetDampProportional("SPH_WAVE", 3, 0, 0, 6.754, 223.54, 0.05, 0.05)

    #Set time integrator
    ret = SapModel.LoadCases.DirHistNonlinear.SetTimeIntegration("SPH_WAVE", int(1), 0.25, 0.5, 0.0, 0.0)

    #Set proportional damping
    ret = SapModel.LoadCases.DirHistNonlinear.SetDampProportional("SPH_WAVE", 3, 0, 0, 6.754, 223.54, 0.05, 0.05)

    #Set time integrator
    ret = SapModel.LoadCases.DirHistNonlinear.SetTimeIntegration("SPH_WAVE", int(1), 0.25, 0.5, 0.0, 0.0)
    
    #Add load cases
    ret = SapModel.LoadCases.DirHistNonlinear.SetLoads("SPH_WAVE",area_num+1,load_type,load_name,th_func,l_scale,t_scale,at,Csys,ang)

    

    
    
    



    ################################################################# "SPH-Wave-withpit-damping"
    #########################################################################

    #Set dynamic cases + Dead Load
    ret = SapModel.LoadCases.DirHistNonlinear.SetCase("SPH_WAVE_noD")
    
    load_type = []
    load_name = []
    th_func = []
    l_scale = []
    t_scale = []
    at = []
    Csys = []
    ang = []
    for i in range(area_num+1):
        load_type.append("Load")
        at.append(0)
        ang.append(0)
        Csys.append("Global")
        if i == 0:
            load_name.append("DEAD")
            th_func.append("UNIFTH")
            l_scale.append(1)
            t_scale.append(1)
        else:
            load_name.append("SPH_wave_{}".format(i))
            th_func.append("element_{}".format(i))
            l_scale.append(0.001) #PA TO KPA
            t_scale.append(1)
    
    
    #Set time integrator, Newmark (gamma=0.5,beta=0.05)
    ret = SapModel.LoadCases.DirHistNonlinear.SetTimeIntegration("SPH_WAVE_noD", int(1), 0.25, 0.5, 0.0, 0.0)

    #Set simulation time and time step
    ret = SapModel.LoadCases.DirHistNonlinear.SetTimeStep("SPH_WAVE_noD",num_dt,time_step)
    
    #Set proportional damping
    ret = SapModel.LoadCases.DirHistNonlinear.SetDampProportional("SPH_WAVE_noD", 1, 0, 0,0,0,0,0)


    #Add load cases
    ret = SapModel.LoadCases.DirHistNonlinear.SetLoads("SPH_WAVE_noD",area_num+1,load_type,load_name,th_func,l_scale,t_scale,at,Csys,ang)
    
    


    
    # Change Dead Load Case to Nonlinear
    ret = SapModel.LoadCases.StaticNonlinear.SetCase("DEAD")
    ret = SapModel.LoadCases.StaticNonlinear.SetLoads("DEAD",1,["Load"],["DEAD"],[1])
    
    #Save the model
    ret = SapModel.File.Save(ModelPath)
    
    return SapModel,area_num,node_num,my_sapobj
 
#%%
def Run_Sap(SapModel):
    ret = SapModel.Analyze.RunAnalysis()
    
#%%
''' This part defines a function that summarizes the critical structural response
For the shell elements:
    1. SPH_wave-dynamic time history analysis (Envelope)
        1. Moment: envelope min and max for M11 and M22, take the one that has the maximum absolute value
        2. Shear: V11,V22
        3. Membrane force: Nc and Nt
    2. Other cases:
        2.1 - Dead: SW
        2.2 - SPH_Static_max
        2.3 - Goda_static
        2.4 - Inundation + Dead load : the difference between 2.4 and 2.2/2.3/1 is the influence of waves
For the support column:
    1. Base Moment
    2. Shear
    3. Axial
'''
def Post_Sap(SapModel,area_num,node_num):
    
    
    #Cases number
    case_num = 6 # 6 cases, sequence from 1 to 2.4 
    
    #Arrays storing results of each cases
    '''
    Shell response
    Rows:
    1: Mmax
    2: Vmax
    3: F_C_MAX
    4:F_T_MAX
    5: Max_disp
    
    
    Column response (node 626)
    Rows:
    1:V-shear
    2:Moment(times height)
    3:F
    '''
    shell_res = np.zeros((5,case_num)) # empty array storing maximum shell response
    col_res = np.zeros((3,case_num)) #empty array storing support column reactions
    case_name = ["SPH_WAVE","SPH_WAVE_noD","Dead","SPH_Static_Max","Goda_Static","Inundation_Static"]
    
    for i in range(case_num):
        ret = SapModel.Results.Setup.DeselectAllCasesAndCombosForOutput() #Deselects cases chosen first
        
        
        ret = SapModel.Results.Setup.SetCaseSelectedForOutput("{}".format(case_name[i])) #Selects the case
        if i == 0: 
            ret = SapModel.Results.Setup.SetOptionDirectHist(1) #For time history analysis, we care for the envelope
        
        
        max_m11 = 0
        max_m22 = 0
        max_v13 = 0
        max_v23 = 0
        max_f11_c =0
        max_f22_c = 0
        max_f11_t = 0
        max_f22_t = 0
        max_disp = 0


        #Iterate over each node:
        for j_node in range(node_num):
            GroupElm = 0
            NumberResults = []
            Obj = []
            Elm = []
            NumberResults = 0
            PointElm = []
            LoadCase = []
            StepType = []
            StepNum = []
            U1 = []
            U2 = []
            U3 = []
            R1 = []
            R2 = []
            R3 = []
            [GroupElm, NumberResults, Obj, Elm, LoadCase, StepType, StepNum, U1, U2, U3, R1, R2, R3] = SapModel.Results.JointDispl("{}".format(j_node+1), GroupElm, NumberResults, Obj, Elm, LoadCase, StepType, StepNum, U1, U2, U3, R1, R2, R3)
            temp_U1= max([abs(max(U1)),abs(min(U1))])
            temp_U2 = max([abs(max(U2)),abs(min(U2))])
            temp_U3 = max([abs(max(U3)),abs(min(U3))])
            temp_disp = np.sqrt(temp_U1**2+temp_U2**2+temp_U3**2)
            if temp_disp > max_disp:
                max_disp = temp_disp

        #Iterate over each element:
        for j in range(area_num):
            
            #Initializing arrays:
            M11 = []
            ObjectElm = 0
            NumberResults = []
            Obj = []
            NumberResults = 0
            PointElm = []
            LoadCase = []
            StepType = []
            F11 = []
            F22 = []
            F12 = []
            FMax = []
            FMin = []
            FAngle = []
            FVM = []
            M11 = []
            M22 = []
            M12 = []
            MMax = [] 
            MMin = []
            MAngle = []
            V13 = [] 
            V23 = [] 
            VMax = [] 
            VAngle= []
            Elm = []
            StepType = []
            StepNum = []
            
        
            [NumberResults, Obj, Elm, PointElm, LoadCase, StepType, StepNum, F11, F22, F12, FMax, FMin, FAngle, FVM, M11, M22, M12, MMax, MMin, MAngle, V13, V23, VMax, VAngle, ret] = SapModel.Results.AreaForceShell("{}".format(j+1), ObjectElm, NumberResults, Obj, Elm, PointElm, LoadCase, StepType, StepNum, F11, F22, F12, FMax, FMin, FAngle, FVM, M11, M22, M12, MMax, MMin, MAngle, V13, V23, VMax, VAngle)
            
            #M11
            temp_m11_max_pos = max(M11)
            temp_m11_max_neg = min(M11)
            temp_m11_max = max(abs(temp_m11_max_pos),abs(temp_m11_max_neg)) #maximum m11 
            if temp_m11_max >= max_m11:
                max_m11 = temp_m11_max
            #M22
            temp_m22_max_pos = max(M22)
            temp_m22_max_neg = min(M22)
            temp_m22_max = max(abs(temp_m22_max_pos),abs(temp_m22_max_neg)) #maximum m22
            if temp_m22_max >= max_m22:
                max_m22 = temp_m22_max
            
            #V13
            temp_v13_max_pos = max(V13)
            temp_v13_max_neg = min(V13)
            temp_v13_max = max(abs(temp_v13_max_pos),abs(temp_v13_max_neg)) #maximum v13
            if temp_v13_max >= max_v13:
                max_v13 = temp_v13_max
            
            #V23
            temp_v23_max_pos = max(V23)
            temp_v23_max_neg = min(V23)
            temp_v23_max = max(abs(temp_v23_max_pos),abs(temp_v23_max_neg)) #maximum v23
            if temp_v23_max >= max_v23:
                max_v23 = temp_v23_max
            
            #F11_C
            temp_f11_c_max = min(F11) #Compression, negative
            if temp_f11_c_max <= max_f11_c:
                max_f11_c = temp_f11_c_max
            
            #F22_C
            temp_f22_c_max = min(F22) #Compression, negative
            if temp_f22_c_max <= max_f22_c:
                max_f22_c = temp_f22_c_max
            
            #F11_T
            temp_f11_t_max = max(F11) #Tension, positive
            if temp_f11_t_max > max_f11_t:
                max_f11_t = temp_f11_t_max
            
            #F22_T
            temp_f22_t_max = max(F22) #Tension, positive
            if temp_f22_t_max > max_f22_t:
                max_f22_t = temp_f22_t_max
    
    
        #Storing it for each case, Shell Response
        #Maximum shell moment
        if max_m11 > max_m22:
            shell_res[0,i] = max_m11
        else:
            shell_res[0,i] = max_m22
            
        #Maximum shell shear force
        if max_v13 > max_v23:
            shell_res[1,i] = max_v13
        else:
            shell_res[1,i] = max_v23
        
        #Maximum shell membrane compression
        if max_f11_c < max_f22_c:
            shell_res[2,i] = max_f11_c
        else:
            shell_res[2,i] = max_f22_c
        
        #Maximum shell membrane tension
        if max_f11_t > max_f22_t:
            shell_res[3,i] = max_f11_t
        else:
            shell_res[3,i] = max_f22_t
        
        #Maximum shell displacement
        shell_res[4,i] = max_disp
         
        ###########################################
        #Storing it for each case, Column Response#
        ###########################################
        
        node_col = node_num+1 # node index of the column node
        x_col = 0
        y_col = 0
        z_col = 0
        [x_col, y_col, z_col,ret] = SapModel.PointObj.GetCoordCartesian(str(node_col), x_col, y_col, z_col)
        #Get the column reaction force from SAP
        ObjectElm = 0
        NumberResults = 0
        Obj = []
        Elm = []
        LoadCase = []
        StepType = []
        StepNum = []
        F1 = []
        F2 = []
        F3 = []
        M1 = []
        M2 = []
        M3 = []

        #
        [NumberResults, Obj, Elm, LoadCase, StepType, StepNum, F1, F2, F3, M1, M2, M3,ret] = SapModel.Results.JointReact("{}".format(int(node_col)), ObjectElm, NumberResults, Obj, Elm, LoadCase, StepType, StepNum, F1, F2, F3, M1, M2, M3)
        
        col_shear = max(max(F1),abs(min(F1)))
        col_moment = col_shear*z_col # Base moment
        if min(F3) >= 0:
            col_axial = -max(F3) # Positive F3 is compression
        else:
            col_axial = -min(F3) # Negative F3 is uplift
            
        col_res[0,i] = col_shear
        col_res[1,i] = col_moment
        col_res[2,i] = col_axial

    return shell_res,col_res



        

    
    
    
    