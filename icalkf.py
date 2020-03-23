# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



# %%
class measurements:
    def __init__(self,file_name):
        self.df = pd.read_csv(file_name,sep="\t")
        # y = df['posx']
        # z = df['posy']
        # x = df['posz']

    # get the measurements in event n
    def get_measurements(self, n=0): 
        df_temp = self.df.query("eid=={0}".format(n)).sort_values(by=["posy"])
        y = df_temp['posx']
        z = df_temp['posy']
        x = df_temp['posz']
        
        KE = np.max(df_temp['tot_E'])
        mom =  np.sqrt(np.max(df_temp['momx']**2 +  df_temp['momy']**2 + df_temp['momz']**2))

        # return a state vector with values initialized to the x,y near the vertex and the slopes from the first two planes, qp is made 0
        sv = state_vector(x[0],y[0],(x[1]-x[0])/(z[1]-z[0]),(y[1]-y[0])/(z[1]-z[0]),0)        
        
        return x,y,z,KE,mom,len(df_temp),sv
    
    # this will return the number of events in the dataframe
    def entries(self): 
        return 1 #self.df.
    
    # get the initial state vector for event n
    def get_initsv(self,n):
        return self.get_measurements(n)[6]

    # print details of event n
    def event_details(self,n):
        KE, mom, nmp = self.get_measurements(n)[3:6] 
        print("=================Event {0}========================================".format(n))
        print("E={0} GeV\t P={1} GeV/c\t Meas.planes={2}".format(np.around(KE,decimals=2),np.around(mom,decimals=2),nmp))
        print("================================================================")

    # the initial state vector 



#%%


class state_vector:
    def __init__(self,x0=0,y0=0,tx0=0,ty0=0,qp0=0):
        self.x = x0
        self.y = y0
        self.tx= tx0
        self.ty= ty0
        self.qp= qp0
        k  = 0.29979  # GeV/c/T/m

    def get_sv(self):
        return self.x,self.y,self.tx,self.ty,self.qp
    def set_sv(self,x0,y0,tx0,ty0,qp0):
        self.x = x0
        self.y = y0
        self.tx= tx0
        self.ty= ty0
        self.qp= qp0
    
    def print_sv(self):
        print(self.x,self.y,self.tx,self.ty,self.qp)
    def h(self):
        return self.k*self.qP*np.sqrt(1+self.tx**2+self.ty**2)


class kf:
    def __init__(self, Bx0=1.5, By0=0, Bz0=0):
        Bx = Bx0
        By = By0
        Bz = Bz0
        dz = dz0
        projector = [[1, 0, 0], [0, 1, 0]]

    def predict_sv(self,sv):
        hh = h(tx,ty,qP)

 
    def predict_x(self,sv,dz):
        return sv.x + sv.tx*dz + sv.h*(sv.tx*sv.ty*self.Sx(dz)- (1+sv.tx**2)*self.Sy(dz))+ sv.h**2*(sv.tx*(3*sv.ty*2+1)*self.Sxx(dz) - sv.ty*(3*sv.tx**2+1)*self.Sxy(dz) - sv.ty*(3*sv.tx**2+1)*self.Syx(dz) + sv.tx*(3*sv.tx**2+3)*self.Syy(dz))

    def predict_y(self,sv,dz):
        return sv.y + sv.ty*dz + sv.h*(-sv.tx*sv.ty*self.Sy(dz)+ (1+sv.ty**2)*self.Sx(sv.dz))+ sv.h**2*(sv.tx*(3*sv.ty*2+1)*self.Sxx(dz) - sv.ty*(3*sv.tx**2+1)*self.Sxy(dz) - sv.ty*(3*sv.tx**2+1)*self.Syx(dz) + sv.tx*(3*sv.tx**2+3)*self.Syy(dz))

    def predict_tx(self, sv, dz):
        return sv.tx + sv.h*(sv.tx*sv.ty*self.Rx(dz) - (1+sv.tx**2)*self.Ry(dz)) + sv.h**2*(sv.tx*(3*sv.ty**2 + 1)*self.Rxx(dz) - sv.ty*(3*sv.tx**2 + 1)*self.Rxy(dz) - sv.ty*(3*sv.tx**2+1)*self.Ryx(dz) + sv.tx(3*sv.tx**2+3)*self.Ryy(dz)) 
    
    def predict_ty(self, sv, dz):
        return sv.ty + sv.h*(-sv.tx*sv.ty*self.Ry(dz) + (1+sv.ty**2)*self.Rx(dz)) + sv.h**2*(sv.ty*(3*sv.ty**2 + 3)*self.Rxx(dz) - sv.tx*(3*sv.ty**2 + 1)*self.Rxy(dz) - sv.tx*(3*sv.ty**2+1)*self.Ryx(dz) + sv.ty(3*sv.tx**2+1)*self.Ryy(dz)) 

    def predict_qp(self, sv, dz):
        return 

    def Sx(self,dz):
        return 0.5*self.Bx*dz**2

    def Rx(self,dz):
        return self.Bx*dz

    def Sy(self,dz):
        return 0.5*self.By*dz**2

    def Ry(self,dz):
        return self.By*dz

    def Sxx(self,dz):
        return (self.Bx**2*dz**3)/6

    def Rxx(self,dz):
        return (self.Bx**2*dz**2)/2

    def Sxy(self,dz):
        return (self.Bx*self.By*dz**3)/6

    def Rxy(self,dz):
        return (self.Bx*self.By*dz**2)/2

    def Syx(self,dz):
        return (self.Bx*self.By*dz**3)/6

    def Ryx(self,dz):
        return (self.Bx*self.By*dz**2)/2

    def Syy(self,dz):
        return (self.By**2*dz**3)/6

    def Ryy(self,dz):
        return (self.By**2*dz**2)/2






def predict_x(x0,tx,ty,qP,dz):
    hh = h(tx,ty,qP)
    return x0 + tx*dz + hh*(tx*ty*Sx(dz)- (1+tx**2)*Sy(dz))+ hh**2*(tx*(3*ty*2+1)*Sxx(dz) - ty*(3*tx**2+1)*Sxy(dz) - ty*(3*tx**2+1)*Syx(dz) + tx*(3*tx**2+3)*Syy(dz))

    

# #print_slopes()56
# for i in range(560):
#     if(i==0):
#         x = predict_x(0.6008,0.00846977093333,-0.014892050833,4.53,0.01)
#     else: 
#         x = predict_x(x,0.00846977093333,-0.014892050833,4.53,0.01)
#     print(x)
#plt.plot(df["posy"],np.sqrt(df["momx"]**2 + df["momy"]**2 + df["momz"]**2))


# %%