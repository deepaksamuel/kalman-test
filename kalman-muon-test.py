# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



# %%
class meas:
    def __init__(self,file_name):
        self.df = pd.read_csv(file_name,sep="\t")
        # y = df['posx']
        # z = df['posy']
        # x = df['posz']

    # get the measurements in event n
    def get_meas(self, n=0): 
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
    
    # print details of event n
    def event_details(self,n):
        KE, mom, nmp = self.get_meas(n)[3:6] 
        print("=================Event {0}========================================".format(n))
        print("E={0} GeV\t P={1} GeV/c\t Meas.planes={2}".format(np.around(KE,decimals=2),np.around(mom,decimals=2),nmp))
        print("================================================================")

    # the initial state vector 



#%%


class state_vector:
    def __init__(self,x0,y0,tx0,ty0,qp0):
        self.x = x0
        self.y = y0
        self.tx= tx0
        self.ty= ty0
        self.qp= qp0
        k  = 0.29979  # GeV/c/T/m

    def print_sv(self):
        print(self.x,self.y,self.tx,self.ty,self.qp)
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
    def __init__(self, Bx0=1.5, By0=0, Bz0=0,dz0=0.01):
        Bx = Bx0
        By = By0
        Bz = Bz0
        dz = dz0
        projector = [[1, 0, 0], [0, 1, 0]]
    
    def Sx(self):
        return 0.5*self.Bx*self.dz**2

    def Rx(self):
        return self.Bx*self.dz

    def Sy(self):
        return 0.5*self.By*self.dz**2

    def Ry(self):
        return self.By*self.dz

    def Sxx(self):
        return (self.Bx**2*self.dz**3)/6

    def Rxx(self):
        return (self.Bx**2*self.dz**2)/2

    def Sxy(self):
        return (self.Bx*self.By*self.dz**3)/6

    def Rxy(self):
        return (self.Bx*self.By*self.dz**2)/2

    def Syx(self):
        return (self.Bx*self.By*self.dz**3)/6

    def Ryx(self):
        return (self.Bx*self.By*self.dz**2)/2

    def Syy(self):
        return (self.By**2*self.dz**3)/6

    def Ryy(self):
        return (self.By**2*self.dz**2)/2










def measurements():
    #the coordinates are transformed to match with KB's Kalman filter
    y = df['posx']
    z = df['posy']
    x = df['posz']
    return x,y,z
def measurement(a):
    x,y,z = measurements()
    return x[a], y[a], z[a]
def print_slopes():
    x,y,z = measurements()
    xx =[]
    sx =[]
    sy =[]
    for i in range(meas_size()-1):
        slpx=(x[i+1]-x[i])/(z[i+1]-z[i])
        slpy=(y[i+1]-y[i])/(z[i+1]-z[i])
        print(i,slpx,slpy)
        xx.append(i)
        sx.append(slpx)
        sy.append(slpy)
    plt.plot(xx,sx)
    plt.plot(xx,sy)

def meas_size():
    return len(df['posx'])



def predict_x(x0,tx,ty,qP,dz):
    hh = h(tx,ty,qP)
    return x0 + tx*dz + hh*(tx*ty*Sx(dz)- (1+tx**2)*Sy(dz))+ hh**2*(tx*(3*ty*2+1)*Sxx(dz) - ty*(3*tx**2+1)*Sxy(dz) - ty*(3*tx**2+1)*Syx(dz) + tx*(3*tx**2+3)*Syy(dz))

    
# def predict_y(y0,tx,ty,qP,dz):
#     hh = h(tx,ty,qP)
#     return y0 + ty*dz + hh*(-tx*ty*Sy(dz)+ (1+ty**2)*Sx(dz))+ hh**2*(tx*(3*ty*2+1)*Sxx(dz) - ty*(3*tx**2+1)*Sxy(dz) - ty*(3*tx**2+1)*Syx(dz) + tx*(3*tx**2+3)*Syy(dz))

# #print_slopes()56
# for i in range(560):
#     if(i==0):
#         x = predict_x(0.6008,0.00846977093333,-0.014892050833,4.53,0.01)
#     else: 
#         x = predict_x(x,0.00846977093333,-0.014892050833,4.53,0.01)
#     print(x)
#plt.plot(df["posy"],np.sqrt(df["momx"]**2 + df["momy"]**2 + df["momz"]**2))


# %%