#%%
import icalkf as f

m = f.measurements("5-Gev-mu-minus.txt")

isv = m.get_initsv(0)
isv.print_sv()
#%%