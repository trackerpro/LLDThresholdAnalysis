import pandas as pd
df = pd.read_csv("lumiByDay.csv")
print(df.columns)
df['integrated_recorded'] = df.iloc[:,1].cumsum();
df['integrated_delivered'] = df.iloc[:,2].cumsum();
df['integrated_recorded'] = df['integrated_recorded']*0.001*0.001*0.001;
df['integrated_delivered'] = df['integrated_delivered']*0.001*0.001*0.001;
for x,y in zip(df.iloc[:,0],df.iloc[:,3]):
    print(x,' ',y)
