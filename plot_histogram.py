import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('csvs/pad_results.csv', sep=',',header=None, index_col =1)

data.plot(kind='bar')
plt.ylabel('Frequency')
plt.xlabel('interference')
plt.title('Title')

plt.show()