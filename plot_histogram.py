import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('csvs/pab_results.csv', header=0, index_col = None)


data['Firefly_Luciferase_interference'].plot(kind='bar')
plt.ylabel('Frequency')
plt.xlabel('interference')
plt.title('Title')

plt.show()

