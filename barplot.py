import numpy as np
import matplotlib.pyplot as plt
 
# set width of bar
barWidth = 0.25
fig = plt.subplots(figsize =(12, 8))
 
# set height of bar
rf = [0.73, 0.66]
mp = [0.99, 0.98]
 
# Set position of bar on X axis
br1 = np.arange(len(rf))
br2 = [x + barWidth for x in br1]

# Make the plot
plt.bar(br1, rf, color ='b', width = barWidth,
        edgecolor ='grey', label ='RF')
plt.bar(br2, mp, color ='g', width = barWidth,
        edgecolor ='grey', label ='MPGNN')

# Adding Xticks
plt.xlabel('Model', fontweight ='bold', fontsize = 15)
# plt.ylabel('', fontweight ='bold', fontsize = 15)
plt.xticks([r + barWidth for r in range(len(rf))],
        ['Balanced Accuracy', 'MCC'])
 
plt.legend()
plt.show()




