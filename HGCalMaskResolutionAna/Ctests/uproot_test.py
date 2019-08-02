import matplotlib.pyplot as plt
import seaborn as sns
import uproot as up

file = up.open("file_after_weights.root")
#trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
tree = file["data"]
geneta = tree.array('geneta')
deltaE_corr = tree.array('deltaE_corr')

#ax = sns.distplot(f1[:,0], bins=100, color='g', kde=False)
ax = sns.jointplot(x=deltaE_corr[:,0], y=geneta, kind='hex', color='k', 
                   joint_kws=dict(gridsize=(100,100)))
ax.ax_marg_x.set_xlim(-1.2,1.2)
ax.ax_marg_y.set_ylim(2.74,3.04)
plt.savefig('pic.png')

#funciona com hexplot
ax = sns.jointplot(x=deltaE_corr[:,0], y=geneta, kind='kde', 
                   joint_kws=dict(gridsize=(100,100)),
                   color='k', cmap='Pastel2_r', shade=True)
ax.ax_marg_x.set_xlim(-1.2,1.2)
ax.ax_marg_y.set_ylim(2.74,3.04)
plt.savefig('pic2.png')
