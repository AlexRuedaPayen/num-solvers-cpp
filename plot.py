import json
import matplotlib.pyplot as plt

data_file=open('./output/CGvsPCG.json','r')
data=json.load(data_file)
data_file.close()
i=0

for d in data :
		i+=1
		ax=plt.subplot(2,2,i%4+1)
    		plt.plot(data[d]['cg']['resvec'],label="Conjugate gradient")
		plt.plot(data[d]['pcg']['resvec'],label="Preconditioned conjugate gradient")
                
		plt.legend(loc='best')

		plt.xlabel('Iteration number')
		plt.ylabel('$\log_{10}(\|r\|/\|b\|)$')

		plt.yscale("log")
                ax.title.set_text('File = '+d)
                
		if (i%4==0) :
                	plt.tight_layout()          
                	plt.savefig('./output/comparison'+str(i/4 +1)+'.png')
			plt.show()
