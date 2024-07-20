from os import read
from networkx import non_edges
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os

def read_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    data = [list(map(float, line.strip().split())) for line in lines]
    return data



if __name__ == '__main__':

    a0 = 1e-2
    ml=1
    mr=1
    fname = "./plots/Dfum_1e-2_ps_1_1_cpp3"

    def plot_phase(fname):
        pts= read_data("./data/FUM2_ps_1e-2_22_22_160_2.txt")
        speed = np.array([])
        phi = np.array([])
        for i in pts:
            speed = np.append(speed, i[1])
            phi = np.append(phi,i[0])
        
        plt.figure(figsize=(5.5,5))
        plt.scatter(phi, speed, s=0.1, color='black')
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$v$')
        plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        plt.title(r'$a_0 = $'+'{:.2e}'.format(a0)+','+r'$m_r/m_l = $'+"{:.3f}".format(mr/ml))

        fname += '_'+str(ml).replace('.', '')+'_'+str(mr).replace('.', '-')
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plot_file_path = os.path.join(script_dir, fname)
        print(plot_file_path)
        #plt.savefig(plot_file_path)
        plt.show()

    #plot_phase(fname)
    data = read_data("./data/FUM2_diff_01k_01k_1k_3.txt")
    #
    def D_n(data):
        plt.figure(figsize=(6.5,5))
        labels = ["1e-1","1e-2","1e-3","1e-4", "1e-5"]
        data = np.array(data)
        data = data[np.argsort(data[:,0])[::-1]]
        print(data)
        for index, row in enumerate(data):
            nonzero_values = np.abs(np.array([val for val in row if val != 0]))
            x_values = np.arange(0, len(nonzero_values)*100, 100)
            y_values = nonzero_values[:len(x_values)]
            x_values[0]=1
            print(y_values)
            plt.plot(x_values, y_values,marker='o',linestyle= '--',
                     label=r'$\epsilon$='+labels[index])
        # nonzero_values = np.abs(np.array([i for i in data if i != 0]))
        # x_values = np.arange(0, len(nonzero_values)*100, 100)
        #plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('N')
        plt.ylabel('D')
        plt.legend()
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        #plt.ylim(1e-11,5e-6)
        plt.title(r'$NP = 10^5$')
        script_dir = os.path.dirname(os.path.abspath(__file__))
        fname = "FUM2_diff_N_2.png"
        plot_file_path = os.path.join(script_dir, fname)

        plt.savefig(plot_file_path)
        plt.show()

    #D_n(data)
    #data = read_data("./data/D_mr_eps_dep.txt")
    def D_m():
        plt.figure(figsize=(6.5,5))
        labels = ["1e-1","1e-2","1e-3","1e-4", "1e-5"]
        #labels=["1e-3", "1e-2","1", "1e2", "1e3"]
        #labels= ["100","200","500","100","200","500"]
        #color = ['green', "orange", "grey" ]

        initial_values = [1e-5,5e-5,1e-4,5e-4,0.001, 0.005, 0.01, 0.05, 0.1, 0.5]

        # Add 1.0 to the list
        x_values = initial_values + [1.0]

        # Add the reversed reciprocals of the initial values
        x_values += list(reversed([1.0 / x for x in initial_values]))
        data = read_data(f"./data/Gfum_01k_01k_Dm16_new.txt")
        data = np.array(data)
        print(data)
     
        #x = np.array([1e-1,1e-2,1e-3,1e-4, 1e-5])[::-1]
        data = np.abs(data[np.argsort(data[:,-1])[::-1]])
        for index, row in enumerate(data):
            nonzero_values = np.abs(np.array([val for val in row if val != 0]))

            y_values = nonzero_values[:len(x_values)]
            print(y_values)
            plt.plot(x_values, y_values,marker='o',linestyle= '--',
                             label=r'$\epsilon$='+ labels[index])
            # if (index>2):
            #     plt.plot(x_values, y_values,marker='o',linestyle= '--',color=color[index-3],
            #             label=r'$N$='+labels[index])
            # else:
            #     plt.plot(x_values, y_values,marker='x',linestyle= '--',color=color[index])
        # nonzero_values = np.abs(np.array([i for i in data if i != 0]))
        #plt.plot(0,0,c="black",marker='x', label=r'$\epsilon =10^{-2}$')
        #plt.plot(0,0,c="black", marker='o',label=r'$\epsilon = 10^{-4}$')

        #x_values = np.arange(0, len(nonzero_values)*100, 100)
        
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$m_l/mr$')
        #plt.xlabel(r'$\epsilon^2$')
        
        plt.ylabel('D')
        plt.legend()
        #plt.ylim(1e-11,5e-6)
        plt.title(r'$NP = 10^4, N=100 $')
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        fname = "Gfum_Dm_16b.png"
        plot_file_path = os.path.join(script_dir, fname)

        plt.savefig(plot_file_path)
        plt.show()

    #D_n(data)
    D_m()

    def Desq(data):
        x_val = np.array([1e-5,1e-4,1e-3,1e-2,1e-1])
        lab = np.array(["1e-4", "1e-2",  "1", "1e2", "1e4"])
        #lab = ["1e-3", "1e-2",  "1", "1e2", "1e3"]
        # Add 1.0 to the list
        x_val = x_val**2 
        data = np.array(data)
        data = np.abs(data[np.argsort(data[:,-1])])
        mrarr = np.array([2,8,10,14,18])
        #mrarr = np.array([0,2,6,10,12])
        #initial_values = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
        print(data[:,2])
        for i in range(len( mrarr)):
            print(data[:, mrarr[i]])
            plt.plot(x_val, data[:, mrarr[i]],marker='o',linestyle= '--',label=r'$m_l/m_r = $'+ lab[i])


        plt.xlabel(r'$\epsilon^2$')
        plt.ylabel('D')
        plt.xscale('log')
        plt.yscale('log')

        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.title(r'$NP=10^5, N=500$')
        plt.legend()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        fname = "D_epssq_mr_dep.png"
        plot_file_path = os.path.join(script_dir, fname)

        plt.savefig(plot_file_path)
        plt.show()
        return None
    
    #data = read_data("./data/Dfum_01k_10k_Dm2_new.txt")
    
    #Desq(data)