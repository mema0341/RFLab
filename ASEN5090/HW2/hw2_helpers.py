import numpy as np
import matplotlib.pyplot as plt

def get_ca(prn, N, shift=0, convert=False):
    # Step 1: Load G1 and G2 with all ones
    shift_reg1 = np.ones(10)
    shift_reg2 = np.ones(10)

    # Step 2: Compute sums
    tapped1 = [2, 9]
    tapped2 = [1, 2, 5, 7, 8, 9]

    ca_code = ""
    ca_code_plot = np.ndarray((N))
    output_plot = np.ndarray((N))

    for i in range(N):
        # Grab Output
        output1 = shift_reg1[len(shift_reg1)-1]
        output2 = shift_reg2[len(shift_reg2)-1]

        # Grab Input
        input1 = np.mod(np.sum(shift_reg1[tapped1]),2)
        input2 = np.mod(np.sum(shift_reg2[tapped2]),2)

        g2i_output = np.mod(np.sum(shift_reg2[prn]),2)

        ca_code += str(int(np.mod(np.sum([output1, g2i_output]),2)))
        ca_code_plot[i] = np.mod(np.sum([output1, g2i_output]),2)

        # Shift
        shift_reg1 = np.roll(shift_reg1,1)
        shift_reg2 = np.roll(shift_reg2,1)

        shift_reg1[0] = input1
        shift_reg2[0] = input2

    # Convert 0 --> 1 and 1 --> -1
    if convert:
        idx0 = np.where(ca_code_plot==0)
        idx1 = np.where(ca_code_plot==1)

        ca_code_plot[idx0] = 1
        ca_code_plot[idx1] = -1
        
    
    return np.roll(ca_code_plot, shift)

def plot16(ca_code, N, prn_val, show_plot=True):
    string_ca_code = ""
    for x in ca_code:
        string_ca_code = string_ca_code + str(int(x))

    fig, (ax1, ax2) = plt.subplots(1,2)

    fig.set_figwidth(fig.get_figwidth()*2)

    # Plot First 16 chips
    ax1.plot(np.arange(1,17,1), ca_code[0:16])
    ax1.set_title('1023-chip C/A-code PRN%s, Epochs 1024-2046\nFirst 16 Chips: %s' % (str(prn_val), str(hex(int(string_ca_code[0:16],2)))))
    ax1.set_xlabel('Chips')
    ax1.set_ylabel('C/A-code PRN%s' % str(prn_val))
    ax1.set_xlim(1,16)
    ax1.set_ylim(-2, 2)
    ax1.grid(True)

    # Plot Last 16 Chips
    ax2.plot(np.arange(1,17,1), ca_code[int(N/2-16-1):int(N/2-1)])
    ax2.set_title('1023-chip C/A-code PRN%s, Epochs 1024-2046\nLast 16 Chips: %s' % (str(prn_val), str(hex(int(string_ca_code[int(N/2-16-1):int(N/2-1)],2)))))
    ax2.set_xlabel('Chips')
    ax2.set_ylabel('C/A-code PRN%s' % str(prn_val))
    ax2.set_xlim(1,16)
    ax2.set_ylim(-2, 2)
    ax2.grid(True)

    if show_plot:
        plt.show()

def get_correlation(x_k, x_i, N):
    R_k = np.zeros(np.shape(x_k))
    # Cross-correlation
    for n in range(1022):
        s = []
        for i in range(1022):
            # If we hit the end of the array we need to loop back to the beginning
            try:
                s.append(x_k[i] * x_i[i+n])
            except: 
                s.append(x_k[i] * x_i[i+n-N])
        R_k[n] = np.sum(s) / 1023

    return R_k

def plot_correlation(R_k, title, show_plot=True):
    # Plot Auto Correlation
    fig, ax = plt.subplots()
    ax.plot(R_k)
    ax.set_title(title)
    ax.set_xlabel('Chip Delay')
    ax.set_ylabel('Correlation')
    ax.set_ylim(-0.2, 1)
    ax.grid(True)

    if show_plot:
        plt.show()

def plot_partf(x1, x2, x3, noise, N, show_plot=True):
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1)

    # fig.set_figwidth(fig.get_figwidth()*2)

    # Plot First 16 chips
    ax1.plot(np.arange(1,N+1,1), x1, color='C0')
    ax1.set_title('x1')
    ax1.set_xlabel('Chip Delay')
    ax1.set_ylabel('x1')
    ax1.set_xlim(1,N)
    ax1.set_ylim(-5, 5)
    ax1.grid(True)

    # Plot Last 16 Chips
    ax2.plot(np.arange(1,N+1,1), x2, color='C1')
    ax2.set_title('x2')
    ax2.set_xlabel('Chip Delay')
    ax2.set_ylabel('x2')
    ax2.set_xlim(1,N)
    ax2.set_ylim(-5, 5)
    ax2.grid(True)

    ax3.plot(np.arange(1,N+1,1), x3, color='C2')
    ax3.set_title('x3')
    ax3.set_xlabel('Chip Delay')
    ax3.set_ylabel('x3')
    ax3.set_xlim(1,N)
    ax3.set_ylim(-5, 5)
    ax3.grid(True)

    ax4.plot(np.arange(1,N+1,1), noise, color='C3')
    ax4.set_title('Noise')
    ax4.set_xlabel('Chip Delay')
    ax4.set_ylabel('Noise')
    ax4.set_xlim(1,N)
    ax4.set_ylim(-5, 5)
    ax4.grid(True)

    plt.tight_layout()

    if show_plot:
        plt.show()
