import numpy as np 

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
    g1 = np.ndarray((N))

    for i in range(N):
        # Grab Output
        output1 = shift_reg1[len(shift_reg1)-1]
        output2 = shift_reg2[len(shift_reg2)-1]

        g1[i] = output1

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

        # Convert g1
        idx0 = np.where(g1==0)
        idx1 = np.where(g1==1)

        g1[idx0] = 1
        g1[idx1] = -1        
        
    
    return np.roll(ca_code_plot, shift), g1