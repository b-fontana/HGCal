import csv
import numpy as np

def dump_tcl_resolution_formula(name, bins, res, bias):
    """
    Creates a Delphes 'ResolutionFormula' in 'name' to be used within a Calorimeter module.
    Both the inner and outer regions are considered in the remianing arguments. A NaN value is used to indicate
        where one region starts and the other finishes: np.array(out1, out2, ..., NaN, in1, in2, ...)
    """
    assert(len(res)==4 and len(bias)==4)
    def find_nan_index(arr):
        v = np.where(np.isnan(arr))
        if v[0].shape != (1,):
            raise ValueError('The array has to contain one and only one NaN.')
        return v[0][0]

    for imask in range(len(bias)):
        with open(name+'_Mask'+str(imask+3)+'.tcl', 'w') as f:
            f.write('set ResolutionFormula {\n')
            f.write('    (abs(eta) > 1.65 && abs(eta) <= 1.75) * sqrt(energy^2*0.006^2 + energy*0.20^2) + \ \n')
            f.write('    (abs(eta) > 1.75 && abs(eta) <= 2.15) * sqrt(energy^2*0.007^2 + energy*0.21^2) + \ \n')
            f.write('    (abs(eta) > 2.15 && abs(eta) <= 2.70) * sqrt(energy^2*0.008^2 + energy*0.24^2) + \ \n')
            for i, (edge1,edge2) in enumerate(zip(bins[:-1],bins[1:])):
                if edge1==1.65 and edge2==2.7: #transition between outer and inner region is not a bin
                    continue
                
                nan_idx = find_nan_index(np.array(res[imask]))
                if nan_idx != find_nan_index(np.array(bias[imask])):
                    raise ValueError('The Nan value must always lie at the same position in the array.')

                sf_in = str(round(res[imask][i]/res[imask][nan_idx+1],3)-1) if res[imask][i]!=0. else "0."
                sf_out = str(round(res[imask][i]/res[imask][nan_idx-1],3)-1) if res[imask][i]!=0. else "0."
                f.write('    (abs(eta) > '+str(edge1)+' && abs(eta) < '+str(edge2)+') ')
                if edge1>=2.7 and edge1<3.:
                    f.write('* ( sqrt(energy^2*(0.008+'+sf_out+')^2 + energy*0.24^2) + '+str(bias[imask][i])+' )')
                elif edge1 >= 3. and edge1<5.:
                    f.write('* ( sqrt(energy^2*(0.010+'+sf_out+')^2 + energy*1.82^2) + '+str(bias[imask][i])+' )')
                elif edge2<=1.65 and edge2>1.5:
                    f.write('* ( sqrt(energy^2*(0.006+'+sf_in+')^2 + energy*0.20^2) + '+str(bias[imask][i])+' )')
                elif edge2<=1.5  and edge2>=1.45:
                    f.write('* ( sqrt(energy^2*(0.009+'+sf_in+')^2 + energy*0.12^2 + 0.45^2) ')
                    f.write('+ '+str(bias[imask][i])+' )')
                else:
                    print(edge1, ' ', edge2)
                    raise ValueError('This cannot be happening!')
                if edge2==bins[-1]:
                    f.write('\n')
                else:
                    f.write('+ \ \n')
            f.write("}")

def save_to_csv(name, bins, res, biases):
    assert(len(res)==4 and len(biases)==4)
    for r in res:
        assert(len(bins)==len(r)+1)
    for bias in biases:
        assert(len(bins)==len(bias)+1)
    with open(name, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(['Edge1', 'Edge2', 
                             'ResolutionMask3', 'ResolutionMask4', 'ResolutionMask5', 'ResolutionMask6',
                             'BiasMask3', 'BiasMask4', 'BiasMask5', 'BiasMask6'])
        for i, (edge1,edge2) in enumerate(zip(bins[:-1],bins[1:])):
            csv_writer.writerow([edge1, edge2, 
          round(res[0][i]/res[0][i],3), round(res[1][i]/res[0][i],3), round(res[2][i]/res[0][i],3),  round(res[3][i]/res[0][i],3),
          round(biases[0][i],3), round(biases[1][i],3), round(biases[2][i],3),  round(biases[3][i],3)])
