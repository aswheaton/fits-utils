import numpy as np

def convert_spectral_type(data):
    """
    Conversion using Jordi et al (2006) equations.
    """
    color_u_g = (0.75 * data[:,1]) + (0.77 * data[:,2]) + 0.72
    color_g_r = (1.646 * data[:,3]) - 0.139
    new_data = [data[:,0], color_u_g, color_g_r]
    new_data = np.column_stack(new_data)
    np.savetxt('SDSS Calibration/spectral_ref.txt', new_data)
    return new_data

def main():
    #: ndarray: Data inside spectral type file.
    spectral_type_data = np.loadtxt('SDSS Calibration/spectral_ref.txt')
    dereddened_sources = np.loadtxt('cat/dr_gr.cat')
    spectral_types = []
    for source in dereddened_sources:
        index = np.where(spectral_type_data[:,2] == source[3]-source[5])[0]
        if np.size(index) != 0:
            spectral_types.append(spectral_type_data[index[0], 0])
        else:
            spectral_types.append(np.nan)
    dereddened_sources = np.concatenate(dereddened_sources, spectral_types,
                                        axis=1)
    np.savetxt('stars_with_type', dereddened_sources)
if __name__ == '__main__':
    main()
