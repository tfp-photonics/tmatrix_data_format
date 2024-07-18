from addatmatrix.data import write_t_matrix_data
from addatmatrix.t_matrix import t_matrix

output_path = "sphere.tmat.h5"
adda_string = "-shape sphere -m 2 0"
tmatrix_data, metadata, indices = t_matrix(adda_string,wavelength=500)
write_t_matrix_data(output_path, tmatrix_data, indices, metadata)