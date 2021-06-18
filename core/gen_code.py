#!/usr/bin/env python3

#%%
import os
from template import *
import shutil

def concatenate_files(src_file, dest_file):
  with open(dest_file, 'wb') as fd:
    print("Start combing file", flush=True)
    for fname in src_file:
      with open(fname, 'rb') as fs:
          shutil.copyfileobj(fs, fd)
    print("File {:s} is ready".format(dest_file), flush=True)


def gen_write_dataset(ftype_name, rank):
  ftype = ftype_name
  h5type = ''
  if ftype_name == 'integer':
    h5type = 'H5T_NATIVE_INTEGER'
  elif ftype_name == 'real':
    ftype = 'real(sp)'
    h5type = 'H5T_NATIVE_REAL'
  elif ftype_name == 'double':
    ftype = 'real(dp)'
    h5type = 'H5T_NATIVE_DOUBLE'
  elif ftype_name == 'complex_double':
    ftype = 'complex(dp)'
    h5type = 'complexd_type_id'
  elif ftype_name == 'character':
    h5type = 'dtype_id'
    ftype = 'character(len=*)'
  else:
    raise ValueError(f'Unkndown ftype {ftype_name}')

  array_comma = ''
  if rank > 0:
    array_comma = ':,' * (rank - 1) + ':'


  if ftype_name == 'complex_double':
    write_part = write_complex_dataset_template.format(h5type=h5type)
  else:
    write_part = write_regular_dataset_template.format(h5type=h5type)

  config = {
    'ftype': ftype,
    'ftype_name': ftype_name,
    'h5type': h5type,
    'rank': rank,
    'write_string': write_part,
    'array_comma': array_comma
  }

  if rank == 0:
    declaration = f'''{ftype}, intent(in) :: array'''
    declaration = ' '*4 + declaration
    declaration = '{:54s}! data to be written'.format(declaration)
    config['declaration'] = declaration

    if ftype_name == 'character':
      return write_char0_template.format(**config)
    else:
      return write_single_value_template.format(**config)
  else:
    declaration = f'    {ftype}, intent(in) :: array({array_comma})'
    declaration = '{:54s}! data to be written'.format(declaration)
    config['declaration'] = declaration

    if ftype_name == 'character':
      length_calculation = '    length=len(array({:s}))'.format('1,'*(rank-1)+'1')
      config['length_calculation'] = length_calculation
      return write_char_array_template.format(**config)
    else:
      return write_array_template.format(**config)

# %%

if __name__ == '__main__':
  with open('new.f90', 'w') as f:
    for ftype_name in ['integer', 'real', 'double', 'character', 'complex_double']:
      f.write(f'''
  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_{ftype_name}--------------------------------
  !!----------------------------------------------------------------------------------------

''')
      if ftype_name == 'character':
        max_dim = 2
      else:
        max_dim = 6
      for rank in range(max_dim+1):
        if rank == 0:
          f.write('''  !  \\brief write a scalar to a hdf5 file\n''')
        else:
          f.write(f'''  !  \\brief write a {rank} - dimension array to a hdf5 file\n''')
        f.write(gen_write_dataset(ftype_name, rank))
        f.write('\n')

  concatenate_files(['hdf5_utils_mpi_head.f90', 'new.f90', 'hdf5_utils_mpi_end.f90'], '../hdf5_utils_mpi.f90')
  os.remove('new.f90')


# %%
