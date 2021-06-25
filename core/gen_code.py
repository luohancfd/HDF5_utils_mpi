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

#%%
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
  elif ftype_name == 'character':
    write_part_1 = write_regular_dataset_template.format(h5type=h5type, data_name='array')
    write_part_2 = write_regular_dataset_template.format(h5type=h5type, data_name='buffer')

    write_part = ['    if (length .eq. length_glob) then'] + ['  '+i for i in write_part_1.split('\n')]
    write_part += ['    else'] + ['  '+i for i in write_part_2.split('\n')] + ['    end if\n']
    write_part = '\n'.join(write_part)
  else:
    write_part = write_regular_dataset_template.format(h5type=h5type, data_name='array')

  config = {
    'ftype': ftype,
    'ftype_name': ftype_name,
    'h5type': h5type,
    'rank': rank,
    'write_string': write_part,
    'array_comma': array_comma,
    'additonal_delcaration': '',
    'dtype_creation': '',
    'additional_close': '',
    'additional_attribute': '',
    'data_name': 'array'
  }

  nspace = 4
  if rank == 0:
    declaration = f'''{ftype}, intent(in) :: array'''
    declaration = ' '*nspace + declaration
    declaration = '{:54s}! data to be written'.format(declaration)
    config['declaration'] = declaration
  else:
    declaration = f'    {ftype}, intent(in) :: array({array_comma})'
    declaration = '{:54s}! data to be written'.format(declaration)
    config['declaration'] = declaration

  dim_string = ','.join([f'dimsm({i})' for i in range(1, rank+1)])

  if ftype_name == 'character':
    config['additonal_delcaration'] = [
      'integer(HSIZE_T) :: length, length_glob',
      'integer(HID_T) :: dtype_id'
    ]

    if rank > 0:
      config['additonal_delcaration'].append('integer :: ' + ', '.join([f'i_{i}' for i in range(1, rank+1)]))

    if rank == 0:
      config['additonal_delcaration'].append(f'character(len=:),allocatable :: buffer')
      buffer_allocation = 'allocate(character(len=length_glob)::buffer)'
    else:
      config['additonal_delcaration'].append(f'character(len=:), dimension({array_comma}), allocatable :: buffer')
      buffer_allocation = f'allocate(character(len=length_glob)::buffer({dim_string}))'

    config['additonal_delcaration'] = '\n'.join([' '*nspace + i for i in config['additonal_delcaration']]) + '\n'

    cps = []
    if rank == 0:
      length_calculation = 'length=len(array)'
      cps = ['buffer = array']
    else:
      length_calculation = 'length=len(array({:s}))'.format('1,'*(rank-1)+'1')
      temp_n = 0
      for i in range(1, rank+1):
        cps.append(temp_n*' ' + f'do i_{i} = 1, dimsm({i})')
        temp_n += 2
      index_string = ','.join([f'i_{i}' for i in range(1, rank+1)])
      cps.append(' '*(temp_n) + f'buffer({index_string}) = array({index_string})')
      temp_n -= 2
      for i in range(1, rank+1):
        cps.append(temp_n*' ' + f'end do')
        temp_n -= 2
    cps = ['    ' + i for i in cps]


    config['dtype_creation'] = [
      length_calculation,
      'length_glob = length',
      'if (processor_write .ne. -1) then',
      '  call MPI_Bcast(length, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)',
      'else',
      '  call MPI_Allreduce(length, length_glob, 1, mpi_hsize_t, MPI_MAX, mpi_comm, mpi_ierr)',
      '  if (length .ne. length_glob) then',
      f'    {buffer_allocation}'] + cps + [
      '  end if',
      'end if',
      'call h5tcopy_f(H5T_FORTRAN_S1, dtype_id, hdferror)',
      'call h5tset_size_f(dtype_id, length_glob, hdferror)'
    ]
    config['dtype_creation'] = '\n'.join([' '*nspace + i for i in config['dtype_creation']]) + '\n'


    config['additional_attribute'] = [
      "call hdf_write_attribute(dset_id, '', 'char_length', int(length_glob, kind=4))"
    ]
    config['additional_attribute'] = '\n'.join([' '*nspace + i for i in config['additional_attribute']]) + '\n'

    config['additional_close'] = [
      'call h5tclose_f(dtype_id, hdferror)',
      'if (length .ne. length_glob .and. processor_write .eq. -1) then',
      '  deallocate(buffer)',
      'end if'
    ]
    config['additional_close'] = '\n'.join([' '*nspace + i for i in config['additional_close']])+'\n'

  if rank == 0:
    return write_single_value_template.format(**config)
  else:
    return write_array_template.format(**config)

#%%
def gen_read_dataset(ftype_name, rank, nspace=4):
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
  buffer_indexing = ''
  if rank > 0:
    array_comma = ':,' * (rank - 1) + ':'
    buffer_indexing = array_comma + ','

  buffer_size = ', '.join([f'dims({i})' for i in range(1, rank+1)])

  if rank == 0:
    declaration = f'{ftype}, intent(out) :: array'
  else:
    declaration = f'{ftype}, intent(out) :: array({array_comma})'
  declaration = ' '*nspace + declaration
  declaration = '{:48s}! data to be read'.format(declaration)

  additional_declaration = ''
  if ftype_name == 'complex_double':
    if rank == 0:
      additional_declaration = ' '*nspace + 'real(dp) :: buffer(2)                       ! buffer to save real and imag part'
    else:
      additional_declaration = ' '*nspace + f'real(dp), allocatable, dimension({array_comma},:) :: buffer ! buffer to save real and imag part'

  set_dims = ' '*nspace + 'dims = shape(array, KIND=HSIZE_T)'
  dims_declaration = ' '*nspace + f'integer(HSIZE_T) :: dims({rank})'
  if rank == 0:
    set_dims = ' '*nspace + 'dims = (/0/)'
    dims_declaration = ' '*nspace + f'integer(HSIZE_T) :: dims(1)'

  if ftype_name == 'complex_double':
    read_string = read_complex_dataset_template.format(buffer_indexing=buffer_indexing)
    if rank != 0:
      read_string = read_string.split('\n')
      read_string = [f'    allocate (buffer({buffer_size}, 2))'] + read_string + [
        '    deallocate (buffer)'
      ]
      read_string = '\n'.join(read_string)
  else:
    read_string = read_regular_dataset_template.format(h5type=h5type)

  config = {
    'ftype_name': ftype_name,
    'declaration': declaration,
    'additional_declaration': additional_declaration,
    'rank': rank,
    'set_dims': set_dims,
    'read_string': read_string,
    'dims_declaration': dims_declaration,
  }

  return read_array_template.format(**config)


# %%

if __name__ == '__main__':
  with open('new.f90', 'w') as f:

    # ====================== write dataset ========================================
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
    # ====================== read dataset ========================================
    for ftype_name in ['integer', 'real', 'double', 'complex_double']:
      max_dim = 6
      f.write(f'''
  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_{ftype_name}--------------------------------
  !!----------------------------------------------------------------------------------------

''')
      for rank in range(max_dim+1):
        if rank == 0:
          f.write('''  !  \\brief read a scalar from a hdf5 file\n''')
        else:
          f.write(f'''  !  \\brief read a {rank} - dimension array from a hdf5 file\n''')
        f.write(gen_read_dataset(ftype_name, rank))
        f.write('\n')

  concatenate_files([
    'hdf5_utils_mpi_head.f90',
    'new.f90',
    'hdf5_utils_mpi_read_character.f90',
    'hdf5_utils_mpi_end.f90'], '../hdf5_utils_mpi.f90')
  os.remove('new.f90')
# %%
