import os
import numpy


data_type_key = {
    'MET_CHAR': numpy.byte,
    'MET_UCHAR': numpy.ubyte,
    'MET_SHORT': numpy.short,
    'MET_USHORT': numpy.ushort,
    'MET_INT': numpy.intc,
    'MET_UINT': numpy.uintc,
    'MET_FLOAT': numpy.single,
    'MET_DOUBLE': numpy.double,
    numpy.byte: 'MET_CHAR',
    numpy.ubyte: 'MET_UCHAR',
    numpy.short: 'MET_SHORT',
    numpy.ushort: 'MET_USHORT',
    numpy.intc: 'MET_INT',
    numpy.uintc: 'MET_UINT',
    numpy.single: 'MET_FLOAT',
    numpy.double: 'MET_DOUBLE'}

# The order for these is important
accepted_tags = ('ObjectType','NDims','BinaryData','BinaryDataByteOrderMSB','CompressedData','CompressedDataSize',
                 'TransformMatrix','Offset','CenterOfRotation','AnatomicalOrientation','ElementSpacing','DimSize',
                 'ElementType','ElementDataFile','Comment','SeriesDescription','AcquisitionDate','AcquisitionTime',
                 'StudyDate','StudyTime')


def load_mhd_header(filename):
    """ Return a dictionary of meta data from meta header file """
    header_file = open(filename, "r")

    meta_dict = {}
    for line in header_file:
        tag, value = line.split(' = ')
        if tag in accepted_tags:
            meta_dict[tag.strip()] = value.strip()
        else:
            print 'Encountered unexpected tag: ' + tag
    header_file.close()

    return meta_dict


def load_mhd(filename):
    meta_dict = load_mhd_header(filename)
    data_dimensions = map(int, meta_dict['DimSize'].split())

    image_dir = os.path.dirname(filename)
    data_filepath = os.path.join(image_dir, meta_dict['ElementDataFile'])
    data_type = data_type_key[meta_dict['ElementType'].upper()]

    image_data = numpy.fromfile(data_filepath, dtype=data_type)
    image_data = numpy.reshape(image_data, data_dimensions)

    return image_data, meta_dict


def write_meta_header(filename, meta_dict):
    header = ''
    # Tag order matters here so I can't just iterate through meta_dict.keys()
    for tag in accepted_tags:
        if tag in meta_dict.keys():
            header += '%s = %s\n'%(tag,meta_dict[tag])

    f = open(filename,'w')
    f.write(header)
    f.close()


def dump_raw_data(filename, data):
    # TODO: THIS
    """ Write the data into a raw format file. Big endian is always used. """
    rawfile = open(filename,'wb')
    a = array.array('f')
    for o in data:
        a.fromlist(list(o))
    #if is_little_endian():
    #    a.byteswap()
    a.tofile(rawfile)
    rawfile.close()


def write_mhd(output_path, image_data, **kwargs):
     metadata = {'ObjectType': 'Image',
                 'BinaryData': 'True',
                 'BinaryDataByteOrderMSB': 'False',
                 'ElementType': data_type_key[image_data.dtype],
                 'NDims': str(image_data.ndim),
                 'DimSize': ' '.join(image_data.shape),
                 'ElementDataFile': os.path.basename(output_path).replace('.mhd','.raw')}

     write_meta_header(output_path, metadata)
     output_dir = os.path.dirname((output_path))
     data_filepath = os.path.join(output_dir, metadata['ElementDataFile'])
     dump_raw_data(data_filepath, image_data)


