import os
import numpy
import array
from itertools import imap

data_type_key = {
    'met_char': numpy.byte,
    'met_uchar': numpy.ubyte,
    'met_short': numpy.short,
    'met_ushort': numpy.ushort,
    'met_int': numpy.intc,
    'met_uint': numpy.uintc,
    'met_float': numpy.single,
    'met_double': numpy.double}


def read_meta_header(filename):
    """ Return a dictionary of meta data from meta header file """
    header_file = open(filename, "r")

    meta_dict = {}
    accepted_tags = ('ObjectType', 'NDims', 'DimSize', 'ElementType', 'ElementDataFile',
                     'BinaryData', 'BinaryDataByteOrderMSB', 'CompressedData', 'CompressedDataSize',
                     'Offset', 'CenterOfRotation', 'AnatomicalOrientation', 'ElementSpacing', 'TransformMatrix',
                     'Comment', 'SeriesDescription', 'AcquisitionDate', 'AcquisitionTime', 'StudyDate', 'StudyTime')

    for line in header_file:
        tag, value = line.split(' = ')
        if tag in accepted_tags:
            meta_dict[tag.strip()] = value.strip()
        else:
            print 'Encountered unexpected tag: ' + tag
    header_file.close()

    return meta_dict


def load_raw_data_with_mhd(filename):
    meta_dict = read_meta_header(filename)
    data_dimensions = tuple(imap(int, meta_dict['DimSize'].split()))

    image_dir = os.path.dirname(filename)
    data_filepath = os.path.join(image_dir, meta_dict['ElementDataFile'])
    print data_filepath
    data_type = data_type_key[meta_dict['ElementType'].lower()]

    image_data = numpy.fromfile(data_filepath, dtype=data_type)
    image_data = numpy.reshape(image_data, data_dimensions)

    return image_data, meta_dict


def write_meta_header(filename, meta_dict):
    header = ''
    # do not use tags = meta_dict.keys() because the order of tags matters
    tags = ['ObjectType','NDims','BinaryData',
       'BinaryDataByteOrderMSB','CompressedData','CompressedDataSize',
       'TransformMatrix','Offset','CenterOfRotation',
       'AnatomicalOrientation',
       'ElementSpacing',
       'DimSize',
       'ElementType',
       'ElementDataFile',
       'Comment','SeriesDescription','AcquisitionDate','AcquisitionTime','StudyDate','StudyTime']
    for tag in tags:
        if tag in meta_dict.keys():
            header += '%s = %s\n'%(tag,meta_dict[tag])
    f = open(filename,'w')
    f.write(header)
    f.close()


def dump_raw_data(filename, data):
    """ Write the data into a raw format file. Big endian is always used. """
    rawfile = open(filename,'wb')
    a = array.array('f')
    for o in data:
        a.fromlist(list(o))
    #if is_little_endian():
    #    a.byteswap()
    a.tofile(rawfile)
    rawfile.close()


def write_mhd_file(mhdfile, data, dsize):
    assert(mhdfile[-4:]=='.mhd')
    meta_dict = {}
    meta_dict['ObjectType'] = 'Image'
    meta_dict['BinaryData'] = 'True'
    meta_dict['BinaryDataByteOrderMSB'] = 'False'
    meta_dict['ElementType'] = 'MET_FLOAT'
    meta_dict['NDims'] = str(len(dsize))
    meta_dict['DimSize'] = ' '.join([str(i) for i in dsize])
    meta_dict['ElementDataFile'] = os.path.split(mhdfile)[1].replace('.mhd','.raw')
    write_meta_header(mhdfile, meta_dict)

    pwd = os.path.split(mhdfile)[0]
    if pwd:
        data_file = pwd +'/' + meta_dict['ElementDataFile']
    else:
        data_file = meta_dict['ElementDataFile']

    dump_raw_data(data_file, data)


