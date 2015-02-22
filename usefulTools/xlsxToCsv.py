#!/usr/bin/python

import xlrd
import csv
import argparse


def generate_parser():
    parser = argparse.ArgumentParser(
        description='Concatenates the sheets from the specified xlsx workbook into a csv file')
    parser.add_argument('-x', '--xlsxPath', help='The path to the xlsx file', required=True)
    parser.add_argument('-c', '--csvPath', help='The path to the csv output file', required=True)
    parser.add_argument('-n', '--numHeaders', type=int, default=0,
                        help='The number of header rows on each sheet of the xlsx file')
    return parser


def get_row_iter(sheet):
    for i in xrange(sheet.nrows):
        yield sheet.row(i)


def skip_rows(irows, n):
    while irows and n > 0:
        irows.next()
        n -= 1


def unpack_cells(cells):
    return [cell.value for cell in cells]


def main():
    parser = generate_parser()
    args = parser.parse_args()
    csv_file = open(args.csvPath, 'w')
    csv_writer = csv.writer(csv_file)

    book = xlrd.open_workbook(args.xlsxPath)
    for sheet in book.sheets():
        irows = get_row_iter(sheet)
        skip_rows(irows, args.numHeaders)
        for cells in irows:
            csv_writer.writerow(unpack_cells(cells))

if __name__ == '__main__':
    main()