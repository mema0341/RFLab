'''
Copyright 2019 Penina Axelrad, Ryan Kain

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from datetime import datetime
import re


blank_expr = re.compile('([0-9]+)X')
value_expr = re.compile('([0-9]+)?([SIFD]{1}[0-9]+)(\.[0-9]+)?')


def parse_rinex_line(line, line_format):
    '''
    Given the RINEX 2 navigation message line and
    line format string, parses and returns list of
    data in that line

    Parameters
    ----------
    line : str
        Line to parse into nav data
        
    line_format: array_like
        Line format used to describe data contents

    Returns
    -------
    list
        Vector of the form [header, nav_data] where header is a list
        containing the parsed header information and nav_data is a
        dictionary containing a list of the navigation records for each
        sat_id.
    '''
    
    line_components = line_format.split(',')
    values = []
    
    # find the letter format specifier
    for expr in line_components:
        blank_match = blank_expr.match(expr)
        if blank_match is not None:
            num_blank = int(blank_match.groups()[0])
            line = line[num_blank:]
            continue
        value_match = value_expr.match(expr)
        if value_match is not None:
            groups = value_match.groups()
            multiple = int(groups[0]) if groups[0] is not None else 1
            val_type = groups[1][0]
            length = int(groups[1][1:])
            precision = int(groups[2][1:]) if groups[2] is not None else 0
            for i in range(multiple):
                x = line[:length].strip()
                if x == '':
                    value = None
                elif val_type is 'S':
                    value = x
                elif val_type is 'I':
                    value = int(x)
                elif val_type is 'F':
                    value = float(x)
                elif val_type is 'D':
                    value = float(x.replace('D', 'E'))
                values.append(value)
                line = line[length:]
    return values

rinex2_nav_record_line_formats = [
    'I2,1X,I2.2,1X,I2,1X,I2,1X,I2,1X,I2,F5.1,3D19.12',
    '3X,4D19.12',
    '3X,4D19.12',
    '3X,4D19.12',
    '3X,4D19.12',
    '3X,4D19.12',
    '3X,4D19.12',
    '3X,4D19.12',
]

rinex2_nav_record_var_names = [
    ['prn', 'yy', 'month', 'day', 'hour', 'minute', 'second', 'a0', 'a1', 'a2'],
    ['iode1', 'c_rs', 'delta_n', 'm_0',],
    ['c_uc', 'e', 'c_us', 'sqrt_a',],
    ['t_oe', 'c_ic', 'omega_0', 'c_is',],
    ['i_0', 'c_rc', 'omega', 'omega_dot',],
    ['i_dot', 'l2_codes', 'week', 'l2p_data',],
    ['accuracy', 'health', 'tgd', 'iodc',],
    ['transmit_time', 'fit_interval'],
]

def parse_RINEX2_nav_records(lines):
    '''
    Given the lines corresponding to navigation data records 
    from a RINEX 2 Nav file, parses the records and returns 
    a list of dictionaries with their contents.

    Parameters
    ----------
    lines : array_like
        List of strings, where each string is a line of nav data

    Returns
    -------
    array_like
        List of records containing data as passed in from `lines`
    '''

    # Assume lines contain complete records
    i = 0
    records = []
    record_length = len(rinex2_nav_record_line_formats)
    while i + record_length <= len(lines):
        record_lines = lines[i:i + record_length]
        record = {}
        for line, line_format, var_names in zip(record_lines, rinex2_nav_record_line_formats, rinex2_nav_record_var_names):
            values = parse_rinex_line(line, line_format)
            for key, val in zip(var_names, values):
                if key:
                    record[key] = val
        records.append(record)
        i += record_length
    return records


def format_RINEX2_nav_records(records, century=2000):
    '''
    Take list of nav file records and sort into a dictionary with 
    PRN as the key and a list of records sorted by epoch as the value

    Parameters
    ----------
    records : array_like
        List of strings, where each string is a line of nav data
    century : int
        Century during which the two-digit year should be interpreted.
        Default 2000.

    Returns
    -------
    array_like
        Dictionary of records containing data as passed in from `records`,
        keyed by PRN.
    '''
    
    ephemerides = {}
    
    # Look through each record and assign it to the respective PRN
    for record in records:
        sat_id = 'G{0:02}'.format(record['prn'])
        ephemeris_list = ephemerides.get(sat_id, [])
        eph = record.copy()
        # Add the epoch time as a datetime object
        eph['epoch'] = datetime(century + record['yy'], *(int(record[k]) for k in ['month', 'day', 'hour', 'minute', 'second']))
        ephemeris_list.append(eph)
        ephemerides[sat_id] = ephemeris_list

    # Sort the data by epoch for each PRN
    for key, ephemeris_list in ephemerides.items():
        ephemerides[sat_id] = sorted(ephemeris_list, key=lambda x: x['epoch'])

    return ephemerides


def parse_rinex(filepath, return_header=True):
    '''
    Given the filepath to a RINEX 2 navigation message file,
    parses and returns header and navigation ephemeris data.
    Automatically removes obsoleted (corrected) data.

    Parameters
    ----------
    filepath : str, file_like
        filepath to or open file object of RINEX 2 navigation file

    Returns
    -------
    array_like
        Vector of the form [header, nav_data] where header is a list
        containing the parsed header information and nav_data is a
        dictionary containing a list of the navigation records for each
        sat_id.
    '''
    
    try:
        with open(filepath, 'r') as f:
            lines = list(f.readlines())
    except TypeError:
        lines = list(filepath.readlines())
    
    # Parse out the header
    for i, line in enumerate(lines):
        if line.find('END OF HEADER') >= 0:
            break
    header_lines = lines[:i + 1]
    
    # Parse through the nav lines
    nav_lines = lines[i + 1:]
    records = parse_RINEX2_nav_records(nav_lines)
    ephemerides = format_RINEX2_nav_records(records)
    
    # Clean up ephemerides
    for prn in ephemerides:
        for ephem in ephemerides[prn]:
            del_index = []
            if ephem['t_oe'] % 3600 != 0:
                # Ephemeris was corrected!
                # Limit to 4 minutes away and remove any updates
                # that happened "later" in time
                for ii, e in enumerate(ephemerides[prn]):
                    toe_diff = e['t_oe'] - ephem['t_oe']
                    # If update is from less than 4 minutes in the future,
                    # it should be superceded by this update
                    if toe_diff < 4*60 and toe_diff > 0:
                        # Keep track of each item to remove
                        del_index.append(ii)
            # Remove all obsoleted data
            for ii in del_index:
                del ephemerides[prn][ii]
    
    if return_header:
        return header_lines, ephemerides
    else:
        return ephemerides
