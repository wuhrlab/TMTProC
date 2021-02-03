
import xml.sax
import os
import re
import numpy
import base64


def _decode_peaks(precision, data):
    size = numpy.float32
    if precision == '64':
            size = numpy.float64
    dtype = numpy.dtype([('mz', size), ('intensity', size)]).newbyteorder('>')
    decoded = base64.b64decode(data.encode('ascii'))
    return numpy.frombuffer(decoded, dtype=dtype)


def _decode_label_data(precision, data):
    size = numpy.float32
    if precision == '64':
            size = numpy.float64
    dtype = numpy.dtype([('mz', size), ('baseline', size), ('noise', size)]).newbyteorder('>')
    decoded = base64.b64decode(data.encode('ascii'))
    return numpy.frombuffer(decoded, dtype=dtype)

class MzXMLStreamHandler(xml.sax.handler.ContentHandler):
    index_offset = None
    index = {}
    last_entry = {}
    scan = {}
    last_name = None

    def startElement(self, name, attrs):
        self.last_name = name
        if name == 'scan':
            self.scan = dict(zip(attrs.keys(), attrs.values()))
        if name != 'root':
            self.last_entry[name] = {'attrs': attrs, 'content': ''}

    def endElement(self, name):
        if name == 'indexOffset':
            self.index_offset = int(self.last_entry['indexOffset']['content'])
        elif name == 'offset':
            offset_id = int(self.last_entry['offset']['attrs'].getValue('id'))
            content = int(self.last_entry['offset']['content'])
            self.index[offset_id] = content
        elif name == 'peaks':
            compression = self.last_entry['peaks']['attrs'].getValue('compressionType')
            if compression.lower() != 'none':
                raise Exception('Compressed data not supported')
            precision = self.last_entry['peaks']['attrs'].getValue('precision')
            self.scan['peaks'] = _decode_peaks(precision, self.last_entry['peaks']['content'])
        elif name == 'labelData':
            compression = self.last_entry['labelData']['attrs'].getValue('compressionType')
            if compression.lower() != 'none':
                raise Exception('Compressed data not supported')
            precision = self.last_entry['labelData']['attrs'].getValue('precision')
            self.scan['noise'] = _decode_label_data(precision, self.last_entry['labelData']['content'])
        elif name == 'precursorMz':
            precursor_mz = self.last_entry['precursorMz']
            self.scan['precursor_charge'] = int(precursor_mz['attrs'].getValue('precursorCharge'))
            self.scan['precursor_mz'] = float(precursor_mz['content'])
            self.scan['precursor_intensity'] = float(precursor_mz['attrs'].getValue('precursorIntensity'))
            self.scan['precursor_scan_num'] = float(precursor_mz['attrs'].getValue('precursorScanNum'))
        elif name == 'root':
            raise StopIteration

    def characters(self, content):
        if self.last_entry:
            self.last_entry[self.last_name]['content'] += content


class MzXML:
    def __init__(self):
        self.handle = None
        self.index = None
        self.parser = xml.sax.make_parser()
        self.size = 0

    def open(self, path):
        self.size = os.stat(path).st_size
        self.handle = open(path)
        self._load_index()

    def close(self):
        if self.handle:
            self.handle.close()

    def get_scan(self, number):
        number = int(number)
        if number not in self.index:
            raise ValueError('Could not find scan in MzXML file: ' + str(number))
        data = MzXMLStreamHandler()
        self.parser.setContentHandler(data)
        self.handle.seek(self.index[number])
        started = False
        buffer = ''
        for line in self.handle:
            if '<scan' in line:
                if started:
                    # Already reached another scan, but we didnt find the end tag
                    # This can happen when ms1 scans nest ms2 scans
                    buffer += '</scan>'
                    break
                started = True
            buffer = buffer + line
            if '</scan>' in line:
                break
        if not buffer:
            return None
        xml.sax.parseString(buffer, data)
        return data.scan

    def _load_index(self):
        index_offset = self._find_index_offset()
        if index_offset is None:
            raise ValueError('Could not read index from mzxml file')
        data = MzXMLStreamHandler()
        self.parser.setContentHandler(data)
        self.handle.seek(index_offset)
        try:
            self.parser.parse(self.handle)
        except xml.sax.SAXParseException:
            # reached EOF with mismatched tags
            pass

        # parser may have tried to close the file.
        if self.handle.closed:
            self.handle = open(self.handle.name)
            
        self.index = data.index

    def _find_index_offset(self):
        self.handle.seek(self.size-200)
        for line in self.handle:
            match = re.search(r'<indexOffset>(\d+)<', line)
            if match:
                return int(match.group(1))
        return None
