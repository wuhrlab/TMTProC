
import configparser
import collections


def match_peak(peaks, target_mz, tolerance, tolerance_units, strategy='nearest'):
    if peaks.size < 1:
        return None
    if strategy == 'nearest':
        window = make_range(target_mz, tolerance, tolerance_units)
        index = lower_bound(peaks, target_mz)
        e = abs(peaks[index][0] - target_mz)
        if index + 1 < peaks.size and abs(peaks[index + 1][0] - target_mz) < e:
            index += 1
        match = index
    elif strategy == 'max':
        window = make_range(target_mz, tolerance, tolerance_units)
        index = lower_bound(peaks, window.low)
        match = None
        max_intensity = 0
        while True:
            index += 1
            if index >= peaks.size:
                break
            if peaks[index][0] < window.low:
                continue
            if peaks[index][0] > window.high:
                break
            if peaks[index][1] > max_intensity:
                match = index
                max_intensity = peaks[index][1]
    else:
        raise ValueError('strategy must be nearest or max')
    
    if window.low < peaks[match][0] < window.high:
        return match
    return None


def lower_bound(peaks, target_mz):
    low = 0
    high = peaks.size
    c = 0
    while True:
        if low == high:
            break
        
        mid = int((high + low) * 0.5)

        if peaks[mid][0] < target_mz:
            low = mid + 1
        else:
            high = mid


    if low > 0 and low >= peaks.size:
        return peaks.size - 1

    if low > 0 and peaks[low][0] > target_mz:
        return low - 1
    return low


def mz_to_mh(mz, charge):
    return (float(charge) * (mz - 1.007276466812)) + 1.007276466812


def mh_to_mz(mh, charge):
    return ((mh - 1.007276466812) + (charge * 1.007276466812)) / charge


def ppm_to_error(mz, ppm):
    return mz * ppm / 1000000


def make_range(mz, tolerance, units):
    e = tolerance
    if units == "'ppm'":
        e = ppm_to_error(mz, tolerance)
    return collections.namedtuple('Range', ['low', 'high'])(mz - e, mz + e)


def calculate_isolation_specificity(scan, target_mz, charge, isolation_window, tolerance=0.04):
    peaks = scan['peaks']
    if peaks.size < 1:
        return 0

    window = make_range(target_mz, isolation_window, 'Da')
    index = lower_bound(peaks, window.low)
    total_intensity = 0
    target_intensity = 0
    while index < peaks.size:
        peak = peaks[index]
        index += 1
        if peak[0] < window.low:
            continue

        if peak[0] > window.high:
            break

        total_intensity += peak[1]
        for i in range(-5, 5):
            mz = target_mz + (i * (1.00286864 / charge))  # Averagine
            if abs(peak[0] - mz) < tolerance:
                target_intensity += peak[1]
                break

    if total_intensity:
        return target_intensity / total_intensity
    return 0


class ReporterIon:
    name = ''
    mz = 0
    charge = 0
    precursor = False
    intensity = 0
    noise = 0
    mz_error = 0


class PeakExtractor:
    reporters = []
    match_tolerance = 0.1
    match_tolerance_units = 'Da'

    def load_options(self, path):
        parser = configparser.RawConfigParser()
        parser.read(path)

        for section in parser.sections():
            if section == 'general':
                self.match_tolerance = self._get_option(parser, section, 'tolerance', 'float', 25.)
                self.match_tolerance_units = self._get_option(parser, section, 'tolerance_unit', 'str', 'ppm')
            else:
                reporter = ReporterIon()
                reporter.name = section
                reporter.mz = self._get_option(parser, section, 'mz', 'float', 0.0)
                reporter.precursor = self._get_option(parser, section, 'precursor', 'bool', False)
                reporter.charge = self._get_option(parser, section, 'charge', 'int', 0)
                self.reporters.append(reporter)

    def channel_names(self):
        output = []
        for reporter in self.reporters:
            output.append(reporter.name)
        return output

    def extract_peaks(self, scan, scan_precursor_mz=None):
        precursor_charge = 0
        precursor_mz = 0
        if 'precursor_mz' in scan:
            precursor_mz = scan['precursor_mz']
        if 'precursor_charge' in scan:
            precursor_charge = scan['precursor_charge']
        if scan_precursor_mz:
            # user has supplied alternate precursor m/z value.
            precursor_mz = scan_precursor_mz

        for reporter in self.reporters:
            reporter.mz_error = 0
            reporter.intensity = 0
            target_mz = reporter.mz
            if reporter.precursor:
                if precursor_charge == 0 or precursor_mz == 0:
                    raise ValueError('Quant relative to precursor but scan has no precursor info')
                m_h = mz_to_mh(precursor_mz, precursor_charge)
                target_charge = precursor_charge + reporter.charge
                target_mz = mh_to_mz(m_h + reporter.mz, target_charge)

            match = match_peak(scan['peaks'],
                               target_mz,
                               self.match_tolerance,
                               self.match_tolerance_units)

            if match is not None:
                peak = scan['peaks'][match]
                reporter.mz_error = abs(peak[0] - target_mz)
                reporter.intensity = peak[1]
                if 'noise' in scan:
                    reporter.noise = scan['noise'][match][2]

    @staticmethod
    def _get_option(parser, section, option, settype='str', default=None):
        """

        :param settype:
            :param parser: RawConfigParser
            :param section: str
            :param option: str
            :param default:
            :return:
            """
        if not parser.has_section(section) or not parser.has_option(section, option):
            return default

        if settype == 'int':
            return parser.getint(section, option)
        elif settype == 'float':
            return parser.getfloat(section, option)
        elif settype == 'bool':
            return parser.getboolean(section, option)
        elif settype == 'str':
            return parser.get(section, option)
        else:
            raise ValueError('Invalid value for type parameter')
