import re


transition_metal_identifiers = [
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'
]

class TddftDataParser():

    """Parser class for TD-DFT output files."""

    def __init__(self, file_path):

        self._id = file_path.split('/')[-1].split('.')[0]

        with open(file_path, 'r') as fh:
            self.lines = fh.read().strip().split('\n')

    def parse(self):

        """Parses the given TD-DFT output file.
        """

        return_dict = {'id': self._id}

        if self._has_failed():
            return_dict['has_failed'] = True
            return return_dict
        else:
            return_dict['has_failed'] = False

        return_dict['homo_lumo_gap'] = self._parse_homo_lumo_gap()
        return_dict['dipole_moment'] = self._parse_dipole_moment()
        return_dict['metal_charge'] = self._parse_metal_charge()

        spec = self._parse_spectrum()
        for i in range(len(spec[0])):
            return_dict['lambda_' + str(i+1)] = spec[0][i]
            return_dict['f_' + str(i+1)] = spec[1][i]

        spec_uv, spec_vis, spec_nir = self._split_spectrum_into_uv_vis_nir(spec)
        lambda_max_uv, f_max_uv, sigma_max_uv = self._get_lambda_max(spec_uv)
        lambda_max_vis, f_max_vis, sigma_max_vis = self._get_lambda_max(spec_vis)
        lambda_max_nir, f_max_nir, sigma_max_nir = self._get_lambda_max(spec_nir)

        return_dict['lambda_max_uv'] = lambda_max_uv
        return_dict['f_max_uv'] = f_max_uv
        return_dict['sigma_uv'] = sigma_max_uv

        return_dict['lambda_max_vis'] = lambda_max_vis
        return_dict['f_max_vis'] = f_max_vis
        return_dict['sigma_vis'] = sigma_max_vis

        return_dict['lambda_max_nir'] = lambda_max_nir
        return_dict['f_max_nir'] = f_max_nir
        return_dict['sigma_nir'] = sigma_max_nir

        return return_dict

    def _has_failed(self):

        """Checks if the job has failed.

        Returns:
            bool: The flag indicating whether the job has failed.
        """

        if 'Normal termination' not in self.lines[-1]:
            return True

        return False

    def _parse_spectrum(self):

        """Parses the spectrum from the output file.

        Returns:
            list[list[float]]: The parsed spectrum.
        """

        nms = []
        os = []
        for line in self.lines:
            if 'Excited State' in line:
                nms.append(float(line.split()[6]))
                os.append(float(line.split()[8].replace('f=', '')))

        spec = [nms, os]

        return spec

    def _parse_metal_charge(self):

        """ Parses the metal charge the output file.

        Returns:
            float: The metal charge.
        """

        for i, line in enumerate(self.lines):

            if 'Mulliken charges:' in line:

                j = i + 1
                while len(self.lines[j+1].split()) == 3:

                    j += 1
                    line_split = self.lines[j].split()
                    if line_split[1] in transition_metal_identifiers:
                        return float(line_split[2])

        print('No metal found.')

    def _parse_homo_lumo_gap(self):

        """ Parses the HOMO-LUMO gap from the output file.

        Returns:
            float: The HOMO-LUMO gap.
        """

        occ = []
        vir = []
        for i, line in enumerate(self.lines):

            if 'Alpha  occ. eigenvalues' in line:
                line_split = re.findall('-{0,1}[0-9]{1,}.[0-9]{1,}', line)
                occ.extend([float(entry) for entry in line_split])

            if 'Alpha virt. eigenvalues' in line:
                line_split = re.findall('-{0,1}[0-9]{1,}.[0-9]{1,}', line)
                vir.extend([float(entry) for entry in line_split])

        return vir[0] - occ[-1]

    def _parse_dipole_moment(self):

        """ Parses the dipole moment from the output file.

        Returns:
            float: The dipole moment.
        """

        for i, line in enumerate(self.lines):

            if 'Dipole moment (field-independent basis, Debye)' in line:
                line_split = self.lines[i + 1].split()
                return float(line_split[7])

    def _split_spectrum_into_uv_vis_nir(self, spec: list):

        """Splits a given spectrum into UV, Vis and nIR sub spectra.

        Arguments:
            spec (list[list[float]]): The spectrum.

        Returns:
            list[list[float]]: The UV spectrum.
            list[list[float]]: The Vis spectrum.
            list[list[float]]: The nIR spectrum.
        """

        spec_uv = [[], []]
        spec_vis = [[], []]
        spec_nir = [[], []]

        for nm, o in zip(spec[0], spec[1]):

            if nm < 350:
                spec_uv[0].append(nm)
                spec_uv[1].append(o)
            elif nm >= 350 and nm <= 825:
                spec_vis[0].append(nm)
                spec_vis[1].append(o)
            elif nm > 825:
                spec_nir[0].append(nm)
                spec_nir[1].append(o)

        return spec_uv, spec_vis, spec_nir

    def _get_lambda_max(self, spec: list):

        """Determines lambda max and band broadness of a given spectrum.

        Arguments:
            spec (list[list[float]]): The spectrum.

        Returns:
            float: The lambda max value.
            float: The corresponding oscillator strength.
            float: The band broadness of the given spectrum.
        """

        # return None if no transitions in range
        if len(spec[0]) == 0:
            return None, None, None

        lambda_max = spec[0][spec[1].index(max(spec[1]))]
        f_max = max(spec[1])
        # peak broadness
        sigma_max = (max(spec[0]) - min(spec[0])) * len(spec[0])

        # return None if maximum excitation is less than 0.01
        if f_max < 0.01:
            return None, None, None

        return lambda_max, f_max, sigma_max
