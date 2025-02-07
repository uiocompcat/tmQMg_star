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

        """Parses the spectrum from the outputfile.

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
