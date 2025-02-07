import re


class NtoDataParser():

    """Parser class for NTO output files."""

    def __init__(self, file_path):

        self._id = file_path.split('/')[-1].split('.')[0]

        with open(file_path, 'r') as fh:
            self.lines = fh.read().strip().split('\n')

    def parse(self):

        """Parses the given NTO output file.
        """

        return_dict = {'id': self._id}

        if self._has_failed():
            return_dict['has_failed'] = True
            return return_dict
        else:
            return_dict['has_failed'] = False

        for i, line in enumerate(self.lines):

            if 'Eigenvalues --' in line:
                nto_data, nto_type = self._extract_nto_data(i)

                if nto_type == 'O':
                    return_dict['occupied_nto'] = nto_data
                elif nto_type == 'V':
                    return_dict['virtual_nto'] = nto_data

        return return_dict

    def _has_failed(self):

        """Checks if the job has failed.

        Returns:
            bool: The flag indicating whether the job has failed.
        """

        if 'Normal termination' not in self.lines[-1]:
            return True

        return False

    def _extract_nto_data(self, start_index):

        """Extracts the NTO data.

        Arguments:
            start_index (int): The line index to start parsing from.

        Returns:
            dict: Dictionary of NTO results.
            str: NTO type, occupied (O) or virtual (V).
        """

        # determine number of columns
        n_columns = len(self.lines[start_index].split()) - 2

        # set the index of the column to extract
        max_column_index = 5

        # determine type
        if 'O' in self.lines[start_index - 1]:
            nto_type = 'O'
        elif 'V' in self.lines[start_index - 1]:
            nto_type = 'V'
        else:
            print('NTO type not recognized.')
            exit()

        # determine eigenvalues
        eigenvalues = [float(_) for _ in self.lines[start_index].split()[2:]]

        # data variable
        nto_data = []
        # start parsing from next line
        i = start_index + 1
        while len(self.lines[i].split()) > n_columns:


            line_split = self.lines[i].split()

            if len(line_split) == 4 + n_columns:

                atom_index = int(line_split[1])
                atom_element = line_split[2]

                ntos = {eigenvalues[j]: {line_split[3]: float(line_split[4 + j])} for j in range(max_column_index)}

                nto_data.append({
                    'atom_index': atom_index,
                    'atom_element': atom_element,
                    'ntos': ntos
                })

            else:

                line_split = re.split(r'\s{5,}', self.lines[i].strip())
                orbital_id = line_split[1]
                line_split = line_split[2].split()

                for j in range(max_column_index):
                    nto_data[-1]['ntos'][eigenvalues[j]][orbital_id] = float(line_split[j])


            i += 1

        return nto_data, nto_type
