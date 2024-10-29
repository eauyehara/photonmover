from photonmover.Interfaces.Experiment import Experiment
from photonmover.utils.plot_utils import plot_graph

# Interfaces/instruments necessary for the experiment
# - You use an Interface if any instrument of that category can be used
# - You use a specific instrument if you can only use that specific model
from photonmover.instruments.Optical_spectrum_analyzers.HP70951B import HP70951B
from photonmover.instruments.Source_meters.Keithley2635A import Keithley2635A
from photonmover.instruments.Lasers.PraeviumVCSEL import PraeviumVCSEL

# General imports
import time
import numpy as np
import csv
from scipy import io


class pcf_opa_spectra(Experiment):

    def __init__(self, instrument_list, visa_lock=None):
        """
        :param instrument_list: list of available instruments. IMPORTANT: WE ASSUME THAT THE INSTRUMENTS HAVE BEEN INITIALIZED ALREADY!
        """
        super().__init__(visa_lock)

        # It is always good practice to initialize variables in the init

        # Instruments
        self.osa = None
        self.sl = None

        self.data = None

        if not self.check_necessary_instruments(instrument_list):
            raise ValueError("The necessary instruments for this experiment are not present!")

    def check_necessary_instruments(self, instrument_list):
        """
        Checks if the instruments necessary to perform the experiment are present.
        :param instrument_list: list of the available instruments
        :return: True if the necessary instruments are present, False otherwise.
        """

        for instr in instrument_list:
            if isinstance(instr, HP70951B):
                self.osa = instr
            if isinstance(instr, PraeviumVCSEL):
                self.sl = instr

        if (self.osa is not None) and (self.sl is not None):
            return True
        else:
            return False

    def get_description(self):
        """
        Returns a string with a brief summary of the experiment.
        """
        return "Sweeps seed wavelength while recording osa trace into single csv file"

    def get_name(self):
        """
        Returns a string with the experiment name
        """
        return "pcf_opa_spectra"

    def perform_experiment(self, params, filename=None):
        """
        Performs the experiment, and saves the relevant data (if there is any)
        to the specified file (if given) - assumes osa set to desired display parameters prior to running experiment
        :param params: dictionary of the parameters necessary for the experiment.
        :param filename: if specified, the data is saved in the specified file.
        :return: [wavelength_list, 2D spectrum_array]
        """

        params = self.check_all_params(params)

        seed_wls = params["seed_wls"]

        meas_volt_list = []
        est_wl_list = []
        pk_wl_arr = []
        [trace_len, [osa_cwl, osa_span, osa_reflev, osa_rbw]] = self.osa.get_osa_parameters()

        spectrum_array = np.array(np.zeros((len(seed_wls), int(trace_len))))
        [wavelength_list, _] = self.osa.read_data()

        # Sweep power supply voltage and get osa trace
        for ind, wl in enumerate(seed_wls):

            print('Setting seed wavelength to %.4f nm...' % wl)
            # Set the voltage
            [est_wl, meas_volt] = self.sl.set_wavelength(wl)

            print('Power Supply voltage set to %0.4f V' % meas_volt)
            est_wl_list.append(est_wl)
            meas_volt_list.append(meas_volt) #[V]

            # Wait [s]
            time.sleep(0.5)

            #Read osa spectrum and store in wavelength_array
            [_, amp] = self.osa.read_data()
            spectrum_array[ind, :] = np.reshape(np.array(amp), (1, int(trace_len)))
            pk_ind = np.argmax(spectrum_array[ind,:])
            pk_wl_arr.append(wavelength_list[pk_ind])

        print(np.shape(spectrum_array))
        print('Finished voltage sweep')
        print('-----------------------------')

        if filename is not None:
            # Save the data in a csv file
            time_tuple = time.localtime()
            complete_filename = "%s-%d-%d-%d_%d-%d-%d.csv" % (filename,
                                                            time_tuple[0],
                                                            time_tuple[1],
                                                            time_tuple[2],
                                                            time_tuple[3],
                                                            time_tuple[4],
                                                            time_tuple[5])

            with open(complete_filename, 'w+') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(wavelength_list)  #[nm]
                for ind in range(len(seed_wls)):
                    writer.writerow(spectrum_array[ind, :])  #[W]

                # Save the parameters in a .mat file
                time_tuple = time.localtime()
                params_filename = "params_%s_%d-%d-%d.mat" % (
                    filename,
                    time_tuple[0],
                    time_tuple[1],
                    time_tuple[2])

                io.savemat(params_filename, {'osa_cwl': osa_cwl,
                                             'osa_span': osa_span,
                                             'osa_reflev': osa_reflev,
                                             'osa_rbw': osa_rbw,
                                             'seed_wls': est_wl_list,
                                             'pump_wl': pump_wl,
                                             'seed_power': seed_power,
                                             'pump_power': pump_power,
                                             'meas_volt_list': meas_volt_list
                                             })

        self.data = [est_wl_list, pk_wl_arr]

        return [wavelength_list, spectrum_array]

    def required_params(self):
        """
        Returns a list with the keys that need to be specified in the params dictionary, in order for
        a measurement to be performed
        """
        return ["seed_wls", "pump_wl", "seed_power", "pump_power"]

    def plot_data(self, canvas_handle, data=None):

        if data is None:
            if self.data is not None:
                data = self.data
            else:
                raise ValueError('plot_data was called before performing the experiment or providing data')

        seed_wl = data[0]
        peak_wl = data[1]

        plot_graph(x_data=seed_wl, y_data=peak_wl, canvas_handle=canvas_handle, xlabel='Seed Wavelength (nm)',
                   ylabel='Peak Wavelength (nm)', title='Signal vs Idler', legend=None)


if __name__ == '__main__':
    # ------------------------------------------------------------
    # SAFETY LIMITS
    i_limit = 100e-9  # current limit

    # OTHER PARAMETERS
    seed_laser = 'Dev1a'

    # pump_laser = 'CW976' #'OEland1038'
    pump_power = 500.0 # [mW] Coupled into PCF
    pump_wl = 1050 #[nm]
    seed_power = 5.0 # [mW] Coupled into PCF

    seed_startwl = 1255 #[nm]
    seed_stopwl = 1260 #[nm]
    stepwl = 1 #[nm]

    FC = 0.5  # fiber coupling efficiency at output
    RBW = 0.5  # nm (OSA)

    # EXPERIMENT PARAMETERS
    wavvolt_file = "wavvolt_Dev1a_25C_CW976_5.0mW"
    # ------------------------------------------------------------

    # INSTRUMENTS
    smu = Keithley2635A(current_compliance=100e-9, voltage_compliance=81) #A, V
    sl = PraeviumVCSEL(smu, wavvolt_file)
    osa = HP70951B()

    # Initialize instruments
    sl.initialize()
    osa.initialize()

    seed_wls = sl.generate_valid_sweep(seed_startwl, seed_stopwl, stepwl)

    file_name = "pcfOPA_%s_%d-%dnm_%dmW_OE%3.2fnm_%dmW" % (
        seed_laser, seed_wls[0], seed_wls[-1], seed_power, pump_wl, pump_power)  # Filename where to save csv data

    # SET UP THE EXPERIMENT
    instr_list = [osa, sl]
    params = {"seed_wls": seed_wls, "pump_wl": pump_wl, "seed_power": seed_power, "pump_power": pump_power}
    exp = pcf_opa_spectra(instr_list)

    # RUN IT
    [wavelength_list, spectrum_array] = exp.perform_experiment(params, filename=file_name)

    # CLOSE INSTRUMENTS
    osa.close()
    sl.close()
