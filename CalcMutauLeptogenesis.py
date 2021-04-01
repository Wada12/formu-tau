import os
import subprocess
import numpy as np
from CalcRHNpara import CalcRHNParameters
from Tools import writeRecord,readInputRecords


class PreProcessing:
    def __init__(self, t12, t23, t13, dm2, Dm2, abLambda, theta, phi):
        self.t12, self.t23, self.t13 = t12, t23, t13
        self.dm2, self.Dm2 = dm2, Dm2
        self.abLambda, self.theta, self.phi = abLambda, theta, phi

    def calcmutauPara(self):
        calcNP = CalcRHNParameters(self.t12, self.t23, self.t13, self.dm2, self.Dm2, self.abLambda, self.theta,
                                   self.phi, "deg")
        calcNP.calcinputTriAngle()
        calcNP.calcallNparameters()
        calcNP.calcRHNMassMatrix()
        calcNP.diagRHNMassMatrix()
        calcNP.calcYukawaabsandphs()
        #calcNP.show_inputall()
        self.RHNmass = calcNP.output["RHNmass[e+20]"] * 10 ** 11  # GeV
        self.Yukawa_abs, self.Yukawa_phs = calcNP.output["Yukawa_abs"], calcNP.output["Yukawa_phs"]

    def inputcards_generator(self, filename, Yukawa_abs, Yukawa_phs, RHNmass):
        # os.makedirs("inputcards", exist_ok=True)
        os.makedirs("inputcards", exist_ok=True)
        os.chdir("inputcards")
        input_file = open(filename, "w", encoding="utf-8")
        line_num = [1, 2, 3]
        for i in line_num:
            for j in line_num:
                Yukawa_abs_ij = Yukawa_abs[i - 1][j - 1]
                line = "Y" + str(i) + str(j) + "_mag" + " " + str(Yukawa_abs_ij) + "\n"
                input_file.write(line)

        for i in line_num:
            for j in line_num:
                Yukawa_phs_ij = Yukawa_phs[i - 1][j - 1]
                line = "Y" + str(i) + str(j) + "_phs" + " " + str(Yukawa_phs_ij) + "\n"
                input_file.write(line)

        for i in line_num:
            log_RHNmass_i = np.log10(RHNmass[i - 1])
            line = "M" + str(i) + "  " + str(log_RHNmass_i) + "\n"
            input_file.write(line)

        input_file.close()
        os.chdir("../")

    def recordInput(self, filename):
        os.makedirs("inputrecords", exist_ok=True)
        os.chdir("inputrecords")
        file_name = os.path.splitext(os.path.basename(filename))[0]
        recordfile = open(file_name + "_record.txt", "w", encoding="utf-8")
        self.record = \
            writeRecord(self.t12, self.t23, self.t13, self.dm2, self.Dm2, self.abLambda, self.theta, self.phi,
                        self.RHNmass[0],self.RHNmass[1],self.RHNmass[2])
        recordfile.write(self.record)
        recordfile.close()
        os.chdir("../")

    def preProcessing(self, filename):
        # filename = "sample.dmg"
        os.chdir("/Users/wadajuntaro/PycharmProjects/mu_tau_Ulysses")
        self.calcmutauPara()
        self.inputcards_generator(filename, self.Yukawa_abs, self.Yukawa_phs, self.RHNmass)
        self.recordInput(filename)


class ControlUlysses:

    def outputfile_generator(self, output_filename, proc, input_records=""):
        os.makedirs("output", exist_ok=True)
        os.chdir("output")
        outputfile = open(output_filename, "w", encoding="utf-8")
        outputfile.write(input_records)
        outputfile.write("--output--\n")
        for line in proc.stdout:
            if "eta_b" in str(line) or "Y_b" in str(line) or "Omega_b h^2" in str(line):
                new_line = str(line).replace("\\n", "\n").replace("b'", "").replace("'", "")
                outputfile.write(new_line)
        outputfile.close()
        os.chdir("../")


    def uls_calc(self, BE, input_filename,input_records=""):
        # input_filename = "sample.dmg"
        os.chdir("/Users/wadajuntaro/PycharmProjects/mu_tau_Ulysses")
        file_name = os.path.splitext(os.path.basename(input_filename))[0]
        proc = subprocess.Popen(["uls-calc", '-m', BE, 'inputcards/' + input_filename]
                                , stdout=subprocess.PIPE)
        # it is independent on choosing the value of m2atm,m2solor
        output_filename = "output_{0}_{1}.txt".format(file_name, BE)
        self.outputfile_generator(output_filename, proc, input_records)



if __name__ == '__main__':
    filename = "test.dmg"
    pre = PreProcessing(33.44, 49.2, 8.57, 7.42 * 10 ** (-5), 2.517 * 10 ** (-3), 1*10**(-10), 60, 30)
    pre.preProcessing(filename)
    input_records = readInputRecords(filename)
    ulysses = ControlUlysses()
    ulysses.uls_calc("1DME", filename,input_records)

