import numpy as np
from CalcRHNpara import CalcRHNParameters


def calcallparameters(t12, t23, t13, dm2, Dm2, abLambda, theta, phi, RorD):

    calcNP = CalcRHNParameters(t12, t23, t13, dm2, Dm2,abLambda, theta, phi, RorD)
    calcNP.calcinputTriAngle()
    calcNP.calcallNparameters()
    calcNP.calcRHNMassMatrix()
    calcNP.diagRHNMassMatrix()
    calcNP.show_inputall()
    #calcNP.show_outputall()
    calcNP.show_output(*["RHNmass[e+20]","U"])
    return calcNP.output


np.set_printoptions(precision=4, suppress=True)

calcallparameters(33.44, 49.2, 8.57, 7.42 * 10 ** (-5), 2.517 * 10 ** (-3), 0.1, 60, 30, "deg")
calcallparameters(33, 40.5, 8.4, 7.42 * 10 ** (-5), 2.517 * 10 ** (-3), 0.1, 60, 30, "deg")
calcallparameters(33, 43, 8.5, 7.42 * 10 ** (-5), 2.517 * 10 ** (-3), 0.1, 60, 30, "deg")

"test"



