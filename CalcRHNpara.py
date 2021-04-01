from scipy import optimize
from scipy import linalg
import numpy as np
import pprint


class CalcNParameters:

    def __init__(self, t12, t23, t13, dm2, Dm2, RorD: str):
        """
        this calculates N's phases such as delta,alpha2,alpha3,and mi
        from input parameter:t12,t13,t23,dm2,Dm2

        :param t12:mixing angle
        :param t13:
        :param t23:
        :param dm2:deltam^2=m2^2-m1^2
        :param Dm2:Deltam^2=m3^2-(m1^2+m2^2)/2
        """

        self.originalInputAngle = {"t12": t12, "t23": t23, "t13": t13}

        if RorD == 'rad':
            self.t12, self.t23, self.t13 = t12, t23, t13
        elif RorD == 'deg':
            self.t12, self.t23, self.t13 = np.deg2rad(t12), np.deg2rad(t23), np.deg2rad(t13)
        else:
            print('Please input str: "rad" or "deg"')

        self.dm2, self.Dm2 = dm2, Dm2

        self.inputAngle = {"t12": self.t12, "t23": self.t23, "t13": self.t13}
        self.ep = self.dm2 / (self.Dm2 + self.dm2 / 2)

        self.c = {}  # cos mixing angle
        self.s = {}  # sin mixing angle
        self.c2 = {}  # cos2t
        self.s2 = {}  # sin2t

        self.inputTriAngle = {}
        self.output = {}

    def calcinputTriAngle(self):

        for Input in self.inputAngle:
            self.c["{}".format(Input)] = np.cos(self.inputAngle["{}".format(Input)])
            self.s["{}".format(Input)] = np.sin(self.inputAngle["{}".format(Input)])
            self.c2["{}".format(Input)] = np.cos(2 * self.inputAngle["{}".format(Input)])
            self.s2["{}".format(Input)] = np.sin(2 * self.inputAngle["{}".format(Input)])

        self.inputTriAngle = {"c": self.c, "s": self.s}

        return self.inputTriAngle


    def calcdelta(self):
        def Cubiceq(x):
            cubiceqep0 \
                = self.s["t13"] ** 2 * (4 * self.s["t13"] ** 2 * self.c2["t12"] ** 2 * self.c2["t23"] ** 2
                                        - self.s["t13"] * 2*self.c2["t12"]*self.s2["t12"] * 2*self.c2["t23"]*self.s2["t23"] * (1 + self.s["t13"] ** 2) * x
                                        + self.s2["t12"] ** 2 * self.s2["t23"] ** 2 * (
                                                self.c["t13"] ** 4 + 4 * self.s["t13"] ** 2 * x ** 2)
                                        ) \
                  * (2 * (2 * self.c2["t12"] * self.c2["t23"] ** 2
                          - self.s["t13"] * self.s2["t12"] * 2*self.c2["t23"]*self.s2["t23"] * x)
                     )

            epterm = -self.ep * (
                    (4 * self.s["t12"] ** 4 * self.c2["t23"] ** 2
                     + self.s["t13"] ** 2 * self.s2["t12"] ** 2 * self.s2["t23"] ** 2
                     + 4 * self.s["t12"] ** 3 * self.c["t12"] * self.s["t13"] * 2*self.c2["t23"]*self.s2["t23"] * x
                     )
                    * (4 * self.c2["t23"] ** 2 * (self.c["t12"] ** 4 * self.c["t13"] ** 4
                                                  - self.s["t13"] ** 4 * self.c2["t12"] ** 2)
                       - self.s["t13"] * 2*self.c2["t23"]*self.s2["t23"]
                       * (4 * self.c["t13"] ** 4 * self.c["t12"] ** 3 * self.s["t12"]
                          - self.s["t13"] ** 2 * 2*self.c2["t12"]*self.s2["t12"] * (1 + self.s["t13"] ** 2)
                          ) * x
                       - 4 * self.s["t13"] ** 4 * self.s2["t12"] ** 2 * self.s2["t23"] ** 2 * x ** 2
                       )
            )
            cubiceq = cubiceqep0 + epterm
            return cubiceq

        sol = optimize.fsolve(Cubiceq, 0)
        delta1 = np.arccos(sol)
        delta2 = delta1 + 2 * (np.pi - delta1)  # up to Pi radian
        self.output["cos_delta"] = sol[0]
        self.output["delta"] = delta2[0]
        return self.output

    def calcmi(self):
        cos_delta = self.output["cos_delta"]

        numerator \
            = 4 * self.s["t12"] ** 4 * self.c2["t23"] ** 2 \
              + 4 * self.s["t12"] ** 3 * self.c["t12"] * self.s["t13"] * 2*self.c2["t23"]*self.s2["t23"] * cos_delta \
              + self.s["t13"] ** 2 * self.s2["t12"] ** 2 * self.s2["t23"] ** 2
        denominator \
            = 2 * (2 * self.c2["t12"] * self.c2["t23"] ** 2
                   - self.s["t13"] * self.s2["t12"] * 2*self.c2["t23"]*self.s2["t23"] * cos_delta)

        m1m1 = self.dm2 * (numerator / denominator)
        m1 = np.sqrt(m1m1)
        m2m2 = self.dm2 + m1m1
        m2 = np.sqrt(m2m2)
        m3m3 = self.Dm2 + (m2m2 + m1m1) / 2
        m3 = np.sqrt(m3m3)

        self.output["m1"] = m1
        self.output["m2"] = m2
        self.output["m3"] = m3

        return self.output

    def calcalpha(self):
        exp_delta = np.exp(1j * self.output["delta"])
        exp_delta2 = exp_delta ** 2
        exp_delta_inv = exp_delta ** (-1)

        R2 \
            = -(2 * self.s["t12"] ** 2 * self.c2["t23"] + self.s2["t12"] * self.s2["t23"] * self.s["t13"] * exp_delta) \
              / (2 * self.c["t12"] ** 2 * self.c2["t23"] - self.s2["t12"] * self.s2["t23"] * self.s["t13"] * exp_delta)

        R3 \
            = -(self.s["t13"] * exp_delta2 * (2 * self.c2["t12"] * self.c2["t23"] * self.s["t13"]
                                              - self.s2["t12"] * self.s2["t23"] * (
                                                      exp_delta_inv + self.s["t13"] ** 2 * exp_delta)
                                              )
                ) / (self.c["t13"] ** 2 * (2 * self.c["t12"] ** 2 * self.c2["t23"]
                                           - self.s2["t12"] * self.s2["t23"] * self.s["t13"] * exp_delta))

        alpha2 = np.angle(R2 * self.output["m2"] / self.output["m1"])
        alpha3 = np.angle(R3 * self.output["m3"] / self.output["m1"])

        self.output["alpha2"] = alpha2
        self.output["alpha3"] = alpha3

        return self.output

    def calcallNparameters(self):
        self.calcdelta()
        self.calcmi()
        self.calcalpha()
        return self.output

    def show_inputall(self):
        print("--Input--")
        print(self.originalInputAngle)
        print("'δm^2'={0},'Δm^2'={1}".format(self.dm2, self.Dm2))
        print("---------")

    def show_outputall(self):
        print("--Output--")
        pprint.pprint(self.output, sort_dicts=False)
        print("---------")




class CalcRHNParameters(CalcNParameters):
    Higgs_vev = 174 * 10 ** 9

    def __init__(self, t12, t23, t13, dm2, Dm2, abLambda, theta, phi, RorD: str):
        super().__init__(t12, t23, t13, dm2, Dm2, RorD)

        self.originalInputLambda = {"abLambda":abLambda,"theta":theta,"phi":phi}
        self.abLambda = abLambda

        if RorD == 'rad':
            self.theta, self.phi = theta, phi
        elif RorD == 'deg':
            self.theta, self.phi = np.deg2rad(theta), np.deg2rad(phi)
        else:
            print('Please input str: "rad" or "deg"')

        self.lambda_e = self.abLambda * np.cos(self.theta)
        self.lambda_mu = self.abLambda * np.sin(self.theta) * np.cos(self.phi)
        self.lambda_tau = self.abLambda * np.sin(self.theta) * np.sin(self.phi)

        self.matrix = {}

    def calcRHNMassMatrix(self):

        exp_delta = np.exp(1j * self.output["delta"])
        exp_alpha2 = np.exp(1j * self.output["alpha2"] / 2)
        exp_alpha3 = np.exp(1j * self.output["alpha3"] / 2)

        Vmatrix = np.array(
            [[self.c["t12"] * self.c["t13"], self.s["t12"] * self.c["t13"], self.s["t13"] * exp_delta ** (-1)],
             [-self.s["t12"] * self.c["t23"] - self.c["t12"] * self.s["t23"] * self.s["t13"] * exp_delta
                 , self.c["t12"] * self.c["t23"] - self.s["t12"] * self.s["t23"] * self.s["t13"] * exp_delta
                 , self.s["t23"] * self.c["t13"]],
             [self.s["t12"] * self.s["t23"] - self.c["t12"] * self.c["t23"] * self.s["t13"] * exp_delta
                 , -self.c["t12"] * self.s["t23"] - self.s["t12"] * self.c["t23"] * self.s["t13"] * exp_delta
                 , self.c["t23"] * self.c["t13"]]])
        Pmatrix = np.array([[1, 0, 0], [0, exp_alpha2, 0], [0, 0, exp_alpha3]])

        Ndiagmatrix_inv \
            = np.array([[1 / self.output["m1"], 0, 0], [0, 1 / self.output["m2"], 0], [0, 0, 1 / self.output["m3"]]])
        Diracmatirix = CalcRHNParameters.Higgs_vev \
                       * np.array([[self.lambda_e, 0, 0], [0, self.lambda_mu, 0], [0, 0, self.lambda_tau]])

        UPMNS = Vmatrix.dot(Pmatrix)
        self.matrix["UPMNS"] = UPMNS

        Nmatrix_inv = UPMNS @ Ndiagmatrix_inv @ UPMNS.T
        MR = - Diracmatirix @ Nmatrix_inv @ Diracmatirix
        self.matrix["MR[e+20]"] = MR * 10 ** (-20)

        return self.matrix

    def diagRHNMassMatrix(self):

        MR = self.matrix["MR[e+20]"]
        MRcon = np.conjugate(self.matrix["MR[e+20]"].T)
        MR2 = MRcon @ MR
        M,omega = linalg.eigh(MR2)
        T = omega.T @ MR @ omega     # T is already diagonalized
        ph = np.diag(np.exp(-1j* np.angle(np.diag(T))/2))
        self.Omega = omega @ ph
        # MRdiag= self.Omega.T @ MR @ self.Omega
        self.output["RHNmass[e+20]"] = np.sqrt(M)
        self.output["Omega"] = self.Omega

        return self.output

    def calcYukawaabsandphs(self):
        Yukawadiag = np.diag([self.lambda_e,self.lambda_mu,self.lambda_tau])
        Yukawa_hat = Yukawadiag @ np.linalg.inv(self.Omega)
        Yukawa_abs = np.abs(Yukawa_hat)
        Yukawa_phs = np.angle(Yukawa_hat)
        self.output["Yukawa_hat"] = Yukawa_hat
        self.output["Yukawa_abs"]= Yukawa_abs
        self.output["Yukawa_phs"]= Yukawa_phs
        #print(Yukawa_hat)
        #print(self.Yukawa_abs)
        #print(self.Yukawa_phs)

        return self.output


    def show_inputall(self):
        super().show_inputall()
        print("'|λ|'={0},'θ'={1},'φ'={2}".format(self.abLambda, self.originalInputLambda["theta"],
                                                 self.originalInputLambda["phi"]))
        print("Higgs vev={}".format(CalcRHNParameters.Higgs_vev))
        print("---------")

    def show_matrix(self):
        print("--Matrix--")
        pprint.pprint(self.matrix, sort_dicts=False)
        print("---------")

    def show_outputall(self):
        super().show_outputall()

    def show_output(self,paramater1,*paramaters):
        print("--Output--")
        print(paramater1)
        print(self.output[paramater1])
        for paramater in paramaters:
            print(paramater)
            print(self.output[paramater])
        print("---------")


if __name__ == '__main__':
    test = CalcRHNParameters(33.44, 49.2, 8.57, 7.42 * 10 ** (-5), 2.517 * 10 ** (-3), 0.1, 60, 30, "deg")
    test.calcinputTriAngle()
    test.calcallNparameters()
    test.calcRHNMassMatrix()
    test.diagRHNMassMatrix()
    test.calcYukawaabsandphs()
    test.show_output(*["RHNmass[e+20]", "Omega"])