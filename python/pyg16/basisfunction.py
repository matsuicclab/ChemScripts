import numpy as np

class GTOBasis:
    def __init__(self, shellcoord, angulars, contractions, exponents):
        # TODO 引数チェック

        self.__shellcoord = np.array(shellcoord)
        #
        self.__angulars = angulars
        l,m,n = angulars
        self.__angularfunc = lambda r : r[:,0]**l * r[:,1]**m * r[:,2]**n
        # 縮約係数
        self.__contractions = np.array(contractions)
        # 指数
        self.__exponents = np.array(exponents)
        # 規格化定数
        self.__normalizes = np.power(2, l+m+n) * \
                                np.power(self.__factorial2(2*l-1)*self.__factorial2(2*m-1)*self.__factorial2(2*n-1), -1/2) * \
                                np.power(2/np.pi, 3/4) * \
                                np.power(exponents, (2*(l+m+n)+3)/4)

    def __factorial2(self, n):
        if n < 0:
            return (-1)**((n-1)/2) * n / self.__factorial2(-n)
        k = n // 2 + n % 2
        if n % 2 == 0:
            # 偶数の場合
            return 2**k * np.math.factorial(k)
        else:
            return np.math.factorial(2*k)/(2**k * np.math.factorial(k))

    def calc(self, r):
        """
        r: 基底関数の値を計算する座標(単位: Bohr): np.ndarray or list: shape: (3,) or shape: (*,3)
        """
        r = np.array(r).reshape(-1,3)
        _r = r - self.__shellcoord
        return self.__angularfunc(_r) * \
                np.sum(
                    self.__normalizes * self.__contractions * \
                        np.exp((-1)*self.__exponents*np.linalg.norm(_r,axis=1)[:,np.newaxis]**2),
                    axis=1
                )
