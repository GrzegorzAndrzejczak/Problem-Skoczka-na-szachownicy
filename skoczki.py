import numpy as np
from numpy import (int8, int16)
from builtins import int
import io
"""
Szachownice 10x10 (lub 6x6) przepisane na Pythona
"""


class Szach:
    M: int = 0
    N: int8 = 10
    E0: np.ndarray
    E: np.ndarray
    save: bool = True
    plik: io.BytesIO

    @staticmethod
    def read_data(file_name=f'skoczki{N}x{N}.bin') -> np.ndarray:
        with io.open(file_name, "rb") as plik:
            pl = plik.read()
        row = Szach.N * Szach.N // 4
        Szach.M = len(pl) // row
        X = np.ndarray((Szach.M,row), dtype=np.int8)
        with io.open(file_name, "rb") as plik:
            plik.readinto(X)
        return X


    @staticmethod
    def nr(x, y: int8) -> int16:
        N = Szach.N
        return (N//2)*(N-1-y)+(x+N-1)//2

    @staticmethod
    def x(nr: int16) -> int8:
        N = Szach.N
        return 2*(nr % N)-(N-1)

    @staticmethod
    def y(nr: int16) -> int8:
        N = Szach.N
        return -2*(nr//N)+(N-1)

    @staticmethod
    def E0_to(n: int16, x1, y1: int8, gdzie: list) -> None:     # zaznacz, jeśli jest OK
        N, nast = Szach.N, gdzie[0]
        # print(f'36: x1={x1}, y1={y1}, n={n}, nast={nast}, N={N}')
        if (abs(x1)<N) & (abs(y1)<N):
            Szach.E0[n][nast] = Szach.nr(x1,y1)
            gdzie[0] += 1

    @staticmethod
    def _del(E0_n: np.ndarray, x, y) -> None:   # wykluczamy kilka połączeń wokół środka
        nr, k = Szach.nr(x, y), 7
        while E0_n[k] < 0:   k -= 1
        i = k
        while E0_n[i]!=nr:   i -= 1
        E0_n[i] = E0_n[k]
        E0_n[k] = -1

    @staticmethod
    def initE0():       # STAŁA tablica możliwych połączeń
        N = Szach.N
        # print("!", N)
        Szach.E0 = np.ndarray((N*N,8), dtype=int16)
        Szach.E0[:, :] = -1
        E0_to = Szach.E0_to
        for n in range(N*N):
            nast = [0]
            x, y = Szach.x(n), Szach.y(n)
            E0_to(n, x + 2, y + 4, nast)
            E0_to(n, x + 2, y - 4, nast)
            E0_to(n, x - 2, y + 4, nast)
            E0_to(n, x - 2, y - 4, nast)
            E0_to(n, x + 4, y + 2, nast)
            E0_to(n, x + 4, y - 2, nast)
            E0_to(n, x - 4, y + 2, nast)
            E0_to(n, x - 4, y - 2, nast)
            # if Szach.e0Nast(x+2, y+4, n, nast):
            #     nast += 1
            # print(f"73: n={n},x={x},y={y},nast={nast[0]}")
            if (abs(x)==1) & (abs(y)==3) | (abs(x)==3) & (abs(y)==1):   # wykluczone
                Szach._del(Szach.E0[n], -y, x)
                Szach._del(Szach.E0[n], y, -x)

    def is_loop(self, p: int16) -> bool:     # długość ewent. pętli z "p" - lub -1
        m = 1
        E = self.E
        n = E[p][0]
        back: int8 = 0
        if E[n][0]!=p:
            back = 1      # E[n][back]--n--E[n][1-back]
        while n!=p:
            n1 = E[n][1-back]
            if n1<0:
                return False
            m += 1
            if E[n1][0]==n: back = 0
            else:           back = 1
            n = n1
        N = Szach.N
        if m == N*N:
            Szach.M += 1
            if Szach.M % 1000==0:
                print(Szach.M)
            if Szach.save:
                seq = np.ndarray((N*N//4,), dtype=np.uint8) # aż do 14*14
                m = 0
                back = 1
                n = 0
                while (n!=N-1) & (n!=N*N-N):
                    n1 = E[n][1-back]
                    seq[m] = n
                    m += 1
                    if E[n1][0]==n: back = 0
                    else:           back = 1
                    n = n1
                Szach.plik.write(seq)
                if m != N*N//4:
                    print("?")
        return True

    def _addEdge(self, p, p1: int16, check: bool) -> bool:  # z obu stron i po obrotach
        x, y = Szach.x(p), Szach.y(p)
        x1, y1 = Szach.x(p1), Szach.y(p1)
        E = self.E
        if E[p][0]<0: E[p][0] = p1
        else:         E[p][1] = p1
        if E[p1][0]<0: E[p1][0] = p
        else:          E[p1][1] = p
        q, q1 = Szach.nr(-y, x), Szach.nr(-y1, x1)
        if E[q][0]<0: E[q][0] = q1
        else:         E[q][1] = q1
        if E[q1][0]<0: E[q1][0] = q
        else:          E[q1][1] = q
        r, r1 = Szach.nr(-x, -y), Szach.nr(-x1, -y1)
        if E[r][0]<0: E[r][0] = r1
        else:         E[r][1] = r1
        if E[r1][0]<0: E[r1][0] = r
        else:          E[r1][1] = r
        s, s1 = Szach.nr(y, -x), Szach.nr(y1, -x1)
        if E[s][0]<0: E[s][0] = s1
        else:         E[s][1] = s1
        if E[s1][0]<0: E[s1][0] = s
        else:          E[s1][1] = s
        if check & (E[p1][1]>=0):         # sprawdzić ewent. pętlę
            if self.is_loop(p):
                return True
        return False

    def _tryComplete(self) -> bool:      # jakieś wymuszone pary krawędzi?
        changed: bool = True
        E, E0 = self.E, Szach.E0
        N = Szach.N
        while changed:
            changed = False
            for p in range(1, N*N):
                if E[p][1]>=0:
                    continue
                nFree = 0
                m1, m2 = -1, -1
                for i in range(8):
                    n1 = E0[p][i]           # dostępne z "p"
                    if (n1<0) | (E[n1][1]>=0):
                        continue
                    if m1<0: m1 = n1
                    else:    m2 = n1        # może się przydadzą
                    nFree +=1
                if nFree==0:
                    return True             # bez loop() - BAD
                if E[p][0]<0:               # potrzebne są 2 krawędzie
                    if nFree==1:
                        return True         # bez loop() - BAD
                    if nFree==2:
                        self._addEdge(p, m1, False)
                        self._addEdge(p, m2, True)
                        changed = True
                elif nFree==1:              # potrzebna jest 1 krawędź
                    self._addEdge(p, m1, True)
                    changed = True
        return False

    def __init__(self, N: int | None = None, save: bool = False):
        if N is not None:
            Szach.save = save
            Szach.M, Szach.N = 0, N
            Szach.initE0()
            self.E = np.ndarray((N*N, 2), dtype=int16)
            self.E[:][:] = -1                   # brak jeszcze krawędzi
            self._addEdge(0, N+2, False)        # w wierzchołkach szachownicy
            self._addEdge(0, 2*N+1, False)

    def duplicate(self):
        kopia = Szach()
        kopia.E = self.E.copy()
        return kopia

    def extend(self, nr, nast: int16):
            # nr ma wolny 1 koniec, i<8 jest wolne i jest propozycją następnej krawędzi
            # .... i WSZYSTKIE konsekwencje
        N = Szach.N
        sNast: Szach = self.duplicate()
        if sNast._addEdge(nr, nast, True):      # powstała jakaś pętla - koniec
            return                              # było w tym sNast.tryComplete()
        if sNast._tryComplete():                # ta wersja skończyła się
            return
            # znaleźć 1! wolny koniec i "extend" o możliwe krawędzie
        for p in range(1, N*N):
            if (sNast.E[p][0]>=0) & (sNast.E[p][1]<0): break
        else:
            return
        for i in range(8):
            n1 = Szach.E0[p][i]                 # nr -- p -- n1?
            if n1<0: break                      # brak wolnych krawędzi
            if n1==sNast.E[p][0]:
                continue                        # krawędź dojścia
            if sNast.E[n1][1]>=0:
                continue                        # wierzchołek jest już zajęty
            sNast.extend(p, n1)

    def extend_to_all(self):
            # od pierwszego wolnego końca
        N = Szach.N
        if Szach.save:
            Szach.plik = io.open(f"skoczki{N}x{N}.bin", "wb+")
        for i in range(8):
            nast = Szach.E0[2*N+1][i]
            if nast<0:
                continue
            if self.E[nast][1]<0:               # ...to także wyklucza cofanie!
                self.extend(2 * N + 1, nast)
        if Szach.save:
            Szach.plik.close()
        print(f'M={Szach.M}')

    def read(self, n: int, dataset=None) -> np.ndarray:
        N = Szach.N
        if dataset is None:
            seq = np.ndarray((N*N//4,), dtype=np.uint8)
            with io.open(f"skoczki{N}x{N}.bin", "rb") as plik:
                plik.seek(n*seq.size)
                plik.readinto(seq)
        else:
            seq = dataset[n]
            if seq.shape[0]!=N*N//4:
                print("Błędny rozmiar danych!!!")
        # print(seq)
        path = np.ndarray((N*N+1,), dtype=np.int16)
        m = N*N//4
        for i in range(m):
            path[i] = seq[i]
            path[2*m+i] = N*N-seq[i]-1
        zx, zy = -1, 1
        if seq[N*N//4-1] > N*N//2:
            zx, zy = 1, -1
        for i in range(m):
            x, y = Szach.x(seq[i]), Szach.y(seq[i])
            path[m+i] = Szach.nr(zy*y, zx*x)
            path[3*m+i] = Szach.nr(-zy*y, -zx*x)
        path[N*N] = path[0]
        return path

def main():
    a6=Szach(6)
    path = a6.read(5)
    xp = Szach.x(path)
    yp = Szach.y(path)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.axis('equal')
    ax.plot(xp,yp, "r")
    # nn = plt.Polygon(np.asarray([xp, yp]).T, fill=True)
    # ax.add_patch(nn)
    plt.show(block=False)


