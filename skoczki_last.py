import matplotlib.pyplot as plt

from skoczki import Szach
import numpy as np
from fractions import Fraction as Fra


import pickle
with open(".\\skoczki10x10.pkl","br") as f:
    Xp = pickle.load(f)


def analyse(xp,yp):
    if len(xp)==37: N=6
    else: N=10

    def pasma(xp, yp):
        if len(xp) == 37:
            N = 6
        else:
            N = 10
        Pasma = [[] for i in range(N - 1)]  # [] * (N-1) jest błędne -> tworzy kopie tego samego
        for j in range(N * N):
            ya, yb = yp[j], yp[j + 1]
            xa, xb = xp[j], xp[j + 1]
            if xa > xb:
                xa, xb, ya, yb = xb, xa, yb, ya
            i = (xa + N) // 2
            if xb - xa == 2:
                Pasma[i].append([ya, yb])
            else:
                yab = (ya + yb) // 2
                Pasma[i].append([ya, yab])
                Pasma[i + 1].append([yab, yb])

        for pas in Pasma:  # lista krawędzi w ramach kolejnych pasm
            pas.sort()
        return [np.asarray(pas) for pas in Pasma]

    Pasma = pasma(xp,yp)
    Poly = [[] for _ in Pasma]
    for k,Pas in enumerate(Pasma):
        rPas, r_index = np.unique(Pas[:, 1], return_index=True)
        lPas, l_index = np.unique(Pas[:, 0], return_index=True)
        poly = Poly[k]
        queue = []
        v_pas = [[] for _ in Pas]
        # v_1 = {(0, 1)}
        for i in range(len(Pas)):
            a, b = Pas[i]
            j = i + 1
            for j in range(i + 1, len(Pas)):
                la, lb = Pas[j]
                da = la - a
                if da == 0:
                    # print(f"{i,j}:0")
                    v_pas[i].append([Fra(0), j])
                    v_pas[j].append([Fra(0), i])
                    # v_1.add((i, j))
                    continue
                db = b - lb
                if db >= 0:
                    fra = Fra(da, da + db)
                    # print(f"{i,j}: {fra}")
                    v_pas[i].append([fra, j])
                    v_pas[j].append([fra, i])
                    # v_1.add((i, j))
                    continue
                if la > b + 4:
                    j = len(Pas)
        for vp in v_pas:
            vp.sort()
            if len(vp) == 0 or vp[-1][0] < 1:
                vp.append([Fra(1), -1])
            if vp[0][0] > 0:
                vp.insert(0, [Fra(0), -1])

        def to_queue(bw, j, ia, dir):  # bw=0/1
            for i in range(len(queue)):
                if queue[i] == (1 - bw, j, ia + dir, -dir):
                    queue.__delitem__(i)  # usunąć, jeśli już było, dir = +/-1
                    return
            queue.append((bw, j, ia, dir))  # j: ia -- ia+dir

        def from_queue():
            return queue.pop(0)  # pobiera i usuwa

        def the_same(u, v):
            if u[0] != v[0]:
                return False
            if type(u[1]) == list:
                n1 = u[1][0]
                if type(v[1]) == list:
                    n2 = v[1][0]
                    return n1 == n2
                else:
                    j = v[1]
                    return n1 == Pas[j, u[0]]
            elif type(v[1]) == list:
                return the_same(v, u)
            else:
                a = u[0]  # = v[0] - Fraction
                u0, u1 = Pas[u[1]]
                v0, v1 = Pas[v[1]]
                y = u0 + a * (u1 - u0)
                return y == v0 + a * (v1 - v0)

        def r_next(j: int, ie: int, bw):  # pas[j][ie-1 -- ie]  ... i co dalej? bw=0/1
            vertex = poly[-1]
            # punkt (a,j) ~ v_pas[j][ie-1] powinien już być ostatni na liście vertex
            b, j1 = v_pas[j][ie]  # raczej v_pas_1
            vertex.append((b, j))
            if the_same((b, j), vertex[1]):
                return -1, b, False
            # to_queue(1-bw, j, ie, -1)   # pominąć, mogło przyjść z kolejki
            if j < j1:  # nawrót - w lewo w górę
                for ix, v_ in enumerate(v_pas[j1]):
                    if v_[1] == j:
                        break
                to_queue(1 - bw, j1, ix - 1, 1)
                vertex.append((b, j1))  # ten sam koniec, inaczej
                if the_same((b, j1), vertex[1]):
                    return -1, b, False  # .........................KONIEC
                return j1, ix - 1, True  # b, l_next...
            if b < 1:  # w prawo w górę
                for ix, v_ in enumerate(v_pas[j1]):
                    if v_[1] == j:
                        break
                to_queue(1 - bw, j1, ix + 1, -1)
                vertex.append((b, j1))
                if the_same((b, j1), vertex[1]):
                    return -1, b, False  # .........................KONIEC
                return j1, ix + 1, False  # b, r_next
            else:  # b=1, w górę, jeśli można
                nr = Pas[j, 1]  # pozycja końca
                vertex.append((1, [nr]))  # ten sam koniec, w pionie
                if nr < rPas[-1]:  # jeszcze jest coś wyżej
                    ve = np.argwhere(rPas == nr)[0, 0]  # od wierzchołka <ve>
                    nr1 = rPas[ve + 1]
                    vertex.append((1, [nr1]))  # kontynuacja w górę
                    if the_same((1, [nr1]), vertex[1]):
                        return -1, 1, False
                    # i dalej w lewo
                    j1 = r_index[ve + 1]
                    to_queue(1 - bw, j1, -2, 1)
                    vertex.append((1, j1))
                    return j1, -2, True  # od brzegu w lewo
                else:
                    if nr < N - 1:
                        vertex.append((1, [N - 1]))
                    vertex.append((0, [N - 1]))  # górna krawędź w poziomie
                    nr = lPas[-1]
                    if nr < N - 1:
                        vertex.append((0, [nr]))
                    if the_same((0, [nr]), vertex[1]):
                        return -1, 0, False
                    j1 = l_index[-1]
                    _, j2 = v_pas[j1][0]
                    j1 = max(j1, j2)  # jeśli są 2 krawędzie - górna
                    to_queue(1 - bw, j1, 1, -1)
                    vertex.append((0, j1))
                    return j1, 1, False

        def l_next(j: int, ie: int, bw):  # pas[j][ie+1 -- ie] ... i co dalej?
            vertex = poly[-1]
            # punkt (a,j)  powinien już być ostatni na liście vertex
            b, j1 = v_pas[j][ie]  # raczej v_pas_1
            vertex.append((b, j))
            if (b, j) == vertex[1] or the_same((b, j), vertex[1]):
                return -1, b, False
            # to_queue(1-bw, j, ie, 1) - dodane przed wywołaniem l_next
            if b > 0:
                for ix, v_ in enumerate(v_pas[j1]):
                    if v_[1] == j:
                        break
                vertex.append((b, j1))
                if the_same((b, j1), vertex[1]):
                    return -1, b, False  # .........................KONIEC
                if j < j1:
                    to_queue(1 - bw, j1, ix + 1, -1)
                    return j1, ix + 1, False  # w prawo
                else:
                    to_queue(1 - bw, j1, ix - 1, 1)
                    return j1, ix - 1, True  # w lewo
            elif 0 <= j1 < j:
                vertex.append((b, j1))
                if the_same((b, j1), vertex[1]):
                    return -1, b, False  # .........................KONIEC
                to_queue(1 - bw, j1, 1, -1)
                return j1, 1, False  # nawrót w prawo
            else:
                nr = Pas[j, 0]  # pozycja początku
                vertex.append((0, [nr]))
                if nr == lPas[0]:  # skrajny dolny
                    if nr > -(N - 1):
                        vertex.append((0, [-(N - 1)]))
                    return -1, 0, False  # lewy dolny .........................KONIEC
                else:  # jeszcze jest coś niżej
                    ve = np.argwhere(lPas == nr)[0, 0]  # od wierzchołka <ve>
                    nr1 = lPas[ve - 1]
                    vertex.append((0, [nr1]))  # kontynuacja w dół i dalej w prawo
                    if the_same((0, [nr1]), vertex[1]):
                        return -1, 0, False  # .........................KONIEC
                    j1 = l_index[ve - 1]
                    _, j2 = v_pas[j1][0]
                    j1 = max(j1, j2)  # jeśli są 2 krawędzie - górna
                    to_queue(1 - bw, j1, 1, -1)
                    return j1, 1, False  # od brzegu w prawo

      # POJEDYNCZE PASMO  np.array: pas
        # !!! vertex=... tworzy nową zmienną
        # print("N=", N)
        vertex = [0, (0, [-N + 1]), (1, [-N + 1])]  # 0=BIAŁE
        poly.append(vertex)     # ... wstawiamy ją na koniec listy ścian
        nr = rPas[0]
        if nr > -N + 1:
            vertex.append((1, [nr]))
        j = r_index[0]  # pierwsza w lewo (od dołu)
        vertex.append((1, j))
        to_queue(1, j, -2, 1)
        left, ie, bw = True, -2, 0
        while j != -1:  # pierwsza biała ściana od dołu
            j, ie, left = l_next(j, ie, bw)
        while len(queue) > 0:       # while queue
            bw, j, ia, dir = from_queue()
            left = dir < 0
            ie = ia + dir
            b, _ = v_pas[j][ia]
            vertex = [bw, (b, j)]  # był błąd
            poly.append(vertex)
            while j != -1:  # ??
                if left:
                    j, ie, left = l_next(j, ie, bw)
                else:
                    j, ie, left = r_next(j, ie, bw)
        # koniec?
        # print(k, len(poly))
    return Poly, Pasma

def arrea_bw(xp, yp):
    N = int(np.sqrt(xp.shape[0]))

    def toArea(poly, pasmo):  # Poly[i], Pasma[i]
        def to_vertex(u):
            if type(u[1]) == list:
                return [2 * u[0], u[1][0]]
            else:
                a = u[0]  # Fraction
                u0, u1 = pasmo[u[1]]
                return [2 * a, u0 + a * (u1 - u0)]

        def to_area(verts):
            x0, y0 = verts[0]
            v, area = verts[1], 0
            v = v[0] - x0, v[1] - y0
            for w in verts[2:]:
                w = w[0] - x0, w[1] - y0
                area += (v[0] * w[1] - v[1] * w[0]) / 2
                v = w
            return area / 4  # dla siatki o kwaqratach 1x1

        area = 0
        for pol in poly:
            bw = pol[0]
            # polys.append(vertices)
            if bw == 1:
                vertices = [to_vertex(pol[1])]
                for u in pol[2:]:
                    v = to_vertex(u)
                    if v != vertices[-1]:
                        vertices.append(v)
                area += to_area(vertices)
        return area

    Poly, Pasma = analyse(xp, yp)
    pole = 0
    # for i, pasmo in enumerate(Pasma):
    #     area = toArea(Poly[i], pasmo)
    #     # print(f"{i}: ", area)
    #     pole += area
    #  # pole /= 4   już podzielone
    for i, pasmo in enumerate(Pasma[: N//2-1]):
        area = toArea(Poly[i], pasmo)
        # print(f"{i}: ", area)
        pole += area
     # pole /= 4   już podzielone
    pole = 2* pole + toArea(Poly[N//2-1], Pasma[N//2-1])
    return pole

# t0 = time_ns()
for nr in [40264,740264,102424,438169]:          #range(200):
    path = Xp[nr]
    xp = Szach.x(path)
    yp = Szach.y(path)
    area = arrea_bw(xp, yp)
    # if nr % 100 == 0:
    #     print(nr,"\n")
    print(f'{nr:3}', area, '\t', int(area*15))
# t1 = time_ns() - t0
# print(f'Czas : {t1 / 1000_000} ms')  # Czas : 1000: 91489.6371 ms


def toPoly(poly, pasmo, pole=False):    # Poly[i], Pasma[i]
    def to_vertex(u):
        if type(u[1]) == list:
            return [2*u[0], u[1][0]]
        else:
            a = u[0]  # Fraction
            u0, u1 = pasmo[u[1]]
            return [2*a, u0 + a*(u1 - u0)]

    def toFloat(u):
        return [float(u[0]),float(u[1])]

    def to_area(verts):
        x0,y0 = verts[0]
        v, area = verts[1], 0
        v = v[0]-x0, v[1]-y0
        for w in verts[2:]:
            w = w[0] - x0, w[1] - y0
            area += (v[0]*w[1]-v[1]*w[0])/2
            v = w
        return area/4   # dla siatki o kwaqratach 1x1

    polys, f_polys, area  = [], [], 0
    for pol in poly:
        bw = pol[0]
        vertices = [to_vertex(pol[1])]
        for u in pol[2:]:
            v = to_vertex(u)
            if v != vertices[-1]:
                vertices.append(v)
        polys.append(vertices)
        if bw==1:
            f_vertices = np.asarray([toFloat(u) for u in vertices])
            f_polys.append(f_vertices)
            if pole:
                area += to_area(vertices)
    if pole:
        return polys, f_polys, area
    return polys, f_polys

def draw_filledSzach(nr, N=10):
    path = Xp[nr]
    xp = Szach.x(path)
    yp = Szach.y(path)
    Poly, Pasma = analyse(xp, yp)
    fig, ax = plt.subplots()
    ax.plot(xp, yp, color="g")  # coś MUSI być wcześniej narysowane
    pole = 0
    for i, pasmo in enumerate(Pasma):
        polys, f_polys, area = toPoly(Poly[i], pasmo, pole=True)
        print(f"{i}: ", area)
        pole += area
        for j in range(len(f_polys)):
            ppp = plt.Polygon(f_polys[j] + [2*i-N+1, 0], facecolor='gold')  # , edgecolor='0.5')
            ax.add_patch(ppp)
    # pole /= 4   już podzielone
    print(pole, float(pole))
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.title(f'nr = {nr}')
    ax.set_xlabel(f'Area = {pole*15}/15')
    plt.show(block=False)

draw_filledSzach(5)

# for j in range(-10, 12, 2):
#     plt.plot([j,j],[-10,10],color="r")
#     plt.plot([-10,10],[j,j],color="r")
#
# for j in [-9,-7,-5]:
#     plt.plot([j,j],[-9,9],color="k")

import pickle
from time import time_ns
with open(".\\areas100K.pkl", "br") as f:
    Ar = pickle.load(f)
print(Ar.min(),Ar.max()) # 340  822
with open(".\\areas10x10.pkl", "br") as f:
    AR = pickle.load(f)
print(AR.min(),AR.max()) # 340  822

A2 = np.hstack([Ar,np.zeros((100000,), dtype=np.int32)])
A3 = np.hstack([A2,np.zeros((100000,), dtype=np.int32)])
A7 = np.hstack([A3,np.zeros((400000,), dtype=np.int32)])
AR = np.hstack([A7,np.zeros((131804,), dtype=np.int32)])
t0 = time_ns()      # 10000: 539123.191 ms (po skróceniu/2)... 10**5: 4915959.4878 ms
                    # 10**5..2*10**5: Czas : 6112.640 s ..3*10**5: Czas : 5554.306 s
Szach.N = 10
for nr in range(1):     # 700000,831804):          #:[40264,740264,102424,438169] (531634 167/3 	 835)
    path = Xp[nr]
    xp = Szach.x(path)
    yp = Szach.y(path)
    area = arrea_bw(xp, yp)
    # if nr % 100 == 0:
    #     print(nr,"\n")
    if area*15 >800 or area*15 <370:
        print(f'{nr:3}', area, '\t', int(area*15))
    AR[nr] = int(area*15)
t1 = time_ns() - t0
print(f'Czas : {t1 / 10**9:6.3f} s')  # Czas : 1000: 91.489 s


# with open(".\\areas10x10.pkl", "bw") as f:
#     pickle.dump(Ar, f)
# print(Ar.min(),Ar.max())

from sklearn.linear_model import LinearRegression
lreg = LinearRegression()
lreg.fit(Xp,Ar)
r2 = lreg.score(Xp,Ar)      # 0.004988885 !!!!!!!!!!!!!!!!!!
