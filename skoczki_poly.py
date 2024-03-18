import matplotlib.pyplot as plt

from skoczki import Szach
import numpy as np
from fractions import Fraction as Fra


# a6=Szach(6,save=True)
def szach(nr, N=6):
    a_ = Szach(N)
    path = a_.read(nr)
    xp = Szach.x(path)
    yp = Szach.y(path)
    return xp,yp

ab = Fra(3,7)
float(ab)

# z list path, xp, yp wydobyć
# krawędzie lub połówki zawarte w pasach:  [-5, -3] i [-3, -1]


def pasma(xp,yp):
    if len(xp)==37: N = 6
    else: N = 10
    Pasma = [[] for i in range(N-1)]  # [] * (N-1) jest błędne -> tworzy kopie tego samego
    for j in range(N*N):
        ya, yb = yp[j], yp[j+1]
        xa, xb = xp[j], xp[j+1]
        if xa > xb:
            xa, xb, ya, yb = xb, xa, yb, ya
        i = (xa+N)//2
        if xb-xa == 2:
            Pasma[i].append([ya,yb])
        else:
            yab = (ya+yb)//2
            Pasma[i].append([ya,yab])
            Pasma[i+1].append([yab,yb])

    for pas in Pasma:     # lista krawędzi w ramach kolejnych pasm
        pas.sort()
    return [np.asarray(pas) for pas in Pasma]

xp, yp = szach(5,6)
Pasma = pasma(xp, yp)
# pas = Pasma[0]
# v_pas = [[] for _ in pas]
# v = {(0,1)}
# for i in range(len(pas)):
#     a, b = pas[i]
#     j = i+1
#     for j in range(i+1, len(pas)):
#         la, lb = pas[j]
#         da = la-a
#         if da == 0:
#             # print(f"{i,j}:0")
#             v_pas[i].append([Fra(0),(i,j)])
#             v_pas[j].append([Fra(0),(i,j)])
#             v.add((i,j))
#             continue
#         db = b-lb
#         if db >= 0:
#             fra = Fra(da,da+db)
#             # print(f"{i,j}: {fra}")
#             v_pas[i].append([fra,(i,j)])
#             v_pas[j].append([fra,(i,j)])
#             v.add((i,j))
#             continue
#         if la > b+4:
#             j = len(pas)
# for vp in v_pas:
#     vp.sort()
# N=6
Pas0 = Pasma[0]
Pas1 = Pasma[1]
# r0 = np.bincount(Pas1[:, 0]+(N-1))  # krotności lewych końców
# r1 = np.bincount(Pas1[:, 1]+(N-1), minlength=2*N-1)  # i ... prawych

# zmiana koncepcji - tylko krawędź przecinająca + [-1] na wolnych końcach 0/1
v_pas_0 = [[] for _ in Pas0]
v_0 = set()     # ewentualnie {(0,1)}
for i in range(len(Pas0)):
    a, b = Pas0[i]
    j = i+1
    for j in range(i+1, len(Pas0)):
        la, lb = Pas0[j]
        da = la-a
        if da == 0:
            # print(f"{i,j}:0")
            v_pas_0[i].append([Fra(0), j])
            v_pas_0[j].append([Fra(0), i])
            v_0.add((i, j))
            continue
        db = b-lb
        if db >= 0:
            fra = Fra(da, da+db)
            # print(f"{i,j}: {fra}")
            v_pas_0[i].append([fra, j])
            v_pas_0[j].append([fra, i])
            v_0.add((i, j))
            continue
        if la > b+4:
            j = len(Pas1)
for vp in v_pas_0:
    vp.sort()
    if len(vp)==0 or vp[-1][0]<1:
        vp.append([Fra(1), -1])
    if vp[0][0]>0:
        vp.insert(0, [Fra(0), -1])

v_pas_1 = [[] for _ in Pas1]
v_1 = {(0,1)}
for i in range(len(Pas1)):
    a, b = Pas1[i]
    j = i+1
    for j in range(i+1, len(Pas1)):
        la, lb = Pas1[j]
        da = la-a
        if da == 0:
            # print(f"{i,j}:0")
            v_pas_1[i].append([Fra(0), j])
            v_pas_1[j].append([Fra(0), i])
            v_1.add((i, j))
            continue
        db = b-lb
        if db >= 0:
            fra = Fra(da, da+db)
            # print(f"{i,j}: {fra}")
            v_pas_1[i].append([fra, j])
            v_pas_1[j].append([fra, i])
            v_1.add((i, j))
            continue
        if la > b+4:
            j = len(Pas1)
for vp in v_pas_1:
    vp.sort()
    if len(vp)==0 or vp[-1][0]<1:
        vp.append([Fra(1), -1])
    if vp[0][0]>0:
        vp.insert(0, [Fra(0), -1])
# krawędź przecinająca + [-1] na wolnych końcach 0/1

# v_pas_1[0][-1]  # [(1/1), (0, 3)] - wspólny koniec 1[0]=1[3]
# v_pas_1[:7] .............
# 0:[[[(0/1),-1], [(1/2), 1], [(1/1), 3]],
# 1: [[(0/1), 2], [(1/2), 0], [(1/1),-1]],
# 2: [[(0/1), 1], [(1/2), 3], [(1/1),-1]],
# 3: [[(0/1), 4], [(1/2), 2], [(1/1), 0]],
# 4: [[(0/1), 3], [(3/5), 5], [(1/1),-1]],
# 5: [[(0/1),-1], [(3/5), 4], [(1/1),-1]],
# 6: [[(0/1),-1], [(1/1),-1]]]
# np.unique(Pas1[:, 1],return_index=True,return_inverse=True,return_counts=True)

# Lpas_1, l_index = np.unique(Pas1[:, 0], return_index=True)
# [-5,-4,-3, 0, 3, 4, 5],
# [ 0, 1, 3, 5, 6, 7, 9]            # początek której krawędzi
# Rpas_1, r_index, r_counts = np.unique(Pas1[:, 1],
#                                       return_index=True,return_counts=True)
# [-5,-4,-3,-1, 1, 2, 3, 4, 5],
# [ 1, 0, 2, 5, 4, 6, 7, 9, 8]      # koniec której krawędzi
# [ 1, 2, 1, 1, 1, 1, 1, 1, 1]      # koniec pojedynczy, czy podwójny
# "na start:"  (0[-1]--1[-1])
# i0, r_i = r_index[0], 0     # -5=Rpas_1[r_i] jest końcem krawędzi i0,
                            # jedynym: r_counts[r_i]=1
# print(v_pas_1[i0][-1])          # koniec:  [(1/1), -1]
# albo  1[-1]==1[i0],
# albo dołożyć  1[-1]--1[i0]    # koniec:  [(1/1), i0]
#      doszliśmy do punktu 1[i0] = R_pas_1[0], nr r_i=0
# a1, i1 = v_pas_1[i0][-2]        # poprzedni punkt na  i0: ((1/2), 0)
# pos = len(v_pas_1[i0])-2
# dołożyć  1[i0]--a1[i0]
# --> "do zbadania nast. ściany": (a1[i0]--1[i0])
# ----------- PĘTLA --------------------------------------------
# jeśli  a1==0:
#   zamknąć ścianę, gdy  0[i0]==0[-1],
#   a jeśli nie: dołożyć  0[i0]--0[-1]
# jeśli  a1>0,  kontynuacja wzdłuż  i1:
# kontrolnie: czy  v_pas_1[i1][pos]==[a1,i0]?
# pos -= 1
# a2, i2 = v_pas_1[i1][pos]       # poprzedni punkt na  i1: ((1/2), 0)
# dołożyć  a1[i1]--a2[i1]
# --> "do zbadania nast. ściany": (a2[i1]--a1[i1])
# a1, i1, i0 = a2, i2, i1
# --- powrót do PĘTLI --------------------------

# r_ix = np.argwhere(Rpas_1==-4)[0,0]     # która krawędź ma koniec -4
# print(r_counts[r_ix])                   # czy jest jedyna?

# Jak pamiętać ścianę i jej wierzchołki?
N = 6
# rPas = Rpas_1
# vertex = [0, (0, [-N+1]), (1, [-N+1])]  # 0=BIAŁE
poly = []  # 1. ściana - do dokończenia / potem kolejne
queue = []
v_pas = v_pas_1   # pasmo 1
Pas = Pas1


rPas, r_index = np.unique(Pas0[:, 1], return_index=True)
lPas, l_index = np.unique(Pas0[:, 0],  return_index=True)
poly = []
queue = []
v_pas = v_pas_0   # pasmo 0
Pas = Pas0


def to_queue(bw, j, ia, dir):           # bw=0/1
    for i in range(len(queue)):
        if queue[i] == (1-bw, j, ia+dir, -dir):
            queue.__delitem__(i)             # usunąć, jeśli już było, dir = +/-1
            # print("removed:",queue)
            return
    # print(queue)
    queue.append((bw, j, ia, dir))      # j: ia -- ia+dir


def from_queue():
    return queue.pop(0)     # pobiera i usuwa


def the_same(u,v):
    # print("u:",u)
    # print("v:",v)
    if u[0] != v[0]:
        return False
    if type(u[1])==list:
        n1 = u[1][0]
        if type(v[1])==list:
            n2 = v[1][0]
            return n1==n2
        else:
            j = v[1]
            return n1==Pas[j,u[0]]
    elif type(v[1])==list:
        return the_same(v,u)
    else:
        a = u[0]    # = v[0] - Fraction
        u0,u1 = Pas[u[1]]
        v0,v1 = Pas[v[1]]
        y = u0 + a*(u1-u0)
        return y == v0 + a*(v1-v0)


def startPas():     # POJEDYNCZE PASMO  np.array: pas
    # vertex.append((0, [-N+1]))... (1, [-N+1]))   !!! vertex=... tworzyło nową zmienną
    vertex = [0, (0, [-N + 1]), (1, [-N + 1])]  # 0=BIAŁE
    poly.append(vertex)
    nr = rPas[0]
    if nr > -N+1:
        vertex.append((1, [nr]))
    j = r_index[0]      # pierwsza w lewo (od dołu)
    vertex.append((1, j))
    to_queue(1, j, -2, 1)
    left, ie, bw = True, -2, 0
    while j!=-1:        # pierwsza biała ściana od dołu
        j, ie, left = l_next(j, ie, bw)
    # print(queue)
    while len(queue)>0:
        bw, j, ia, dir = from_queue()
        # print(bw, j, ia, dir)
        left = dir<0
        ie = ia+dir
        b,_ = v_pas[j][ia]
        vertex = [bw, (b, j)]       # był błąd
        poly.append(vertex)
        # print(poly)
        while j!=-1:    # ??
            if left:
                j, ie, left = l_next(j, ie, bw)
            else:
                j, ie, left = r_next(j, ie, bw)
            # if len(poly)>6:
            #     return
    # koniec?


def r_next(j: int, ie: int, bw):    # pas[j][ie-1 -- ie]  ... i co dalej? bw=0/1
    vertex = poly[-1]
    # punkt (a,j) ~ v_pas[j][ie-1] powinien już być ostatni na liście vertex
    b, j1 = v_pas[j][ie]    # raczej v_pas_1
    vertex.append((b,j))
    if the_same((b,j),vertex[1]):
        return -1, b, False
    # to_queue(1-bw, j, ie, -1)   # pominąć, mogło przyjść z kolejki
    if j<j1:    # nawrót - w lewo w górę
        for ix, v_ in enumerate(v_pas[j1]):
            if v_[1]==j:
                break
        to_queue(1-bw, j1, ix-1, 1)
        vertex.append((b,j1))       # ten sam koniec, inaczej
        if the_same((b,j1), vertex[1]):
            return -1, b, False  # .........................KONIEC
        return j1, ix-1, True       # b, l_next...
    if b<1:     # w prawo w górę
        for ix, v_ in enumerate(v_pas[j1]):
            if v_[1]==j:
                break
        to_queue(1-bw, j1, ix+1, -1)
        vertex.append((b,j1))
        if the_same((b,j1), vertex[1]):
            return -1, b, False  # .........................KONIEC
        return j1, ix+1, False      # b, r_next
    else:       # b=1, w górę, jeśli można
        nr = Pas[j,1]       # pozycja końca
        vertex.append((1,[nr]))     # ten sam koniec, w pionie
        if nr < rPas[-1]:   # jeszcze jest coś wyżej
            ve = np.argwhere(rPas==nr)[0,0]     # od wierzchołka <ve>
            nr1 = rPas[ve+1]
            vertex.append((1,[nr1]))            # kontynuacja w górę
            if the_same((1,[nr1]), vertex[1]):
                return -1, 1, False
            # i dalej w lewo
            j1 = r_index[ve+1]
            to_queue(1-bw, j1, -2, 1)
            vertex.append((1, j1))
            return j1, -2, True     # od brzegu w lewo
        else:
            if nr < N-1:
                vertex.append((1, [N-1]))
            vertex.append((0, [N-1]))       # górna krawędź w poziomie
            nr = lPas[-1]
            if nr < N-1:
                vertex.append((0, [nr]))
            if the_same((0,[nr]), vertex[1]):
                return -1, 0, False
            j1 = l_index[-1]
            _, j2 = v_pas[j1][0]
            j1 = max(j1,j2)         # jeśli są 2 krawędzie - górna
            to_queue(1-bw, j1, 1, -1)
            vertex.append((0, j1))
            return j1, 1, False


def l_next(j: int, ie: int, bw):    # pas[j][ie+1 -- ie] ... i co dalej?
    vertex = poly[-1]
    # punkt (a,j)  powinien już być ostatni na liście vertex
    b, j1 = v_pas[j][ie]    # raczej v_pas_1
    vertex.append((b,j))
    if (b,j)==vertex[1] or the_same((b,j),vertex[1]):
        return -1, b, False
    # to_queue(1-bw, j, ie, 1) - dodane przed wywołaniem l_next
    if b>0:
        for ix, v_ in enumerate(v_pas[j1]):
            if v_[1] == j:
                break
        vertex.append((b, j1))
        if the_same((b,j1), vertex[1]):
            return -1, b, False  # .........................KONIEC
        if j < j1:
            to_queue(1-bw, j1, ix+1, -1)
            return j1, ix+1, False  # w prawo
        else:
            to_queue(1-bw, j1, ix-1, 1)
            return j1, ix-1, True   # w lewo
    elif 0<=j1<j:
        vertex.append((b, j1))
        if the_same((b,j1), vertex[1]):
            return -1, b, False  # .........................KONIEC
        to_queue(1-bw, j1, 1, -1)
        return j1, 1, False             # nawrót w prawo
    else:
        nr = Pas[j, 0]                  # pozycja początku
        vertex.append((0, [nr]))
        if nr==lPas[0]:                 # skrajny dolny
            if nr> -(N-1):
                vertex.append((0, [-(N-1)]))
            return -1, 0, False         # lewy dolny .........................KONIEC
        else:               # jeszcze jest coś niżej
            ve = np.argwhere(lPas==nr)[0,0]   # od wierzchołka <ve>
            nr1 = lPas[ve-1]
            vertex.append((0, [nr1]))   # kontynuacja w dół i dalej w prawo
            if the_same((0, [nr1]),vertex[1]):
                return -1, 0, False     #  .........................KONIEC
            j1 = l_index[ve-1]
            _, j2 = v_pas[j1][0]
            j1 = max(j1,j2)             # jeśli są 2 krawędzie - górna
            to_queue(1-bw, j1, 1, -1)
            return j1, 1, False         # od brzegu w prawo


# startPas()


# [[0, (0, [-5]), (1, [-5]), (1, 1)], [0, (0, [-5]), (1, [-5]), (1, 1)], [0, (0, [-5]), (1, [-5]), (1, 1)], [1, (1, -2)]]

def analyse(xp,yp):
    if len(xp)==37: N=6
    else: N=10

    def pasma(xp, yp):
        if len(xp) == 37:
            N = 6
        else:
            N = 10
        Pasma = [[] for i in range(N - 1)]  # [] * (N-1) jest błędne -> tworzy kopie tego samego
        # print("xp=",xp)
        # print("N=",N)
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
        while len(queue) > 0:
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

xp, yp = szach(5,6)
Poly, Pasma = analyse(xp, yp)

def toVertex(u, pasmo):
    if type(u[1]) == list:
        return [2*u[0], u[1][0]]
    else:
        a = u[0]  # Fraction
        u0, u1 = pasmo[u[1]]
        return [2*a, u0 + a*(u1 - u0)]


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

# i=1
# poly, pasmo = Poly[i], Pasma[i]
# N = len(Pasma)+1
# # polys, f_polys = toPoly(poly, pasmo)
# # f_polys[:3]
# # f_polys[1]+[1,0]
# fig, ax = plt.subplots()
# ax.plot(xp,yp,color="g") # coś MUSI być wcześniej narysowane
# for i, pasmo in enumerate(Pasma):
#     polys, f_polys = toPoly(Poly[i], pasmo)
#     for j in range(len(f_polys)):
#         ppp = plt.Polygon(f_polys[j]+[2*i-N+1,0], facecolor='gold')  # , edgecolor='0.5')
#         ax.add_patch(ppp)
# ax.axis('equal')
# ax.set_xticks([])
# ax.set_yticks([])
# plt.show(block=False)

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
    for i, pasmo in enumerate(Pasma[: N//2-1]):
        area = toArea(Poly[i], pasmo)
        # print(f"{i}: ", area)
        pole += area
     # pole /= 4   już podzielone
    pole = 2* pole + toArea(Poly[N//2-1], Pasma[N//2-1])
    return pole

def draw_filledSzach(nr, N=6):
    xp, yp = szach(nr, N)
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

draw_filledSzach(5,6)
draw_filledSzach(6,6)
draw_filledSzach(8,6)

draw_filledSzach(118,10)        # 766/15   766
draw_filledSzach(40264,10)      # 382/15 	 382
draw_filledSzach( 4300,10)      # 398/15 	 398
draw_filledSzach( 8182,10)      # 382/15 	 382
draw_filledSzach( 14817,10)     #  163/3 	 815
draw_filledSzach(26614, 10)     # 361/15 	 361
draw_filledSzach(30195, 10)     # 371/3 	 355
draw_filledSzach(40263,10)      # 68/3 	 340 = sym(102424)
draw_filledSzach(63507,10)      # 274/5 	 822
draw_filledSzach(125759,10)      # 55 	 825
draw_filledSzach(125989,10)      # 55 	 825
draw_filledSzach(126462,10)      # 55 	 825
draw_filledSzach(740264,10)
draw_filledSzach(102424,10)     # 68/3 	 340
draw_filledSzach(438169,10)     # 167/3 	 835
draw_filledSzach(531634,10)     # 167/3 	 835


M6, M10 = 10, 831804

# sz = Szach(10)
sz = Szach(10)
X = Szach.read_data()
from time import time_ns


def save_paths():
    N=10
    Xp = np.ndarray((M10, N*N+1), dtype=np.int16)

    t0 = time_ns()
    for nr in range(M10):
        seq = X[nr]
        path = Xp[nr]
        m = N * N // 4
        for i in range(m):
            path[i] = seq[i]
            path[2 * m + i] = N * N - seq[i] - 1
        zx, zy = -1, 1
        if seq[N * N // 4 - 1] > N * N // 2:
            zx, zy = 1, -1
        for i in range(m):
            x, y = Szach.x(seq[i]), Szach.y(seq[i])
            path[m + i] = Szach.nr(zy * y, zx * x)
            path[3 * m + i] = Szach.nr(-zy * y, -zx * x)
        path[N * N] = path[0]

    t1 = time_ns() - t0
    print(f'Czas : {t1 / 1000_000} ms')

    import pickle
    with open(".\\skoczki10x10.pkl","bw") as f:
        pickle.dump(Xp,f)

import pickle
with open(".\\skoczki10x10.pkl","br") as f:
    Xp = pickle.load(f)

Ar = np.zeros((100000,), dtype=np.int32)
A2 = np.hstack([Ar,np.zeros((100000,), dtype=np.int32)])
t0 = time_ns()      # 10000: 539123.191 ms (po skróceniu/2)... 10**5: 4915959.4878 ms
Szach.N = 10
for nr in range(100000,200000):          #:[40264,740264,102424,438169]
    path = Xp[nr]
    xp = Szach.x(path)
    yp = Szach.y(path)
    area = arrea_bw(xp, yp)
    # if nr % 100 == 0:
    #     print(nr,"\n")
    if area*15 >800 or area*15 <400:
        print(f'{nr:3}', area, '\t', int(area*15))
    A2[nr] = int(area*15)
t1 = time_ns() - t0
print(f'Czas : {t1 / 1000_000} ms')  # Czas : 1000: 91489.6371 ms

#   4300 398/15 	 398
#   8182 398/15 	 398
#   8464 26 	     390


# 1363 262/5 	 786
# 1654 83/3 	 415
# 4000 791/15 	 791
# 4279 139/5 	 417
# 4300 398/15 	 398
# 6111 139/5 	 417
# 6417 412/15 	 412
# 7716 419/15 	 419
# 8001 416/15 	 416
# 8182 398/15 	 398
# 8256 80/3 	 400
# 8464 26 	     390
# 8497 406/15 	 406
# 10491 406/15 	 406
# 10519 418/15 	 418
# 10522 416/15 	 416
# 10534 413/15 	 413
# 10541 418/15 	 418
# 10800 133/5 	 399
# 12941 80/3 	 400
# 14122 416/15 	 416
# 14155 80/3 	 400
# 14655 157/3 	 785
# 14678 163/3 	 815
# 14679 811/15 	 811
# 14793 811/15 	 811
# 14817 163/3 	 815
# 14819 157/3 	 785
# 15155 137/5 	 411
# 15253 139/5 	 417
# 15338 398/15 	 398
# 15341 401/15 	 401
# 15351 137/5 	 411
# 15406 133/5 	 399
# 15416 419/15 	 419
# 15424 412/15 	 412
# 15597 419/15 	 419
# 15676 389/15 	 389
# 15807 137/5 	 411
# 16121 401/15 	 401
# 16248 409/15 	 409
# 16403 137/5 	 411
# 16404 82/3 	 410
# 16849 419/15 	 419
# 17039 83/3 	 415
# 17040 76/3 	 380
# 17384 403/15 	 403
# 17402 401/15 	 401
# 17639 79/3 	 395
# 18523 137/5 	 411
# 18557 419/15 	 419
# 18560 391/15 	 391
# 18566 394/15 	 394
# 19182 419/15 	 419
# 19437 138/5 	 414
# 19523 137/5 	 411
# 19584 137/5 	 411
# 19679 409/15 	 409
# 19780 262/5 	 786
# 19952 781/15 	 781

#   0 172/5 	 516
#   1 179/5 	 537
#   2 574/15 	 574
#   3 106/3 	 530
#   4 118/3 	 590
#   5 602/15 	 602
#   6 629/15 	 629
#   7 218/5 	 654
#   8 643/15 	 643
#   9 683/15 	 683
#  10 581/15 	 581
#  11 592/15 	 592
#  12 583/15 	 583
#  13 508/15 	 508
#  14 184/5 	 552
#  15 526/15 	 526
#  16 177/5 	 531
#  17 517/15 	 517
#  18 216/5 	 648
#  19 44 	     660
#  20 229/5 	 687
#  21 712/15 	 712
#  22 206/5 	 618
#  23 216/5 	 648
#  24 193/5 	 579
#  25 202/5 	 606
#  26 622/15 	 622
#  27 199/5 	 597
#  28 199/5 	 597
#  29 37 	     555
#  30 623/15 	 623
#  31 130/3 	 650
#  32 677/15 	 677
#  33 228/5 	 684
#  34 653/15 	 653
#  35 481/15 	 481
#  36 103/3 	 515
#  37 35 	     525
#  38 478/15 	 478
#  39 159/5 	 477
#  40 532/15 	 532
#  41 568/15 	 568
#  42 107/3 	 535
#  43 551/15 	 551
#  44 169/5 	 507
#  45 547/15 	 547
#  46 512/15 	 512
#  47 169/5 	 507
#  48 598/15 	 598
#  49 529/15 	 529
#  50 203/5 	 609
#  51 484/15 	 484
#  52 673/15 	 673
#  53 553/15 	 553
#  54 589/15 	 589
#  55 611/15 	 611
#  56 556/15 	 556
#  57 572/15 	 572
#  58 176/5 	 528
#  59 106/3 	 530
#  60 198/5 	 594
#  61 208/5 	 624
#  62 613/15 	 613
#  63 662/15 	 662
#  64 643/15 	 643
#  65 538/15 	 538
#  66 593/15 	 593
#  67 199/5 	 597
#  68 209/5 	 627
#  69 613/15 	 613
#  70 623/15 	 623
#  71 202/5 	 606
#  72 202/5 	 606
#  73 208/5 	 624
#  74 629/15 	 629
#  75 628/15 	 628
#  76 37 	     555
#  77 191/5 	 573
#  78 557/15 	 557
#  79 554/15 	 554
#  80 589/15 	 589
#  81 142/3 	 710
#  82 216/5 	 648
#  83 187/5 	 561
#  84 224/5 	 672
#  85 186/5 	 558
#  86 136/3 	 680
#  87 593/15 	 593
#  88 629/15 	 629
#  89 689/15 	 689
#  90 136/3 	 680
#  91 692/15 	 692
#  92 719/15 	 719
#  93 136/3 	 680
#  94 127/3 	 635
#  95 197/5 	 591
#  96 204/5 	 612
#  97 707/15 	 707
#  98 628/15 	 628
#  99 229/5 	 687
# 100 248/5 	 744
# 101 41 	     615
# 102 212/5 	 636
# 103 698/15 	 698
# 104 136/3 	 680
# 105 43 	     645
# 106 556/15 	 556
# 107 189/5 	 567
# 108 119/3 	 595
# 109 40 	     600
# 110 596/15 	 596
# 111 608/15 	 608
# 112 644/15 	 644
# 113 112/3 	 560
# 114 208/5 	 624
# 115 584/15 	 584
# 116 679/15 	 679
# 117 127/3 	 635
# 118 766/15 	 766
# 119 704/15 	 704
# 120 623/15 	 623
# 121 219/5 	 657
# 122 622/15 	 622
# 123 193/5 	 579
# 124 583/15 	 583
# 125 196/5 	 588
# 126 131/3 	 655
# 127 127/3 	 635
# 128 649/15 	 649
# 129 45 	     675
# 130 143/3 	 715
# 131 631/15 	 631
# 132 542/15 	 542
# 133 649/15 	 649
# 134 44 	     660
# 135 526/15 	 526
# 136 613/15 	 613
# 137 613/15 	 613
# 138 107/3 	 535
# 139 37 	     555
# 140 556/15 	 556
# 141 592/15 	 592
# 142 176/5 	 528
# 143 572/15 	 572
# 144 556/15 	 556
# 145 192/5 	 576
# 146 214/5 	 642
# 147 569/15 	 569
# 148 202/5 	 606
# 149 562/15 	 562
# 150 662/15 	 662
# 151 679/15 	 679
# 152 41 	     615
# 153 127/3 	 635
# 154 133/3 	 665
# 155 207/5 	 621
# 156 613/15 	 613
# 157 208/5 	 624
# 158 529/15 	 529
# 159 183/5 	 549
# 160 206/5 	 618
# 161 218/5 	 654
# 162 629/15 	 629
# 163 207/5 	 621
# 164 637/15 	 637
# 165 593/15 	 593
# 166 116/3 	 580
# 167 229/5 	 687
# 168 41 	     615
# 169 36 	     540
# 170 724/15 	 724
# 171 134/3 	 670
# 172 662/15 	 662
# 173 587/15 	 587
# 174 38 	     570
# 175 614/15 	 614
# 176 202/5 	 606
# 177 571/15 	 571
# 178 578/15 	 578
# 179 188/5 	 564
# 180 557/15 	 557
# 181 171/5 	 513
# 182 174/5 	 522
# 183 539/15 	 539
# 184 110/3 	 550
# 185 197/5 	 591
# 186 547/15 	 547
# 187 581/15 	 581
# 188 608/15 	 608
# 189 208/5 	 624
# 190 182/5 	 546
# 191 601/15 	 601
# 192 526/15 	 526
# 193 206/5 	 618
# 194 158/5 	 474
# 195 197/5 	 591
# 196 547/15 	 547
# 197 569/15 	 569
# 198 623/15 	 623
# 199 611/15 	 611

import io
N=10
import os
size = os.path.getsize(f"skoczki{N}x{N}.bin")

with io.open(f"skoczki{N}x{N}.bin", "rb") as plik:
    plik.seek(0, os.SEEK_END)
    ssizz = plik.tell()
    # res = plik.read()
print(size)

import os
with  open(r'C:\mydata\testdata2.txt') as size:
    size.seek(0, os.SEEK_END)
    print('File size:', size.tell(), 'bytes')

import urllib.request
file=urllib.request.urlopen("https://sample-videos.com/xls/Sample-Spreadsheet-10-rows.xls")
print("File size:", file.length, "bytes")

import pickle

with open(".\\areas100K.pkl", "bw") as f:
    pickle.dump(Ar, f)
A2.min() # 340

A2.max() # 822

with open(".\\areas10x10.pkl", "br") as f:
    AR = pickle.load(f)
print(AR.min(),AR.max()) # 340  822

a_num, a_where, a_count = np.unique(AR,return_index=True,return_counts=True)

fig, ax = plt.subplots()
ax.plot(a_num, a_count)     # wygląda na rozkład normalny
plt.show(block=False)

for nr in a_where[:10]:
    draw_filledSzach(nr, 10)    # 503032:   362
for nr in a_where[-10:]:
    draw_filledSzach(nr, 10)    # 79966:    818

fig, ax = plt.subplots()
ax.plot(Ar[1234:2345], linewidth=.2)   # wygląda losowo,
plt.show(block=False)

from sklearn.linear_model import LinearRegression
lreg = LinearRegression()
lreg.fit(Xp,Ar)
r2 = lreg.score(Xp,Ar)      # 0.004988885 !!!!!!!!!!!!!!!!!! - TO JEST LOSOWOŚĆ
m = Ar.mean()               # 594.74170
std = Ar.std()              #  62.35921
print(m, std)
(Ar<=m-3*std-2).sum()       # Out[57]: 524
(Ar>=m+3*std-2).sum()       # Out[58]: 544

from scipy.stats import normaltest, skewtest, kurtosistest
k2, pval = normaltest(Ar[:10000])   # (63.473491606578094, 1.647807305913591e-14) już NIE!
sk, pval = skewtest(Ar)
kt, pval = kurtosistest(Ar)

from scipy.special import boxcox1p, boxcox, tests

par = [0,.3,.6,.9,1.3]
Res = boxcox(Ar[:,None], par)
for i in range(5):
    print(normaltest(Res[:1000,i]))

