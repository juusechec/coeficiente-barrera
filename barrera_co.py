# Cristian David Aguazaco Rodriguez - 20121005083
# Cristian David Vargas Sierra - 20121005103

# Coeficientes de transmision y reflexion para una barrera de potencial

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import sys
from scipy.integrate import *
from scipy.constants import *
from matplotlib.widgets import Slider, Button, RadioButtons

sem = 0
num_sam = 1999
Vo = 3.0
a = 0.3e-18
fig, ax = plt.subplots(figsize=(20, 10))

lastKey = None
def graficar():
    global Vo, a, sem, num_sam, point, point1, fig, ax

    font1 = {'family': 'century', 'color': 'black',
             'weight': 'normal', 'size': 14, }
    font2 = {'family': 'century', 'color': '#138706',
             'weight': 'normal', 'size': 16, }
    font3 = {'family': 'century', 'color': '#138706',
             'weight': 'normal', 'size': 18, }


    E = np.linspace(0, 5 * Vo, 2000)

    k1a = np.sqrt(2 * m_e * E[0:400]) / hbar
    k2a = np.sqrt(2 * m_e * (Vo - E[0:400])) / hbar
    p1 = ((k2a**2) + (k1a**2)) / (2 * k1a * k2a)
    q1 = (np.exp(k2a * a) - np.exp(-k2a * a)) / 2
    R1 = ((p1**2) * (q1**2)) / (1 + ((p1**2) * (q1**2)))
    T1 = 1 - R1

    k1b = np.sqrt(2 * m_e * E[400:2000]) / hbar
    k2b = np.sqrt(2 * m_e * (E[400:2000] - Vo)) / hbar
    p2 = ((k2b**2) - (k1b**2))
    q2 = np.sin(k2b * a)
    R2 = ((p2**2) * (q2**2)) / ((4 * (k1b**2) * (k2b**2)) + ((p2**2) * (q2**2)))
    T2 = 1 - R2

    R = list(R1) + list(R2)
    T = list(T1) + list(T2)


    plt.axvline(Vo, color='k')
    plt.axvline(Vo + a, color='k')
    plt.plot(E, R, 'b')
    plt.plot(E, T, 'g')
    plt.title('Barrera de potencial', fontdict=font2, color='b')
    plt.xlabel('Potencial (eV)', fontdict=font2)
    plt.xlim(0, 5 * Vo)
    plt.ylim(-0.2, 1.2)
    plt.text(6, 0.6, 'Presione "x" para mover el cursor hacia adelante.')
    plt.text(6, 0.55, 'Presione "z" para mover el cursor hacia atras.')
    plt.text(6, 0.50, 'Presione "2","3" o "4" para cambiar el Vo.')

    plt.text(6, 0.9, 'Coeficiente de Transmision')
    plt.text(6, 0.05, 'Coeficiente de Reflexion')
    plt.grid(True)

    msg = r'$E = %f eV,$ $R = %f $' % (E[sem], R[sem])
    msg1 = r'$E = %f eV,$ $T = %f $' % (E[sem], T[sem])
    point, = plt.plot(E[0], R[0], "ro", color='b', marker='o')
    point1, = plt.plot(E[0], T[0], "ro", color='g', marker='o')

    plt.legend([point, point1], [msg, msg1])

    def eventoClickB(num):
        global Vo, lastKey
        if lastKey != num:
            lastKey = num
            Vo = num
            print('hola ana')
            plt.clf()
            #plt.cla()
            #plt.gcf().clear()
            #plt.close()
            graficar()

    def limites(tecla):
        global sem
        if tecla == 'x':
            if sem == num_sam - 1:
                return sem
            sem += 2
        if tecla == 'z':
            if sem == 0:
                return sem
            sem -= 2
        if tecla == '2':
            eventoClickB(2)
        if tecla == '3':
            eventoClickB(3)
        if tecla == '4':
            eventoClickB(4)
        return sem


    def press(event):
        global point, point1
        sys.stdout.flush()
        numero = limites(event.key)

        msg = r'$E = %f eV,$ $R = %f $' % (E[numero], R[numero])
        msg1 = r'$E = %f eV,$ $T = %f $' % (E[numero], T[numero])

        point.set_ydata(R[numero])
        point.set_xdata(E[numero])

        point1.set_ydata(T[numero])
        point1.set_xdata(E[numero])

        fig.canvas.draw_idle()

        plt.legend([point, point1], [msg, msg1])

    fig.canvas.mpl_connect('key_press_event', press)


if __name__=='__main__':
    graficar()
    plt.show()
