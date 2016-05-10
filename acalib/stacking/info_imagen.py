#This file is part of ChiVO, the Chilean Virtual Observatory
#A project sponsored by FONDEF (D11I1060)
#Copyright (C) 2015 Universidad Tecnica Federico Santa Maria Mauricio Solar
#                                                            Marcelo Mendoza
#                   Universidad de Chile                     Diego Mardones
#                   Pontificia Universidad Catolica          Karim Pichara
#                   Universidad de Concepcion                Ricardo Contreras
#                   Universidad de Santiago                  Victor Parada
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import math

def distancia(x1,y1,x2,y2):
    # Funcion que calcula la distancia entre dos puntos mediante pitagoras
    return math.sqrt(math.pow((x2-x1),2) + math.pow((y2-y1),2))

def centro_geometrico(dis,x1,y1,x2,y2,angulo):
    # Corresponde al centro geometrico de la figura (mitad eje mayor)

    if (angulo < 90):
        punto_central = [round(x1 + 0.5*math.cos(angulo*math.pi/180)*dis ), round(y1 + 0.5*math.sin(angulo*math.pi/180)*dis)]
    if (angulo > 90):
        punto_central = [round(x2 - 0.5*math.cos(angulo*math.pi/180)*dis ), round(y2 - 0.5*math.sin(angulo*math.pi/180)*dis)]
    if (angulo == 90):
        punto_central = [abs(x1-x2)/2,y1]

    return punto_central

def angulo_a_rotar(x1,y1,x2,y2):
    # Indica el angulo a rotar para dejar el eje mayor paralelo al eje x
    if (((y2 - y1)*1.0) == 0):
        return 90
    else:
        return 180 - math.degrees(math.atan(abs(((y2 - y1)*1.0)/((x2 - x1)*1.0))))


def eje_mayor(borde):

    # info_elipse = [Dimension,x1,y1,x2,y2]
    info_elipse = [0,0,0,0,0]

    for i in range (0, len(borde)-2):
        if ( i % 2 != 0):
            continue
        for j in range (i+2, len(borde)-1):
            if ( j % 2 != 0):
                continue
            dis = distancia(borde[i],borde[i+1],borde[j],borde[j+1])
            if dis > info_elipse[0]:
                info_elipse[0] = dis          # Distancia del eje mayor
                info_elipse[1] = borde[i]     # (X1,
                info_elipse[2] = borde[i+1]   #  Y1) Eje mayor
                info_elipse[3] = borde[j]     # (X2,
                info_elipse[4] = borde[j+1]   #  Y2) Eje mayor

    return info_elipse

if __name__ == "__main__":
    vec = [0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    eje_mayor(vec)


